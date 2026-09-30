module SpectraFromScratch

using Statistics
using Distributions
using FFTW
using OffsetArrays

export FourierTransform
export RegularTimeseries
export centered_fft
export centered_ifft
export time_average
export band_average
export confid, total_spectral_energy
export spectral_power_law, spectral_basis
export convolve
export periodogram
export expand
export phase

import Base: /

"""
    FourierTransform{T}

Store a Fourier transform in an efficient way such that
- not all frequencies must be stored but they can be retrieved with `x.freq`
- be sure that the mean and Nyquist frequencies have real coefficients `x.coeff`
- use an `OffsetArray` so that all positive and negative modes can be directly accessed

# Arguments
- `coeff::OffsetArray`: complex coefficients
- `df<: Number`: fundamental frequency
"""
struct FourierTransform{T<:Number, R<:Number}
    coeff::OffsetVector{T}
    df::R
end

Base.propertynames(x::FourierTransform, private::Bool=false) =
      private ? (:freq,  fieldnames(typeof(x))...) : fieldnames(typeof(x))

function Base.getproperty(x::FourierTransform, d::Symbol)
    if d === :freq
        # reconstruct fourier frequencies
        ind = first(axes(x.coeff))
        return OffsetArray(x.df*ind, ind)
    else
        return getfield(x, d)
    end
end

"""
    RegularTimeseries{T<:Number, R <:Number}

Store a uniformly-sampled timeseries in an efficient way such that
- not all times must be stored but they can be retrieved with `x.time`
- be sure that the mean and Nyquist frequencies have real coefficients `x.x`
- use an `OffsetArray` in symmetry with `FourierTransform`

# Fields
- `x::OffsetVector`: timeseries values
- `dt<: Number`: temporal spacing (fixed)
"""
struct RegularTimeseries{T <: Number, R <: Number}
    x::OffsetVector{T}
    dt::R
end
Base.length(y::RegularTimeseries) = length(y.x)

Base.propertynames(x::RegularTimeseries, private::Bool=false) =
      private ? (:time,  fieldnames(typeof(x))...) : fieldnames(typeof(x))

function Base.getproperty(x::RegularTimeseries, d::Symbol)
    if d === :time
        # reconstruct times
        # should this be an OffsetArray?
        ind = first(axes(x.x))
        # return x.dt* ind
        return OffsetArray(x.dt* ind, ind)
        # range(start=x.dt*0, step=x.dt, length=length(x.x)) 
    else
        return getfield(x, d)
    end
end

"""
    function RegularTimeseries(x::AbstractVector, t::AbstractVector)
"""
function RegularTimeseries(x::AbstractVector, t::AbstractVector{R}) where R <: Number
    length(x) != length(t) && error("lengths do not match")
    
    # if all( abs.(diff(diff(t))) .< 1e-12*oneunit(eltype(t)))
    if all( abs.(diff(diff(t))) .< 1e-12*oneunit(R))

        # minimize machine error
        dt = (last(t)-first(t))/(length(t)-1)

        # times an integer multiple of sampling time?
        ind0 = t./dt
        offset = Integer(first(t)./dt)-1

        # force all timeseries to be OffsetArrays in symmetry with `FourierTransform`
        x_offset = OffsetArray(x, offset)
        return RegularTimeseries{eltype(x), eltype(t)}(x_offset, dt)
    else
        error("not evenly spaced")
    end
end


"""
    FrequencySpectrum{T}

One-sided frequency spectrum.

# Fields
- `psd::AbstractVector`: power spectral density
- `freq::AbstractVector`: frequencies
"""
struct FrequencySpectrum{T}
    psd::AbstractVector
    freq::AbstractVector
    function FrequencySpectrum(psd, f)
        isnegative = (x -> x < zero(x))
        if any(isnegative.(f))
            error("one sided spectrum where frequency must be positive ")
        else
            new{eltype(real.(psd))}(real.(psd), f)
        end
    end
end

"""
    fourier_modes(y::RegularTimeseries)
    fourier_modes(y::Number)
"""
fourier_modes(N::Number) = iseven(N) ?
	                   (m = (-convert(Int,N/2):convert(Int,(N/2)-1))) :
	                   (m = (-convert(Int,(N-1)/2):convert(Int,((N-1)/2))))

fourier_modes(y::RegularTimeseries) = fourier_modes(length(y))

fourier_modes(Ψ::FrequencySpectrum; even=true) =
    (even = true) ?
    fourier_modes(2length(Ψ.psd))  :
    fourier_modes(2length(Ψ.psd) + 1)
    
"""
    fourier_frequencies(m, T)
"""
fourier_frequencies(m, T) = OffsetArray(m/T, m)

function fourier_frequencies(y::RegularTimeseries)
    m = fourier_modes(y)
    T = record_length(y)
    #the dimensional frequency scale, this is an "iterator", not a vector, in julia
    return fourier_frequencies(m, T)
end

"""
    sampling_resolution(y::RegularTimeseries) 
"""
sampling_resolution(y::RegularTimeseries) = y.dt

"""
    record_length(y::RegularTimeseries) 
"""
record_length(y::RegularTimeseries) = length(y) * sampling_resolution(y)

"""
    function centered_fft(y::RegularTimeseries)

Computes FFT, with zero frequency in the center, and returns 
dimensional frequency vector.

Adapted from a function written by Quan Quach of blinkdagger.com 
Modified by Tom Farrar, 2016, jfarrar@whoi.edu. Julia version,
<G Jake Gebbie, ggebbie@whoi.edu>, 2021.

# Arguments
- `y::RegularTimeseries`

# Output
- `x̂`::FourierTransform
"""
function centered_fft(y::RegularTimeseries)
    m = fourier_modes(y)
    T = record_length(y) 

    #the dimensional frequency scale, this is an "iterator", not a vector, in julia
    f = fourier_frequencies(m, T)

    # df = fundamental frequency
    df = f[1]
    
    #=swaps the halves of the FFT vector so that 
    the zero frequency is in the center.
    If you are going to compute an IFFT, 
    first use X=ifftshift(X) to undo the shift =#
    x̂ = fftshift(fft(OffsetArrays.no_offset_view(y.x)))
    return FourierTransform(OffsetArray(x̂, m), df)
end

"""
    function centered_ifft(beta::FourierTransform)

Computes inverse FFT

# Output
- `x`::RegularTimeseries
"""
function centered_ifft(beta::FourierTransform)
    y = ifft(ifftshift(OffsetArrays.no_offset_view(beta.coeff)))
    f_nyquist = -beta.df*first(eachindex(beta.coeff))
    dt = 1 / (2*f_nyquist)

    # assume indices start at zero
    return RegularTimeseries( OffsetArray(real.(y), -1), dt)
end

function FourierTransform_manual(y::RegularTimeseries)
    #the dimensional frequency scale, this is an "iterator", not a vector, in julia
    m = fourier_modes(y)
    f = fourier_frequencies(y)
    dt = -1 / (2*f[begin])

    # make a β coefficient for every value of m
    ft_type = eltype(first(y.x)*im)
    # β = OffsetArray(zero(Vector{ComplexF64}(undef, length(y))), m)
    β = OffsetArray(zero(Vector{ft_type}(undef, length(y))), m)

    for m in eachindex(f)
        # check that eachindex correctly pulls indices
        for n in eachindex(y.x)
            # here assuming (n-1) is ok
            β[m] += exp(-2π*im*f[m]*dt*(n-1)) * y.x[n]
        end
    end
    offset = first(eachindex(f)) - 1
    return FourierTransform(OffsetArray(β, offset), f[1])
end

function FourierTransform(y::RegularTimeseries; alg=:centered_fft)
    if alg==:centered_fft
        return centered_fft(y)
    elseif alg==:manual
        return FourierTransform_manual(y)
    else
        error("SpectraFromScratch.jl: incorrect keyword")
    end
end

function RegularTimeseries(y::FourierTransform; alg=:centered_ifft)
    if alg==:centered_ifft
        return centered_ifft(y)
    elseif alg==:manual
        return RegularTimeseries_manual(y)
    else
        error("SpectraFromScratch.jl: incorrect keyword")
    end
end

Base.length(x::FourierTransform) = 1

"""
    expand(t, beta::FourierTransform)

Expand the complex exponentials at time,
t, the time elapsed from record start, t=0.
"""
function expand(t, beta::FourierTransform{C, T}) where {C, T}
    N = length(beta.coeff) # number of observations
    # y = 0 * real(first(beta.coeff))
    y = zero(C)
    for n in eachindex(beta.coeff)
        # assume time starts at zero
        y += expand(t, n, beta)
        # y += real.(beta.coeff[j] * exp(2π*im*beta.df*j*t))
    end
    abs(imag(y)) > 1e-10*oneunit(eltype(imag(y))) && println("note: imaginary =", imag(y))
    return real(y) 
end

expand(t::Number, n::Number, beta::FourierTransform) =
    beta.coeff[n] * exp(2π*im*beta.df*n*t) / length(beta.coeff)

"""
derivative(t, beta::FourierTransform)

t is the time elapsed from record start, t=0
"""
function derivative(t, beta::FourierTransform{C, T}) where {C, T}
    N = length(beta.coeff) # number of observations
    y = 0 * real(first(beta.coeff)*beta.df)
    for n in eachindex(beta.coeff)
        # assume time starts at zero
        y += derivative(t, n, beta)
    end
    return real(y) 
end

derivative(t::Number, n::Number, beta::FourierTransform) =
    2π*im*n*beta.df*beta.coeff[n] * exp(2π*im*beta.df*n*t) / length(beta.coeff)

function RegularTimeseries_manual(beta::FourierTransform)
    N = length(beta.coeff) # number of observations
    f_nyquist = -beta.df*first(eachindex(beta.coeff))
    dt = 1 / (2*f_nyquist)

    # dumb to do a calculation just to get the type
    y_eltype = eltype(expand(dt, beta))
    y = zeros(y_eltype, 0:N-1) # an OffsetArray
    
    # assume ok to start at index 0
    for  i in eachindex(y)
        y[i] = expand(dt*i, beta)
    end
    # again assume that indices start at zero
    return RegularTimeseries( y, dt)
end

function convolve(w::RegularTimeseries,y::RegularTimeseries)
    # require time sampling to be equal
    w.dt != y.dt && error("time sampling required to be consistent")

    i0 = 0 # by construction with OffsetArrays
    h = zero(y.x) # output
    nmin = minimum(eachindex(y.x))
    nmax = maximum(eachindex(y.x))
    for n in eachindex(y.x)
	for m in eachindex(w.x)
	    if (nmin <= (n-m+i0) <= nmax) # check bounds
		h[n] += w.x[m] * y.x[n-m+i0]
            elseif (n-m+i0) < nmin
                # assume equilibrium at start
                h[n] += w.x[m] * y.x[nmin]
            elseif (n-m+i0) < nmax
                # assume equilibrium at end
                h[n] += w.x[m] * y.x[nmax]
	    end
	end
    end
    return RegularTimeseries(h, y.dt)
end

function Base.:(/)(h::FourierTransform, x::FourierTransform)
    (h.freq != x.freq) && error("frequencies do not match")
    return FourierTransform(h.coeff ./ x.coeff, h.freq)
end

periodogram(y::RegularTimeseries) = periodogram(FourierTransform(y))   

#  is this function updated?
function periodogram(ŷ::FourierTransform)
    T = 1 / ŷ.freq[1] #SpectraFromScratch.record_length(y)
    N = length(ŷ.coeff) #length(y.x)
    psd = zeros(eltype(abs(first(ŷ.coeff))^2), maximum(abs.(eachindex(ŷ.coeff))))
    f = zeros(eltype(first(ŷ.freq)), maximum(abs.(eachindex(ŷ.coeff))))
    for m in eachindex(ŷ.coeff)
        if m < 0
            psd[-m] += abs(ŷ.coeff[m])^2
            f[-m] = abs(ŷ.freq[m])
        elseif m > 0
            psd[m] += abs(ŷ.coeff[m])^2
            f[m] = ŷ.freq[m] # overwrite just to be sure
        end
    end
    return FrequencySpectrum((T/N^2)*psd, f)     
end

"""
    function band_avg(yy,num,dimension)

Inputs:
yy, quantity to be averaged (must be vector or matrix)

num, number of bands to average
dimension (optional), dimension to average along; if specified, must be 1 or 2

Tom Farrar, 2016, jfarrar@whoi.edu
Ported to Julia, Jake Gebbie, 2021, jgebbie@whoi.edu =#
"""
function band_average(yy::AbstractVector{T}, num; dim=missing) where T <: Number
    numdims = ndims(yy)
    nyy = size(yy)

    if (numdims > 2) error("Dimension must be equal to 1 or 2 for band_avg") end

    # shortcut execution
    if numdims == 1
        # initialize yy_avg
        yy_avg = fill(zero(T),floor(Integer,nyy[1]/num))
        for n = 1:num
            yy_avg += yy[n:num:end-(num-n)]
        end
        
    elseif numdims == 2
        if ismissing(dim) 
            greaterthanone = x -> x>1
            if count(greaterthanone,yy) > 1
                error("Dimension must be specified for 2D input to band_avg")
            else
                dim = findfirst(greaterthanone,yy) 
            end
        end
        if dim==1
            # initialize yy_avg
            nyy_avg = (floor(Integer,nyy[1]/num),nyy[2])
            yy_avg = fill(zero(T),nyy_avg)
            for n=1:num
                yy_avg += yy[n:num:end-(num-n),:]
            end
        elseif dim==2
            #initialize yy_avg
            nyy_avg = (nyy[1],floor(Integer,nyy[2]/num))
            yy_avg = fill(0,nyy)
            for n=1:num
                yy_avg += yy[:,n:num:end-(num-n)]
            end
        end
    end

    # take the average
    return yy_avg./num
end

function band_average(psi::FrequencySpectrum, num; dim=missing)
    yy_avg = band_average(psi.psd, num, dim=dim)
    f_avg = band_average(psi.freq, num, dim=dim)
    return FrequencySpectrum(yy_avg, f_avg)
end

"""
function confid(α,ν)

Help with computing confidence intervals

should be sigma^2/S^2 confidence bounds where sigma^2 is true variance
check value (J&W) is alpha =.05, nu=19, lower bound is .58
upper bound is 2.11

"""
function confid(α,ν)

    upperv = quantile(Chisq(ν),1-α)
    lowerv = quantile(Chisq(ν),α)
    lower=ν/upperv;
    upper=ν/lowerv;

    return lower, upper

end

"""
function total_spectral_energy(Φ,f)

# Arguments
- `Φ`: power spectral density
- `f`: Fourier frequencies
# Output
- `e`: total energy
"""
function total_spectral_energy(Ψ,f)
    !iszero(first(f)) ? (Δf = first(f)) : (Δf = f[2])
    return e = sum(Ψ)*Δf
end
function total_spectral_energy(Ψ::FrequencySpectrum)
    f = Ψ.freq
    psd = Ψ.psd
    return total_spectral_energy(psd, f)
end
function total_spectral_energy(x::FourierTransform)
    N = length(x.coeff)
    e = zero(eltype((abs(first(x.coeff))^2)))
    for m in eachindex(x.coeff)
        # do not include energy in mean
        if m ≠ 0
            e += abs(x.coeff[m])^2
        end
    end
    return e/N^2
end

"""
    spectral_power_law(β,f)

Create a `FrequencySpectrum` according to f.^-β.  
For units, type of `f` requires uniform vector.

# Arguments
- `f`: frequencies
- `β`: power law coefficient, low frequencies
- `e`: total energy
# Optional Arguments
- `βhi`: power law coefficient, high frequencies
- `fbreak`:: break in power lawrence
# Output
- `Φ`::`FrequencySpectrum`
"""
function spectral_power_law(f, βlo, σ2=1.0; βhi=nothing, fbreak=nothing)
    nf = length(f)
    fnondim = f ./ first(f)
    T = 1 / first(f)
    Ψnondim = fnondim.^-βlo 
    if !isnothing(βhi)
        # high-low frequency break point, add to arguments
        fbreak_nondim =  fbreak ./ first(f)
        scale = fbreak_nondim^(βlo - βhi)
        Ψnondim .+= (1/scale)*fnondim.^-βhi
    end

    σ2nondim = sum(Ψnondim)/T # nf^2
    return FrequencySpectrum((σ2/σ2nondim) .* Ψnondim, f)
end

"""
function spectralbasis(t,f)

basis function to reconstruct mean ocean temperature (Θ̄)
on the t temporal grid

# Arguments
- `t`: times of interest
- `f`: Fourier frequencies
- `includemean=false::Bool`: include the mean value in the basis set?, 
# Output
- `A::Matrix`: each column is an independent basis function,
first (nt-1)/2 columns are sine coefficients
second (nt-1)/2 columns are cosine coefficients
last column represents the mean value
"""
function spectral_basis(t,f,includemean=false)
    
    Acos = Matrix{Float64}(undef,length(t),length(f))
    Asin = Matrix{Float64}(undef,length(t),length(f))
    for (ii,ff) in enumerate(f)
        Acos[:,ii] = cos.(2π*ff.*t)
        Asin[:,ii] = sin.(2π*ff.*t)
    end

    if includemean
        # add a column for the mean.
        return hcat(Acos,Asin,ones(length(t)))
    else
        return hcat(Acos,Asin)
    end
    
end

function time_average(tstart::Number, tend::Number, x::FourierTransform)
    (tend ≤ tstart) && error("times out of order") 
    # summation of frequencies
    b = zero(eltype(real(first(x.coeff)))) # allocate a complex number
    N = length(x.coeff)
    # iterate over frequencies
    # why first? the first dimension is Frequency
    for m in eachindex(x.coeff)
        b += time_average(tstart, tend, x, m) 
    end
    return b
end

time_average(tstart, tend, x::FourierTransform, m::Number) =
    ( m ≠ 0) ?
    (real(integrate(tstart, tend, x, m)) / (length(x.coeff) * (tend - tstart))) :
    (real(x.coeff[0])/length(x.coeff))
                                                             
function integrate(tstart, tend, x::FourierTransform, m::Number)
    # N = length(x.coeff)
    # A = x.coeff[m] / (2π * im * x.freq[m]) # amplitude of wave 
    A = x.coeff[m] / (2π * im * m * x.df ) # amplitude of wave 
    # limit1 = exp(2π*im*x.freq[m]*tstart)
    # limit2 = exp(2π*im*x.freq[m]*tend)
    limit1 = exp(2π*im*m*x.df*tstart)
    limit2 = exp(2π*im*m*x.df*tend)
    return   A * (limit2 - limit1)
end

function phase(x::FourierTransform)
    ind = first(axes(x.coeff))
    phi = OffsetArray(zeros(length(ind)), ind)
    for n in eachindex(x.coeff)
        phi[n] = atan(imag(x.coeff[n])/real(x.coeff[n]))        
    end
    return phi
end

function FourierTransform(Ψ::FrequencySpectrum)
    ## get timeseries that goes with frequency spectrum
    # NOTE: This function is not deterministic.
    # It is one timeseries realization of the spectrum, but there are others.
    modes = fourier_modes(Ψ)
    nf = maximum(abs.(eachindex(Ψ.psd)))

    # careful here about N even or odd
    ϕnyquist = rand() > 0.5 ? 0.0 : π # Nyquist must be 0° or 180° phase
    
    ϕ = vcat(2π *rand(nf-1) .- π, ϕnyquist) 
    N = length(modes) #2nf
    df = first(Ψ.freq)

    amps = amplitudes(Ψ)
    unt = first(amps)*exp(im)
    # xhat = OffsetArray(zeros(ComplexF64,N),modes)
    xhat = OffsetArray(zeros(eltype(unt),N),modes)
    for m in modes
        if m > 0
            xhat[m] = amps[m]*exp(im*ϕ[m]) # apply phase
        elseif m < 0
            xhat[m] = conj(amps[m]*exp(im*ϕ[-m])) # apply phase
        else
            xhat[0] = zero(eltype(xhat))
        end
    end
    return FourierTransform(xhat, df)
end

"""
    function amplitudes(Ψ::FrequencySpectrum; even = true)

Retrieve amplitudes for individual positive + negative frequency waves.        
"""
function amplitudes(Ψ::FrequencySpectrum; even = true)
    nf = length(Ψ.psd)
    T = 1/first(Ψ.freq)
    modes  = fourier_modes(Ψ)
    N = length(modes)
    physvar = sqrt(first(Ψ.psd)*first(Ψ.freq))
    amp = OffsetArray(zeros(eltype(physvar),N), modes)
    for m in modes
        if even && m == -nf # solo Nyquist frequency
            # keep all energy in one wave
            amp[m] = √(N^2*Ψ.psd[-m]/T)
        elseif !iszero(m)
            # split energy evenly betwwen + and - frequencies
            amp[m] = √(N^2*Ψ.psd[abs(m)]/(2T))
        end
    end
    return amp
end


RegularTimeseries(Ψ::FrequencySpectrum) = RegularTimeseries(FourierTransform(Ψ))
    
end
