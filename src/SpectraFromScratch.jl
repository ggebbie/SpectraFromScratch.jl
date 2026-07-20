module SpectraFromScratch

using Statistics
using Distributions
using FFTW
using OffsetArrays

export FourierTransform
export RegularTimeseries
export centered_fft
export centered_ifft
export band_average
export confid, total_spectral_energy
export spectral_power_law, spectral_basis
#, observationalmatrix
export convolve
export periodogram

import Base: /

struct FourierTransform{T<:Number, R<:Number}
    # coeff::AbstractVector{T}
    coeff::OffsetVector{T}
    df::R
end

Base.propertynames(x::FourierTransform, private::Bool=false) =
      private ? (:freq,  fieldnames(typeof(x))...) : fieldnames(typeof(x))

function Base.getproperty(x::FourierTransform, d::Symbol)
    if d === :freq
        # reconstruct fourier frequencies
        # should this be an OffsetArray?
        ind = first(axes(x.coeff))
        return OffsetArray(x.df*ind, ind)
        #range(start=x.dt*0, step=x.dt, length=length(x.x))
    else
        return getfield(x, d)
    end
end

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

function RegularTimeseries(x::AbstractVector, t::AbstractVector)
    length(x) != length(t) && error("lengths do not match")
    # if all( abs.(diff(diff(t))) .< 1e-8*one(eltype(t)))
    if all( isapprox.(diff(diff(t)), zero(eltype(t))))

        # minimize machine error
        dt = (last(t)-first(t))/(length(t)-1)

        # times an integer multiple of sampling time?
        ind0 = t./dt
        # if all(isinteger.(ind0))
        #     ind = Integer.(ind0)
            offset = Integer(first(t)./dt)-1
        # else
        #     # revert to starting indices at 0
        #     error("SpectraFromScratch: time not an integer multiple")
        #     # ind = range(start=0, step=1, length=length(x))
        # end
        # println("ind", ind)
        
        # force all timeseries to be OffsetArrays in symmetry with `FourierTransform`
        # x_offset = OffsetArray(x, ind)
        x_offset = OffsetArray(x, offset)
        return RegularTimeseries{eltype(x), eltype(t)}(x_offset, dt)
    else
        error("not evenly spaced")
    end
end

struct FrequencySpectrum{T}
    psi::AbstractVector
    freq::AbstractVector
    function FrequencySpectrum(psi, f)
        isnegative = (x -> x < zero(x))
        if any(isnegative.(f))
            error("one sided spectrum where frequency must be positive ")
        else
            new{eltype(real.(psi))}(real.(psi), f)
        end
    end
end

fourier_modes(N::Number) = iseven(N) ?
	                       (m = (-convert(Int,N/2):convert(Int,(N/2)-1))) :
	                       (m = (-convert(Int,(N-1)/2):convert(Int,((N-1)/2))))

fourier_modes(y::RegularTimeseries) = fourier_modes(length(y))

fourier_frequencies(m, T) = OffsetArray(m/T, m)
function fourier_frequencies(y::RegularTimeseries)
    m = fourier_modes(y)
    T = record_length(y)
    #the dimensional frequency scale, this is an "iterator", not a vector, in julia
    return fourier_frequencies(m, T)
end

sampling_resolution(y::RegularTimeseries) = y.dt

record_length(y::RegularTimeseries) = length(y) * sampling_resolution(y)

# function fourier_basis(y::RegularTimeseries{Tin}) where Tin
#     m = fourier_modes(y)
#     T = record_length(y)
#     f = fourier_frequencies(m, T)
#     U = Vector{Vector{C}}(undef, length(y))
#     for i in eachindex(f)
#         U[i] = exp.(2π*im*f[i].*y.dt)
#     end
# end

"""
 function centered_fft(x,Δt)

 Computes FFT, with zero frequency in the center, and returns 
  dimensional frequency vector.

- Adapted from a function written by Quan Quach of blinkdagger.com 
- Tom Farrar, 2016, jfarrar@whoi.edu
- Julia version, G Jake Gebbie, 2021, ggebbie@whoi.edu

# Arguments
- `x::RegularTimeseries`

# Output
- `x̂`: centered discrete Fourier transform
- `f`: dimensional frequency scale
"""
function centered_fft(y::RegularTimeseries)
    m = fourier_modes(y)
    T = record_length(y) 

    # println("list of Fourier indices: ",m)
    # println("record length:",T)

    #the dimensional frequency scale, this is an "iterator", not a vector, in julia
    f = fourier_frequencies(m, T)
    # println("Fourier frequencies:", f)

    # df = fundamental frequency
    df = f[1]
    
    #=swaps the halves of the FFT vector so that 
    the zero frequency is in the center.
    If you are going to compute an IFFT, 
    first use X=ifftshift(X) to undo the shift =#
    # x̂ = fftshift(x̂)
    x̂ = fftshift(fft(OffsetArrays.no_offset_view(y.x)))
    return FourierTransform(OffsetArray(x̂, m), df)
end

# function centered_ifft(beta::FourierTransform, t::AbstractVector)
#     # assume that time axis needs shifting for beta from the FFT.
#     # tshift = t .- first(t)
#     y = ifft(ifftshift(OffsetArrays.no_offset_view(beta.coeff)))
#     # println("largest complex component is ", maximum(abs.(real.(im.*y))))
#     println("largest complex component is ", maximum(abs.(imag.(y))))
#     return RegularTimeseries(real.(y), t)
# end

function centered_ifft(beta::FourierTransform)
    y = ifft(ifftshift(OffsetArrays.no_offset_view(beta.coeff)))
    # println("largest complex component is ", maximum(abs.(real.(im.*y))))
    # println("largest complex component is ", maximum(abs.(imag.(y))))
    f_nyquist = -beta.df*first(eachindex(beta.coeff))
    # f_nyquist = -beta.freq[begin]
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
    β = OffsetArray(zero(Vector{ComplexF64}(undef, length(y))), m)

    for m in eachindex(f)
        # println(m)

        # check that eachindex correctly pulls indices
        for n in eachindex(y.x)
        # for n in eachindex(y.time)
            # println(n)
            # println(exp(-2π*im*f[m]*y.t[n]) * y.x[n])
            # println(β[m])
            # β[m] += exp(-2π*im*f[m]*y.time[n]) * y.x[n]

            # here assuming (n-1) is ok
            β[m] += exp(-2π*im*f[m]*dt*(n-1)) * y.x[n]
        end
    end
    offset = first(eachindex(f)) - 1
    return FourierTransform(OffsetArray(β, offset), f[1])
end

# default version: fast FFT
FourierTransform(y::RegularTimeseries) = FourierTransform(y; alg=:centered_fft)

function FourierTransform(y; alg=:centered_fft)
    if alg==:centered_fft
        return centered_fft(y)
    elseif alg==:manual
        return FourierTransform_manual(y)
    else
        error("SpectraFromScratch.jl: incorrect keyword")
    end
end

function RegularTimeseries(y; alg=:centered_ifft)
    if alg==:centered_ifft
        return centered_ifft(y)
    elseif alg==:manual
        return RegularTimeseries_manual(y)
    else
        error("SpectraFromScratch.jl: incorrect keyword")
    end
end

Base.length(x::FourierTransform) = 1

# function RegularTimeseries(beta::FourierTransform, t::AbstractVector)
function RegularTimeseries_manual(beta::FourierTransform)
    N = length(beta.coeff) # number of observations
    y = zeros(0:N-1) # an OffsetArray
    # y = zeros(Float64, N)
    # f_nyquist = -beta.freq[begin]
    f_nyquist = -beta.df*first(eachindex(beta.coeff))
    dt = 1 / (2*f_nyquist)
    
    # assume ok to start at index 0
    for  i in eachindex(y)
        for j in eachindex(beta.coeff)
            # assume time starts at zero
            # y[i] += real.(beta.coeff[j] * exp(2π*im*beta.freq[j]*dt*(i-1)))
            y[i] += real.(beta.coeff[j] * exp(2π*im*beta.df*j*dt*i))
        end
    end
    # again assume that indices start at zero
    return RegularTimeseries( y ./N, dt)
end

function convolve(w::RegularTimeseries,y::RegularTimeseries)
    # require time sampling to be equal
    (first(diff(w.time)) != first(diff(y.time))) &&
    error("time sampling required to be consistent")

    # w required to have a zero time for reference
    i0 = findfirst(iszero, w.time)
            
    if isempty(i0)
        # if no zero, could add code to extrapolate off end of time grid
        error("time grid not consistent")
    end
            
    h = zero(y.x) # output
    nmin = minimum(eachindex(y.x))
    nmax = maximum(eachindex(y.x))
    for n in eachindex(y.x)
        # println(n)
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
    return RegularTimeseries(h, y.time)
end

function Base.:(/)(h::FourierTransform, x::FourierTransform)
    (h.freq != x.freq) && error("frequencies do not match")
    return FourierTransform(h.coeff ./ x.coeff, h.freq)
end

function periodogram(y::RegularTimeseries)
        ŷ = FourierTransform(y)
    return periodogram(ŷ)   
    #     # compute spectrum
    #     ispositive = x -> x > 0
    #     ff = findall(ispositive, ŷ.f)
    #     Y = ŷ.coeff[ff]
    #     freq_i = ŷ.f[ff]
    #     T = SpectraFromScratch.record_length(y)
    # N = length(y.x)
    # # check that "2" is appropriate for zero-frequency coefficient
    # return FrequencySpectrum((2*T/N^2).*Y.*conj(Y), freq_i)
end

#  is this function updated?
function periodogram(ŷ::FourierTransform)
    # # compute spectrum
    # ispositive = x -> x > zero(x)
    # ff = findall(ispositive, ŷ.f)
    # freq_i = ŷ.f[ff]
    T = 1 / ŷ.freq[1] #SpectraFromScratch.record_length(y)
    N = length(ŷ.coeff) #length(y.x)
    psi = zeros(eltype(abs(first(ŷ.coeff))^2), maximum(abs.(eachindex(ŷ.coeff))))
    f = zeros(eltype(first(ŷ.freq)), maximum(abs.(eachindex(ŷ.coeff))))
    for m in eachindex(ŷ.coeff)
        if m < 0
            psi[-m] += abs(ŷ.coeff[m])^2
            f[-m] = abs(ŷ.freq[m])
        elseif m > 0
            psi[m] += abs(ŷ.coeff[m])^2
            f[m] = ŷ.freq[m] # overwrite just to be sure
        end
    end
    return FrequencySpectrum((T/N^2)*psi, f)     
end

"""
 function band_avg.jl   Block averages for band averaging
 [yy_avg]=band_avg(yy,num,dimension)

 Inputs:
	yy, quantity to be averaged (must be vector or matrix)

	num, number of bands to average
	dimension (optional), dimension to average along; if specified, must be 1 or 2

 Tom Farrar, 2016, jfarrar@whoi.edu
 Ported to Julia, Jake Gebbie, 2021, jgebbie@whoi.edu =#
"""
function band_average(yy, num; dim=missing)
    numdims = ndims(yy)
    nyy = size(yy)

    if (numdims > 2) error("Dimension must be equal to 1 or 2 for band_avg") end

    # shortcut execution
    if numdims == 1
        # initialize yy_avg
        yy_avg = fill(0,floor(Integer,nyy[1]/num))
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
            yy_avg = fill(0,nyy_avg)
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
    yy_avg = band_average(psi.psi, num, dim=dim)
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
    psi = Ψ.psi
    return total_spectral_energy(psi, f)
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

# """
#     spectral_power_law(β,f) = f.^-β  

# # Arguments
# - `f`: frequencies
# - `β`: power law coefficient, low frequencies
# - `e`: total energy
# - `βhi`: power law coefficient, high frequencies
# # Output
# - `Φ`: spectral energy density
# """
# function spectral_power_law(f,βlo,e=1.0,βhi=0.0)
#     nf = length(f)
#     Ψ = f.^-βlo  

#     if !iszero(βhi)
#         scale = 0.01^(βlo - βhi)
#         println("scale ",scale)
#         Ψ .+= (1/scale)*f.^βhi
#     end

#     e₀ = 2sum(Ψ)/nf^2
#     Ψ .*= e/e₀

#     return Ψ
# end

# for units
# type of `f` requires uniform vector
#function spectral_power_law(f::StepRangeLen{<:Quantity{<:Number}},βlo,σ2=1.0,βhi=0.0)
function spectral_power_law(f, βlo, σ2=1.0; βhi=nothing, fbreak=nothing)
    nf = length(f)
    fnondim = f ./ first(f)
    T = 1 / first(f)
    Ψnondim = fnondim.^-βlo 
    if !isnothing(βhi)
        # high-low frequency break point, add to arguments
        fbreak_nondim =  fbreak ./ first(f)
        scale = fbreak_nondim^(βlo - βhi)
        #println("scale ",scale)
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

end
