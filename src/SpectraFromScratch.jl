module SpectraFromScratch

using Statistics
using Distributions
using FFTW
using OffsetArrays

export FourierTransform
export EvenlySampledTimeseries
export centered_fft
export centered_ifft
export band_average
export confid, totalspectralenergy
export spectralpowerlaw, spectralbasis, observationalmatrix
export convolve
export periodogram

import Base: /

struct FourierTransform{T<:Number,C<:Number}
    xhat::AbstractVector{C}
    f::AbstractVector{T}
end

struct FrequencySpectrum{T}
    psi::AbstractVector
    f::AbstractVector
    function FrequencySpectrum(psi, f)
        isnegative = (x -> x < zero(x))
        if any(isnegative.(f))
            error("one sided spectrum where frequency must be positive ")
        else
            new{eltype(real.(psi))}(real.(psi), f)
        end
    end
end

struct EvenlySampledTimeseries{T <: Number}
    x::AbstractVector{T}
    t::AbstractVector
    function EvenlySampledTimeseries(x::AbstractVector, t::AbstractVector)
        length(x) != length(t) && error("lengths do not match")
        if all(iszero(diff(diff(t))))
            return new{eltype(x)}(x, t)
        else
            error("not evenly spaced")
        end
    end
end
Base.length(y::EvenlySampledTimeseries) = length(y.x)

fourier_modes(N::Number) = iseven(N) ?
	                       (m = (-convert(Int,N/2):convert(Int,(N/2)-1))) :
	                       (m = (-convert(Int,(N-1)/2):convert(Int,((N-1)/2))))
fourier_modes(y::EvenlySampledTimeseries) = fourier_modes(length(y))

fourier_frequencies(m, T) = OffsetArray(m/T, m)
function fourier_frequencies(y::EvenlySampledTimeseries)
    m = fourier_modes(y)
    T = record_length(y)
    #the dimensional frequency scale, this is an "iterator", not a vector, in julia
    return fourier_frequencies(m, T)
end

sampling_resolution(y::EvenlySampledTimeseries) = -first(diff(y.t))

record_length(y::EvenlySampledTimeseries) = length(y) * sampling_resolution(y)

# function fourier_basis(y::EvenlySampledTimeseries{Tin}) where Tin
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
- `x::EvenlySampledTimeseries`

# Output
- `x̂`: centered discrete Fourier transform
- `f`: dimensional frequency scale
"""
function centered_fft(y::EvenlySampledTimeseries)
    m = fourier_modes(y)
    T = record_length(y) 
    println(m)
    println(T)
    #the dimensional frequency scale, this is an "iterator", not a vector, in julia
    f = fourier_frequencies(m, T)
    println(f)
    #=swaps the halves of the FFT vector so that 
    the zero frequency is in the center.
    If you are going to compute an IFFT, 
    first use X=ifftshift(X) to undo the shift =#
    # x̂ = fftshift(x̂)
    x̂ = fftshift(fft(OffsetArrays.no_offset_view(y.x)))
    return FourierTransform(OffsetArray(x̂, m), f)
end

function centered_ifft(beta::FourierTransform, t::AbstractVector)
    # assume that time axis needs shifting for beta from the FFT.
    tshift = t .- first(t)
    y = ifft(ifftshift(OffsetArrays.no_offset_view(beta.xhat)))
    println("largest complex component is ", maximum(abs.(real.(im.*y))))
    return EvenlySampledTimeseries(real.(y), t)
end

function FourierTransform(y::EvenlySampledTimeseries)
    #the dimensional frequency scale, this is an "iterator", not a vector, in julia
    m = fourier_modes(y)
    f = fourier_frequencies(y)

    # make a β coefficient for every value of m
    β = OffsetArray(zero(Vector{ComplexF64}(undef, length(y))), m)

    for m in eachindex(f)
        # println(m)
        for n in eachindex(y.t)
            # println(n)
            # println(exp(-2π*im*f[m]*y.t[n]) * y.x[n])
            # println(β[m])
            β[m] += exp(-2π*im*f[m]*y.t[n]) * y.x[n]
        end
    end
    return FourierTransform(β, f)
end

Base.length(x::FourierTransform) = 1

function EvenlySampledTimeseries(beta::FourierTransform, t::AbstractVector)
    N = length(beta.f) # number of observations
    y = zeros(Float64, N)
    for  i in eachindex(t) 
        for j in eachindex(beta.xhat)
            y[i] += real.(beta.xhat[j] * exp(2π*im*beta.f[j]*t[i]))
        end
    end
    return EvenlySampledTimeseries(y./N, t)
end

function convolve(w::EvenlySampledTimeseries,y::EvenlySampledTimeseries)
    # require time sampling to be equal
    (first(diff(w.t)) != first(diff(y.t))) &&
    error("time sampling required to be consistent")

    # w required to have a zero time for reference
    i0 = findfirst(iszero, w.t)
            
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
    return EvenlySampledTimeseries(h, y.t)
end

function Base.:(/)(h::FourierTransform, x::FourierTransform)
    (h.f != x.f) && error("frequencies do not match")
    return FourierTransform(h.xhat ./ x.xhat, h.f)
end

function periodogram(y::EvenlySampledTimeseries)
        ŷ = FourierTransform(y)
    return periodogram(ŷ)   
    #     # compute spectrum
    #     ispositive = x -> x > 0
    #     ff = findall(ispositive, ŷ.f)
    #     Y = ŷ.xhat[ff]
    #     freq_i = ŷ.f[ff]
    #     T = SpectraFromScratch.record_length(y)
    # N = length(y.x)
    # # check that "2" is appropriate for zero-frequency coefficient
    # return FrequencySpectrum((2*T/N^2).*Y.*conj(Y), freq_i)
end

function periodogram(ŷ::FourierTransform)
    # # compute spectrum
    ispositive = x -> x > zero(x)
    ff = findall(ispositive, ŷ.f)
    # # println(ff)
    # # Y = ŷ.xhat[ff]
    # # println(Y)
    freq_i = ŷ.f[ff]
    T = 1 / ŷ.f[1] #SpectraFromScratch.record_length(y)
    println("T ",T)
    N = length(ŷ.xhat) #length(y.x)
    println(N)
    # # # check that "2" is appropriate for zero-frequency coefficient
    # # # return FrequencySpectrum((2*T/N^2).*Y.*conj(Y), freq_i)
    # # return FrequencySpectrum((2*T/N^2)*(abs.(Y).^2), freq_i)
    # # e
    psi = zeros(eltype(abs(first(ŷ.xhat))^2), length(ff))
    for m in eachindex(ŷ.xhat)
        println(m)
        if m < 0
            psi[-m] += abs(ŷ.xhat[m])^2
        elseif m > 0
            psi[m] += abs(ŷ.xhat[m])^2
        end
    end
    return FrequencySpectrum((T/N^2)*psi, freq_i)     
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
    f_avg = band_average(psi.f, num, dim=dim)
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
    function totalspectralenergy(Φ,f)

# Arguments
- `Φ`: power spectral density
- `f`: Fourier frequencies
# Output
- `e`: total energy
"""
function totalspectralenergy(Ψ,f)
    #nf = length(f)
    !iszero(first(f)) ? (Δf = first(f)) : (Δf = f[2])
    return e = 2sum(Ψ)*Δf
    #return e = 2sum(Ψ)/nf^2
end
function totalspectralenergy(Ψ::FrequencySpectrum)
    f = Ψ.f
    psi = Ψ.psi
    #nf = length(f)
    !iszero(first(f)) ? (Δf = first(f)) : (Δf = f[2])
    return e = 2sum(psi)*Δf
    #return e = 2sum(Ψ)/nf^2
end

function totalspectralenergy(x::FourierTransform)
    N = length(x.xhat)
    e = zero(eltype((abs(first(x.xhat))^2)))
    for m in eachindex(x.xhat)
        # do not include energy in mean
        if m ≠ 0
            e += abs(x.xhat[m])^2
        end
    end
    return e/N^2
    # #nf = length(f)
    # !iszero(first(f)) ? (Δf = first(f)) : (Δf = f[2])
    # return e = 2sum(Ψ)*Δf
    #return e = 2sum(Ψ)/nf^2
end


# """
#     spectralpowerlaw(β,f) = f.^-β  

# # Arguments
# - `f`: frequencies
# - `β`: power law coefficient, low frequencies
# - `e`: total energy
# - `βhi`: power law coefficient, high frequencies
# # Output
# - `Φ`: spectral energy density
# """
# function spectralpowerlaw(f,βlo,e=1.0,βhi=0.0)
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
#function spectralpowerlaw(f::StepRangeLen{<:Quantity{<:Number}},βlo,σ2=1.0,βhi=0.0)
function spectralpowerlaw(f, βlo, σ2=1.0; βhi=nothing, fbreak=nothing)
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

    σ2nondim = 2sum(Ψnondim)/T # nf^2
    return (σ2/σ2nondim) .* Ψnondim
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
function spectralbasis(t,f,includemean=false)
    
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
