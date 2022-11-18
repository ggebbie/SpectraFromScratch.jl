module SpectraFromScratch

using Statistics, Distributions, FFTW

export centeredFFT, band_avg, confid, totalspectralenergy,
    spectralpowerlaw, spectralbasis, observationalmatrix

"""
 function centeredFFT(x,Δt)

 Computes FFT, with zero frequency in the center, and returns 
  dimensional frequency vector.

- Adapted from a function written by Quan Quach of blinkdagger.com 
- Tom Farrar, 2016, jfarrar@whoi.edu
- Julia version, G Jake Gebbie, 2021, ggebbie@whoi.edu

# Arguments
- `x`: vector to be transformed
- `Δt`: time increment

# Output
- `x̂`: centered discrete Fourier transform
- `f`: dimensional frequency scale
"""
function centeredFFT(x,Δt)
    
    n=length(x)

    #Generate frequency index
    if n%2 == 0
        m=-n/2:n/2-1 # N even
    else
        m=-(n-1)/2:(n-1)/2 # N odd
    end

    #the dimensional frequency scale, this is an "iterator", not a vector, in julia
    f = m/(n*Δt)  
    x̂ = fft(x)
    
    #=swaps the halves of the FFT vector so that 
    the zero frequency is in the center.
    If you are going to compute an IFFT, 
    first use X=ifftshift(X) to undo the shift =#
    x̂ = fftshift(x̂)
    return x̂,f
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
function band_avg(yy,num,dim=missing)
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
    yy_avg=yy_avg./num
    
    return yy_avg
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
function totalspectralenergy(Φ,f)
    nf = length(f)
    return e = 2sum(Φ)/nf^2
end

"""
    spectralpowerlaw(β,f) = f.^-β  

# Arguments
- `f`: frequencies
- `β`: power law coefficient, low frequencies
- `e`: total energy
- `βhi`: power law coefficient, high frequencies
# Output
- `Φ`: spectral energy density
"""
function spectralpowerlaw(f,βlo,e=1.0,βhi=0.0)
    nf = length(f)
    Ψ = f.^-βlo  

    if !iszero(βhi)
        scale = 0.01^(βlo - βhi)
        println("scale ",scale)
        Ψ .+= (1/scale)*f.^βhi
    end

    e₀ = 2sum(Ψ)/nf^2
    Ψ .*= e/e₀

    return Ψ
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
