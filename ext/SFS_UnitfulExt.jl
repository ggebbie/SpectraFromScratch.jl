module SFS_UnitfulExt

# Load main package and extension dependencies
using SpectraFromScratch, Unitful, FFTW, OffsetArrays

# const SFS = SpectraFromScratch

import SpectraFromScratch: centered_fft
import SpectraFromScratch: RegularTimeseries

# Extend functionality in main package with types from the extension dependencies
# MyPackage.func(x::ExtDep.SomeStruct) = ...
# function SpectraFromScratch.FourierTransform(y::RegularTimeseries{T,R}; alg=:centered_fft) where T <: Quantity where R <: Quantity

#     # function FourierTransform(y::RegularTimeseries; alg=:centered_fft)

    
#     if alg==:centered_fft
#         return centered_fft(y)
#     elseif alg==:manual
#         return FourierTransform_manual(y)
#     else
#         error("SpectraFromScratch.jl: incorrect keyword")
#     end
# end


# end

function SpectraFromScratch.centered_fft(y::SpectraFromScratch.RegularTimeseries{<:Quantity,<:Quantity})
    m = SpectraFromScratch.fourier_modes(y)
    T = SpectraFromScratch.record_length(y) 

    #the dimensional frequency scale, this is an "iterator", not a vector, in julia
    f = SpectraFromScratch.fourier_frequencies(m, T)

    # df = fundamental frequency
    df = f[1]
    
    #=swaps the halves of the FFT vector so that 
    the zero frequency is in the center.
    If you are going to compute an IFFT, 
    first use X=ifftshift(X) to undo the shift =#

    # assumes input is uniform (handled by struct definition?)
    unt = Unitful.unit(first(y.x))
    x̂ = unt.*SpectraFromScratch.fftshift(fft(ustrip.(OffsetArrays.no_offset_view(y.x))))
    return SpectraFromScratch.FourierTransform(OffsetArray(x̂, m), df)
end

function SpectraFromScratch.centered_ifft(beta::FourierTransform{<:Quantity})
    unt = Unitful.unit(real(first(beta.coeff)))
    y = unt.*ifft(ifftshift(ustrip.(OffsetArrays.no_offset_view(beta.coeff))))
    f_nyquist = -beta.df*first(eachindex(beta.coeff))
    dt = 1 / (2*f_nyquist)

    # assume indices start at zero
    return RegularTimeseries( OffsetArray(real.(y), -1), dt)
end

end
