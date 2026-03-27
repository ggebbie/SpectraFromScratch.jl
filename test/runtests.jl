using Revise
using SpectraFromScratch
using Distributions
using Test
using Statistics

@testset "SpectraFromScratch.jl" begin
    N  = 2000 # underscore just for visual appearance
    Δt = 1      # could make \Delta in Julia REPL but not in notebook
    t  = Δt:Δt:N*Δt
    f  = 20/((N-1)*Δt)
    noise_val = 0.2 # desired noise std deviation
    yb = 1 .+ noise_val.*randn(N) .+ 0.75 .* sin.(2π*f*t)

    #plot(t,yb, leg=false)
    #title!("A sinusoid plus noise")
    #xlabel!("Time")

    @testset "FourierTransform struct" begin
        yy = yb # just renaming yb to mimic matlab code
        N = length(yy)
        T = N * Δt
        yy .-= mean(yy) # remove the mean

        # Compute the FFT of the entire tapered record.
        # Y,freq_i = centered_fft(yy,Δt)
        # ŷ_orig = FourierTransform(Y,freq_i)

        y = EvenlySampledTimeseries(yy, t)
        @time ŷ_orig = centered_fft(y)
        @time ŷ = FourierTransform(y)

        # Inverse Fourier Transform
        @time ỹ = EvenlySampledTimeseries(ŷ, y.t)
        @test isapprox(y.x, ỹ.x)
        @test isapprox(y.t, ỹ.t)

        # needs time correction, fft assumes t=0 at first obs
        @time ỹ_orig = EvenlySampledTimeseries(ŷ_orig, y.t.-first(y.t))
        @test isapprox(y.x, ỹ_orig.x) # now passes
        @test isapprox(y.t, ỹ_orig.t .+ first(y.t)) 

        @time ỹ_orig2 = centered_ifft(ŷ_orig, y.t)
        @test isapprox(y.x, ỹ_orig2.x) # now passes
        @test isapprox(y.t, ỹ_orig2.t) 
    end

    @testset "bin averaging" begin
        navg = 20
        y_avg = band_avg(yb,navg)
        t_avg = band_avg(t,navg)

        @test isapprox(length(y_avg),N/navg)

        #plot(t_avg,y_avg,leg = false)
        #title!("Bin averaged version of the timeseries above")
        #xlabel!("Time")

    end

    @testset "spectrum" begin

        # OK, now let's try doing a spectrum. First, though, let's recall what we are wanting to do overall. We want to make a function to estimate the band-averaged spectrum and plot it with a 95% confidence interval. To do that, we will follow the steps from Section 4.7.1 of the class notes on spectral analysis
        yy = yb # just renaming yb to mimic matlab code
        N = length(yy)
        T = N * Δt
        yy .-= mean(yy) # remove the mean

        # Compute the FFT of the entire tapered record.
        # Y,freq_i = centeredFFT(yy,Δt)
        y = EvenlySampledTimeseries(yy, t)
        ŷq = FourierTransform(y)

        Ψraw = SpectraFromScratch.periodogram(y)
        @test all(Ψraw.psi .> zero(first(Ψraw.psi)))

        # Band average the raw spectrum over 𝑛𝑑 frequency bands-- this could be done by an algorithm like equation ??? or by computing a running average and subsampling. Generate the new frequency vector, either by subsampling the Fourier frequencies at the interval of 𝑛𝑑/𝑇 or by band averaging the frequency vector. --> We will use our band-averaging function on both the spectrum and the frequency vector
        nbands = 11
        Ψavg = band_average(Ψraw, nbands)

        @test length(Ψraw.psi) > length(Ψavg.psi)
        
        #plot(freq,real(Ψavg),leg=false)
        #plot!(freq,imag(Ψavg),leg=false)
        #title!("Band-averaged spectral estimate")

        # We were expecting the imaginary part to be zero. It appears to be nonzero. This may be a numerical artifact, but let's check this out.
        #plot(freq,imag(Ψavg),leg=false)
        #title!("Imaginary part of spectral estimate")

        #         OK, just as I suspected, the imaginary part is a numerical artifact, close to floating point precision. The best way to handle this in the function we will build would probably be to check if the inputs are real, and, if they are, then just throw away the imaginary part of our spectral estimate. (If the input time series is complex, we should expect some imaginary part in the spectrum.)

        #plot(freq,real(Ψavg), xaxis=:log, yaxis=:log)
        #title!("Band-averaged spectral estimate")

        @testset "confidence intervals" begin

            #=
            In MATLAB, chi2inv(alpha/2,nu) returns 𝜒2𝜈;𝛼/2 . To compute confidence intervals, you simply need to compute the multiplicative factors 𝜈/𝜒2𝜈;𝛼/2 and 𝜈/𝜒2𝜈;1−𝛼/2. Check values (e.g., Jenkins and Watts, 1968, p. 81) are 𝜈/𝜒219;0.025=0.58 and 𝜈/𝜒219;0.975=2.11

(I get 0.5783 and 2.1333 in MATLAB).

            It seems that the analogous function in python's scipy.stats module is the "Inverse survival function" (stats.chi2.isf), which is defined as 11−cdf , where cdf is the cumulative distribution function of the chi-square distribution: https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.chi2.html#scipy.stats.chi2 =#

            ν = 19 # unicode nu looks like "v" on my computer
            α = 0.025
            χ = cquantile(Chisq(ν), α) 
            println(ν/χ)

            #=    should be sigma^2/S^2 confidence bounds where sigma^2 is true variance
    check value (J&W) is alpha =.05, nu=19, lower bound is .58
            upper bound is 2.11 =#
            lower,upper = confid(α,ν)

            @test 0.59 > lower > 0.57
            @test 2.14 > upper > 2.08 # should be able to narrow this range and still pass
            println(lower)
            println(upper)

        end
    
    end
    @testset "convolution" begin
        function rectangle(T,τ) 
            M = length(τ) # length(rectangle(Trectangle,τ))
            Mmax = convert(Int,(M-1)/2) 
	    w = zeros(length(τ))
	    for i in eachindex(w)
		if abs(τ[i]) <= T
		    w[i] += 1
		end
	    end
            w /= sum(w)
            return EvenlySampledTimeseries(w, τ)
            # return OffsetArray(w, -Mmax:Mmax)
        end
        
        τ = range(-10,10,step=1)
        Trectangle = 0
        M = length(rectangle(Trectangle,τ))
        Mmax = convert(Int,(M-1)/2) 
        w = rectangle(Trectangle,τ)
        N_convolve = 51
        t_convolve = 1:N_convolve
        x = EvenlySampledTimeseries( randn(N_convolve), t_convolve)
        h = convolve(w,x)
        @test h.x == x.x
        @test h.t == x.t

        # test the convolution theorem
        @time ĥ  = centered_fft(h)
        @time x̂  = centered_fft(x)
        @time ŵ  = centered_fft(w)
        
        ŵ_residual = ĥ / x̂
       ŵ_residual = FourierTransform(ĥ.xhat ./ x̂.xhat, ĥ.f)

        @test maximum(abs.(ŵ_residual.xhat)) < 1.1
        @test minimum(abs.(ŵ_residual.xhat)) > 0.9

        @test maximum(abs.(ŵ.xhat)) < 1.1
        @test minimum(abs.(ŵ.xhat)) > 0.9

        N_padded = convert(Int, floor(N_convolve/2))
        τ_padded = range(-N_padded, N_padded, step=1)
        w_padded = rectangle(Trectangle,τ_padded)
        @time ŵ_padded  = centered_fft(w_padded)

        @test maximum(abs.(ŵ_padded.xhat)) < 1.1
        @test minimum(abs.(ŵ_padded.xhat)) > 0.9

    end    
end
