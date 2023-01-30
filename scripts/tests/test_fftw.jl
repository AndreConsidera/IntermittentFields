using FFTW
using Statistics
using Distributions

N = 2^12;
x = rand(Normal(0,1), N);

P = plan_rfft(x; flags = FFTW.MEASURE, timelimit = Inf)
y = P * x

Pinv = plan_irfft(y, length(x); flags = FFTW.ESTIMATE, timelimit = Inf)
xinv = Pinv * y

# check inversion
x ≈ xinv
x == xinv

# testing speed
@time for i in 1:2000
    rfft(x)
end
@time P = plan_rfft(x; flags = FFTW.MEASURE, timelimit = Inf)
@time for i in 1:2000
    P * x
end

# inverse
@time for i in 1:2000
    irfft(y, length(x))
end
@time Pinv = plan_irfft(y, length(x); flags = FFTW.MEASURE, timelimit = Inf)
@time for i in 1:2000
    Pinv * y
end
