module CPS

using LinearAlgebra

author = Dict{Symbol, String}(
    :index => "415104",
    :name  => "Maciej Kucharski",
    :email => "mackuchar@student.agh.edu.pl",
    :group => "4",
)

# SygnaĹy ciÄgĹe
cw_rectangular(t::Real; T=1.0)::Real = abs(t) < T/2 ? 1 : 0
cw_triangle(t::Real; T=1.0)::Real = abs(t) < T ? 1 - 1/T*abs(t) : 0
cw_literka_M(t::Real; T=1.0)::Real = abs(t) < T ? (t < 0 ? -t + 1 : t + 1) : 0
cw_literka_U(t::Real; T=1.0)::Real = abs(t) < T/2 ? t^2 : 0

ramp_wave(t::Real)::Real = t - floor(t)
sawtooth_wave(t::Real)::Real = -2*(t-floor(t+1/2))
triangular_wave(t::Real)::Real = 4*abs((t+1/4)-floor((t+1/4)+1/2)) - 1
square_wave(t::Real)::Real =  sign(sin(2*pi*t))
pulse_wave(t::Real, ρ::Real=0.2)::Real = (t - floor(t)) <= ρ ? 1 : 0
impulse_repeater(g::Function, t1::Real, t2::Real)::Function = x -> g(mod(x - t1, t2 - t1))

function ramp_wave_bl(t; A=1.0, T=1.0, band=20.0)
    out = 0
    N = 1
    while (arg = (2*pi*N)/T) <= band * 2*pi
        out += (-1)^N * sin(arg * t) / N
        N += 1
    end
    return -2 * A / pi * out
end


function sawtooth_wave_bl(t; A=1.0, T=1.0, band=20.0)
    out = 0
    N = 1
    while (arg = (2pi*N)/T) <= band * 2pi
        out += (-1)^N * sin(arg * t) / N
        N += 1
    end
    return 2 * A / pi * out
end

function triangular_wave_bl(t; A=1.0, T=1.0, band=20.0)
    out = 0
    N = 1
    while (arg = (2pi*(2*N-1))) <= band * 2pi
        out += ((-1)^N)/(2N-1)^2 * sin(arg*t)
        N+=1
    end
    
    return -8 / pi^2 * out
end

function square_wave_bl(t; A=1.0, T=1.0, band=20.0)
    out = 0
    N = 1
    while (arg = (2pi*(2*N-1))) <= band * 2pi
        out += 1/(2*N-1) * sin(arg*t)
        N+=1
    end

    return 4/pi * out
end

function pulse_wave_bl(t; ρ=0.2, A=1.0, T=1.0, band=20.0)
    out = 0
    N = 1
    while (arg = 2*pi*N/T) <= band * 2pi
        out += 1/N * sin(pi*N*ρ/T)*cos(arg*t)
        N+=1
    end

    return A*ρ/T + 2*A/pi * out

end


function impulse_repeater_bl(g::Function, t1::Real, t2::Real, band::Real)::Function
    T = t2 - t1
    omega= (2pi / T)
    n= div(band * 2pi, omega)

    N = 1000
    t = range(t1, t2, length=N + 1)
    dt = (t2 - t1) / N
end

function rand_signal_bl(f1::Real, f2::Real)::Function
    f = f1 .+ rand(1000) .* (f2 - f1)
    fi = -pi .+ rand(1000) * 2pi
    A = randn(1000)
    A = A ./ sqrt(0.5 * sum(A .^ 2))
    return t -> sum(A .* sin.(2pi * f .* t .+ fi))
end


# SygnaĹy dyskretne
kronecker(n::Integer)::Real = n == 0 ? 1 : 0
heaviside(n::Integer)::Real = n >= 0 ? 1 : 0

# Okna
rect(N::Integer)::AbstractVector{<:Real} = ones(N)
triang(N::Integer)::AbstractVector{<:Real} = [1 - (2abs(n - ((N - 1) / 2))) / (N - 1) for n = 0:N-1]
hanning(N::Integer)::AbstractVector{<:Real} = [0.5*(1 - cos((2pi*n)/(N-1))) for n = 0:N-1]
hamming(N::Integer)::AbstractVector{<:Real} = [0.54 - 0.46*cos((2pi*n)/(N-1)) for n = 0:N-1]
blackman(N::Integer)::AbstractVector{<:Real} = [0.42 - 0.5*cos((2pi*n)/(N-1)) + 0.08*cos((4pi*n)/(N-1)) for n = 0:N-1]

# Parametry sygnaĹĂłw
mean(x::AbstractVector)::Number = sum(x) / length(x)
peak2peak(x::AbstractVector)::Real = abs(maximum(x) - minimum(x))
energy(x::AbstractVector)::Real = sum(x.^2)
power(x::AbstractVector)::Real = energy(x) / length(x)
rms(x::AbstractVector)::Real = sqrt(power(x))

function running_mean(x::AbstractVector, M::Integer)::Vector
    N = -M/2:M/2
    out = zeros(length(x))
    for n in 1:length(x)
        for m in N
            if n+m > 0 && n+m <= length(x)
                out[n]+=x[n+m]
            end
        end
    end
    return out/length(out)
end

function running_energy(x::AbstractVector, M::Integer)::Vector
    N = -M/2:M/2
    out = zeros(length(x))
    for n in 1:length(x)
        for m in N
            if n+m > 0 && n+m <= length(x)
                out[n]+=x[n+m]^2
            end
        end
    end
    return out
end

function running_power(x::AbstractVector, M::Integer)::Vector
    N = -M/2:M/2
    out = zeros(length(x))
    for n in 1:length(x)
        for m in N
            if n+m > 0 && n+m <= length(x)
                out[n]+=x[n+m]^2
            end
        end
    end
    return out/length(out)
end



# PrĂłbkowanie
function interpolate(
    m::AbstractVector,
    s::AbstractVector,
    kernel::Function=sinc
)
    return x -> begin
    out = 0
    T = m[2] - m[1]
    for i in eachindex(m)
        out += s[i] * kernel((x - m[i]) / T)
    end
    return out
end


# Kwantyzacja
quantize(L::AbstractVector)::Function = x -> L[argmin(abs.(x.-L))]
SQNR(N::Integer)::Real = 6.02*N
SNR(Psignal, Pnoise)::Real = 10*log10(Psignal/Pnoise)


# Obliczanie DFT
function dtft(freq::Real; signal::AbstractVector, fs::Real)
    out = 0.0
    for N in eachindex(signal)
        out += signal[n] * exp(-2im*pi* freq * N / fs)
    end
    return out
end
function dft(x::AbstractVector)::Vector
    N = length(x)
    X = zeros(ComplexF64,N)
    for k in 1:N
        for n in 1:N
            X[k]+=x[n]* exp(-2im*pi*(n-1)*(k-1)/N)
        end
    end
    return X/N    
end

function idft(X::AbstractVector)::Vector
    N = length(X)
    k = OffsetArray([ exp(-2im*pi* n / N) for n in 0:(N-1)],0:(N-1))
    [
        (1 / N) * sum((
            X[n+1] * k[(n*f)%N]
            for n in 0:(N-1)
        ))
        for f in 0:(N-1)
    ]
end
function rdft(x::AbstractVector)::Vector
    N = length(x)
    k = OffsetArray([exp(-2im*pi* n / N) for n in 0:(N-1)],0:(N-1))
    [
        sum((x[n+1] * k[(n*f)%N] for n in 0:(N-1))) for f in 0:(N÷2)
    ]
end

function irdft(X::AbstractVector, N::Integer)::Vector
    L = length(X)
    X1 = [n <= L ? X[n] : conj(X[2S-n+(N % 2 == 0 ? 0 : 1)]) for n in 1:N]
    real.(idft(X1))
end

function fft_radix2_dit_r(x::AbstractVector)::Vector
   missing
end

function ifft_radix2_dif_r(X::AbstractVector)::Vector
   missing
end

function fft(x::AbstractVector)::Vector
    dft(x) # # Może da rade lepiej?
end

function ifft(X::AbstractVector)::Vector
    idft(X) # # Może da rade lepiej?
end


fftfreq(N::Integer, fs::Real)::Vector = [n * N / fs for n in 0:(N-1)]
rfftfreq(N::Integer, fs::Real)::Vector = [n * N / fs for n in 0:(N÷2)]
amplitude_spectrum(x::AbstractVector, w::AbstractVector=rect(length(x)))::Vector = abs.(fft(x .* w)) / (length(x) * mean(w))
power_spectrum(x::AbstractVector, w::AbstractVector=rect(length(x)))::Vector = amplitude_spectrum(x, w) .^ 2
psd(x::AbstractVector, w::AbstractVector=rect(length(x)), fs::Real=1.0)::Vector = (fft(x .* w)).^2 / (sum(w.^2) * fs)


function periodogram(
    x::AbstractVector,
    w::AbstractVector=rect(length(x)),
    L::Integer = 0,
    fs::Real=1.0)::Vector
    N = length(x)
    K = length(w)
    if L == 0
        L = K
    end
    M = div(N - K, L) + 1
    out = zeros(K)
    for n in 0:M-1
        s1 = n * L + 1
        if s1 + K - 1 > N
            break
        end
        lop = x[s1:s1+K-1] .* w
        X = fft(lop)
        out += X.^2 / (sum(w.^2) * fs)
    end
    return out
end



function stft(x::AbstractVector, w::AbstractVector, L::Integer)::Matrix
    missing
end


function istft(X::AbstractMatrix{<:Complex}, w::AbstractVector{<:Real}, L::Integer)::AbstractVector{<:Real}
    missing
end

function conv(f::AbstractVector, g::AbstractVector)::Vector
    n = length(f)
    m = length(g)
    out = zeros(n + m - 1)
    for i in 1:n
        for j in 1:m
            out[i+j-1] += f[i] * g[j]
        end
    end
    return out
end

function fast_conv(f::Vector, g::Vector)::Vector
    N = length(f) + length(g) - 1
    f1 = vcat(f, zeros(N - length(f)))
    g1 = vcat(g, zeros(N - length(g)))
    F = fft(f1)
    G = fft(g1)
    Y = F .* G
    out = Real(ifft(Y))

    return out
end

function overlap_add(x::Vector, h::Vector, L::Integer)::Vector
    missing
end

function overlap_save(x::Vector, h::Vector, L::Integer)::Vector
    missing
end

function lti_filter(b::Vector, a::Vector, x::Vector)::Vector
    N = length(x)
    M = length(b) - 1
    K = length(a) - 1
    out = zeros(Float64, N)
    for n in 1:N
        for k in 0:M
            if n - k > 0
                out[n] += b[k+1] * x[n-k]
            end
        end
        for k in 1:K
            if n - k > 0
                out[n] -= a[k+1] * out[n-k]
            end
        end
    end
    return out
end

function filtfilt(b::Vector, a::Vector, x::Vector)::Vector
    y_fwd_filt = lti_filter(b, a, x)
    y_rvs = reverse(y_fwd_filt)
    y_rvs_filt = lti_filter(b, a, y_rvs)
    return reverse(y_rvs_filt)
end

function lti_amp(f::Real, b::Vector, a::Vector)::Real
    M = length(b)
    K = length(a)
    a = sum(b[m+1] *  exp(-2im*pi* * f * m) for m in 0:M-1)
    b = sum(a[k+1] * exp(-2im*pi * f * k) for k in 0:K-1)
    H = a / b
    return abs(H)
end

function lti_phase(f::Real, b::Vector, a::Vector)::Real
    M = length(b)
    K = length(a)
    a = sum(b[m+1] * exp(-2im*pi * f * m) for m in 0:M-1)
    b = sum(a[k+1] * exp(-2im*pi * f * k) for k in 0:K-1)
    H = a / b
    return angle(H)
end


function firwin_lp_I(order, F0)
    for n in -order/2:order/2
    H =2F0*sinc(2*F0*n)
    return H
    end
end

function firwin_hp_I(order, F0)
    for n in -order/2:order/2
    H = kronecker(n) - 2F0*sinc(2F0*n)
    return H
    end
end

function firwin_bp_I(order, F1, F2)
    for n in -order/2:order/2
    H = 2F2*sinc(2F2*n) - 2F1*sinc(2F1*n)
    return H
    end
end

function firwin_bs_I(order, F1, F2)
    for n in -order/2:order/2
    H = kronecker(n) - (2F2*sinc(2F2*n) - 2F1*sinc(2F1*n))
    return H
    end
end

function firwin_lp_II(N, F0)
    n = range(-N/2,N/2,length=N)
    for i in eachindex(n)
    H =2F0*sinc(2*F0*n[i])
    return H
    end
end

function firwin_bp_II(N, F1, F2)
    n = range(-N/2,N/2,length=N)
    for i in eachindex(n)
    H = 2F2*sinc(2F2*n[i]) - 2F1*sinc(2F1*n[i])
    return H
    end
end

function firwin_diff(N::Int)
    result = []
    for n in -order:1:order
        if n == 0
            push!(result, 0)
        else
            push!(result, cospi(n) / n)
        end
    end
    return result
end

function resample(x::Vector, M::Int, N::Int)
    org = 1:length(x)
    l_temp = div(length(x) * M, N)
    i_temp = range(1, length(x), l_temp)
    i_func = interpolate(org, x)
    out = [i_func(t) for t in i_temp]
    return out
end



end