
##
function rozwiazanie(;
    order::Int = 12,
    fp::Float64 = 193.0,
    f0::Float64 = 77.2,
    z::Vector{Int} = [2, 7, 5, 1, 8, 3],
)
    #0.15693388554399187
    delta(n) = n == 0 ? 1 : 0
    hanning_window(N) = [0.5 + 0.5cospi(2*n/(2N+1)) for n in -N:N]
    n = -order/2:order/2
    h = zeros(Float64,length(n))
    F0 = f0/fp
    for i in 1:length(n)
        h[i] = delta(n[i]) - 2F0*sinc(2*F0*n[i])
    end

    w = hanning_window(order/2)
    hw = h.*w
    out = 0
    for i in z
        out += hw[i]
    end
    return out

end
rozwiazanie()

##
function rozwiazanie(;
    fp::Int = 1806,
    x::Vector{ComplexF64} = ComplexF64[0.11 + 0.33im, -0.26 + 0.18im, -0.01 - 0.66im, -1.42 + 0.89im, -0.79 + 0.64im, 0.23 - 0.78im, -0.23 + 0.29im, -0.59 + 1.19im, -0.13 + 1.01im, -0.84 - 0.07im, 1.1 + 0.43im, 1.12 + 0.05im, 1.58 + 0.24im, -0.96 + 0.74im, 0.55 + 0.53im, -1.51 - 0.23im, -0.6 - 0.62im, -0.03 - 0.09im, -1.16 + 0.25im, 0.28 + 0.05im, -1.25 - 0.5im, 0.34 + 1.17im, 0.35 - 0.19im, 1.25 + 0.2im, 0.62 - 0.41im, 0.29 - 0.27im, -0.68 + 0.13im, 1.52 + 0.36im, 0.57 - 0.5im, -1.09 - 1.05im, -1.01 + 0.76im, 1.37 - 0.34im, 0.28 + 1.01im, -0.48 - 1.21im, -0.77 - 0.69im, -0.39 - 0.54im, 0.24 + 0.98im, 0.14 + 0.1im, -0.25 + 0.44im, 0.26 + 0.86im, -0.8 + 0.23im, 0.34 + 0.39im, -0.01 + 0.19im],
    f::Vector{Int} = [-672, 882, -714, -756, 462],
)
    #-1.7167917797544863
    function dft(x)
    N = length(x)
    X = zeros(ComplexF64,N)
        for k in 1:N
            for n in 1:N
                X[k]+=x[n]*cispi(-2*(n-1)*(k-1)/N)
            end
        end
        return X/N
    end

    N = length(x)
    X = dft(x) 
    delta_f = fp/N
    F = zeros(Float64,N)
    F[Int(floor(N/2)+2)] = -(floor((N-1)/2))*delta_f

    for i in 1:N-1
        if i+1 != floor(N/2)+2
            F[i+1] = F[i]+delta_f
        end
    end

    out = 0

    for i in 1:N
        if F[i] in f
            out += angle(X[i])
        end
    end

    return out

end
rozwiazanie()
## bipolar prostokat / moc

function rozwiazanie_1(;
    fp::Float64 = 397.38,
    t1::Float64 = 8.22,
    N::Int = 926,
)
    #0.010000000000000005
    g(t)=sign(sin(2*π*t))
    y(t)=0.1*g(2.7*t-4.3)
    x = range(start=t1,step=1/fp,length=N)
    out = zeros(length(x))
    for i in 1:length(x)
        out[i] = y(x[i])
    end
    return sum(out.^2)/length(out)
end
rozwiazanie_1()
## trojkat
g(t)= 4*abs(t+1/4-floor(t+1/4+1/2)) - 1
##  fala piłokształtna opadający / energia dyskretna
using CairoMakie
function rozwiazanie_1(;
    fp::Float64 = 481.83,
    t1::Float64 = 0.15,
    N::Int = 91,
)
    #54.711181331623294
    g(t)=-2*(t-floor(t+1/2))
    y(t)=1.3*g(4.8*t-1.5)
    x = range(start=t1,step=1/fp,length=N)
    out = zeros(length(x))
    for i in 1:length(x)
        out[i] = y(x[i])
    end
    return sum(out.^2)
end
rozwiazanie_1()

## trojkatna energia ez

function rozwiazanie_1(;
    fp::Float64 = 339.14,
    t1::Float64 = 7.3,
    N::Int = 435,
)
    #2943.062264350335
    g(t)= 4*abs(t+1/4-floor(t+1/4+1/2)) - 1
    y(t)=4.5*g(3.1*t -1)
    x = range(start=t1,step=1/fp,length=N)
    out = zeros(length(x))
    for i in 1:length(x)
        out[i] = y(x[i])
    end
    return sum(out.^2)
end
rozwiazanie_1()

## whittaker
function rozwiazanie_2(;
    m::Vector{Float64} = [1.3, 1.3093, 1.3186, 1.3279, 1.3372, 1.3465, 1.3558, 1.3651, 1.3744, 1.3837, 1.393, 1.4023, 1.4116, 1.4209, 1.4302, 1.4395, 1.4488, 1.4581, 1.4674, 1.4767, 1.486, 1.4953, 1.5046, 1.5139, 1.5232, 1.5325, 1.5418, 1.5511, 1.5604, 1.5697, 1.579, 1.5883, 1.5976, 1.6069, 1.6162, 1.6255, 1.6348, 1.6441, 1.6534, 1.6627, 1.672, 1.6813, 1.6906, 1.6999, 1.7092, 1.7185, 1.7278, 1.7371, 1.7464, 1.7557, 1.765, 1.7743],
    s::Vector{Float64} = [0.468, 0.5912, 0.3033, 0.8138, 0.4577, 0.6599, 0.5161, 0.1859, 0.7704, 0.1431, 0.0421, 0.2117, 0.2661, 0.1776, 0.6362, 0.1396, 0.6982, 0.5936, 0.0051, 0.6195, 0.2461, 0.5844, 0.7542, 0.6688, 0.0486, 0.1126, 0.6487, 0.5459, 0.543, 0.5742, 0.4165, 0.6827, 0.8477, 0.2168, 0.9653, 0.2773, 0.7595, 0.9532, 0.423, 0.6763, 0.4473, 0.5563, 0.9264, 0.5081, 0.6717, 0.7113, 0.2017, 0.6533, 0.0463, 0.5605, 0.693, 0.9104],
    t::Vector{Float64} = [1.35394, 1.70362, 1.5976, 1.68967, 1.49437, 1.77337, 1.6813],
)
    #4.811456472683071
    T = m[2] - m[1]
    out = zeros(length(t))
    for i in 1:length(t)
        for j in 1:length(s)
            out[i]+= s[j]*sinc((t[i]-m[j])/T)
        end
    end
    return sum(out)
end
rozwiazanie_2()

## trojkatna skuteczna
function rozwiazanie_1(
    fp::Float64 = 422.39,
    t1::Float64 = 5.9,
    N::Int = 37,
)
    #1.5823603047137287
    g(t)= 4*abs(t+1/4-floor(t+1/4+1/2)) - 1
    y(t)=2.8*g(0.7*t-0.3)
    x = range(start=t1,step=1/fp,length=N)
    out = zeros(length(x))
    for i in 1:length(x)
        out[i]=y(x[i])
    end
    return sqrt(sum(out.^2)/length(out))
end
rozwiazanie_1()
##
