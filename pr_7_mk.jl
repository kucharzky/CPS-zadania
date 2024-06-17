# TRANSMITACJA

## średnie wzmocnienie
function rozwiazanie_7(;
    b::Vector{Float64} = [5.979578037000322e-5, 0.0002989789018500161, 0.0005979578037000322, 0.0005979578037000322, 0.0002989789018500161, 5.979578037000322e-5],
    a::Vector{Float64} = [1.0, -3.984543119612337, 6.434867090275871, -5.253615170352271, 2.165132909724133, -0.35992824506355625],
    F::Vector{Float64} = [0.09, 0.11, 0.39],
)
    M = length(b)
    K = length(a)
    N = length(F)
    out = zeros(N)
    for i in 1:N
        mianownik = 0
        licznik = 0
        for m in 1:M
            licznik+=b[m]*exp(-2im*pi*F[i])^-(m-1)
        end
        for k in 1:K
            mianownik+=a[k]*exp(-2im*pi*F[i])^-(k-1)
        end
        out[i] = abs(licznik/mianownik)
    end
    return sum(out)/length(out)
end
rozwiazanie_7()