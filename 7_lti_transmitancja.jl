## srednie wzmocnienie
function rozwiazanie(;
    b::Vector{Float64} = [0.5597939753636629, -2.2391759014546517, 3.3587638521819776, -2.2391759014546517, 0.5597939753636629],
    a::Vector{Float64} = [1.0, -2.8548662700170597, 3.176160357993334, -1.6123057300173769, 0.31337124779083303],
    F::Vector{Float64} = [0.0, 0.02, 0.05, 0.11, 0.12],
)
    # 0.4469452099071517
    M = length(b)
    L = length(a)
    out = zeros(Float64,length(F))

    for n in 1:length(F)
        licznik = 0
        mianownik = 0
        for m in 1:M
            licznik += b[m]*cispi(2F[n])^-(m-1)
        end

        for l in 1:L
            mianownik += a[l]*cispi(2F[n])^-(l-1)
        end

        out[n] = abs(licznik/mianownik)
    end

    return sum(out)/length(out)
end

rozwiazanie()

## srednie przesuniecie fazowe

          
          
function rozwiazanie(;
    zz::Vector{ComplexF64} = ComplexF64[0.8007730344245322 + 0.5989679017597914im, 0.8007730344245322 - 0.5989679017597914im, 0.9627435618567148 + 0.2704160389167881im, 0.9627435618567148 - 0.2704160389167881im],
    pp::Vector{ComplexF64} = ComplexF64[0.17275551359950622 + 0.682065150719473im, 0.17275551359950622 - 0.682065150719473im, 0.1634593060409503 + 0.21486496675147843im, 0.1634593060409503 - 0.21486496675147843im],
    k::Float64 = 0.18223710611141924,
    F::Vector{Float64} = [0.27, 0.32, 0.5, 0.5],
)
    # 0.7787987427554584
    out = zeros(Float64,length(F))
    for n in 1:length(F)
        licznik = 1
        mianownik = 1
        for l in 1:length(zz)
            licznik *= (cispi(2*F[n])-zz[l])
        end
        for m in 1:length(pp)
            mianownik *= (cispi(2*F[n])-pp[m])
        end
        out[n] = angle(k*licznik/mianownik)
    end

    return sum(out)/length(out)
end

rozwiazanie()

## srednie przesuniecie fazowe

          
function rozwiazanie(;
    zz::Vector{ComplexF64} = ComplexF64[0.3087713233400143 - 0.9511363045761929im, 0.3087713233400143 + 0.9511363045761929im, -0.4340048224389007 - 0.9009105472241837im, -0.4340048224389007 + 0.9009105472241837im],
    pp::Vector{ComplexF64} = ComplexF64[0.6132771641136482 - 0.6909940723386615im, 0.6132771641136482 + 0.6909940723386615im, 0.6603221644043624 - 0.30811201595669935im, 0.6603221644043624 + 0.30811201595669935im],
    k::Float64 = 0.02964301816148452,
    F::Vector{Float64} = [0.04, 0.1, 0.1, 0.28, 0.37, 0.39],
)
    # -1.241568389827319
    out = zeros(Float64,length(F))
    for n in 1:length(F)
        licznik = 1
        mianownik = 1
        for l in 1:length(zz)
            licznik*=(cispi(2*F[n])-zz[l])
        end
        for m in 1:length(pp)
            mianownik*=(cispi(2*F[n])-pp[m])
        end
        out[n] = angle(k*licznik/mianownik)
    end

    return sum(out)/length(out)
end

rozwiazanie()


## srednie przesuniecie fazowe

          
function rozwiazanie(;
    zz::Vector{ComplexF64} = ComplexF64[1.0 + 0.0im, 1.0 + 0.0im, 1.0 + 0.0im, 1.0 + 0.0im, 1.0 + 0.0im, 1.0 + 0.0im],
    pp::Vector{ComplexF64} = ComplexF64[0.8956223942668234 - 0.25133617864363955im, 0.8956223942668234 + 0.25133617864363955im, 0.802065089637981 - 0.16477104702341286im, 0.802065089637981 + 0.16477104702341286im, 0.7564436771955586 - 0.05687993784022431im, 0.7564436771955586 + 0.05687993784022431im],
    k::Float64 = 0.5777931270964739,
    F::Vector{Float64} = [0.06, 0.1, 0.21, 0.41, 0.45],
)
    # -0.08531785202541793
    out = zeros(Float64,length(F))
    for n in 1:length(F)
        licznik = 1
        mianownik = 1
        for l in 1:length(zz)
            licznik *= cispi(2*F[n])-zz[l]
        end
        for m in 1:length(pp)
            mianownik *= cispi(2*F[n])-pp[m]
        end

        out[n] = angle(k*licznik/mianownik)
    end

    return sum(out)/length(out)
end

rozwiazanie()

        


        

        
        