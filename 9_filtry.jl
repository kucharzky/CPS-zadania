
blackman_window(M) = [0.42 + 0.5cospi(2*n/(2M+1)) + 0.08cospi(4*n/(2M+1)) for n in -M:M]
hamming_window(M) = [0.54 + 0.46*cospi(2*n/(2M+1)) for n in -M:M]
hanning_window(M) = [0.5 + 0.5cospi(2*n/(2M+1)) for n in -M:M]
triangle_window(M) = [1-abs(n)/(M+1) for n in -M:M]


## pasmowozaporowy, okno Blackmana
  
function rozwiazanie(;
    order::Int = 98,
    fp::Float64 = 172.0,
    f1::Float64 = 15.48,
    f2::Float64 = 32.68,
    z::Vector{Int} = [74, 75, 65, 9, 43, 25],
)
    # -0.021036425845695254
    delta(n) = n == 0 ? 1 : 0
    blackman_window(N) = [0.42 + 0.5cospi(2*n/(2N+1)) + 0.08cospi(4*n/(2N+1)) for n in -N:N]

    F1 = f1/fp
    F2 = f2/fp
    n = -order/2:order/2
    h = zeros(Float64,length(n))
    for i in 1:length(n)
        h[i] = delta(n[i]) - (2*F2*sinc(2*F2*n[i]) - 2*F1*sinc(2*F1*n[i]))
    end

    w = blackman_window(order/2)
    hw = h.*w
    out = 0
    for i in z
        out += hw[i]
    end
    return out
end

rozwiazanie()

## pasmowoprzepustowy, okno Hamminga
          
function rozwiazanie(;
    order::Int = 46,
    fp::Float64 = 182.0,
    f1::Float64 = 42.77,
    f2::Float64 = 89.18,
    z::Vector{Int} = [31, 21, 28, 2, 2],
)
    # 0.16636771256903723
    hamming_window(N) = [0.54 + 0.46*cospi(2*n/(2N+1)) for n in -N:N]

    n = -order/2:order/2
    F1 = f1/fp
    F2 = f2/fp
    h = zeros(Float64,length(n))
    for i in 1:length(n)
        h[i] = 2*F2*sinc(2*F2*n[i]) - 2*F1*sinc(2*F1*n[i])
    end

    w = hamming_window(order/2)
    hw = h .* w
    out = 0
    for i in z
        out += hw[i]
    end

    return out
end

rozwiazanie()

## gornoprzepustowy, okno Hamminga

          
function rozwiazanie(;
    order::Int = 84,
    fp::Float64 = 134.0,
    f0::Float64 = 25.46,
    z::Vector{Int} = [19, 67, 22, 83],
)
    # 0.00529904694004262
    delta(n) = n == 0 ? 1 : 0
    hamming_window(N) = [0.54 + 0.46*cospi(2n/(2N+1)) for n in -N:N]

    F0 = f0/fp
    n = -order/2:order/2
    h = zeros(Float64,length(n))
    for i in 1:length(n)
        h[i] = delta(n[i]) - 2F0*sinc(2F0*n[i])
    end
    w = hamming_window(order/2)
    hw = h.*w
    out = 0
    
    for i in z
        out += hw[i]
    end
    return out
end

rozwiazanie()

## gornoprzepustowy, okno Hanninga    
          
function rozwiazanie(;
    order::Int = 86,
    fp::Float64 = 149.0,
    f0::Float64 = 16.39,
    z::Vector{Int} = [8, 16, 39, 44],
)
    # 0.7976517876318747
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

## gornoprzepustowy, okno trojkatne
          
function rozwiazanie(;
    order::Int = 100,
    fp::Float64 = 142.0,
    f0::Float64 = 51.12,
    z::Vector{Int} = [3, 16, 83, 96, 34, 86],
)
    # -0.005901056067326823
    delta(n) = n == 0 ? 1 : 0
    triangle_window(M) = [1-abs(n)/(M+1) for n in -M:M]
    n = -order/2:order/2
    h = zeros(Float64,length(n))
    F0 = f0/fp
    for i in 1:length(n)
        h[i] = delta(n[i]) - 2F0*sinc(2F0*n[i])
    end
    w = triangle_window(order/2)
    hw = h.*w
    out = 0
    for i in z
        out += hw[i]
    end

    return out
end

rozwiazanie()


## pasmowoprzepustowy, okno Blackmana

          
function rozwiazanie(;
    order::Int = 66,
    fp::Float64 = 169.0,
    f1::Float64 = 45.63,
    f2::Float64 = 65.91,
    z::Vector{Int} = [35, 9, 39, 31],
)
    # 0.007996269374684395
    blackman_window(M) = [0.42 + 0.5cospi(2*n/(2M+1)) + 0.08cospi(4*n/(2M+1)) for n in -M:M]
    F1 = f1/fp
    F2 = f2/fp
    n = -order/2:order/2
    h = zeros(Float64,length(n))
    
    for i in 1:length(n)
        h[i] = 2*F2*sinc(2*F2*n[i]) - 2*F1*sinc(2*F1*n[i])
    end

    w = blackman_window(order/2)
    hw = h.*w
    out = 0
    for i in z
        out += hw[i]
    end

    return out
end

rozwiazanie()

## pasmowoprzepustowy, okno Blackmana

          
function rozwiazanie(;
    order::Int = 68,
    fp::Float64 = 156.0,
    f1::Float64 = 52.26,
    f2::Float64 = 75.66,
    z::Vector{Int} = [1, 43, 36],
)
    # -0.23616888968610308
    blackman_window(M) = [0.42+0.5cospi(2*n/(2M+1))+0.08*cospi(4*n/(2M+1)) for n in -M:M]
    n = -order/2:order/2
    F1 = f1/fp
    F2 = f2/fp
    h = zeros(Float64,length(n))
    for i in 1:length(n)
        h[i] = 2*F2*sinc(2*F2*n[i]) - 2*F1*sinc(2*F1*n[i])
    end

    w = blackman_window(order/2)
    hw = h.*w
    out = 0

    for i in z
        out += hw[i]
    end
    
    return out
end

rozwiazanie()

        

        

        
        