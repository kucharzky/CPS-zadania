## 36-punktowe dft, suma amplitud
          
function rozwiazanie(;
    fp::Int = 1116,
    x::Vector{ComplexF64} = ComplexF64[0.73 - 0.43im, 0.82 + 0.15im, 1.25 + 0.45im, 1.52 - 0.03im, -0.16 + 0.43im, 0.11 - 0.37im, -0.77 + 1.38im, -0.64 - 2.39im, 1.47 - 0.83im, 0.24 + 0.59im, -0.1 - 0.8im, 0.2 - 0.32im, 0.6 + 1.53im, 0.06 + 0.79im, 0.01 + 1.5im, 1.9 - 0.32im, -0.63 + 1.61im, 0.08 - 1.08im, 0.56 - 0.18im, -0.8 + 0.9im, -0.2 - 0.1im, 0.33 - 0.42im, 1.56 + 0.13im, 0.51 - 1.01im, 0.82 + 0.46im, -0.09 + 0.22im, 0.4 + 0.04im, -1.69 + 0.47im, -0.09 + 2.15im, 0.4 - 0.56im, -0.64 + 1.04im, 1.43 - 0.75im, 0.87 + 0.8im, 0.07 - 1.3im, -1.33 - 0.52im, -1.3 - 0.85im],
    f::Vector{Int} = [-31, -217, -403],
)
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
    delta_f = fp/N
    X = dft(x)
    F = zeros(Float64,N)
    F[Int(floor(N/2)+2)] = -(floor((N-1)/2))*delta_f
    out = 0

    for i in 1:N-1
        if i+1 != floor(N/2)+2
            F[i+1] = F[i]+delta_f
        end
    end
    

    for i in 1:N
        if F[i] in f
            out += abs(X[i])
        end
    end

    return out
    # 0.3683319476378634
end

rozwiazanie()

## 46-punktowe dft, suma faz
          
function rozwiazanie(;
    fp::Int = 1656,
    x::Vector{ComplexF64} = ComplexF64[0.39 + 0.89im, 0.39 - 0.68im, -1.0 + 0.64im, 1.43 - 0.36im, 0.51 + 0.7im, -1.41 - 0.02im, 0.63 + 0.06im, 0.8 - 0.08im, 1.05 - 0.11im, -0.34 + 0.66im, 0.32 - 0.09im, 0.97 + 0.13im, -0.4 - 0.97im, 0.89 + 0.61im, -0.29 + 0.63im, 0.02 + 0.15im, -0.52 + 0.39im, -0.64 - 0.19im, 0.6 + 0.27im, 1.31 - 0.18im, 2.62 + 1.5im, -0.11 + 0.61im, 0.77 - 0.78im, -0.92 + 0.77im, -0.14 + 0.14im, -0.27 - 0.51im, -0.44 + 0.66im, -0.33 + 0.05im, 0.2 + 1.19im, 1.04 + 1.04im, -0.59 + 0.41im, 0.52 + 0.28im, -1.11 - 0.39im, -0.65 + 0.9im, -0.05 + 0.97im, -0.55 + 0.21im, 0.15 + 0.01im, -0.48 - 0.88im, -0.83 + 0.11im, -1.3 - 1.43im, -0.33 - 0.07im, -0.29 + 1.34im, -0.29 - 0.29im, -0.13 + 1.73im, -0.56 + 1.4im, -0.29 - 0.03im],
    f::Vector{Int} = [-684, -540, 360, 720, -468, 576],
)
    # 4.85177043749461
    function dft(x)
        N = length(x)
        X = zeros(ComplexF64,N)
        for k in 1:N
            for n in 1:N
                X[k]+=x[n]*cispi(-2*(n-1)*(k-1)/N)
            end
        end
        return X        
    end

    N = length(x)
    X = dft(x)
    delta_f = fp/N
    F = zeros(Float64, N)
    out = 0
    F[Int(floor(N/2)+2)] = -(floor((N-1)/2))*delta_f
    for i in 1:N-1
        if i+1 != floor(N/2)+2
            F[i+1] = F[i]+delta_f
        end
    end

    for i in 1:N
        if F[i] in f
            out += angle(X[i])
        end
    end

    return out
end

rozwiazanie()

## 31-puntkowe dft, suma faz

          
function rozwiazanie(;
    fp::Int = 527,
    x::Vector{ComplexF64} = ComplexF64[0.82 - 1.61im, 0.46 - 1.06im, -0.82 + 1.19im, 0.86 + 0.7im, 1.21 - 1.63im, 0.46 - 0.03im, 0.1 - 0.3im, -0.27 + 0.7im, -0.51 - 0.92im, -0.15 - 0.26im, -0.15 - 0.42im, 0.52 + 0.55im, 0.26 - 0.26im, 0.54 + 1.5im, -1.26 - 0.42im, -0.07 + 0.95im, 1.33 + 0.59im, 0.32 - 1.04im, -0.46 + 0.21im, 0.06 - 0.41im, 0.13 + 0.11im, 0.59 - 0.68im, 0.92 - 0.2im, -0.06 - 0.64im, 0.45 - 0.72im, -0.59 + 0.72im, 0.26 + 0.43im, 0.01 - 0.3im, -1.72 - 0.44im, 0.44 - 0.27im, 0.02 - 1.05im],
    f::Vector{Int} = [85, 0, -204, 238],
)
    # 2.0330200024563085
    function dft(x)
        N = length(x)
        X = zeros(ComplexF64,N)
        for k in 1:N
            for n in 1:N
                X[k]+=x[n]*cispi(-2*(n-1)*(k-1)/N)
            end
        end
        return X        
    end

    N = length(x)
    X = dft(x)
    delta_f = fp/N
    F = zeros(Float64,N)
    F[Int(floor(N/2)+2)] = -(floor((N-1)/2))*delta_f
    out = 0
    for i in 1:N-1
        if i+1 != floor(N/2)+2
            F[i+1] = F[i]+delta_f
        end
    end

    for i in 1:N
        if F[i] in f
            out += angle(X[i])
        end
    end

    return out
end

rozwiazanie()

## 30-puntkowe dft, suma amplitud
          
function rozwiazanie(;
    fp::Int = 720,
    x::Vector{ComplexF64} = ComplexF64[-0.48 + 0.23im, -0.64 - 0.32im, 0.59 - 0.78im, 0.14 - 0.37im, 0.84 - 0.72im, 1.38 - 0.03im, 0.41 - 0.65im, -1.18 - 0.05im, 1.17 + 0.31im, -0.44 - 0.5im, -0.81 - 0.44im, 0.6 + 0.6im, -1.16 - 0.82im, -0.65 - 0.26im, 1.13 + 0.59im, -1.6 - 0.72im, -1.32 + 0.27im, 0.29 - 1.45im, -0.47 - 0.25im, 0.66 - 0.17im, 0.68 + 0.39im, -0.64 - 0.71im, -0.2 - 0.47im, 0.03 - 0.18im, 0.06 + 0.52im, -0.02 - 2.16im, -0.39 - 0.06im, 1.78 - 0.15im, 1.1 - 0.3im, -0.32 + 0.71im],
    f::Vector{Int} = [-120, 312, 192],
)
    # 0.12269746799699506
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
    out = 0

    for i in 1:N-1
        if i+1 != floor(N/2)+2
            F[i+1] = F[i]+delta_f
        end
    end

    for i in 1:N
        if F[i] in f
            out += abs(X[i])
        end
    end

    return out
end

rozwiazanie()

## 30-puntkowe dft, suma amplitud
          
function rozwiazanie(;
    fp::Int = 210,
    x::Vector{ComplexF64} = ComplexF64[-0.93 - 1.15im, 0.18 - 0.19im, 0.58 + 0.34im, 1.01 - 1.44im, -0.6 + 1.16im, 0.1 + 0.76im, -0.38 + 0.69im, 0.45 - 0.51im, -0.97 + 0.92im, -0.97 + 0.37im, -0.06 - 1.31im, -0.53 + 0.08im, 0.7 - 0.07im, -2.44 - 0.23im, -0.87 + 0.21im, -0.59 + 0.8im, 0.53 + 0.13im, 1.42 - 0.79im, 0.7 - 1.1im, 1.42 + 0.14im, -0.61 - 0.14im, -0.34 + 1.4im, 0.11 + 0.04im, -0.17 + 0.33im, -1.05 - 0.97im, -0.42 + 1.12im, 0.97 - 1.21im, -0.19 + 0.11im, 0.28 + 0.79im, 0.67 - 0.03im],
    f::Vector{Int} = [21, -7, -84],
)
    # 0.3870796117698453
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
    F = zeros(Float64,N)
    delta_f = fp/N
    out = 0
    F[Int(floor(N/2)+2)] = -(floor((N-1)/2))*delta_f

    for i in 1:N-1
        if i+1 != floor(N/2)+2
            F[i+1] = F[i]+delta_f
        end
    end

    for i in 1:N
        if F[i] in f
            out += abs(X[i])
        end
    end
    return out
end

rozwiazanie()

## 46-puntkowe dft, suma faz
          
function rozwiazanie(;
    fp::Int = 1748,
    x::Vector{ComplexF64} = ComplexF64[-0.95 - 0.22im, -0.42 + 0.89im, -0.86 - 0.57im, -0.64 - 0.6im, -0.38 + 0.28im, -0.13 + 0.31im, -0.86 - 1.47im, 0.07 - 0.01im, -0.13 - 0.41im, 0.05 - 0.7im, 0.74 - 0.66im, 0.09 - 0.33im, 1.0 + 1.42im, 0.64 + 0.59im, -0.65 + 0.61im, -0.12 + 0.66im, -0.63 + 0.05im, 0.56 - 0.01im, 0.18 + 0.92im, 0.43 + 0.32im, 0.72 + 0.16im, 0.57 + 0.42im, 0.57 + 0.63im, 0.45 - 0.97im, -0.38 + 0.34im, 0.33 - 1.24im, -1.01 + 0.43im, -0.2 - 0.12im, -0.41 - 1.2im, -0.71 - 0.21im, 0.24 - 0.25im, 0.81 + 0.94im, 0.95 - 0.43im, 1.13 + 0.1im, -0.36 - 0.76im, 0.12 - 0.83im, -0.95 - 0.51im, 2.26 - 0.18im, 1.11 - 0.71im, -0.49 - 0.46im, -0.38 + 0.25im, -0.19 + 1.51im, 0.5 - 0.41im, 0.59 + 0.71im, -0.28 + 0.26im, -0.45 + 0.07im],
    f::Vector{Int} = [-304, 76, 114, -152, 646, -190],
)
    # 1.2146918200169898
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

## 31-puntkowe dft, suma faz

function rozwiazanie(;
    fp::Int = 310,
    x::Vector{ComplexF64} = ComplexF64[-0.31 + 0.23im, 0.56 - 0.58im, -1.02 + 1.03im, -0.55 - 1.26im, -0.11 + 0.07im, 0.78 + 0.23im, -0.73 - 0.33im, 0.39 + 0.2im, 1.34 + 0.07im, 0.2 - 0.68im, 0.52 + 0.79im, 0.39 - 1.97im, 1.07 + 0.42im, -1.79 + 1.04im, -0.58 + 0.45im, -0.22 - 0.14im, 0.63 + 0.65im, -0.01 - 0.67im, -0.13 - 0.56im, 0.2 - 0.34im, 0.77 + 0.57im, -1.4 + 1.17im, 1.41 + 0.78im, -0.9 + 0.32im, 0.43 + 0.52im, 1.43 - 0.38im, 0.7 - 0.65im, 0.61 - 1.41im, 0.69 + 1.28im, -0.01 + 0.51im, -0.97 + 0.27im],
    f::Vector{Int} = [20, -50, 0],
)
    # 0.9996448945865202
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
    out = 0
    for i in 1:N-1
        if i+1 != floor(N/2)+2
            F[i+1] = F[i] + delta_f
        end
    end

    for i in 1:N
        if F[i] in f
            out += angle(X[i])
        end
    end

    return out
end

rozwiazanie()

        

        

        
        

        

        