## dane wspolczynniki rownania roznicowego i sygnal wejsciowy, podac rms sygnalu wyjsciowego
          
function rozwiazanie(;
    b::Vector{Float64} = [0.5789906005334758, -3.48520178334795, 10.183079411110358, -18.34821771261318, 22.17749899873266, -18.34821771261318, 10.183079411110352, -3.485201783347948, 0.5789906005334755],
    a::Vector{Float64} = [1.0, -5.690228452549177, 15.750920281139942, -26.977146113349587, 31.10913066813113, -24.656521766438736, 13.167183691526109, -4.357153198749529, 0.7029712393573617],
    x::Vector{Float64} = [0.21, -0.94, 0.72, -0.16, 0.18, -0.78, -0.57, 0.6, -0.2, 0.6, -0.28, 0.95, 0.73, 0.11, 0.41, -0.22, 0.26, -0.75, -1.0, 0.05, -0.94, -0.27, 0.11, 0.7, 0.32, 0.64, 0.23, -0.31, -0.05, -0.47, -0.07, 0.33, 0.64, 0.22, 0.98, 0.54, -0.7, 0.81, -0.93, -0.64],
    L::Int = 56,
)
    # 0.35100789442238917
    N = length(x)
    M = length(b)
    K = length(a)
    y = zeros(Float64,L)
    for n in 1:L
        for m in 1:M
            if n-m+1 > 0 && n-m+1 <=N
                y[n]+=b[m]*x[n-m+1]
            end
        end
        for k in 2:K
            if n-k+1 > 0 && n-k+1 <=L
                y[n]-=a[k]*y[n-k+1]
            end
        end
    end

    return sqrt(sum(y.^2)/L)
end

rozwiazanie()

## energia sygnalu wyjsciowego
          
function rozwiazanie(;
    b::Vector{Float64} = [2.1396152038135912e-5, 0.00010698076019067957, 0.00021396152038135913, 0.00021396152038135913, 0.00010698076019067957, 2.1396152038135912e-5],
    a::Vector{Float64} = [1.0, -4.187300047864399, 7.069722752792468, -6.009958148187328, 2.5704293025241, -0.44220918239961987],
    x::Vector{Float64} = [0.89, -0.21, 0.66, -0.42, -0.45, 0.5, -0.34, 0.9, 0.53, -0.64],
    L::Int = 57,
)
    # 0.13589043516400656
    y = zeros(L)
    M = length(b)
    K = length(a)
    N = length(x)

    for n in 1:L
        for m in 1:M
            if n-m+1 > 0 && n-m+1 <= N
                y[n]+=b[m]*x[n-m+1]
            end
        end

        for k in 2:K
            if n-k+1 > 0 && n-k+1 <= L
                y[n]-=a[k]*y[n-k+1]
            end
        end
    end

    return sum(y.^2)

end

rozwiazanie()

## wartosc skuteczna

          
function rozwiazanie(;
    b::Vector{Float64} = [1.3442585342965897e-8, 8.065551205779538e-8, 2.0163878014448846e-7, 2.688517068593179e-7, 2.0163878014448846e-7, 8.065551205779538e-8, 1.3442585342965897e-8],
    a::Vector{Float64} = [1.0, -5.817838325182689, 14.14056508392353, -18.378425794462224, 13.470873201877469, -5.279507504324549, 0.8643343034695056],
    x::Vector{Float64} = [0.95, -0.39, -0.69, -0.37, 0.85, 0.55, 0.96, 0.68, 0.9, 0.17, -0.49],
    L::Int = 47,
)
    # 0.0733507680956055
    y = zeros(Float64,L)
    M = length(b)
    K = length(a)
    N = length(x)

    for n in 1:L
        for m in 1:M
            if n-m+1>0 && n-m+1<=N
                y[n]+=b[m]*x[n-m+1]
            end
        end

        for k in 2:K
            if n-k+1>0 && n-k+1<=L
                y[n]-=a[k]*y[n-k+1]
            end
        end
    end

    return sqrt(sum(y.^2)/length(y))
end

rozwiazanie()


## wartosc skuteczna

          
function rozwiazanie(;
    b::Vector{Float64} = [0.018006595843758914, 0.10803957506255349, 0.2700989376563837, 0.36013191687517826, 0.2700989376563837, 0.10803957506255349, 0.018006595843758914],
    a::Vector{Float64} = [1.0, -0.5935005660740473, 0.9092169601382036, -0.29435024395670467, 0.1490450096291632, -0.020465432500182384, 0.0024764067641382887],
    x::Vector{Float64} = [-0.34, 0.69, 0.48, -0.59, -0.61, 0.18, -0.47, -0.78, -0.4, -0.86, 0.09],
    L::Int = 59,
)
    # 0.17902849829292963
    M = length(b)
    K = length(a)
    N = length(x)
    y = zeros(Float64,L)

    for n in 1:L
        for m in 1:M
            if n-m+1>0 && n-m+1<=N
                y[n]+=b[m]*x[n-m+1]
            end
        end
        for k in 2:K
            if n-k+1>0 && n-k+1<=L
                y[n]-=a[k]*y[n-k+1]
            end
        end
    end

    return sqrt(sum(y.^2)/length(y))
end

rozwiazanie()

## wartosc skuteczna

          
function rozwiazanie(;
    b::Vector{Float64} = [0.009830446855539143, -0.0757909558568951, 0.297305855056573, -0.7689652308440064, 1.4496864120195627, -2.0889018743668784, 2.3540871515319663, -2.0889018743668784, 1.4496864120195625, -0.7689652308440065, 0.297305855056573, -0.07579095585689512, 0.009830446855539143],
    a::Vector{Float64} = [1.0, -8.068057927287136, 32.84941301282283, -87.26548885206223, 167.07064336155634, -241.58986038576, 269.9246106641525, -234.5694518691555, 157.49876909869332, -79.87085173391544, 29.189406275962227, -6.959985857994676, 0.8375814860752814],
    x::Vector{Float64} = [-0.75, 0.16, 0.95, -0.44, -0.75, 0.66, 0.57, 0.52, -0.72, -0.24, 0.92, -0.59, -0.97, 0.9, -0.2, 0.85, 0.34, 0.93, -0.64, -0.14, -0.89, 0.93, -0.65, 0.32, -0.73, 0.08, 0.35, 0.64, -0.92, 0.27, 0.73, 0.74, 0.8, 0.68, -0.33, 0.39, -0.83, -0.07, 0.56, -0.73, -0.03, 0.37, -0.45, 0.7, 0.86, 0.92, 0.12, -0.19, -0.01, -0.41],
    L::Int = 77,
)
    # 0.09212895190423866
    M = length(b)
    K = length(a)
    N = length(x)
    y = zeros(Float64,L)

    for n in 1:L
        for m in 1:M
            if n-m+1 > 0 && n-m+1 <= N
                y[n]+=b[m]*x[n-m+1]
            end
        end
        for k in 2:K
            if n-k+1 > 0 && n-k+1 <= L
                y[n]-=a[k]y[n-k+1]
            end
        end
    end

    return sqrt(sum(y.^2)/length(y))
end

rozwiazanie()

## srednia 

          
function rozwiazanie(;
    b::Vector{Float64} = [0.4584720039431156, -2.5247478557188705, 6.007120062727516, -7.880077860290833, 6.007120062727517, -2.524747855718871, 0.4584720039431158],
    a::Vector{Float64} = [1.0, -4.267098280963875, 7.7575809555090585, -7.618916729022581, 4.1152994269255245, -1.0435585617421164, 0.05830375090667995],
    x::Vector{Float64} = [-0.62, -0.71, -0.53, 0.26, 0.85, -0.26, 0.71, 0.81, 0.71, -0.71, -0.65, 0.4, 0.37, -0.75, -0.11, 0.14, -0.52, -0.63, 0.33, 0.25, -0.65, -0.38, 0.32, -0.97, 0.58, 0.8, 0.43, 0.67, -0.56, -0.39, 0.44, -0.29, -0.78, 0.73, -0.27, -0.65, -0.4, 0.57, 0.57, 0.46, 0.17, -0.03, -0.8],
    L::Int = 70,
)
    # -0.015252108662856947
    M = length(b)
    K = length(a)
    N = length(x)
    y = zeros(Float64,L)
    for n in 1:L
        for m in 1:M
            if n-m+1 >0 && n-m+1 <= N
                y[n]+=b[m]*x[n-m+1]
            end
        end
        for k in 2:K
            if n-k+1 > 0 && n-k+1 <= L
                y[n]-=a[k]*y[n-k+1]
            end
        end
    end

    return sum(y)/length(y)
end

rozwiazanie()

        

        

        