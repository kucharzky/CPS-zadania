## fala piloksztaltna o narastajacym zboczu, podac srednia sygnalu
          
function rozwiazanie(;
    fp::Float64 = 424.31,
    t1::Float64 = -1.28,
    N::Int = 161,
)
    # -0.07303364834482504
    g(t) = 2*(t-floor(t+1/2))
    x = range(start = t1, step = 1/fp, length = N)
    y = [g(4.8*t-0.4) for t in x]

    return sum(y)/length(y)
end

rozwiazanie()

## fala piloksztaltna o opadajacym zboczu, podac wartosc skuteczna
          
function rozwiazanie(;
    fp::Float64 = 187.59,
    t1::Float64 = 0.4,
    N::Int = 774,
)
    # 1.6494416842989619
    g(t) = -2*(t-floor(t+1/2))
    x = range(start=t1,step=1/fp,length=N)
    y = [2.9*g(t-2.5) for t in x]

    return sqrt(sum(y.^2)/length(y))
end

rozwiazanie()

## fala piloksztaltna o narastajacym zboczu, podac srednia
          
function rozwiazanie(;
    fp::Float64 = 372.33,
    t1::Float64 = -0.06,
    N::Int = 934,
)
    # -0.0003953651670124931
    g(t) = 2*(t-floor(t+1/2))
    x = range(start=t1,step = 1/fp, length=N)
    y = [0.9*g(2.1*t-1) for t in x]
    return sum(y)/length(y)
end

rozwiazanie()

## fala piloksztaltna o opadajacym zboczu, wartosc skuteczna
          
function rozwiazanie(;
    fp::Float64 = 344.23,
    t1::Float64 = 1.97,
    N::Int = 908,
)
    # 1.0968833460255591
    g(t) = -2*(t-floor(t+1/2))
    x = range(start = t1, step = 1/fp, length = N)
    y = [1.9*g(4.9*t - 4.9) for t in x]
    return sqrt(sum(y.^2)/length(y))
end

rozwiazanie()

## fala trojkatna, wartosc skuteczna
          
function rozwiazanie(;
    fp::Float64 = 156.17,
    t1::Float64 = 8.16,
    N::Int = 586,
)
    # 1.332055312709623
    g(t) = 2*abs(2*((t+1/4) - floor((t+1/4)+1/2)))-1
    x = range(start = t1, step = 1/fp, length = N)
    y = [2.3*g(4.9*t-3.9) for t in x]
    return sqrt(sum(y.^2)/length(y))
end

rozwiazanie()

## bipolarna fala prostokątna, wartosc skuteczna
          
function rozwiazanie(;
    fp::Float64 = 165.29,
    t1::Float64 = 8.1,
    N::Int = 93,
)
    # 2.899999999999999
    g(t) = sign(sin(2pi*t))
    x = range(start = t1, step =1/fp, length = N)
    y = [2.9*g(0.9*t-2) for t in x]
    
    return sqrt(sum(y.^2)/length(y))
end

rozwiazanie()

## pila o narastajacym zboczu, moc

          
function rozwiazanie(;
    fp::Float64 = 205.35,
    t1::Float64 = 8.04,
    N::Int = 284,
)
    # 3.0175106925275874
    g(t) = 2*(t-floor(t+1/2))
    x = range(start = t1, step = 1/fp, length = N)
    y = [2.8*g(1.7*t-1.3) for t in x]
    return sum(y.^2)/length(y)
end

rozwiazanie()

        
        


