#SYGNAŁY

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
## prostokatna energia
function rozwiazanie_1(;
    fp::Float64 = 246.81,
    t1::Float64 = -5.79,
    N::Int = 46,
)
    #90.15999999999991
    g(t) = sign(sin(2*pi*t))
    y(t) = 1.4*g(2.2*t-3.3)
    x = range(start=t1,step=1/fp,length=N)
    out = zeros(length(x))
    for i in 1:length(x)
        out[i]=y(x[i])
    end
    return sum(out.^2)
end
rozwiazanie_1()
## narastajaca piła średnia
function rozwiazanie_1(;
    fp::Float64 = 179.95,
    t1::Float64 = -1.05,
    N::Int = 357,
)
    #-0.05778173675694153
    g(t)=2*(t-floor(t+1/2))
    y(t)= 4.2*g(3.0*t-1.9)
    x = range(start=t1,step=1/fp,length=N)
    #y = [4.2*g(3.0*t-1.9)] for t in x]
    out=zeros(length(x))
    for i in 1:length(x)
        out[i]=y(x[i])
    end
    return sum(out)/length(out)
end
rozwiazanie_1()
## bipolar prostokat moc
function rozwiazanie_1(;
    fp::Float64 = 219.45,
    t1::Float64 = 2.49,
    N::Int = 327,
)
    #0.16000000000000003
    g(t)=sign(sin(2*pi*t))
    y(t)=0.4*g(3.0*t - 1.6)
    x = range(start=t1,step=1/fp,length=N)
    out = zeros(length(x))
    for i in 1:length(x)
        out[i]=y(x[i])
    end
    return sum(out.^2)/length(out)
end
rozwiazanie_1()
## trojkat energia
function rozwiazanie_1(;
    fp::Float64 = 119.25,
    t1::Float64 = -1.62,
    N::Int = 596,
)
    #7.9393628544974995
    g(t) = 4*abs(t+1/4-floor(t+1/4+1/2)) - 1
    y(t) = 0.2*g(4.7*t-2.6)
    x = range(start=t1,step=1/fp,length=N)
    out = zeros(length(x))
    for i in 1:length(x)
        out[i]=y(x[i])
    end
    return sum(out.^2)
end
rozwiazanie_1()


## opadająca piła średnia - coś nie tak
function rozwiazanie_1(;
    fp::Float64 = 221.43,
    t1::Float64 = 9.89,
    N::Int = 950,
)
    #1.8808523410399326e-5
    g(t) = -2*(t-floor(t+1/2))
    y(t) = 2.8 * g(1.4*t - 2.1)
    x = range(start=t1,step=1/fp,length= N)
    out = zeros(length(x))
    for i in 1:length(x)
        out[i]=y(x[i])
    end
    return sum(out)/length(out)
end
rozwiazanie_1()

## opadajaca piła energia
function rozwiazanie_1(;
    fp::Float64 = 259.13,
    t1::Float64 = 6.52,
    N::Int = 292,
)
    #1218.1069429240088
    g(t) = -2*(t-floor(t+1/2))
    y(t) = 3.5*g(2.6*t -2.0)
    x = range(start=t1,step=1/fp,length=N)
    out = zeros(length(x))
    for i in 1:length(x)
        out[i]=y(x[i])
    end
    return sum(out.^2)
end
rozwiazanie_1()

## opadajaca piła energia
function rozwiazanie_1(;
    fp::Float64 = 358.87,
    t1::Float64 = 2.39,
    N::Int = 533,
)
    #4400.716868603447
    g(t) = -2*(t-floor(t+1/2))
    y(t) = 4.0*g(0.4*t - 4.7)
    x = range(start=t1,step=1/fp,length=N)
    out = zeros(length(x))
    for i in 1:length(x)
        out[i]=y(x[i])
    end
    return sum(out.^2)
end
rozwiazanie_1()

## prostokat moc
function rozwiazanie_1(;
    fp::Float64 = 411.02,
    t1::Float64 = 0.17,
    N::Int = 841,
)
    #4.0
    g(t) = sign(sin(2*pi*t))
    y(t) = 2.0*g(2.6*t-0.1)
    x = range(start=t1,step=1/fp,length=N)
    out = zeros(length(x))
    for i in 1:length(x)
        out[i] = y(x[i])
    end
    return sum(out.^2)/length(out)
end
rozwiazanie_1()

## trojkat moc też coś zle wychodz
function rozwiazanie_1(;
    fp::Float64 = 487.67,
    t1::Float64 = -8.17,
    N::Int = 112,
)
    #1.2740002440325982
    g(t) = 4*abs(t+1/4-floor(t-3/4)) - 1
    y(t) = 2.0 * g(4.0*t-4.1)
    x = range(start=t1,step=1/fp,length=N)
    out=zeros(length(x))
    for i in 1:length(x)
        out[i] = y(x[i])
    end
    return sum(out.^2)/length(out)
end
rozwiazanie_1()

## trojkat skuteczna tez zjebane
function rozwiazanie_1(;
    fp::Float64 = 113.89,
    t1::Float64 = -0.22,
    N::Int = 931,
)
    #0.057028693748095355
    g(t) = 4*abs(t+1/4 - floor(t-3/4)) - 1
    y(t) = 0.1*g(0.6*t - 2.6)
    x = range(start=t1,step=1/fp,length=N)
    out = zeros(length(x))
    for i in 1:length(x)
        out[i]=y(x[i])
    end
    return sqrt(sum(out.^2)/length(out))
end
rozwiazanie_1()