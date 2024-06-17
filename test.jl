##
#1.1
function silnia_rek(n::Int64)
    if n == 0 
        return 1
    else
        return n*(silnia_rek(n-1))
    end
end

silnia_rek(20)

##
#1.1
function silnia_iter(n::Int64)
    if n == 0
        return 1
    else
        x = 1
        for i in 1:n
            x *= i
        end
        return x
    end
end
silnia_iter(10)

##
#1.3
function is_parzysta(a)
    if a%2 != 0
        return false
    else
        return true
    end
end
is_parzysta(672)

##
#1.4
function is_pierwsza(a)
    if a < 2
        return false
    end
    for i in 2:sqrt(a)
        if a%i == 0
            return false
        end
    end
    return true
end
is_pierwsza(2)

##
#1.5
function odwroc_ciag(s::String)
    return reverse(s)
end
odwroc_ciag("maciek")

##
#1.6
function is_palindrom(s::String)
    x = odwroc_ciag(s)
    if s == x
        return true
    else
        false
    end
end
is_palindrom("kajak1")

##
#1.7
function pole_sierp(n)
    if n == 0
        return 1
    else
        return (3/4)^n
    end
end
pole_sierp(2)

##
#1.8
function newton(a, epsilon = 1e-10)
    if a < 0
        return false
    end
    
    x = a 
    while abs(x * x - a) > epsilon
        x = 0.5 * (x + a / x) 
    end
    
    return x
end
newton(5)
  
##
#1.10
function is_przesunieta_pierwsza(n::Int)
    str_n = string(n)
    for _ in 1:length(str_n)
        if !is_pierwsza(parse(Int, str_n))
            return false
        end
        str_n = circshift(str_n, 1)
    end
    return true
end


liczba_przesunieta_pierwsza = 0
for num in 2:999_999
    if is_przesunieta_pierwsza(num)
        liczba_przesunieta_pierwsza += 1
    end
end

println(liczba_przesunieta_pierwsza)
##