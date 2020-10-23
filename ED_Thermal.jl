"""Thermaldynamic quanties"""

function partitian(T::Float64, eng::Array{Float64,1}) 
    num = length(eng)
    z = 0.0
    for i = 1 : num
        z += exp(-(eng[i]-eng[1])/T)
    end
    return z
    end;

function energy(T::Float64, eng::Array{Float64,1}) 
    num = length(eng)
    e = 0.0
    z = partitian(T, eng)
    for i = 1 : num
        e += (eng[i]-eng[1]) * exp(-(eng[i]-eng[1])/T)
    end
    e = e / z
    return e
    end;

function specific_heat(T::Float64, eng::Array{Float64,1}) 
    num = length(eng)
    c = 0.0
    for i = 1 : num
        c += (eng[i]-eng[1]) * (eng[i]-eng[1])  * exp(-(eng[i]-eng[1])/T)
    end
    c = c / ( T * T * partitian(T, eng))
    c -= (energy(T,eng) / T)^2
    return c
    end;

function entropy(T::Float64, E::Float64, Z::Float64) 
    s = energy(T,eng) / T + log( partitian(T, eng) )
    return s
    end;
