function FTLM(A::AbstractMatrix, M::Integer, R::Integer, Op = nothing)
    """Partition Function by the Finite Temperature Lanczos Method
       Input: A := Hamiltonian Matrix
              M := The number of Lanczos step
              R := The number of random sampling
              temp := Temperature
              Op := A general operator
        Output: V := [E(rj),  <v psi>*<psi v>, <v psi>*<psi O v>]
                dim/R
    """
    dim = size(A)[1]; 
    if Op == nothing
        n = 2
    else
        n = 3
    end
    
    V = zeros(R, M, n)
    V0 = zeros(dim, R)
    for r = 1:R
        T, Q = itFOLM(A, nev = M)
        V0[:,r] = Q[:,1]
        vals, vecs = eigen(T) 
        emin = minimum(vals)
        for j = 1:M
            V[r,j,1] = vals[j] - emin
            if abs(1 -vecs[:,j]'* vecs[:,j]) > 10^-5
                println("neq")
            end
            V[r,j,2] = vecs[1,j] * vecs[1,j]'
            if Op != nothing
                V[r,j,3] = vecs[1,j] * (vecs[:,j]' * Q' * Op * v0)
            end
        end
    end
    return V, dim/R, V0      
end

function FTLM_partition(V::AbstractArray, fac::Number, t::Number)
    Z = 0.
    R, M, n = size(V);
    for r = 1:R, j = 1:M
        Z += exp(-V[r,j,1]/t) * V[r,j,2]
    end
    return Z * fac
end

function FTLM_EandC(V::AbstractArray, fac::Number, t::Number, return_c = true)
    E = 0. ;  
    R, M, n = size(V); Z = FTLM_partition(V, fac, t);
    for r = 1:R, j = 1:M
        E += V[r,j,1] * exp(-V[r,j,1]/t)* V[r,j,2]
    end
    E = E * fac / Z
    
    if return_c
        C = 0.
        for r = 1:R, j = 1:M
            C += V[r,j,1] * V[r,j,1] * exp(-V[r,j,1]/t) * V[r,j,2]
        end
        C = C * fac /(Z * t * t)
        C -= E*E /(t * t)
        return E, C
    else
        return E
    end

end
