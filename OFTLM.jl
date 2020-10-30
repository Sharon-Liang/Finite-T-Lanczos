function OFTLM(A::AbstractMatrix; R = 50, M = 90, Ne =10, Op = nothing)
    """Partition Function by the Finite Temperature Lanczos Method
       Input: A := Hamiltonian Matrix
              M := The number of Lanczos step
              R := The number of random sampling
              Ne := Number of Exact eigenstates
              temp := Temperature
              Op := A general operator
        Output: V := [E(rj),  <v psi>*<psi v>, <v psi>*<psi O v>]
                dim/R
    """
    Ee, Ve = eigs(A, nev = Ne, which =:SR)
    
    dim = size(A)[1]; 
    if Op == nothing
        n = 2
    else
        n = 3
    end
    
    V = zeros(R + Ne, M, n)
    
    fac = (dim - Ne)/R
    for r = 1: R
        T, Q = itFOLM(A, nev = M, lb = Ve)
        vals, vecs = eigen(T) 
        for j = 1:M
            V[r,j,1] = vals[j] - Ee[1]
            V[r,j,2] = vecs[1,j] * vecs[1,j]' * fac
            if Op != nothing
                V[r,j,3] = vecs[1,j] * (vecs[:,j]' * Q' * Op * Q[:, 1])
            end
        end
    end
       
    fac = 1
    for r = 1 : Ne
        V[R + r,1,1] = Ee[r] - Ee[1]
        V[R + r,1,2] = fac
        if Op != nothing
            V[R + r,1,3] = Ve[:,r]' * Op * Ve[:, r]
        end
    end
    
    return V
end


