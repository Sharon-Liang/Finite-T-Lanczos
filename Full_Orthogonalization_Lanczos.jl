include("package.jl")

function icgs(u::AbstractVector, Q::AbstractMatrix, return_norm = true)
    """Iterative Classical Gram-Schmidt Algorithm.
       Input: u := the vector to be orthogonalized
              Q := the orthonormal basis
       Output: u := the orthogonalized vector
    """
    a = 0.5; itmax = 3; 
    r0 = norm(u); r1 = r0
    for it = 1: itmax
        u = u - Q * (Q' * u)
        r1 = norm(u)
        if r1 > a * r0
            break
        end
        it += 1; r0 = r1
    end
    if r1 <= a * r0
        println("Warning: Loss of Orthogonality!")
    end
    if return_norm 
        return u, norm(u) 
    else 
        return u
    end
end

function random_init(N::Int) 
    vr = rand(Float64, N)
    vr = vr/norm(vr)
    return vr
end

function itFOLM(A::AbstractMatrix, nev = 50, return_basis = true)
    """Iterative Full Orthogonalized Lanczos Method
       Input: A:= Symmetric Matrix
              nev:= number of Lanczos steps
       Output: T := tridiagonal matrix
               Q := Orthonormal basis of Krylov space
    """
    dim = size(A)[1];
    ncv = min(nev, dim); v0 = random_init(dim); # random initiation vector
    T, Q = zeros(ncv, ncv), zeros(dim, ncv);
    w = v0; r = zeros(dim); k = 0;
    while k < ncv
        r += A * w; k += 1;
        Q[:, k] = w; 
        T[k,k] = w' * r; r -= T[k,k] * w
        r, b = icgs(r, Q, true)
        for i = 1:dim
            t = w[i]; w[i] = r[i] / b; r[i] = -b * t
        end
        if k < ncv
            T[k, k+1] = b; T[k+1, k] = b
        end
    end
    if return_basis
        return T, Q
    else
        return T
    end
end

function itFOLM_new(A::AbstractMatrix, nev = 50, return_basis = true)
    """Iterative Full Orthogonalized Lanczos Method
       Input: A:= Symmetric Matrix
              nev:= number of Lanczos steps
       Output: T := tridiagonal matrix
               Q := Orthonormal basis of Krylov space
    """
    dim = size(A)[1];
    ncv = min(nev, dim); v0 = random_init(dim); # random initiation vector
    T, Q = zeros(ncv, ncv), zeros(dim, ncv);
    w = v0; r = zeros(dim); k = 0;
    for k =1: ncv
        Q[:, k] = w; 
        r = A * w;
        T[k,k] = w' * r;
        r, b = icgs(r, Q, true)
        w = r/b;
        if k < ncv
            T[k, k+1] = b; T[k+1, k] = b
        end
    end
    
    #T = Q' * A * Q;
    if return_basis
        return T, Q
    else
        return T
    end
end


