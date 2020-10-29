include("package.jl")

function icgs(u::AbstractArray, Q::AbstractArray)
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
    return u
end

function random_init(N::Int) 
    vr = rand(Float64, N) .- 0.5
    vr = vr/norm(vr)
    return vr
end

function itFOLM(A::AbstractMatrix; nev = 50, return_basis = true)
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
        r = icgs(r, Q)
        b = norm(r)
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

function BLM(A::AbstractMatrix, nev = 50, return_basis = true)
    """Block Lanczos Method
       Input: A:= Symmetric Matrix
              nev:= number of Lanczos steps
       Output: T := block tridiagonal matrix
               Q := Orthonormal basis of Krylov space
    """
    dim = size(A)[1]; 
    ncv = min(nev, dim); 
    s = 10
    while mod(ncv,s) !=0
        s -= 1
    end
    p = Int(ncv / s);
    T, Q = zeros(ncv, ncv), zeros(dim, ncv); 
    X = zeros(dim,p)
    u = random_init(dim);
    X[:,1] = u
    for i = 2:p
        u = random_init(dim);
        u = icgs(u, X)
        u = u / norm(u)
        X[:,i] = u
    end
    
    #X = X0;
    for j = 1:s
        pos = (j-1)*p
        Q[:, pos+1 : pos+p] = X
        R = A * X
        T[pos+1 : pos+p, pos+1 : pos+p] = X' * R;
        
        R = icgs(R, Q)
        F = qr(R)
        X = Matrix(F.Q)
        B = Matrix(F.R)
        if j < s
            T[pos+1 : pos+p, j*p+1 : j*p + p] = B'
            T[j*p+1 : j*p + p, pos+1 : pos+p] = B
        end
    end   
    return T, Q        
end

function BLM_old(A::AbstractMatrix, nev = 50, return_basis = true)
    """Block Lanczos Method
       Input: A:= Symmetric Matrix
              nev:= number of Lanczos steps
       Output: T := block tridiagonal matrix
               Q := Orthonormal basis of Krylov space
    """
    dim = size(A)[1]; 
    ncv = min(nev, dim); 
    s = 10
    while mod(ncv,s) !=0
        s -= 1
    end
    p = Int(ncv / s);
    T, Q = zeros(ncv, ncv), zeros(dim, ncv);
    M1, X0 = itFOLM(A, p);
    X = X0;
    for j = 1:s
        pos = (j-1)*p
        Q[:, pos+1 : pos+p] = X
        R = A * X
        T[pos+1 : pos+p, pos+1 : pos+p] = X' * R;
        
        R = icgs(R, Q)
        F = qr(R)
        X = Matrix(F.Q)
        B = Matrix(F.R)
        if j < s
            T[pos+1 : pos+p, j*p+1 : j*p + p] = B'
            T[j*p+1 : j*p + p, pos+1 : pos+p] = B
        end
    end   
    return T, Q        
end


