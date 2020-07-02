"""Tool Functions"""
function bits(i::Integer, num::Integer)
    # Chenck the i-th element of the binary representation of a dicimal number
    # the bits count from right to left
    mask = 2^(i - 1)
    if num & mask == mask
        return 1
    else
        return 0
    end
    end;

function flip(i::Int, j::Int, tag::Int)
    #=Flip the spin on i,j site.
    Inputs: tag: tag of a state
            i,j: position of spins that are flipped
    Output: The tag of new state, type: int =#
    mask = 2^(i-1) + 2^(j-1)
    return xor(tag, mask)
    end;

function append_data(j::Int, colptr::Int, count::Int, row::Array{Int,1})
    point = colptr - 1
    for pos = colptr: count
        point += 1
        if row[pos] == j
            return pos, count
            break
        end
    end
            
    if point == count
        return count, count+1
    end
    end;

"""Generate Hamiltonian"""
function apply_Kitaev(ind::Int, bond::Int, incell::Int, si::Int, K::Array{Float64,1}, la::Lattice)
    # apply the ind-th Kitaev term on bond to si state
    pos = findposition(ind, la)
    if bond != pos.atom + 1
        nind = index(pos, incell, bond, la)
        val = 0.0
        if bond == 1
            sf = flip(ind, nind, si)
            val = K[1] /4
        elseif bond == 2
            sf = flip(ind, nind, si)
            val = -K[2] * (bits(ind, si) - 0.5) * (bits(nind, si) - 0.5)
        elseif bond == 3
            sf = si
            val = K[3] * (bits(ind, si) - 0.5) * (bits(nind, si) - 0.5)
        end
    elseif bond == pos.atom + 1
        sf = si
        val = 0.0
    end
    return sf, val
    end;

function apply_J(ind::Int, bond::Int, incell::Int, isz::Int, si::Int,J::Array{Float64,1}, la::Lattice)
    # Nearest neighbor Heisenberg Hamiltonian
    pos = findposition(ind, la)
    if bond != pos.atom + 1
        nind = index(pos, incell, bond, la)
        val = 0.0
        if isz == 1
            sf = si
            val = J[3] * (bits(ind, si) - 0.5) * (bits(nind, si) - 0.5)
        else
            sf = flip(ind, nind, si)
            val = J[1]/4 - J[2] * (bits(ind, si) - 0.5) * (bits(nind, si) - 0.5)
        end
    elseif bond == pos.atom + 1
        sf = si
        val = 0.0
    end
    return sf, val
    end;

function Kitaev(K::Array{Float64,1}, la::Lattice)
    dim = dimension(la)
    num = sitenum(la)
    max = dim * num * 4
    error = 1.e-8
    
    col, row, data = zeros(Int, dim+1), zeros(Int, max), zeros(max)
    colptr, count = 1, 1
    
    K = K/2
    
    for si = 0: (dim-1)
        col[si+1] = colptr
        for ind = 1:num, bond = 1:3, incell = 0:1
            sf, val = apply_Kitaev(ind, bond, incell, si, K, la)
            if abs(val) > error
                pos, count = append_data(sf + 1, colptr, count, row)
                row[pos] = sf + 1
                data[pos] += val
            end
        end
        colptr = count
    end
    col[dim + 1] = count
    H = SparseMatrixCSC(dim, dim, col, row[1:count-1], data[1:count-1])
    H = sparse(H')
    return (1/2)*(H + H')
end

function Heisenberg(J::Array{Float64,1}, la::Lattice)
    dim = dimension(la)
    num = sitenum(la)
    max = dim * num * 2
    error = 1.e-8
    
    col, row, data = zeros(Int, dim+1), zeros(Int, max), zeros(max)
    colptr, count = 1, 1
    
    J = J/2
    
    for si = 0: (dim - 1)
        col[si + 1] = colptr
        for ind = 1:num, bond = 1:3, isz = 0:1, incell = 0:1
            sf, val = apply_J(ind, bond, incell, isz, si, J, la)
            if abs(val) > error
                pos, count = append_data(sf + 1, colptr, count, row)
                row[pos] = sf + 1
                data[pos] += val
            end
        end
        colptr = count
    end
    col[dim + 1] = count
    H = SparseMatrixCSC(dim, dim, col, row[1:count-1], data[1:count-1])
    H = sparse(H')
    return (1/2)*(H + H')
    end

"""Direct product"""
function eye(dim::Int)
    # make sparse identity matrix
    return sparse(Matrix{Float64}(I, dim, dim))
end

function KH_Hamiltonian(J::Array{Float64,1}, K::Array{Float64,1}, la::Lattice)
    sx = 0.5 * sparse([0 1;1 0])
    sy = 0.5 * 1im * sparse([0 -1;1 0])
    sz = 0.5 * sparse([1 0;0 -1])
    
    J = J/2; K = K/2
    num = sitenum(la)
    dim = dimension(la)
    H = spzeros(Float64, dim, dim)

    for ind = 1:num, bond = 1:3, incell = 0:1
        fac = zeros(3)
        pos = findposition(ind, la)
        if bond != pos.atom + 1
            fac[bond] = 1.0
            nind = index(pos, incell, bond, la)
            l = max(ind, nind); s = min(ind, nind)
            di = 2^(s-1); dm = 2^(l-s-1); df = 2^(num-l)
            H += (J[1] + fac[1] * K[1])* kron( eye(di), kron( sx, kron( eye(dm), kron(sx, eye(df)))))
            H += (J[2] + fac[2] * K[2])* kron( eye(di), kron( sy, kron( eye(dm), kron(sy, eye(df)))))
            H += (J[3] + fac[3] * K[3])* kron( eye(di), kron( sz, kron( eye(dm), kron(sz, eye(df)))))
        else
        end
    end
    return real(H)
    end
