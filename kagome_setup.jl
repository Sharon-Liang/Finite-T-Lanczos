struct Lattice{T <: Integer}
    # Lattice size set up
    N1::T
    N2::T
    end;

struct Position{T <: Integer}
    # position of an atom in the lattice
    row::T
    col::T
    atom::T
    end;

dimension(la::Lattice) = 2^(la.N1 * la.N2 * 3)
sitenum(la::Lattice) = la.N1 * la.N2 * 3

function index(pos::Position, incell::Int, which::Int, la::Lattice)
    #= Calculate the index of the atoms in the (r,c) unit cell.
       c, l start from 0.
       A-sublattice atom = 0; B-subkattice: atom = 1; C-subkattice: atom = 2.
       incell = 0,1 corresponds atoms out or inside of the unitcell
       which = 0, 1, 2, 3 corresponds to the original atom and 
        the one linked to it via x,y,z bonds =#
    sgn = 0
    if which == pos.atom + 1
        return(println("ERROR: No such bond!"))
    elseif which == 0
        if incell == 1
            a = pos.atom
            c = pos.col % la.N1
            r = pos.row % la.N2
        elseif incell == 0
            return println("ERROR: directions should be 1,2,3.")
        end
    elseif which != pos.atom + 1
        a = 5 - (pos.atom + 1) - which
        if incell == 1
            c = pos.col % la.N1
            r = pos.row % la.N2
        elseif incell == 0
            if which == 1
                sgn = pos.atom - a
                c = (pos.col + sgn)% la.N1
                r = pos.row % la.N2
            elseif which == 2
                sgn = pos.atom + 1 - which
                c = pos.col% la.N1
                r = (pos.row + sgn)% la.N2
            elseif which == 3
                sgn = pos.atom - a
                c = (pos.col + sgn)% la.N1
                r = (pos.row - sgn)% la.N2
            end
        end
    end
   
    c < 0 ? c += la.N1 : c += 0
    r < 0 ? r += la.N2 : r += 0
    
    #println(r, c, a)
    n = r * la.N1 + c
    return 3 * n + 1 + a
    end;

function findposition(ind::Int, la::Lattice)
    #= find the position of ind in the lattice
       ind starts from 1 =#
    atom = (ind -1)%3 
    n = div((ind -1 -atom),3)
    r, c = divrem(n, la.N1)
    return Position(r, c, atom)
    end;

function abs_coo(ind::Int, la::Lattice)
    r = zeros(3,2)
    
    #n1 = [1.,0.]; n2 = [1/2, sqrt(3)/2]
    #r[1,:] = [1., sqrt(3)/5]
    #r[2,:] = [3/4, 5*sqrt(3)/12]
    #r[3,:] = [5/4, 5*sqrt(3)/12]
    
    #n1 = [1/2,sqrt(3)/2]; n2 = [1/2, -sqrt(3)/2]
    #r[1,:] = [1/4, 7*sqrt(3)/12]
    #r[2,:] = [-1/4, -7*sqrt(3)/12]
    #r[3,:] = [0, 5*sqrt(3)/6]
    
    n1 = [1/2,sqrt(3)/2]; n2 = [1/2, -sqrt(3)/2]
    r[1,:] = [1/4, sqrt(3)/4]
    r[2,:] = [-1/4, sqrt(3)/4]
    r[3,:] = [0, sqrt(3)/2]
    
    pos = findposition(ind,la)
    coo = zeros(2)
    coo += pos.row * n2 + pos.col * n1
    coo += r[pos.atom+1, :]
    return coo
end
