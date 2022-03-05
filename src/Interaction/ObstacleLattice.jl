"""
Define the geometry of obstacle lattices 
"""


export
    ObstacleLattice,
    FreeSpace,
    SquareLattice,
    TriangularLattice




"""
Free FreeSpace
"""
struct FreeSpace <: ObstacleLattice end

function fold(pos::SV, ob::FreeSpace)
    return pos
    
end

"""
Square Lattice 
consider a square with a obstacle located at center
"""
struct SquareLattice <: ObstacleLattice
    d::Float64 # latice constant
    r::Float64 # radius of obstacle
end


function fold(pos::SV, ob::SquareLattice)
    return mod.(pos, ob.d)
end

function getCent_lst(ob::SquareLattice)
    d = ob.d
    return (SV(0, 0), SV(d, 0), SV(d,d), SV(0,d))
end

"""
Triangle Lattice 
consider parallelogram with 4 obstacle located at corner 
"""
struct TriangularLattice <: ObstacleLattice
    d::Float64
    r::Float64
end


function fold(pos::SV, ob::TriangularLattice)
    d = ob.d
    x, y = pos
    δ = y * sind(60)

    x = mod(x - δ, d) + δ
    y = mod(y, d)
    return SV(x, y)
end

function getCent_lst(ob::TriangularLattice)
    d = ob.d
    return (SV(0, 0), SV(d, 0), SV(3d / 2, d * sind(60)), SV(d / 2, d * sind(60)))
end