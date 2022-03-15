
export
    Chemotaxis3D,
    diffusion
"""

"""
# @with_kw mutable struct Chemotaxis{T} <: Interaction
#     field::T = zeros(3,3)
# end
function Chemotaxis3D(ndx::Int = 5, ndy::Int = 5, ndz::Int = 5)
   return  Chemotaxis(zeros(ndx, ndy, ndz))
end


function getChemotaxisForce(p::ChemoDroplet3D, inter::Chemotaxis, para::ParaChemoDroplet, du)
    inter.field = diffusion(du, inter.field, p, para)
    force = getForce(inter.field, p.pos, para) * para.α
    # return SV(0, 0)
    return force
end



function diffusion(du::Array{Float64,3}, u, p, para)
    @unpack D, dx, dy, dz, nx, ny, nz, ddt, dnt = para
    @unpack pos, src, srctype = p

    coord, ratio = kernel(pos, para)
    spread!(u, coord, ratio, src)
    if srctype == "const_src"
        for _ = 1:dnt
            gridUpdate3D!(du, u, pos, ddt, D, dx, dy, dz, nx, ny, nz)
            spread!(du, coord, ratio, src)
            u, du = du, u
        end
        return u

    else
        srctype == "free"
        for _ = 1:dnt
            gridUpdate3D!(du, u, pos, ddt, D, dx, dy, dz, nx, ny, nz)
            u, du = du, u
        end
        # spread!(u, coord, ratio, src)
        return u
    end
end

function gridUpdate3D!(du, u, pos, dt, D, dx, dy, dz, nx, ny, nz)
    _dx2, _dy2, _dz2 = 1 / dx^2, 1 / dy^2, 1 / dz^2
    @tturbo for k in 2:nx-1
        for i in 2:ny-1
            for j in 2:nz-1
                du[k, i, j] = u[k, i, j] + dt * D * ((u[k, i+1, j] - 2 * u[k, i, j] + u[k, i-1, j]) * _dx2
                                                        +
                                                        (u[k, i, j+1] - 2 * u[k, i, j] + u[k, i, j-1]) * _dy2
                                                        +
                                                        (u[k+1, i, j] - 2 * u[k, i, j] + u[k-1, i, j]) * _dz2)
            end
        end
    end
end




function getForce(field, pos::SV3, para)
    coord, ratio = kernel(pos, para)
    force = SV3(0, 0, 0)
    for i in eachindex(coord)
        force = force + ratio[i] .* ∇(coord[i], field, para)
        # force = force +  ∇(coord[i], field, dx, dy)
    end
    # return SV(0,0)
    return force
end


function ∇(coord::Tuple{Int64,Int64,Int64}, field, para)
    x, y, z = coord
    dx = (field[x+1, y, z] - field[x-1, y, z]) / (2para.dx)
    dy = (field[x, y+1, z] - field[x, y-1, z]) / (2para.dy)
    dz = (field[x, y, z+1] - field[z, y, z+1]) / (2para.dz)
    return SV3(dx, dy, dz)
end


function kernel(pos::SV3, para)
    dx = para.dx
    dy = para.dy
    dz = para.dz
    x, y, z = pos
    rx, ry, rz = rem(x, dx), rem(y, dy), rem(z, dz)
    _rx = dx - rx
    _ry = dy - ry
    _rz = dz - rz

    Ai = SA[(rx*ry)*rz,
        _rx*ry*rz,
        _rx*_ry*rz,
        rx*_ry*rz,
        (rx*ry)*_rz,
        _rx*ry*_rz,
        _rx*_ry*_rz,
        rx*_ry*_rz]

    A = dx * dy * dz
    ratio = @SVector fill(A, 8)
    ratio .-= Ai
    ratio = ratio ./ (7 * A)

    # @show ratio
    i, j, k = Int(div(x, dx)) + 1, Int(div(y, dy)) + 1, Int(div(z, dz)) + 1
    coord = SA[(i, j, k), (i + 1, j, k), (i + 1, j + 1, k), (i, j + 1, k),
        (i, j, k + 1), (i + 1, j, k + 1), (i + 1, j + 1, k + 1), (i, j + 1, k + 1)]

    return coord, ratio
end



function spread!(field, coord, ratio::SVector{8,Float64}, src)
    r = ratio .* src
    for l in eachindex(r)
        if r[l] != 0
            i, j, k = coord[l]
            field[i, j, k] = r[l]
        end
    end
end









