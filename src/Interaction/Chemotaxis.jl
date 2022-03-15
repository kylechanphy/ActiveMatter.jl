
export
    Chemotaxis,
    diffusion,
    spread!,
    kernel,
    kernel2

"""

"""
# mutable struct Chemotaxis{T} <: Interaction
#     field::T
# end
mutable struct Chemotaxis{T} <: Interaction
    field::T
end
Chemotaxis(ndx::Int = 5, ndy::Int = 5) = Chemotaxis(zeros(ndx, ndy))



function getChemotaxisForce(p::ChemoDroplet, inter::Chemotaxis, para::ParaChemoDroplet, du)
    inter.field = diffusion(du, inter.field, p, para)
    force = getForce(inter.field, p.pos, para) * para.α
    # return SV(0, 0)
    return force
end



function diffusion(du::Matrix{Float64}, u, p, para)
    @unpack D, dx, dy, nx, ny, ddt, dnt = para
    @unpack pos, src, srctype = p

    coord, ratio = kernel(pos, para)
    spread!(u, coord, ratio, src)
    if srctype == "const_src"
        for _ = 1:dnt
            gridUpdate!(du, u, pos, ddt, D, dx, dy, nx, ny)
            spread!(du, coord, ratio, src)
            u, du = du, u
        end
        return u

    else
        srctype == "free"
        for _ = 1:dnt
            gridUpdate!(du, u, pos, ddt, D, dx, dy, nx, ny)
            u, du = du, u
        end
        # spread!(u, coord, ratio, src)
        return u
    end
end

function gridUpdate!(du, u, pos, dt, D, dx, dy, nx, ny)
    _dx2, _dy2 = 1 / dx^2, 1 / dy^2
    @tturbo for i in 2:nx-1
        # for i in 2:nx-1
        # @turbo for i in 2:nx-1
        for j in 2:ny-1
            du[i, j] = u[i, j] + dt * D * ((u[i+1, j] - 2 * u[i, j] + u[i-1, j]) * _dx2
                                           +
                                           (u[i, j+1] - 2 * u[i, j] + u[i, j-1]) * _dy2)
        end
    end
end




function getForce(field, pos::SV, para)
    coord, ratio = kernel(pos, para)
    force = SV(0, 0)
    for i in eachindex(coord)
        force = force + ratio[i] * ∇(coord[i], field, para)
        # force = force +  ∇(coord[i], field, dx, dy)
    end
    # return SV(0,0)
    return force
end


function ∇(coord::Tuple{Int64,Int64}, field, para)
    x, y = coord
    dx = (field[x+1, y] - field[x-1, y]) / (2para.dx)
    dy = (field[x, y+1] - field[x, y-1]) / (2para.dy)

    return SV(dx, dy)
end


# function kernel(pos::SV, para)
#     dx = para.dx
#     dy = para.dy
#     x, y = pos
#     rx, ry = rem(x, dx), rem(y, dy)
#     ratio = SA[rx*ry,
#         (dx-rx)*ry,
#         (dx-rx)*(dy-ry),
#         rx*(dy-ry)]

#     ratio = ratio ./ ((dx * dy))
#     # @show sum(ratio)
#     # @show ratio
#     ratio = 1 ./ ratio
#     for r in ratio
#         # @show r
#         if isinf(r)
#             ratio = SA[1.0, 0.0, 0.0, 0.0]
#             break
#         end
#     end
#     ratio = ratio ./ sum(ratio)
#     # @show ratio
#     i, j = Int(div(x, dx)) + 1, Int(div(y, dy)) + 1
#     coord = SA[(i, j), (i + 1, j), (i + 1, j + 1), (i, j + 1)]

#     return coord, ratio
# end



function kernel(pos::SV, para)
    dx = para.dx
    dy = para.dy
    x, y = pos
    rx, ry = rem(x, dx), rem(y, dy)
    Ai = SA[rx*ry,
        (dx-rx)*ry,
        (dx-rx)*(dy-ry),
        rx*(dy-ry)]

    A = dx * dy
    ratio = SA[A, A, A, A] .- Ai
    ratio = ratio ./ (3 * A)

    # @show ratio
    i, j = Int(div(x, dx)) + 1, Int(div(y, dy)) + 1
    coord = SA[(i, j), (i + 1, j), (i + 1, j + 1), (i, j + 1)]

    return coord, ratio
end

function spread!(field, coord, ratio::SVector{4,Float64}, src)

    r = ratio .* src
    for k in eachindex(r)
        if r[k] != 0
            i, j = coord[k]
            field[i, j] = r[k]
        end
    end
end







# function diffusion(du, u, μ, dx, dy, ddt, dnt, nx, ny)
#     for _ = 1:dnt
#         for j = 2:ny-1
#             uj = view(u, :, j-1:j+1)
#             duj = view(du, :, j-1:j+1)
#             for i = 2:nx-1
#                 # u[i, j] = @. (u[i, j] + ddt * μ * ((u[i-1, j] - 2 * u[i, j] + u[i+1, j]) / dx^2
#                 #                                   +(u[i, j-1] - 2 * u[i, j] + u[i, j+1]) / dy^2))
#                 duj[i, 2] = @. uj[i, 2] + ddt * μ * ((uj[i-1, 2] - 2 * uj[i, 2] + uj[i+1, 2]) / dx^2
#                                                     +
#                                                     (uj[i, 1] - 2 * uj[i, 2] + uj[i, 3]) / dy^2)
#             end
#         end
#         u, du = du, u
#     end
#     return u
# end