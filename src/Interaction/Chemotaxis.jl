import Base.Threads.@threads, Base.Threads.@sync, Base.Threads.@spawn
tid = Threads.threadid
nth = Threads.nthreads()

using LoopVectorization
mutable struct Chemotaxis{T} <: Interaction
    field::T
end

Chemotaxis(ndx::Int = 5, ndy::Int = 5) = Chemotaxis(zeros(ndx, ndy))



function getChemotaxisForce(p::ChemoDroplet, inter::Chemotaxis, para::ParaChemoDroplet, du)
    inter.field = diffusion(du, inter.field, p.pos, para)
    force = getForce!(inter.field, p.pos, para)
    return force
end



function diffusion(du, u, pos, para)
    @unpack D, dx, dy, nx, ny, ddt, dnt = para

    coord, ratio = kernel(pos, para)
    spread!(u, coord, ratio, para)
    for _ = 1:dnt
        # coord, ratio = kernel(pos, para)
        # spread!(u, coord, ratio, para)
        gridUpdate!(du, u, pos, ddt, D, dx, dy, nx, ny)
        # gridUpdate!(du, u, pos, para)
        # @show pos
        # coord, ratio = kernel(pos, para)
        spread!(du, coord, ratio, para)
        # du[pos[1],pos[2]] = 1
        u, du = du, u
    end
    return u
end

function gridUpdate!(du, u, pos, dt, D, dx, dy, nx, ny)
# function gridUpdate!(du, u, pos, para)
    @tturbo for i in 2:nx-1
        # for i in 2:nx-1
        # @turbo for i in 2:nx-1
        for j in 2:ny-1
            du[i, j] = u[i, j] + dt * D * ((u[i+1, j] - 2 * u[i, j] + u[i-1, j]) / dx^2
                                                     +
                                                     (u[i, j+1] - 2 * u[i, j] + u[i, j-1]) / dy^2)
        end
    end
end


function getForce(field, pos, para)
    coord, ratio = spread!(field, pos, para)
    force = SV(0, 0)
    for i in eachindex(coord)
        force = force + ratio[i] * ∇(coord[i], field, para)
        # force = force +  ∇(coord[i], field, dx, dy)
    end
    return force
end


function ∇(coord, field, para)
    x, y = coord
    dx = (field[x+1,y] - field[x-1,y]) / (2para.dx)
    dy = (field[x,y+1] - field[x,y-1]) / (2para.dy)

    return SV(dx, dy)
end


function kernel(pos, para)
    dx = para.dx
    dy = para.dy
    x, y = pos
    rx, ry = rem(x, dx), rem(y, dy)
    ratio = [rx * ry,
        (dx - rx) * ry,
        (dx - rx) * (dy - ry),
        rx * (dy - ry)]

    ratio = (ratio ./ ((dx * dy)))
    # @show ratio
    for i in eachindex(ratio)
        r = 1 / ratio[i]
        # @show r
        if isinf(r)
            ratio = [1.0, 0.0, 0.0, 0.0]
            break
            # elseif r == 0
            #     ratio[i] = 1
        else
            ratio[i] = r
        end
    end
    ratio = ratio ./ sum(ratio)

    i, j = Int(div(x, dx)), Int(div(y, dy))
    coord = SA[(i, j), (i + 1, j), (i + 1, j + 1), (i, j + 1)]

    return coord, ratio
end

function spread!(field, coord, ratio, para)

    r = ratio * para.src
    for k in eachindex(r)
        if r[k] != 0
            i,j = coord[k]
            field[i,j] = r[k]
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