import Base.Threads.@threads, Base.Threads.@sync, Base.Threads.@spawn
tid = Threads.threadid
nth = Threads.nthreads()

using LoopVectorization
mutable struct Chemotaxis{T} <: Interaction
    field::T
end

Chemotaxis(ndx::Int = 5, ndy::Int = 5) = Chemotaxis(zeros(ndx, ndy))




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

macro Δx2(i, j)
    esc(:(@inbounds (u[$i+1, $j] - 2 * u[$i, $j] + u[$i-1, $j]) / dx^2))
end
macro Δy2(i, j)
    esc(:(@inbounds (u[$i, $j+1] - 2 * u[$i, $j] + u[$i, $j-1]) / dy^2))
end


function gridUpdate!(du, μ, u, nx, ny, dx, dy, ddt)
    @tturbo for i in 2:nx-1
        # for i in 2:nx-1
        # for i in 2:nx-1
        for j in 2:ny-1
            du[i, j] =  u[i,j] +  ddt* μ * ((u[i+1, j] - 2 * u[i, j] + u[i-1, j]) / dx^2
                            +
                            (u[i, j+1] - 2 * u[i, j] + u[i, j-1]) / dy^2)
        end
    end
    return
end

function diffusion(du, u, μ, dx, dy, ddt, dnt, nx, ny)

    for _ = 1:dnt
        gridUpdate!(du, μ, u, nx, ny, dx, dy, ddt)
        u, du = du, u
    end
    return u
end



function ∇(coord, field,dx,dy)
    x, y = coord
    dx = (field[x+1,y] - field[x-1,y]) / (2dx)
    dy = (field[x,y+1] - field[x,y-1]) / (2dy)

    return SV(dx, dy)
end

function spread!(field, src, pos, dx, dy)
    x, y = pos
    rx, ry = rem(x, dx), rem(y, dy)
    ratio = [rx * ry,
        (dx - rx) * ry,
        (dx - rx) * (dy - ry),
        rx * (dy - ry)]

    ratio = ratio ./ ((dx * dy))
    # @show A
    for i in eachindex(ratio)
        r = 1 / ratio[i]
        if isinf(r)
            ratio[i] = 0
        else
            ratio[i] = r
        end
    end
    ratio = ratio ./ sum(ratio)
    # A = (1 ./ A) ./ sum(1 ./ A)
    # @show A
    r = ratio * src
    # @show ratio
    i, j = Int(div(x, dx)) + 1, Int(div(y, dy)) + 1
    # field[i, j] = ratio[1]
    field[i+1, j] = ratio[2]
    field[i+1, j+1] = ratio[3]
    field[i, j+1] = ratio[4]
    coord = SA[(i, j), (i + 1, j), (i + 1, j + 1), (i, j + 1)]
    return coord, ratio

    # return ratio
end




function getChemoForce(field, pos, src, dx, dy)
    coord, ratio = spread!(field, src, pos, dx, dy)
    force = SV(0,0)
    for i in eachindex(coord)
        force = force + ratio[i] * ∇(coord[i], field, dx, dy)
        # force = force +  ∇(coord[i], field, dx, dy)
    end
    return force
end

