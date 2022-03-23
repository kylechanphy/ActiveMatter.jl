using Plots
using BenchmarkTools
import Base.Threads.@threads, Base.Threads.@sync, Base.Threads.@spawn
tid = Threads.threadid
nth = Threads.nthreads()


dims = 1
μ, src, dx, dy, ddt, dnt = 1, 1, 0.5, 0.5, 0.001, 50
D = μ

r= -80:80
r = r * dx
nx = length(r)
ny = nx

dims = 2
pd = μ, src, dx, dy, ddt, dnt, nx, ny
D = μ

# u1d = zeros(nx)
# du1d = copy(u1d)
# u1d[div(nx, 2)+1] = 1

u0 = zeros(nx, ny)
du = copy(u0)
u0[div(nx, 2)+1, div(nx, 2)+1] = 1




# function diffusion!(du, u, pd)
#     μ, src, dx, dy, ddt, dnt, nx, ny = pd
#     for _ = 1:dnt
#         # @spawn begin
#             for j = 2:ny-1
#                 uj = view(u, :, j-1:j+1)
#                 duj = view(du, :, j-1:j+1)
#                 for i = 2:nx-1
#                     # u[i, j] = @. (u[i, j] + ddt * μ * ((u[i-1, j] - 2 * u[i, j] + u[i+1, j]) / dx^2
#                     #                                   +(u[i, j-1] - 2 * u[i, j] + u[i, j+1]) / dy^2))
#                     duj[i, 2] = @. uj[i, 2] + ddt * μ * ((uj[i-1, 2] - 2 * uj[i, 2] + uj[i+1, 2]) / dx^2
#                                                          +
#                                                          (uj[i, 1] - 2 * uj[i, 2] + uj[i, 3]) / dy^2)
#                 end
#             end
#         # end
#         # du[grid_x, grid_y] = src
#         u, du = du, u
#     end
#     return u
# end

function diffusion!(du, u, μ, dx, dy, ddt, dnt, nx, ny)
    # μ, src, dx, dy, ddt, dnt, nx, ny = pd
    for _ = 1:dnt
        # @spawn begin
        for j = 2:ny-1
            uj = view(u, :, j-1:j+1)
            duj = view(du, :, j-1:j+1)
            for i = 2:nx-1
                # u[i, j] = @. (u[i, j] + ddt * μ * ((u[i-1, j] - 2 * u[i, j] + u[i+1, j]) / dx^2
                #                                   +(u[i, j-1] - 2 * u[i, j] + u[i, j+1]) / dy^2))
                duj[i, 2] = @. uj[i, 2] + ddt * μ * ((uj[i-1, 2] - 2 * uj[i, 2] + uj[i+1, 2]) / dx^2
                                                     +
                                                     (uj[i, 1] - 2 * uj[i, 2] + uj[i, 3]) / dy^2)
            end
        end
        # end
        # du[grid_x, grid_y] = src
        u, du = du, u
    end
    return u
end


function diffusion2!(du, u, pd)
    μ, src, dx, dy, ddt, dnt, nx, ny = pd
    for _ = 1:dnt
        for i in 2:nx-1
            for j in 2:ny-1
                # @show i,j
                du[i, j] = u[i, j] + ddt * μ * ((u[i+1, j] - 2u[i, j] + u[i-1, j]) / dx^2
                                                +
                                                (u[i, j+1] - 2u[i, j] + u[i, j-1]) / dy^2)
            end
        end
        u, du = du, u
    end
    return u
end


function diffusion3!(du, u, pd)
    μ, src, dx, dy, ddt, dnt, nx, ny = pd
    for _ = 1:dnt
        @views du[2:end-1, 2:end-1] .= @. ((u0[1:end-2, 2:end-1] - 2u0[2:end-1, 2:end-1] + u0[3:end, 2:end-1]) / dx^2
                                           +
                                           (u0[2:end-1, 1:end-2] - 2u0[2:end-1, 2:end-1] + u0[2:end-1, 3:end]) / dy^2)
        u, du = du, u
    end
    return u
end

function diffusion1D!(du, u, pd)
    μ, src, dx, dy, ddt, dnt, nx, ny = pd
    for _ = 1:dnt
        for i in 2:nx-1
            # @show i,j
            du[i] = u[i] + ddt * μ * ((u[i+1] - 2u[i] + u[i-1]) / dx^2)
        end
        u, du = du, u
    end
    return u
end



# t = ddt * dnt


# @btime M = diffusion!($du, $u0, $pd)
# @btime M3 = diffusion2!($du, $u0, $pd)
@btime M1d = diffusion!(du, u0, μ, dx, dy, ddt, dnt, nx, ny)
# # # hm = heatmap(r, r, M)
# # pic = plot(r, M3[div(nx, 2)+1, :])
# # # pic = plot(r, M1d)

# G(r, t) = @. exp((-r^2) / (4D * t)) / (4π * D * t)^(dims / 2)
# # plot!(pic, r, G(r, t)*dx*dy, c = :green)



