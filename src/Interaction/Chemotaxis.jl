import Base.Threads.@threads, Base.Threads.@sync, Base.Threads.@spawn
tid = Threads.threadid
nth = Threads.nthreads()


mutable struct Chemotaxis{T} <: Interaction
    field::T
end

Chemotaxis(ndx::Int = 5, ndy::Int = 5) = Chemotaxis(zeros(ndx,ndy))




function diffusion(du, u, μ, dx, dy, ddt, dnt, nx, ny)
    for _ = 1:dnt
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
        u, du = du, u
    end
    return u
end

