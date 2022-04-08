
export
    Chemotaxis,
    diffusion,
    spread!,
    kernel,
    getForce,
    ∇
    


"""

"""

mutable struct Chemotaxis{T} <: Interaction
    field::T
end
Chemotaxis(ndx::Int = 5, ndy::Int = 5) = Chemotaxis(zeros(ndx, ndy))
# function Chemotaxis(ndx::Int = 5, ndy::Int = 5)
#      h = 



# function getChemotaxisForce(p::ChemoDroplet, inter::Chemotaxis, para::ParaChemoDroplet, du)
#     inter.field = diffusion(du, inter.field, p, para)
#     force = getForce(inter.field, p.pos, para) * para.α
#     # return SV(0, 0)
#     return force
# end

""" test move in constanst field """
# function getChemotaxisForce(p::ChemoDroplet, inter::Chemotaxis, para::ParaChemoDroplet, du)
#     u = inter.field
#     len = size(u)[1]
#     rang = 0:len
#     for i in 1:len
#         for j in 1:len
#             u[(j-1)*len + i] = rang[i]
#         end
#      end

#     # inter.field = diffusion(du, inter.field, p, para)
#     force = getForce(u, p.pos, para) * para.α
#     # return SV(0, 0)
#     return force
# end
"""--------------------------------"""

function getChemotaxisForce(p::ChemoDroplet, inter::Chemotaxis, para::ParaChemoDroplet, du)
    inter.field = diffusion3(du, inter.field, p, para)
    force = getForce2(inter.field, p.pos, para) * para.α
    # return SV(0, 0)
    return force
end

function diffusion3(du::Matrix{Float64}, u, p, para)
    @unpack D, dx, dy, nx, ny, ddt, dnt = para
    @unpack pos, vel, src, srctype = p

    ff = copy(du)
    # coord, ratio = kernel(pos, para)
    ii,jj = Int.(round.(pos ./ SA[para.dx, para.dy])) .+ 1
    du[ii,jj] = src
    pos = SA[ii,jj]
    # spread!(u, coord, ratio, src)
    if srctype == "const_src"
        ff = farfield(ff, pos)
        for _ = 1:dnt
            gridUpdateAdvection!(du, u, pos, vel, ddt, D, dx, dy, nx, ny, ff)
            # spread!(du, coord, ratio, src)
            du[ii, jj] = src
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

function diffusion2(du::Matrix{Float64}, u, p, para)
    @unpack D, dx, dy, nx, ny, ddt, dnt = para
    @unpack pos, src, srctype = p

    # coord, ratio = kernel(pos, para)
    ii, jj = Int.(round.(pos ./ SA[para.dx, para.dy])) .+ 1
    du[ii, jj] = src
    # spread!(u, coord, ratio, src)
    if srctype == "const_src"
        for _ = 1:dnt
            gridUpdate!(du, u, pos, ddt, D, dx, dy, nx, ny)
            # spread!(du, coord, ratio, src)
            du[ii, jj] = src
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

function diffusion(du::Matrix{Float64}, u, p, para)
    @unpack D, dx, dy, nx, ny, ddt, dnt = para
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

function gridUpdateAdvection!(du, u, pos, vel, dt, D, dx, dy, nx, ny, ff)
    _dx2, _dy2 = 1 / dx^2, 1 / dy^2
    _dx, _dy = 1 / dx, 1 / dy
    @tturbo for i in 2:nx-1
#    Threads.@threads for i in 2:nx-1
       for j in 2:ny-1
           du[i, j] = u[i, j] + dt * (D * (u[i+1, j] - 2 * u[i, j] + u[i-1, j]) * _dx2 - (u[i+1, j] - u[i-1, j]) * (-vel[1] * _dx) * ff[i, j]
                                      +
                                      D * (u[i, j+1] - 2 * u[i, j] + u[i, j-1]) * _dy2 - (u[i, j+1] - u[i, j-1]) * (-vel[2] * _dy) * ff[i, j])
           # # diffusion
           # new = D * ((u[i+1, j] - 2 * u[i, j] + u[i-1, j]) * _dx2
           #            +
           #            (u[i, j+1] - 2 * u[i, j] + u[i, j-1]) * _dy2)
   
           # # advection
           # new = @. new - (ff[i, j]) * ((u[i+1, j] - u[i-1, j]) * _dx * (-vel[1])
           #                              +
           #                              (u[i, j+1] - u[i, j-1]) * _dy * (-vel[2]))
   
           # du[i, j] = u[i, j] + new * dt
           # # du[ii, jj] = src
   
   
           # diffusion
           # du[i, j] = u[i, j] + dt * D * ((u[i+1, j] - 2 * u[i, j] + u[i-1, j]) * _dx2
           #                           +
           #                           (u[i, j+1] - 2 * u[i, j] + u[i, j-1]) * _dy2)
   
           #                     -  
           #                     # advection               
           #                     ff[i,j] * dt * ((u[i+1, j] - u[i-1, j]) * _dx * (-vel[1])
           #                         +
           #                         (u[i, j+1] - u[i, j-1]) * _dy * (-vel[2]))
   
           # du[i, j] = u[i, j] + new * dt
           # du[ii, jj] = src
   
       end
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


function getForce2(field, pos::SV, para)
    ii, jj = Int.(round.(pos ./ SA[para.dx, para.dy])) .+ 1
    force = SV(0, 0)
    # for i in eachindex(coord)
    #     # force = force + ratio[i] * ∇(coord[i], field, para)
    #     force = force .+ 0.25 .* ∇(coord[i], field, para)
    #     # force = force +  ∇(coord[i], field, dx, dy)
    # end
    force = ∇((ii, jj), field, para)
    # return SV(0,0)
    return force
end


function getForce(field, pos::SV, para)
    coord, ratio = kernel(pos, para)
    force = SV(0, 0)
    for i in eachindex(coord)
        # force = force + ratio[i] * ∇(coord[i], field, para)
        force = force .+  0.25.* ∇(coord[i], field, para)
        # force = force +  ∇(coord[i], field, dx, dy)
    end
    # return SV(0,0)
    return force
end

# function getForce(field, coords, para)
#     forces = SV(0,0)
#     for i in eachindex(coords)
#         forces = @. forces +  
# end
function ∇(coord::Tuple{Int64,Int64}, field, para)
    step = 1
    x, y = coord
    dx = (field[x + step, y] - field[x - step, y]) / (step*para.dx)
    dy = (field[x, y + step] - field[x, y - step]) / (step*para.dy)

    return SV(dx, dy)
end





function kernel(pos::SV, para)
    dx = para.dx
    dy = para.dy
    x, y = pos
    rx, ry = rem(x, dx), rem(y, dy)
    Ai = SA[rx*ry,
            ry*(dx-rx),
            (dy-ry)*(dx-rx),
            (dy-ry)*rx]


    A = dx * dy
    ratio = SA[A, A, A, A] .- Ai
    ratio = ratio ./ (3* A)

    # @show ratio
    # @show A, sum(Ai) 
    # @show sum(ratio)

    i, j = Int(div(x, dx)) + 1, Int(div(y, dy)) + 1
    coord = SA[(i, j), (i+1, j), (i + 1, j + 1), (i, j+1)]
    # @show coord

    return coord, ratio
end

function spread!(field, coord, ratio::SVector{4,Float64}, src)

    r = ratio .* src
    for k in eachindex(r)
            i, j = coord[k]
            field[i, j] = r[k]
    end
end


# function farfield(ii,jj, i, j)
#     pos = SA[ii,jj]
#     fpos = SA[i,j]

#     return (norm(pos .- fpos))^-2
 
# end


# function farfield(ii, jj, i, j)
#     pos = SA[ii, jj]
#     fpos = SA[i, j]
#     if fpos == pos
#         return 1
#     else
#         return (norm(pos .- fpos))^-2
#     end
# end



function surface(pos)
    i,j = pos
    coords = SA[(i-1, j-1), (i-1, j), (i-1, j+1),
                (i, j-1), (i, j+1),
                (i+1, j-1), (i+1, j), (i+1, j+1)]
    return coords
end

function farfield(field, pos)
    nx, ny = size(field)
    for i in 2:nx-1
        for j in 2:ny-1
            field[i, j] = (norm(pos .- SA[i,j]))^-2
        end
    end
    field[pos[1], pos[2]] = 1
    return field
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