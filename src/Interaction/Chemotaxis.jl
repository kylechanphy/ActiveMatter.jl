
export
    Chemotaxis,
    diffusion,
    spread!,
    kernel,
    getForce,
    ∇,
    periodicbound,
    dipoleimage,
    farfield!,
    nearestImage,
    surface2,
    getForce2





"""

"""

""" define interaction """
mutable struct Chemotaxis{T} <: Interaction
    field::T
    flow::Vector{Vector{SV}}
end
function Chemotaxis(ndx::Int=5, ndy::Int=5)


    return Chemotaxis(zeros(ndx, ndy), [[SV(0.0, 0.0) for _ in 1:ndx] for _ in 1:ndy])
end



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
#     force = getForce2(u, p.pos_fold, para) * para.α
#     # return SV(0, 0)
#     return force, inter.flow
# end
"""--------------------------------"""

function getChemotaxisForce(p::ChemoDroplet, inter::Chemotaxis, para::ParaChemoDroplet, dfield, du_fold)
    inter.field, inter.flow = diffusion3(dfield, inter, p, para, du_fold)
    force = getForce2(inter.field, p.pos_fold, para) * para.α
    # return SV(0, 0)
    return force, inter.flow
end

function diffusion3(du::Matrix{Float64}, inter, p, para, pos_old)
    @unpack D, dx, dy, nx, ny, ddt, dnt = para
    @unpack pos_fold, vel, ω, src, srctype = p
    pos = pos_fold
    # @show pos
    u = inter.field
    du = copy(u)
    ff = [[SV(0, 0) for _ in 1:nx] for _ in 1:ny]
    ii, jj = Int.(round.(pos ./ SA[para.dx, para.dy])) .+ 1 ### julia array start from 1
    i_old, j_old = Int.(round.(pos_old ./ SA[para.dx, para.dy])) .+ 1

    ### current position
    pos = SA[ii, jj]
    ff = farfield!(ff, pos, para, p)

    _dx2, _dy2 = 1 / dx^2, 1 / dy^2
    _dx, _dy = 1 / dx, 1 / dy
    if srctype == "const_flux"
        flux = 1
    
        x0, y0 = pos
        # @show pos
        x1, y1 = periodicbound((x0 + 1, y0), nx, ny)
        x2, y2 = periodicbound((x0 - 1, y0), nx, ny)
        x3, y3 = periodicbound((x0, y0 + 1), nx, ny)
        x4, y4 = periodicbound((x0, y0 - 1), nx, ny)
        ls = SA[SA[x0, y0], SA[x1, y1], SA[x2, y2], SA[x3, y3], SA[x4, y4]]
    
    
        for _ in 1:dnt
            u[x0, y0] = src
            # u[i_old, j_old] = 0
    
            point_src!(du, u, ddt, D, dx, dy, nx, ny, ff, flux, ls)
            updateGrid_advection_diffusion!(du, u, pos, vel, ddt, D, dx, dy, nx, ny, ff, flux, ls)
            periodicbound!(du, u, D, ff, ddt, nx, ny, _dx, _dx2, _dy, _dy2)
            du[x0, y0] = src
    
            u, du = du, u
        end
    elseif srctype == "const_src"
        # @show ii, jj, pos
        flux = 0
        if para.dx < 1 || para.dy < 1
            x0, y0 = pos
            # @show pos
            x1, y1 = periodicbound((x0 + 1, y0), para)
            x2, y2 = periodicbound((x0 - 1, y0), para)
            x3, y3 = periodicbound((x0, y0 + 1), para)
            x4, y4 = periodicbound((x0, y0 - 1), para)
            x5, y5 = periodicbound((x0 + 1, y0 + 1), para)
            x6, y6 = periodicbound((x0 - 1, y0 - 1), para)
            x7, y7 = periodicbound((x0 + 1, y0 - 1), para)
            x8, y8 = periodicbound((x0 - 1, y0 + 1), para)
    
            ls = SA[SA[x0, y0], SA[x1, y1], SA[x2, y2], SA[x3, y3], SA[x4, y4],
                SA[x5, y5], SA[x6, y6], SA[x7, y7], SA[x8, y8]]
    
        else
            ls = SA[SA[ii, jj]]
        end
        # u[i_old, j_old] = 0
        for _ in 1:dnt
            for (x, y) in ls
                u[x, y] = src
            end
    
            updateGrid_advection_diffusion!(du, u, pos, vel, ddt, D, dx, dy, nx, ny, ff, flux, ls)
            periodicbound!(du, u, D, ff, ddt, nx, ny, _dx, _dx2, _dy, _dy2)
            #  DirichletBoundary!(du, nx, ny, 0)
            for (x, y) in ls
                du[x, y] = src
                # @show du[x, y] 
            end
    
            u, du = du, u
        end
    
    end
    return u, ff
end


function point_src!(du, u, dt, D, dx, dy, nx, ny, ff, flux, ls)
    _dx2, _dy2 = 1 / dx^2, 1 / dy^2
    _dx, _dy = 1 / dx, 1 / dy

    flux = flux / 4

    x, y = ls[1]
    du[x, y] = u[x, y] + dt * (2D * ((u[x+1, y] - u[x, y]) * _dx2 - flux * _dx) - flux * ff[x][y][1]
                               +
                               D * (u[x, y+1] - 2 * u[x, y] + u[x, y-1]) * _dy2 - (u[x, y+1] - u[x, y-1]) * 2 * _dy * ff[x][y][2])

    x, y = ls[2]
    du[x, y] = u[x, y] + dt * (2D * ((u[x-1, y] - u[x, y]) * _dx2 - flux * _dx) - flux * ff[x][y][1]
                               +
                               D * (u[x, y+1] - 2 * u[x, y] + u[x, y-1]) * _dy2 - (u[x, y+1] - u[x, y-1]) * 2 * _dy * ff[x][y][2])

    x, y = ls[3]
    du[x, y] = u[x, y] + dt * (D * (u[x+1, y] - 2 * u[x, y] + u[x-1, y]) * _dx2 - (u[x+1, y] - u[x-1, y]) * 2 * _dx * ff[x][y][1]
                               +
                               2D * ((u[x, y+1] - u[x, y]) * _dy2 - flux * _dy) - flux * ff[x][y][2])

    x, y = ls[4]
    du[x, y] = u[x, y] + dt * (D * (u[x+1, y] - 2 * u[x, y] + u[x-1, y]) * _dx2 - (u[x+1, y] - u[x-1, y]) * 2 * _dx * ff[x][y][1]
                               +
                               2D * ((u[x, y-1] - u[x, y]) * _dy2 - flux * _dy) - flux * ff[x][y][2])

end

function updateGrid_advection_diffusion!(du, u, pos, vel, dt, D, dx, dy, nx, ny, ff, flux, ls)
    _dx2, _dy2 = 1 / dx^2, 1 / dy^2
    _2dx, _2dy = 1 / 2*dx, 1 / 2*dy
    x, y = pos
    ### 2d 
    # for i in 2:nx-1
    Threads.@threads for i in 2:nx-1
        for j in 2:ny-1
            if !(SA[i, j] in ls)
                # if (i, j) != (x, y)
                du[i, j] = u[i, j] + dt * (D * (u[i+1, j] - 2 * u[i, j] + u[i-1, j]) * _dx2 - (u[i+1, j] - u[i-1, j]) * _2dx * ff[i][j][1]
                                           +
                                           D * (u[i, j+1] - 2 * u[i, j] + u[i, j-1]) * _dy2 - (u[i, j+1] - u[i, j-1]) * _2dy * ff[i][j][2])

                # du[i, j] = u[i, j] + dt * (D * (u[i+1, j] - 2 * u[i, j] + u[i-1, j]) * _dx2
                #                            +
                #                            D * (u[i, j+1] - 2 * u[i, j] + u[i, j-1]) * _dy2)

            end
        end
    end
    # @show du[x,y]
end

function periodicbound!(du, u, D, ff, dt, nx, ny, _dx, _dx2, _dy, _dy2)
    i, j = 1, 1
    du[i, j] = u[i, j] + dt * (D * (u[i+1, j] - 2 * u[i, j] + u[nx, j]) * _dx2 - (u[i+1, j] - u[nx, j]) * 2 * _dx * ff[i][j][1]
                               +
                               D * (u[i, j+1] - 2 * u[i, j] + u[i, ny]) * _dy2 - (u[i, j+1] - u[i, ny]) * 2 * _dy * ff[1][1][2])

    i, j = 1, ny
    du[i, j] = u[i, j] + dt * (D * (u[i+1, j] - 2 * u[i, j] + u[nx, j]) * _dx2 - (u[i+1, j] - u[nx, j]) * 2 * _dx * ff[i][j][1]
                               +
                               D * (u[i, 1] - 2 * u[i, j] + u[i, j-1]) * _dy2 - (u[i, 1] - u[i, j-1]) * 2 * _dy * ff[i][j][2])

    i, j = nx, 1
    du[i, j] = u[i, j] + dt * (D * (u[1, j] - 2 * u[i, j] + u[i-1, j]) * _dx2 - (u[1, j] - u[i-1, j]) * 2 * _dx * ff[i][j][1]
                               +
                               D * (u[i, j+1] - 2 * u[i, j] + u[i, ny]) * _dy2 - (u[i, j+1] - u[i, ny]) * 2 * _dy * ff[i][j][2])

    j, j = nx, ny
    du[i, j] = u[i, j] + dt * (D * (u[1, j] - 2 * u[i, j] + u[i-1, j]) * _dx2 - (u[1, j] - u[i-1, j]) * 2 * _dx * ff[i][j][1]
                               +
                               D * (u[i, 1] - 2 * u[i, j] + u[i, j-1]) * _dy2 - (u[i, 1] - u[i, j-1]) * 2 * _dy * ff[i][j][2])

    for i in 2:nx-1
        j = 1
        du[i, j] = u[i, j] + dt * (D * (u[i+1, j] - 2 * u[i, j] + u[i-1, j]) * _dx2 - (u[i+1, j] - u[i-1, j]) * 2 * _dx * ff[i][j][1]
                                   +
                                   D * (u[i, j+1] - 2 * u[i, j] + u[i, ny]) * _dy2 - (u[i, j+1] - u[i, ny]) * 2 * _dy * ff[i][j][2])

        j = ny
        du[i, j] = u[i, j] + dt * (D * (u[i+1, j] - 2 * u[i, j] + u[i-1, j]) * _dx2 - (u[i+1, j] - u[i-1, j]) * 2 * _dx * ff[i][j][1]
                                   +
                                   D * (u[i, 1] - 2 * u[i, j] + u[i, j-1]) * _dy2 - (u[i, 1] - u[i, j-1]) * 2 * _dy * ff[i][j][2])
    end

    for j in 2:ny-1
        i = 1
        du[i, j] = u[i, j] + dt * (D * (u[i+1, j] - 2 * u[i, j] + u[nx, j]) * _dx2 - (u[i+1, j] - u[nx, j]) * 2 * _dx * ff[i][j][1]
                                   +
                                   D * (u[i, j+1] - 2 * u[i, j] + u[i, j-1]) * _dy2 - (u[i, j+1] - u[i, j-1]) * 2 * _dy * ff[i][j][2])

        i = nx
        du[i, j] = u[i, j] + dt * (D * (u[1, j] - 2 * u[i, j] + u[i-1, j]) * _dx2 - (u[1, j] - u[i-1, j]) * 2 * _dx * ff[i][j][1]
                                   +
                                   D * (u[i, j+1] - 2 * u[i, j] + u[i, j-1]) * _dy2 - (u[i, j+1] - u[i, j-1]) * 2 * _dy * ff[i][j][2])
    end
end

function DirichletBoundary!(du, nx, ny, value)
    for i in 1:nx
        du[i, 1] = value
        du[i, ny] = value
    end
    for j in 1:ny
        du[1, j] = value
        du[nx, j] = value
    end
end



function getForce2(field, pos::SV, para)
    ii, jj = Int.(round.(pos ./ SA[para.dx, para.dy])) .+ 1
    force = SV(0, 0)
    # ls = urface2((ii, jj), para)
    if para.dx < 1
        step = 2
    else
        step = 1
    end

    for (i, j) in surface2((ii, jj), step, para)
        # @show (i,j)
        dx,dy = nearestImage(SA[i,j], SA[ii,jj], para.nx, para.ny).*(para.dx,para.dy)
        vec = SA[dx,dy]
        vhat = vec ./ norm(vec)
        # @show vhat
        force += (∇((i, j), field, step, para) ⋅ vhat) * vhat
        # force += ∇((i, j), field, step, para) 
    end
    # force += ∇((ii, jj), field, para)
    return force
end



function surface2(grid_id, step, para)
    ii, jj = grid_id
    nx, ny = para.nx, para.ny
    # step = 2
    pt1, pt2, pt3, pt4 = (ii + step, jj), (ii - step, jj), (ii, jj + step), (ii, jj - step)

    pt1 = periodicbound(pt1, para)
    pt2 = periodicbound(pt2, para)
    pt3 = periodicbound(pt3, para)
    pt4 = periodicbound(pt4, para)


    # # @show (pt1, pt2, pt3, pt4)
    return (pt1, pt2, pt3, pt4)
end

# function getForce8(field, pos::SV, para)
#     ii, jj = Int.(round.(pos ./ SA[para.dx, para.dy])) .+ 1
#     force = SV(0, 0)
#     surface = ((ii - 1, jj - 1), (ii - 1, jj), (ii - 1, jj + 1),
#         (ii, jj - 1), (ii, jj + 1),
#         (ii + 1, jj - 1), (ii + 1, jj), (ii + 1, jj + 1))
#     for (ii, jj) in surface
#         force += ∇((ii, jj), field, para)
#     end
#     # return SV(0,0)
#     return force
# end

# function getForce3(field, pos::SV, para)
#     ii, jj = Int.(round.(pos ./ SA[para.dx, para.dy])) .+ 1
#     force = SV(0, 0)
#     # surface = ((ii + 1, jj), (ii - 1, jj), (ii, jj + 1), (ii, jj - 1))
#     # for (ii, jj) in surface
#     force += ∇((ii, jj), field, para)
#     # end
#     # return SV(0,0)
#     return force
# end

# function getForce(field, pos::SV, para)
#     coord, ratio = kernel(pos, para)
#     force = SV(0, 0)
#     for i in eachindex(coord)
#         # force = force + ratio[i] * ∇(coord[i], field, para)
#         force = force .+ 0.25 .* ∇(coord[i], field, para)
#         # force = force +  ∇(coord[i], field, dx, dy)
#     end
#     # return SV(0,0)
#     return force
# end

# function ∇(coord::Tuple{Int64,Int64}, field, para)
#     step::Int = 1
#     @unpack nx, ny = para

#     # @show coord
#     (x, y) = coord
#     # @show typeof((x + step, y))
#     x1, y1 = periodicbound((x + step, y), nx, ny)
#     x2, y2 = periodicbound((x - step, y), nx, ny)
#     x3, y3 = periodicbound((x, y + step), nx, ny)
#     x4, y4 = periodicbound((x, y - step), nx, ny)
#     dx = (field[x1, y1] - field[x2, y2]) / (2 * step * para.dx)
#     dy = (field[x3, y3] - field[x4, y4]) / (2 * step * para.dy)

#     # @show field[x+step, y]
#     return SV(dx, dy)
# end

function ∇(coord::Tuple{Int64,Int64}, field, step, para)
    # step::Int = 2
    @unpack nx, ny = para

    # @show coord
    (x, y) = coord
    # @show typeof((x + step, y))
    x1, y1 = periodicbound((x + step, y), nx, ny)
    x2, y2 = periodicbound((x - step, y), nx, ny)
    x3, y3 = periodicbound((x, y + step), nx, ny)
    x4, y4 = periodicbound((x, y - step), nx, ny)
    dx = (field[x1, y1] - field[x2, y2]) / (2 * step * para.dx)
    dy = (field[x3, y3] - field[x4, y4]) / (2 * step * para.dy)

    # @show x1, y1
    # @show x2, y2
    # @show x3, y3
    # @show x4, y4

    # @show field[x+step, y]
    return SV(dx, dy)
end


function farfield!(field, pos, para, part)
    v = part.vel
    ω = part.ω
    nx = para.nx
    ny = para.ny
    # field = [[SV(0, 0) for _ in 1:nx] for _ in 1:ny]
    ii, jj = pos
    # @show pos
    # for pos0 in dipoleimage(pos, nx, ny)
    # images = dipoleimage(pos, nx, ny)
    Threads.@threads for i in 1:nx
        # for i in 1:nx
        for j in 1:ny
            p = SA[i, j]
            # field[i][j] += dipole2D(v, pos, p, para) + rotlet(ω, pos, p, para) + f_dipole(v, pos, p, para)
            # field[i][j] += dipole2D(v, pos, p, nx, ny)
            field[i][j] += dipole2D(v, pos, p, para) + rotlet(ω, pos, p, para)

            # field[i][j] = SV(0,0)
        end
    end
    # end
    # @show field[ii][jj]
    # @show ii, jj
    field[ii][jj] = SV(0.0, 0.0)
    # field[ii][jj] = v
    # @show field
    return field

end

function dipoleimage(pos, nx, ny)
    x0, y0 = pos
    pos1 = copy(pos)

    pos1 = SA[nx, ny]
    pos1 = pos .- SA[nx, 0]
    pos2 = pos .+ SA[nx, 0]
    pos3 = pos .- SA[0, ny]
    pos4 = pos .+ SA[0, ny]

    pos5 = pos .- SA[nx, ny]
    pos6 = pos .+ SA[nx, ny]
    pos7 = pos .+ SA[nx, -ny]
    pos8 = pos .+ SA[-nx, ny]

    return SA[pos0, pos1, pos2, pos3, pos4,
        pos5, pos6, pos7, pos8]
    # return SA[pos0, pos1, pos2, pos3, pos4]

end

function nearestImage(pos0, pos, nx, ny)
    dx, dy = pos - pos0
    nx = nx-1
    ny = ny-1 
    if dx > nx * 0.5
        dx = dx - nx - 1
        # @show pos0, pos
        # @show dx, dy
    
    elseif dx <= -nx * 0.5
        dx = dx + nx + 1
        # @show pos0, pos
        # @show dx, dy
    end

    if dy > ny * 0.5
        dy = dy - ny - 1
        # @show pos0, pos
        # @show dx, dy
    elseif dy <= -ny * 0.5
        dy = dy + ny + 1
        # @show pos0, pos
        # @show dx, dy
    end

    return dx, dy
end


""" 
source dipole

v: instantiate velocity 
pos0: grid id  of dipole 
pos: grid of space 
"""
function dipole2D(v, pos0, pos, para)
    nx, ny = para.nx, para.ny
    dx, dy = para.dx, para.dy
    stokesflow = SA[0.0, 0.0]
    D = norm(v)
    # D = 0
    e0 = atand(v[2], v[1])
    # r = norm(pos - pos0)
    # x, y = pos - pos0
    x, y = nearestImage(pos, pos0, nx, ny) .* (dx, dy)
    r = norm(SA[x, y])
    e = atand(y, x)

    vr = D * cosd(e - e0) / (2π * r^2)
    ve = D * sind(e - e0) / (2π * r^2)

    vx = vr * cosd(e) - ve * sind(e)
    vy = vr * sind(e) + ve * cosd(e)
    stokesflow += SA[vx, vy]

    return stokesflow
end


function f_dipole(v, pos0, pos, para)
    nx, ny = para.nx, para.ny
    dx, dy = para.dx, para.dy
    x, y = nearestImage(pos0, pos, nx, ny) .* (dx, dy)
    r = SA[x, y]
    # r = pos - pos0
    # @show rv
    # x, y = pos - pos0
    # e = atand(y, x)
    hate = v ./ norm(v)
    # hatr = r ./ norm(r)
    # e0 = atand(v[2], v[1])
    flow = (-r / norm(r)^2 + 2(dot(hate, r)^2) * r / norm(r)^4)
    # flow = norm(v) .* (hate + dot(hatr, hate) .* hatr) ./ (4π .* norm(r))
    # flow = norm(v) .* (hate + dot(hatr, hate) .* hatr) ./ (4norm(r))

    # return norm(v) * flow / (4π)
    return 1 * flow / (4π)
end


function rotlet(ω, pos0, pos, para)
    nx, ny = para.nx, para.ny
    dx, dy = para.dx, para.dy
    # x, y = pos - pos0
    x, y = nearestImage(pos0, pos, nx, ny) .* (dx, dy)
    r = norm(SA[x, y])

    rvec = SA[x, y, 0.0]
    flow = cross(SA[0.0, 0.0, 1.0], rvec) ./ r^2
    flow = ω .* flow ./ (4π)
    flow = SA[flow[1], flow[2]]
    return flow
end


function periodicbound(id::Tuple, para)
    nx = para.nx
    ny = para.ny
    x0, y0 = id



    if x0 > nx
        x0 =  (x0 - nx)
        # @show x0
    elseif x0 < 1
        x0 = nx - (1-x0) + 1
        # @show x0
    end

    if y0 > ny
        y0 =  (y0 - ny)
        # @show y0
    elseif y0 < 1
        y0 = ny - (1-y0) + 1
        # @show y0
    end

    return (x0, y0)
end

function periodicbound(id::Tuple, nx, ny)
    xlim = nx
    ylim = ny

    x0, y0 = id

    if x0 > nx
        x0 =  (x0 - xlim)
        # @show x0
    elseif x0 < 1
        x0 = nx - (1 - x0) +1
        # @show x0
    end

    if y0 > ny
        y0 = (y0 - ylim)
        # @show y0
    elseif y0 < 1
        y0 = ny - (1 - y0) + 1
        # @show y0
    end
    return (x0, y0)
end




""" 
# function kernel(pos::SV, para)
#     dx = para.dx
#     dy = para.dy
#     x, y = pos
#     rx, ry = rem(x, dx), rem(y, dy)
#     Ai = SA[rx*ry,
#         ry*(dx-rx),
#         (dy-ry)*(dx-rx),
#         (dy-ry)*rx]


#     A = dx * dy
#     ratio = SA[A, A, A, A] .- Ai
#     ratio = ratio ./ (3 * A)

#     # @show ratio
#     # @show A, sum(Ai) 
#     # @show sum(ratio)

#     i, j = Int(div(x, dx)) + 1, Int(div(y, dy)) + 1
#     coord = SA[(i, j), (i + 1, j), (i + 1, j + 1), (i, j + 1)]
#     # @show coord

#     return coord, ratio
# end

# function spread!(field, coord, ratio::SVector{4,Float64}, src)

#     r = ratio .* src
#     for k in eachindex(r)
#         i, j = coord[k]
#         field[i, j] = r[k]
#     end
# end


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



# function surface(pos)
#     i, j = pos
#     coords = SA[(i - 1, j - 1), (i - 1, j), (i - 1, j + 1),
#         (i, j - 1), (i, j + 1),
#         (i + 1, j - 1), (i + 1, j), (i + 1, j + 1)]
#     return coords
# end

# function farfield(field, pos, dx, dy)
#     nx, ny = size(field)
#     for i in 2:nx-1
#         for j in 2:ny-1
#             field[i, j] = norm(pos - SA[i,j])^-1
#         end
#     end
#     field[pos[1], pos[2]] = 1
#     return field
# end


# function farfield2(field, pos, v, dx, dy)
#     nx = size(field)[1]
#     ny = size(field[1])[1]
#     pos0 = (pos .- 1) .* SA[dx, dy]

#     for i in 2:nx-1
#         for j in 2:ny-1
#             p = (SA[i-1, j-1]) .* SA[dx, dy]
#             # field[i][j] = dipole2D(v, pos0, p) + stokeslet(-v, pos0, p)
#             field[i][j] =  stokeslet(-v, pos0, p)
#             # field[i][j] = SV(0,0)
#         end
#     end
#     # field[pos[1]][pos[2]] = SV(0,0)
#     # @show field
#     return field

# end


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

"""