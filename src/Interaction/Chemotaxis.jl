
export
    Chemotaxis,
    diffusion,
    spread!,
    kernel,
    getForce,
    ∇,
    periodicbound



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

function getChemotaxisForce(p::ChemoDroplet, inter::Chemotaxis, para::ParaChemoDroplet, dfield)
    inter.field, inter.flow = diffusion3(dfield, inter, p, para)
    force = getForce2(inter.field, p.pos_fold, para) * para.α
    # return SV(0, 0)
    return force, inter.flow
end

function diffusion3(du::Matrix{Float64}, inter, p, para)
    @unpack D, dx, dy, nx, ny, ddt, dnt = para
    @unpack pos_fold, vel, ω, src, srctype = p
    pos = pos_fold
    u = inter.field
    # ff =  copy(inter.flow)
    # ff = copy(flow)
    # ff = copy(du)
    ff = [[SV(0, 0) for _ in 1:nx] for _ in 1:ny]
    # coord, ratio = kernel(pos, para)

    ii, jj = Int.(round.(pos ./ SA[para.dx, para.dy])) .+ 1 ### julia array start from 1
    # ii, jj = periodicbound((ii, jj), para)
    # @show ii, jj
    # if ii == nx || jj == ny
    #     @show x,y
    # end
    #  if ii== 0 || jj == 0
    #     @show x,y
    # end  

    # du[ii, jj] = src

    ### current position
    pos = SA[ii, jj]
    # spread!(u, coord, ratio, src)
    ff = farfield2(ff, pos, vel, ω, dx, dy)
    # inter.flow = copy(ff)
    # if srctype == "const_src"
    # ff = farfield(ff, pos, dx,dy)
    for _ = 1:dnt
        gridUpdateAdvection!(du, u, pos, vel, ddt, D, dx, dy, nx, ny, ff)
        # spread!(du, coord, ratio, src)
        # du[ii, jj] = src
        u, du = du, u
    end
    return u, ff

    # else
    #     srctype == "free"
    #     for _ = 1:dnt
    #         gridUpdate!(du, u, pos, ddt, D, dx, dy, nx, ny)
    #         u, du = du, u
    #     end
    #     # spread!(u, coord, ratio, src)
    #     return u
    # end
end

function gridUpdateAdvection!(du, u, pos, vel, dt, D, dx, dy, nx, ny, ff)
    _dx2, _dy2 = 1 / dx^2, 1 / dy^2
    _dx, _dy = 1 / dx, 1 / dy
    x, y = pos
    flux = 1
    # 2d 
    flux = flux / 2
    # @tturbo for i in 2:nx-1
    # if x == nx || y == ny
    #     @show x,y
    # end
    #  if x == 0 || y == 0
    #     @show x,y
    # end  
    x1, y1 = periodicbound((x + 1, y), nx, ny)
    # x2, y2 = periodicbound((x - 1, y), para)
    x3, y3 = periodicbound((x, y + 1), nx, ny)
    # x4, y4 = periodicbound((x, y - 1), para)

    du[x, y] = u[x, y] + dt * (D * (2 * flux * _dx) + 2 * D * (u[x1, y1] - u[x, y]) * _dx2 + flux * ff[x][y][1]
                               +
                               D * (2 * flux * _dy) + 2 * D * (u[x3, y3] - u[x, y]) * _dy2 + flux * ff[x][y][2])
    Threads.@threads for i in 2:nx-1
        for j in 2:ny-1
            if (i, j) != (x, y)
                du[i, j] = u[i, j] + dt * (D * (u[i+1, j] - 2 * u[i, j] + u[i-1, j]) * _dx2 - (u[i+1, j] - u[i-1, j]) * 2 * _dx * ff[i][j][1]
                                           +
                                           D * (u[i, j+1] - 2 * u[i, j] + u[i, j-1]) * _dy2 - (u[i, j+1] - u[i, j-1]) * 2 * _dy * ff[i][j][2])

            end
        end
    end
    periodicbound!(du, u, D, ff, dt, nx, ny, _dx, _dx2, _dy, _dy2)

end

function periodicbound!(du, u, D, ff, dt, nx, ny, _dx, _dx2, _dy, _dy2)

    du[1, 1] = u[1, 1] + dt * (D * (u[2, 1] - 2 * u[1, 1] + u[nx, 1]) * _dx2 - (u[2, 1] - u[nx, 1]) * 2 * _dx * ff[1][1][1]
                               +
                               D * (u[1, 2] - 2 * u[1, 1] + u[1, ny]) * _dy2 - (u[1, 2] - u[1, ny]) * 2 * _dy * ff[1][1][2])

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



function getForce2(field, pos::SV, para)
    ii, jj = Int.(round.(pos ./ SA[para.dx, para.dy])) .+ 1
    # ii, jj = periodicbound((ii, jj), para)

    force = SV(0, 0)
    # surface = ((ii + 1, jj), (ii - 1, jj), (ii, jj + 1), (ii, jj - 1))
    for (ii, jj) in surface2((ii, jj), para)
        # @show ii, jj
        force += ∇((ii, jj), field, para)
    end
    # return SV(0,0)
    return force
end
function surface2(grid_id, para)
    ii, jj = grid_id
    nx, ny = para.nx, para.ny
    pt1, pt2, pt3, pt4 = (ii + 1, jj), (ii - 1, jj), (ii, jj + 1), (ii, jj - 1)
    if ii == nx
        pt1 = (1, jj)
    elseif ii == 1
        pt2 = (nx, jj)
    elseif jj == ny
        pt3 = (ii, 1)
    elseif jj == 1
        pt4 = (ii, ny)
    end

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

function ∇(coord::Tuple{Int64,Int64}, field, para)
    step::Int = 1
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

    # @show field[x+step, y]
    return SV(dx, dy)
end


function farfield2(field, pos, v, ω, dx, dy)
    nx = size(field)[1]
    ny = size(field[1])[1]
    ii, jj = pos
    # @show ii,jj
    # pos = (pos .- 1) .* SA[dx, dy]
    # pos = periodicbound(pos, dx, dy, nx, ny)

    Threads.@threads for i in 1:nx
        for j in 1:ny
            # p = (SA[i, j] .- 1) .* SA[dx, dy]
            # p = periodicbound(p, dx, dy, nx, ny)

            p = SA[i, j]

            # field[i][j] = dipole2D(v, pos, p) + f_dipole(v, pos, p) + rotlet(v, pos, p)
            # field[i][j] = dipole2D(v, pos, p) + f_dipole(v, pos, p)
            field[i][j] = dipole2D(v, pos, p) + rotlet(ω, pos, p)
            # field[i][j] = rotlet(v, pos, p)
            # field[i][j] = dipole2D(v, pos, p)
            # field[i][j] =  f_dipole(v, pos, p)

        end
    end
    # @show field[ii][jj]
    field[ii][jj] = SV(0.0, 0.0)
    # field[ii][jj] = v
    # @show field
    return field

end

# source dipole
function dipole2D(v, pos0, pos)
    r = norm(pos - pos0)
    D = norm(v)
    # D = 0.01
    x, y = pos - pos0
    e0 = atand(v[2], v[1])
    e = atand(y, x)

    vr = D * cosd(e - e0) / (2π * r^2)
    ve = D * sind(e - e0) / (2π * r^2)

    vx = vr * cosd(e) - ve * sind(e)
    vy = vr * sind(e) + ve * cosd(e)

    vx = vx
    vy = vy
    # @show SA[vx,vy]
    return SA[vx, vy]
end


function f_dipole(v, pos0, pos)
    r = pos - pos0
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


function rotlet(ω, pos0, pos)
    r = pos - pos0
    # @show rv
    # x, y = pos - pos0
    # e = atand(y, x)
    # hate = v ./ norm(v)
    # hatr = r ./ norm(r)
    # e0 = atand(v[2], v[1])
    flow = SA[r[2], -r[1]] ./ norm(r)^2
    # flow = norm(v) .* (hate + dot(hatr, hate) .* hatr) ./ (4π .* norm(r))
    # flow = norm(v) .* (hate + dot(hatr, hate) .* hatr) ./ (4norm(r))

    return -ω * flow / (4π)
end


function periodicbound(id::Tuple, para)
    xlim = para.nx
    ylim = para.ny

    x0, y0 = id
    if x0 > xlim
        x0 = x0 - xlim
    elseif x0 < 1
        x0 = xlim + x0
    end

    if y0 > ylim
        y0 = y0 - ylim
    elseif y0 < 1
        y0 = ylim + y0
    end

    return (x0, y0)
end

function periodicbound(id::Tuple, nx, ny)
    xlim = nx
    ylim = ny

    x0, y0 = id
    if x0 > xlim
        x0 = x0 - xlim
    elseif x0 < 1
        x0 = xlim + x0
    end

    if y0 > ylim
        y0 = y0 - ylim
    elseif y0 < 1
        y0 = ylim + y0
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