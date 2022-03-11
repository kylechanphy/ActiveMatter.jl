using BenchmarkTools
using Pkg
Pkg.activate(".")
include("../src/ActiveMatter.jl")
gr();
default(show = true)


# D, src, dx, dy, ddt, dnt = 1, 1, 0.5, 0.5, 0.001, 100


function clfCondi(dx)
    dims = length(dx)

    return prod(dx .^ 2) / ((dims + 0.2) * sum(dx .^ 2))
end
# ddt = clfCondi([dx, dy])
# # ddt = dx^2 * dy^2 / ((2*dims + 0.1) * (dx^2 + dy^2))
# rx = 0:dx:dx*(nx-1)
# ry = 0:dy:dy*(ny-1)

dims = 2
# dx = 0.5;
# dy = 0.5;
# nx = 100;
# ny = 100;
# ddt = clfCondi([0.5, 0.5])
@with_kw mutable struct Para
    src::Float64 = 1
    D::Float64 = 1

    dx::Float64 = 0.5
    dy::Float64 = 0.5
    nx::Int = 10
    ny::Int = 10

    ddt::Float64 = clfCondi([dx,dy])
    dnt::Int = 0
end
para = Para()

@unpack dx, dy, nx, ny, ddt, dnt, D = para

rx = 0:dx:dx*(nx-1)
ry = 0:dy:dy*(ny-1)

u0 = zeros(nx, ny)
du = copy(u0)
pos = [div(nx, 2) + 1, div(ny, 2) + 1] # coord 
# u0[pos[1], pos[2]] = para.src
pos = pos.*[dx,dy] + randn(2)*dx # physical


println("simulation started")
# @time M = diffusion(du, u0, μ, dx, dy, ddt, dnt, nx, ny)
@time M = diffusion(du, u0, pos, para)

"""Matrix in julia is colum based !!!  """
# M = transpose(M)
r = rx

# G(r, t) = @. exp((-r^2) / (4D * t)) / (4π * D * t)^(dims / 2)
# t = ddt*dnt
pic = plot(r, M[div(nx, 2)+1, :], label = "simulation")
# plot!(pic, r, G(r .- r[div(nx, 2)+1], t) * dx * dy, c = :green, label = "theory")
# gui(pic)
display(pic)
