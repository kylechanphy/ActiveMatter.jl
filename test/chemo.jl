using BenchmarkTools
using Pkg
Pkg.activate(".")
include("../src/ActiveMatter.jl")


dims = 2
D, src, dx, dy, ddt, dnt = 1, 1, 0.5, 0.5, 0.001, 100


function clfCondi(dx)
    dims = length(dx)

    return prod(dx .^ 2) / ((2dims + 0.1) * sum(dx .^ 2))
end
ddt = clfCondi([dx,dy])
# ddt = dx^2 * dy^2 / ((2*dims + 0.1) * (dx^2 + dy^2))
nx = 50
ny = 50
rx = 0:dx:dx*(nx-1)
ry = 0:dy:dy*(ny-1)

pd = D, src, dx, dy, ddt, dnt, nx, ny


u0 = zeros(nx, ny)
du = copy(u0)
u0[div(nx, 2)+1, div(ny, 2)+1] = 1

print("simulation started")
@time M = diffusion(du, u0, μ, dx, dy, ddt, dnt, nx, ny)

"""Matrix in julia is colum based !!!  """
# M = transpose(M)
r = rx

G(r, t) = @. exp((-r^2) / (4D * t)) / (4π * D * t)^(dims / 2)
t = ddt * dnt
pic = plot(r, M[div(nx, 2)+1, :], label = "simulation")
plot!(pic, r, G(r.-r[div(nx, 2)+1], t) * dx * dy, c = :green, label = "theory")
display(pic)
