using BenchmarkTools
using Pkg
Pkg.activate(".")
include("../src/ActiveMatter.jl")


dims = 2
μ, src, dx, dy, ddt, dnt = 1, 1, 0.5, 0.5, 0.001, 1000
D = μ

r = -200:200
r = r * dx
nx = length(r)
ny = nx

dims = 2
pd = μ, src, dx, dy, ddt, dnt, nx, ny
D = μ


u0 = zeros(nx, ny)
du = copy(u0)
u0[div(nx, 2)+1, div(nx, 2)+1] = 1

print("simulation started")
@time M = diffusion(du, u0, μ, dx, dy, ddt, dnt, nx, ny)

G(r, t) = @. exp((-r^2) / (4D * t)) / (4π * D * t)^(dims / 2)
t = ddt * dnt
pic = plot(r, M[div(nx, 2)+1, :], label = "simulation")
plot!(pic, r, G(r, t) * dx * dy, c = :green, label = "theory")
display(pic)