include("../src/ActiveMatter.jl")

circle(R, x, y) = (θ = LinRange(0, 2π, 30);
(x .+ R .* cos.(θ), y .+ R .* sin.(θ)))
R = 0.4

c = SV(1, 1)

u0 = SV(0.3,0.8 )

F = SV(0.4, 0.4)
N = ((c - u0) / norm(c - u0))
du = u0 + F
dN = u0 + N

CF = (F ⋅ N) * N
ddu = (F - CF) + u0
ddd2 = CF + u0
function bg()
    plt = plot(circle(0.5, 1, 1), c = :red, alpha = 0.3,
        axis = nothing,
        aspect_ratio = :equal, label = "")

    scatter!(plt, [1], [1])
    scatter!(plt, [u0[1]], [u0[2]])
    scatter!(plt, [du[1]], [du[2]])
    scatter!(plt, [dN[1]], [dN[2]])
    scatter!(plt, [ddu[1]], [ddu[2]])
    scatter!(plt, [ddd2[1]], [ddd2[2]])
    return plt
end

plt = bg()