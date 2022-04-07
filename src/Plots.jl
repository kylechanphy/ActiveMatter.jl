export
    plot
    # makeMoive

###
circle(R, x, y) = (θ = LinRange(0, 2π, 30);
(x .+ R .* cos.(θ), y .+ R .* sin.(θ)))
###            


function RecipesBase.plot(ob::SquareLattice, N = 4)
    d = ob.d
    R = ob.r
    # N = N+1
    # ob = collect.(Iterators.product(-N*d+d/2:d:N*d-d/2, -N*d+d/2:d:N*d-d/2))
    ob = collect.(Iterators.product(-N*d:d:N*d, -N*d:d:N*d))
    x = [x[1] for x in ob[:]]
    y = [y[2] for y in ob[:]]

    x = reshape(x, (1, length(x)))
    y = reshape(y, (1, length(y)))
    plt = plot(circle(R, x, y), c = :red, alpha = 0.3,
        xlim = (-N * d + d / 2, N * d + d / 2), ylim = (-N * d + d / 2, N * d + d / 2),
        axis = nothing,
        aspect_ratio = :equal, label = "")


    return plt
end


function RecipesBase.plot(ob::SquareLattice, N, traj::AbstractArray)
    plt = plot(ob, N)

    x = [v[1] for v in traj]
    y = [v[2] for v in traj]

    plot!(plt, x, y, label = "")
    # scatter!(plt, [x], [y], label = "")
    # scatter!(plt, [1.0592437756635953], [1.0592437756635953], label = "")

    return plt
end


function RecipesBase.plot(traj::Vector{SV})
    x = [v[1] for v in traj]
    y = [v[2] for v in traj]

    plt = plot(x, y, label = "")

    return plt
end

function RecipesBase.plot(traj::Vector{SV3})
    x = [v[1] for v in traj]
    y = [v[2] for v in traj]
    z = [v[2] for v in traj]

    plt = plot(x, y, z, label="")

    return plt
end

function RecipesBase.plot(field::Matrix{Float64}, para::Parameter)
    hmx = (0:para.nx-1) * (para.dx)
    hmy = (0:para.ny-1) * (para.dy)
    hm = heatmap(hmx, hmy, transpose(field), aspect_ratio=1)


    return hm
end


function RecipesBase.plot(field::Matrix{Float64}, para::Parameter, traj::Vector{SV})
    # hmx = (0:para.nx) * (para.dx)
    # hmy = (0:para.ny) * (para.dy)
    hm = plot(field, para)
    # hm = heatmap(hmx, hmy, transpose(field), aspect_ratio=1)
    # hm = heatmap(hmx, hmy, field, aspect_ratio=1)
    traj = traj[1:end-1]
    x = [v[1] for v in traj]
    y = [v[2] for v in traj]
    plt = plot!(hm, x, y, label="", c=:white)
    scatter!(plt, [traj[1][1]], [traj[1][2]])

    return plt
end



# function makeMoive(inter::Chemotaxis, para, traj::Vector{SV3})
#     set_theme!(theme_black())

#     fig, ax, l = lines(points, color = colors,
#         colormap = :inferno, transparency = true,
#         axis = (; type = Axis3, protrusions = (0, 0, 0, 0),
#             viewmode = :fit, limits = (-30, 30, -30, 30, 0, 50)))

#     record(fig, "test.mp4", 1:120) do frame
#         for i in 1:50
#             # push!(points[], step!(attractor))
#             # push!(colors[], frame)
#         end
#         # ax.azimuth[] = 1.7pi + 0.3 * sin(2pi * frame / 120)
#         notify.(traj)
#         l.colorrange = (0, frame)
#     end
# end