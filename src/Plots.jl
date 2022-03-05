export
    plot

###
circle(R, x, y) = (θ = LinRange(0, 2π, 30);
                   (x .+ R .* cos.(θ), y .+ R .* sin.(θ)))
###            


function RecipesBase.plot(ob::SquareLattice, N=4)
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
        xlim = (-N * d + d/2, N * d+ d/2 ), ylim = (-N * d + d/2,  N * d + d/2),
        axis = nothing,
        aspect_ratio = :equal, label = "")


    return plt
end

function RecipesBase.plot(ob::SquareLattice, N, traj::TrajLogger)
    plt = plot(ob, N)

    x = [v[1] for v in traj.coord]
    y = [v[2] for v in traj.coord]

    plot!(plt, x, y, label = "")
    scatter!(plt, [x[1]], [y[1]], label="")

    return plt
end


function RecipesBase.plot(ob::SquareLattice, N, traj::PosLogger)
    plt = plot(ob, N)

    x = [v[1] for v in traj.coord]
    y = [v[2] for v in traj.coord]

    plot!(plt, x, y, label = "")

    return plt
end

function RecipesBase.plot(ob::SquareLattice, N, traj::AbstractArray)
    plt = plot(ob, N)

    x = [v[1] for v in traj]
    y = [v[2] for v in traj]

    plot!(plt, x, y, label = "")

    return plt
end




function RecipesBase.plot(traj::Vector{SV})
    x = [v[1] for v in traj]
    y = [v[2] for v in traj]

    plt = plot( x, y, label = "")

    return plt
end