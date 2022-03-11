export
    ParaCAP
###


@with_kw mutable struct ParaCAP <: Parameter
    # p::Particle = Particle(0.0, 0.0)
    # ob::ObstacleLattice = FreeSpace()

    Dr::Float64 = 1
    v0::Float64 = 1
    ω0::Float64 = 1

    flow::Float64 = 0
    flow_dir::Float64 = 0

    dt::Float64 = 0.01
    n_step::Int = 100_000
end
# ParaCAP(Dr = 1, v0 = 1, ω0 = 1, n_step = 100, dt = 0.005, flow)



function clfCondi(dx)
    dims = length(dx)

    return prod(dx .^ 2) / ((dims + 0.1) * sum(dx .^ 2))
end
@with_kw mutable struct ParaChemoDroplet <: Parameter
    Dr::Float64 = 1
    v0::Float64 = 1
    ω0::Float64 = 0

    flow::Float64 = 0
    flow_dir::Float64 = 0

    dt::Float64 = 0.01
    nt::Int = 1
    n_step::Int = 100_000

    α::Float64 = 1

    ### field parameter
    D::Float64 = 1
    dx::Float64 = v*dt
    dy::Float64 = v*dt
    nx::Int = 100
    ny::Int = 100
    ddt::Float64 = clfCondi([dx,dx])
    dnt::Float64 = dt/ddt

end