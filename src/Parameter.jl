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


@with_kw mutable struct ParaChemoDroplet <: Parameter
    Dr::Float64 = 1
    v0::Float64 = 1
    ω0::Float64 = 0

    flow::Float64 = 0
    flow_dir::Float64 = 0

    dt::Float64 = 0.01
    n_step::Int = 100_000

    α::Float64 = 1
    
end