include("src/ActiveMatter.jl")

function setPara(pack, lp, r, ob)
    latticeConst(R, pack, ob::SquareLattice) = SquareLattice(sqrt(π * (R)^2 / pack), R)
    latticeConst(R, pack, ob::TriangularLattice) = TriangularLattice(sqrt(π * R^2 / (pack * sin(π / 3))), R)
    latticeConst(R, pack, ob::FreeSpace) = FreeSpace()

    ### some defult values
    R = 1
    ob = latticeConst(R, pack, ob)
    pos0 = [0.5, 0.5] .+ randn(2) * ob.d / 4
    ϕ0 = 2π * rand()
    v0 = 1
    vg = 0

    Dr = v0 / lp
    ω0 = v0 / r

    p = Particle(pos0, v0, ϕ0)

    # n_step = 100_000

    dt = 0.005 * R / (v0 + vg)
    n_step = 100 / Dr / dt

    inters = ObstacleCollision(ob)
    return p, inters, ParaCAP(Dr = Dr, v0 = v0, ω0 = ω0, n_step = n_step, dt = dt), ob
end

