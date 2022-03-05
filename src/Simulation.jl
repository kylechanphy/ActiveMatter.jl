"""
Implement of simulation 
"""

export
    runSim



###
function runSim(s::System)
    p = s.particles
    para = s.parameter
    inter = s.interactions
    logger = s.loggers
    
    Langevin!(p, para, inter, logger)

end


# single particle system / non- interacting particles system
function Langevin!(p::Particle, para::Parameter, inter::ObstacleCollision, logger)
    @unpack v0, ω0, flow, flow_dir, Dr, dt, n_step = para
    ob = inter.ob
    vg = flow * SV(cosd(flow_dir), sind(flow_dir))

    ### initialize unit cell positon
    p.cell_pos = p.pos
    cell_du = copy(p.cell_pos)

    ## check if pos0 inside the obstacle
    p.collide, _ = getCollsionForoce(p, cell_du, SV(0, 0), ob)
    while p.collide == 1
        p.pos = p.pos +  ob.d/4 * randn(2)
        p.cell_pos = fold(p.pos, ob)
        p.collide, _ = getCollsionForoce(p, p.cell_pos, SV(0, 0), ob)
        @show p.collide
        if p.collide == 0
            break
        end
    end

    du = copy(p.pos)
    ϕ = atan(p.vel[2], p.vel[1])
    getHead(ϕ) = SV(cos(ϕ), sin(ϕ))

    ### initialize and pre-allocate logger size 
    setLogger!(logger, para)

    for i in 1:n_step
    
        hat_p = getHead(ϕ)
        forces = v0 .* hat_p .+ vg
    
        cell_du = p.cell_pos .+ forces .* dt
        cell_du = fold(cell_du, ob)
        iscollided, collided_F = getCollsionForoce(p, cell_du, forces, ob)
        p.collide = iscollided
    
        total_forces = forces .+ collided_F
        # @show total_forces
        du = p.pos .+ total_forces * dt
    
        p.vel = total_forces
        p.pos, du = du, p.pos
    
        cell_du = p.cell_pos .+ total_forces * dt
        p.cell_pos = fold(cell_du, ob)

        ϕ += ω0*dt +  sqrt(2 * Dr * dt) * randn()
    
        runLogger!(logger, p, i, para::Parameter)
    end
end


# function Langevin!(p::Particle, para::Parameter, inter::ObstacleCollision, logger)
#     @unpack v0, flow, flow_dir, Dr, dt, n_step = para
#     ob = inter.ob
#     vg = flow * SV(cosd(flow_dir), sind(flow_dir))

#     du = p.pos
#     u0 = copy(du)
#     ϕ = randn()*2π
#     for i in 1:n_step
#         du = u0 + v0*SV(cos(ϕ), sin(ϕ))
#         u0, du = du, u0 
#         k, pos =  getCollsionForoce(p,du,u0,ob)
#         # @show pos
#          p.pos = pos
#         ϕ += rand()
#     end
# end

