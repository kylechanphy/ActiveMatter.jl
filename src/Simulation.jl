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
    cent_lst = getCent_lst(ob)

    ## check if pos0 inside the obstacle
    p.collide, _ = getCollsionForoce(p.pos, fold(p.pos, ob), SV(0, 0), ob)
    while p.collide == 1
        @show p.collide
        p.pos = p.pos + ob.d / 4 * randn(2)
        p.collide, _ = getCollsionForoce(p.pos, fold(p.pos, ob), SV(0, 0), ob)
        if p.collide == 0
            break
        end
    end

    u0 = p.pos
    du = copy(u0)

    cell_du = copy(du)
    cell_u0 = copy(u0)

    ϕ = atan(p.vel[2], p.vel[1])
    getHead(ϕ) = SV(cos(ϕ), sin(ϕ))

    ### initialize and pre-allocate logger size 
    setLogger!(logger, para)

    for i in 1:n_step
        # if norm(p.vel) > 1
        #     @show norm(p.vel), p.collide
        # end
        hat_p = getHead(ϕ)
        forces = v0 .* hat_p .+ vg
    
        cell_du = cell_u0 .+ forces .* dt
        iscollided, collided_F = getCollsionForoce(cell_u0, cell_du, forces, cent_lst, ob)
        p.collide = iscollided
    
        total_forces = forces .+ collided_F
    
        cell_du = cell_u0 .+ total_forces * dt
        cell_du = fold(cell_du, ob)
        cell_u0, cell_du = cell_du, cell_u0
        p.cell_pos = cell_u0
    
        du = u0 + total_forces * dt
        u0, du = du, u0
        p.pos = du
    
       """ check correct collision
        for c in cent_lst
            dis = norm(c - p.cell_pos)
            if dis < ob.r
                @show dis, c
            end
        end
        """
    
        ϕ += ω0 * dt + sqrt(2Dr * dt)randn()
    
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

