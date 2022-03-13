export
    ObstacleCollision
    
"""
Implement of collision between particles and obstacle.

"""


###
struct ObstacleCollision <: Interaction
    ob::ObstacleLattice
end



function collide(cell_u0::SV, cell_du::SV, forces::SV, cent_lst, ob::ObstacleLattice)
    iscollided = 0
    for c in cent_lst
        dis = norm(cell_du .- c)
        if dis <= ob.r
            iscollided = 1
            # @show c, dis
            N = c .- cell_u0
            # @show c, cell_du, p.cell_pos, cell_u0
            N = N / norm(N)
            forces = -(forces ⋅ N) * N # noramal to collison surface

            # return iscollided, forces
            break
        end
    end
    if iscollided == 1
        return iscollided, forces
    else
        return iscollided, SV(0, 0)
    end
end




function getCollsionForoce(p::Particle, cell_du::SV, forces::SV, ob::TriangularLattice)

    # check collision 
    cent_lst = getCent_lst(ob)
    iscollided, forces = collide(p.cell_pos, cell_du, forces, cent_lst, ob)
    return iscollided, forces
end

# function getCollsionForoce(p::Particle, cell_du::SV, forces::SV, ob::SquareLattice)

#     # check collision 
#     cent_lst = getCent_lst(ob)
#     iscollided, forces = collide(cell_du, forces, cent_lst, ob)
#     return iscollided, forces
# end


function getCollsionForoce(cell_u0::SV, cell_du::SV, forces::SV, cent_lst, ob::SquareLattice)

    # check collision 
    iscollided, forces = collide(cell_u0, cell_du, forces, cent_lst, ob)
    return iscollided, forces
end


function getCollsionForoce(p::Particle, du::SV, forces::SV, cent_lst, ob::FreeSpace)
    return 0, forces
end



