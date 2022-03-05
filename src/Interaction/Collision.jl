"""
Implement of collision between particles and obstacle.

"""


###
struct ObstacleCollision <: Interaction
    ob::ObstacleLattice
end



function collied(cell_u0::SV, cell_du::SV, forces::SV, cent_lst, ob::ObstacleLattice)
    iscollided = 0
    for c in cent_lst
        dis = norm(cell_du .- c)
        if dis <= ob.r
            iscollided = 1
            # @show c, cell_du

            N = (cell_u0 - c) ./ norm(cell_u0 - c)
            forces = -(dot(forces, N)) .* N # noramal to collison surface

            break
        end
    end
    return iscollided, SV(forces)
end




function getCollsionForoce(p::Particle, cell_du::SV, forces::SV, ob::TriangularLattice)

    # check collision 
    cent_lst = getCent_lst(ob)
    iscollided, forces = collied(p.cell_pos, cell_du, forces, cent_lst, ob)
    return iscollided, forces
end

function getCollsionForoce(p::Particle, cell_du::SV, forces::SV, ob::SquareLattice)

    # check collision 
    cent_lst = getCent_lst(ob)
    iscollided, forces = collied(p.cell_pos, cell_du, forces, cent_lst, ob)
    return iscollided, forces
end


function getCollsionForoce(p::Particle, du::SV, forces::SV,  ob::FreeSpace)
    return 0, forces
end



