
# export
#     Chemotaxis3D,
#     diffusion
# """

# """

# # @with_kw mutable struct Chemotaxis{T} <: Interaction
# #     field::T = zeros(3,3)
# # end
# Chemotaxis3D(ndx::Int=5, ndy::Int=5, ndz::Int=5) = Chemotaxis(zeros(ndx,ndy,ndz))
# # function Chemotaxis3D(ndx::Int=5, ndy::Int=5, ndz::Int=5)
# #     h = SharedArray{Float64,3}((ndx, ndy, ndz))
# #     return Chemotaxis(h)
# # end

# function getChemotaxisForce(p::ChemoDroplet3D, inter::Chemotaxis, para::ParaChemoDroplet, du, ff)
#     inter.field = diffusion2(du, inter.field, p, para, ff)
#     force = getForce3D(inter.field, p.pos, para) * para.α
#     # return SV(0, 0)
#     return force
#     # return 0
# end



# function diffusion2(du::Array{Float64,3}, u, p, para, ff)
#     @unpack D, dx, dy, dz, nx, ny, nz, ddt, dnt = para
#     @unpack pos, vel, src, srctype = p

#     ii, jj, kk = Int.(round.(pos ./ SA[para.dx, para.dy, para.dz])) .+ 1
#     du[ii, jj, kk] = src
#     pos = SA[ii, jj, kk]
#     farfield3D!(ff, pos)
#     # @show pos
#     if srctype == "const_src"
#         for _ = 1:dnt
#             gridUpdateAdvection3D!(du, u, pos, vel, ddt, D, dx, dy, dz, nx, ny, nz, ff)
#             du[ii, jj, kk] = src
#             u, du = du, u
#         end
#         ff = nothing
#         return u
    
#     else
#         srctype == "free"
#         for _ = 1:dnt
            
#            gridUpdate3D!(du, u, pos, ddt, D, dx, dy, dz, nx, ny, nz)
#             u, du = du, u
#         end
#         # spread!(u, coord, ratio, src)
#         ff = nothing
#         return u
#     end

# end

# # function gridUpdate3D!(du, u, pos, dt, D, dx, dy, dz, nx, ny, nz)
# #     _dx2, _dy2, _dz2 = 1 / dx^2, 1 / dy^2, 1 / dz^2
# #     @tturbo for k in 2:nx-1
# #         for i in 2:ny-1
# #             for j in 2:nz-1
# #                 du[k, i, j] = u[k, i, j] + dt * D * ((u[k, i+1, j] - 2 * u[k, i, j] + u[k, i-1, j]) * _dx2
# #                                                         +
# #                                                         (u[k, i, j+1] - 2 * u[k, i, j] + u[k, i, j-1]) * _dy2
# #                                                         +
# #                                                         (u[k+1, i, j] - 2 * u[k, i, j] + u[k-1, i, j]) * _dz2)
# #             end
# #         end
# #     end
# # end

# # function gridUpdateAdvection3D!(du, u, pos, vel, dt, D, dx, dy, dz, nx, ny, nz, ff)
# #     _dx2, _dy2, _dz2 = 1 / dx^2, 1 / dy^2, 1 / dz^2
# #     _dx, _dy, _dz = 1 / dx, 1 / dy, 1 / dz
# #     @tturbo for i in 2:nx-1
# #     for j in 2:ny-1
# #         for k in 2:nz-1
# #             du[i, j, k] = u[i, j, k] + dt * (D * (u[i+1, j, k] - 2 * u[i, j, k] + u[i-1, j, k]) * _dx2 - (u[i+1, j, k] - u[i-1, j, k]) * (-vel[1] * _dx) * ff[i, j, k]
# #                                              +
# #                                              D * (u[i, j+1, k] - 2 * u[i, j, k] + u[i, j-1, k]) * _dy2 - (u[i, j+1, k] - u[i, j-1, k]) * (-vel[2] * _dy) * ff[i, j, k]
# #                                              +
# #                                              D * (u[i, j, k+1] - 2 * u[i, j, k] + u[i, j, k-1]) * _dz2 - (u[i, j, k+1] - u[i, j, k-1]) * (-vel[3] * _dz) * ff[i, j, k])


# #         end
# #     end
# # end
# # end


# function gridUpdateAdvection3D!(du, u, pos, vel, dt, D, dx, dy, dz, nx, ny, nz, ff)
#     _dx2, _dy2, _dz2 = 1 / dx^2, 1 / dy^2, 1 / dz^2
#     _dx, _dy, _dz = 1 / dx, 1 / dy, 1 / dz
#         for i in 2:nx-1
#             @tturbo for j in 2:ny-1
#                 for k in 2:nz-1
#                     du[i, j, k] = u[i, j, k] + dt * (D * (u[i+1, j, k] - 2 * u[i, j, k] + u[i-1, j, k]) * _dx2 - (u[i+1, j, k] - u[i-1, j, k]) * (-vel[1] * _dx) * ff[i, j, k]
#                                                      +
#                                                      D * (u[i, j+1, k] - 2 * u[i, j, k] + u[i, j-1, k]) * _dy2 - (u[i, j+1, k] - u[i, j-1, k]) * (-vel[2] * _dy) * ff[i, j, k]
#                                                      +
#                                                      D * (u[i, j, k+1] - 2 * u[i, j, k] + u[i, j, k-1]) * _dz2 - (u[i, j, k+1] - u[i, j, k-1]) * (-vel[3] * _dz) * ff[i, j, k])
        
        
#                 end
#             end
#         end
# end

# # function farfield3D(field, pos)
# #     nx, ny, nz = size(field)
# #     for i in 2:nx-1
# #         for j in 2:ny-1
# #             for k in 2:nz-1
# #                 field[i, j, k] = (norm(pos .- SA[i, j, k]))^-2
# #             end
# #         end
# #     end
# #     field[pos[1], pos[2], pos[3]] = 1
# #     return field
# # end

# function farfield3D!(field, pos)
#     nx, ny, nz = size(field)
#     for i in 2:nx-1
#         for j in 2:ny-1
#             for k in 2:nz-1
#                 field[i, j, k] = (norm(pos .- SA[i, j, k]))^-2
#             end
#         end
#     end
#     field[pos[1], pos[2], pos[3]] = 1
#     # return field
# end

# # function ff3D(pos, i,j,k)
# #     if pos == SA[i,j,k]
# #         return 1
# #     else
# #         return((norm(pos .- SA[i, j, k]))^-2)
# #     end
# # end



# function getForce3D(field, pos::SV3, para)
#     ii, jj, kk = Int.(round.(pos ./ SA[para.dx, para.dy, para.dz])) .+ 1
#     force = SV3(0, 0, 0)
#     # for i in eachindex(coord)
#     #     # force = force + ratio[i] * ∇(coord[i], field, para)
#     #     force = force .+ 0.25 .* ∇(coord[i], field, para)
#     #     # force = force +  ∇(coord[i], field, dx, dy)
#     # end
#     force = ∇((ii, jj, kk), field, para)
#     # return SV(0,0)
#     return force
# end


# function ∇(coord::Tuple{Int64,Int64,Int64}, field, para)
#     x, y, z = coord
#     dx = (field[x+1, y, z] - field[x-1, y, z]) / (2para.dx)
#     dy = (field[x, y+1, z] - field[x, y-1, z]) / (2para.dy)
#     dz = (field[x, y, z+1] - field[x, y, z+1]) / (2para.dz)
#     return SV3(dx, dy, dz)
# end


# # function kernel(pos::SV3, para)
# #     dx = para.dx
# #     dy = para.dy
# #     dz = para.dz
# #     x, y, z = pos
# #     rx, ry, rz = rem(x, dx), rem(y, dy), rem(z, dz)
# #     _rx = dx - rx
# #     _ry = dy - ry
# #     _rz = dz - rz

# #     Ai = SA[(rx*ry)*rz,
# #         _rx*ry*rz,
# #         _rx*_ry*rz,
# #         rx*_ry*rz,
# #         (rx*ry)*_rz,
# #         _rx*ry*_rz,
# #         _rx*_ry*_rz,
# #         rx*_ry*_rz]

# #     A = dx * dy * dz
# #     ratio = @SVector fill(A, 8)
# #     ratio .-= Ai
# #     ratio = ratio ./ (7 * A)

# #     # @show ratio
# #     i, j, k = Int(div(x, dx)) + 1, Int(div(y, dy)) + 1, Int(div(z, dz)) + 1
# #     coord = SA[(i, j, k), (i + 1, j, k), (i + 1, j + 1, k), (i, j + 1, k),
# #         (i, j, k + 1), (i + 1, j, k + 1), (i + 1, j + 1, k + 1), (i, j + 1, k + 1)]

# #     return coord, ratio
# # end



# # function spread!(field, coord, ratio::SVector{8,Float64}, src)
# #     r = ratio .* src
# #     for l in eachindex(r)
# #         if r[l] != 0
# #             i, j, k = coord[l]
# #             field[i, j, k] = r[l]
# #         end
# #     end
# # end









