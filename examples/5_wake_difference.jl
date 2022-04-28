module_path = splitdir(@__FILE__)[1]           

# example 5 - difference between two flow fields

using JLD
using LinearAlgebra
# Load WFVisual module
import WFVisual
using Statistics
wfv=WFVisual

# Load GeometricTools: https://github.com/byuflowlab/GeometricTools.jl
import GeometricTools
gt=GeometricTools
                                              # Case of Julia v0.7
data_path = joinpath(module_path, "../datav07/")

# Create save path where to store vtks
save_path = "example_5/"         # Save path of this example
gt.create_path(save_path, true)


# --------------------- WIND FARM LAYOUT ---------------------------------------
turbine_x = [  0.0, 1000.0, -1000.0 ]

turbine_y = [  0.0, 0.0, 0.0]
turbine_z = [ 0.,  0.,  0.]

hub_height = [ 100.,  100.,  100.]
rotor_diameter = [ 100.,  100.,  100.]
nBlades = [3, 3, 3]

nturbs = length(turbine_x)
hub_height = zeros(nturbs) .+ 100.0
rotor_diameter = zeros(nturbs) .+ 100.0
nBlades = zeros(Int8,nturbs) .+ 3

wind_direction = 270.0

yaw = [ 0.,  0.,  0.] .- wind_direction


# --------------------- PERIMETER AND FLUID DOMAIN -----------------------------
NDIVSx = 200              # Cells in the parametric x-direction
NDIVSy = 200            # Cells in the parametric y-direction
NDIVSz = 50              # Cells in the geometric z-direction

# Dummy perimeter
Rper = 1500.0
perimeter_points = Rper.*[ [cos(a), sin(a), 0] for a in range(0, 2*pi, 179)]


# --------------------- GENERATE FLOW FIELD ------------------------------------

_zmin = 0
_zmax = maximum(hub_height) + 1.25*maximum(rotor_diameter)/2
z_off = 0.0
spl_s=0.01
fdom = wfv.generate_perimetergrid(perimeter_points,
                                  NDIVSx, NDIVSy, NDIVSz;
                                  z_min=_zmin+z_off, z_max=_zmax+z_off,
                                  verify_spline=false,
                                  spl_s=spl_s,
                                  save_path=nothing
                                )
# Guassian wake function
magVinf = 12.0                          # (m/s) wind speed
Vinf_Oaxis = gt.rotation_matrix2(0, 0, wind_direction) # Free-stream coordinate system
println(Vinf_Oaxis)
invVinf_Oaxis = Vinf_Oaxis'             # Inverse transformation
refH = mean(hub_height)                 # (m) reference height for shear layer


function wake(X)

  k = 0.0325 # we can let this vary
  CT = 8.0/9.0
  loss1 = 0
  alpha = 0.12

  loss2 = 0
  alpha_jensen = 0.1

  # Calculate loss in freestream direction
  for i = 1:length(turbine_x)     # Iterate over turbines

      # Distance from this turbine in global coordinates
      DX = X - [turbine_x[i], turbine_y[i], turbine_z[i]+hub_height[i]]

      # Distance in free-stream coordinates
      dX = invVinf_Oaxis*DX
      dx = dX[1]
      dy = dX[2]
      dz = dX[3]
      if dx > 0.
          sigma = k*dx + rotor_diameter[i]/sqrt(8.0)
          this_loss = (1.0-sqrt(1.0-CT/(8*sigma^2/(rotor_diameter[i]^2))))*
                          exp(-1.0/2.0*(dy/sigma)^2)*exp(-1.0/2.0*(dz/sigma)^2)
          loss1 += this_loss^2
      end

      r = rotor_diameter[i]/2.0 + alpha_jensen*dx
      if dx > 0. && abs(dy) <= r && abs(dz) <= r
        this_loss = 2.0/3.0 * (rotor_diameter[i]/2.0 / r)^2
        loss2 += this_loss^2
      end
  end
  loss1 = sqrt(loss1)
  loss2 = sqrt(loss2)

  windspeed= magVinf
  # Uncomment this line adds the shear layer
  # windspeed *= (X[3]/refH)^alpha

  # Return velocity vector in global coordinate system
  # return Vinf_Oaxis*[ windspeed*(1-loss), 0, 0 ]
  speed1 = windspeed*(1-loss1)
  speed2 = windspeed*(1-loss2)

  return Vinf_Oaxis*[ speed1 - speed2, 0, 0 ]
  # return Vinf_Oaxis*[ speed2, 0, 0 ]
end

gt.calculate_field(fdom, wake, "wake", "vector", "node")

# --------------------- GENERATE WIND FARM -------------------------------------
wfv.generate_windfarm(rotor_diameter, hub_height, nBlades,
                              turbine_x, turbine_y, turbine_z,
                              yaw,
                              perimeter_points, fdom;
                              NDIVSx=NDIVSx, NDIVSy=NDIVSy, NDIVSz=NDIVSz,
                              save_path=save_path, spl_s=0.01,
                              data_path=data_path, paraview=false);
nothing
