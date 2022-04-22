module_path = splitdir(@__FILE__)[1]            # Path to this file

# NOTE: In Julia v0.7, JLD throws the error "Cannot `convert` an object of type
# getfield(JLD, Symbol("##JLD.AssociativeWrapper" unless, imported outside of
# WFVisual through the following line. I'm not sure, but I think it's because
# JLD must be imported before GeometricTools to be able to read Grid JLDs.
# Also, notice that `generate_windfarm(...)` takes an abnormal long time to run
# in v0.7.

using JLD
using LinearAlgebra
# Load WFVisual module
import WFVisual
using Statistics
import floris
wfv=WFVisual


# Load GeometricTools: https://github.com/byuflowlab/GeometricTools.jl
import GeometricTools
gt=GeometricTools

# Data path with geometry JLDs
#if Int(VERSION.major)==0 && Int(VERSION.minor)==6   # Case of Julia v0.6
 #   data_path = joinpath(module_path, "../data/")
#else                                                # Case of Julia v0.7
data_path = joinpath(module_path, "../datav07/")
#end

# Create save path where to store vtks
save_path = "FLORIS/"         # Save path of this example
gt.create_path(save_path, true)



# --------------------- WIND FARM LAYOUT ---------------------------------------
turbine_x = [ 405.20435818,  633.89371231,  928.20332558, 1252.10240706,
1884.38199485, 1204.15409295, 1777.25274162,  228.19760716,
1072.91841174, 1419.05981665] .- 1000

turbine_y = [4.00826899e+02, 1.65573597e+03, 1.17733070e+03, 1.91209078e+03,
1.13953462e+03, 1.88146156e+00, 2.87335724e+02, 1.32424102e+03,
5.49657496e+02, 1.41876751e+03] .- 1000

turbine_z = [ 0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0., 0. ]

hub_height = [90.0, 90.0, 90.0, 90.0, 90.0, 90.0, 90.0, 90.0, 90.0, 90.0]

rotor_diameter = [126.0, 126.0, 126.0, 126.0, 126.0, 126.0, 126.0, 126.0, 126.0, 126.0]

nBlades = [3, 3, 3, 3, 3, 3, 3, 3, 3, 3]

wind_direction = 270.0

yaw = [ 0.,  0.,  0., 0., 0., 0., 0., 0., 0., 0. ] .- wind_direction #direction against wind


# --------------------- PERIMETER AND FLUID DOMAIN -----------------------------
NDIVSx = 100              # Cells in the parametric x-direction
NDIVSy = 100            # Cells in the parametric y-direction
NDIVSz = 25               # Cells in the geometric z-direction

# Dummy perimeter
Rper = 1000.0*sqrt(2)
println("creating perimeter points")
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
invVinf_Oaxis = Vinf_Oaxis'             # Inverse transformation
refH = mean(hub_height)                 # (m) reference height for shear layer

println("generating wake")
function wake(X)

  k = 0.0325
  CT = 8.0/9.0
  loss = 0
  alpha = 0.12

  # Calculate loss in freestream direction
  for i = 1:length(turbine_x)     # Iterate over turbines

      # Distance from this turbine in global coordinates
      DX = X - [turbine_x[i], turbine_y[i], turbine_z[i]+hub_height[i]]

      # Distance in free-stream coordinates
      dX = invVinf_Oaxis*DX
      dx = dX[1]

      if dx > 0.
          dy = dX[2]
          dz = dX[3]
          sigma = k*dx + rotor_diameter[i]/sqrt(8.0)
          this_loss = (1.0-sqrt(1.0-CT/(8*sigma^2/(rotor_diameter[i]^2))))*
                          exp(-1.0/2.0*(dy/sigma)^2)*exp(-1.0/2.0*(dz/sigma)^2)
          loss += this_loss^2
      end
  end
  loss = sqrt(loss)

  windspeed= magVinf
  # Uncomment this line adds the shear layer
  # windspeed *= (X[3]/refH)^alpha

  # Return velocity vector in global coordinate system
  return Vinf_Oaxis*[ windspeed*(1-loss), 0, 0 ]
end


gt.calculate_field(fdom, wake, "wake", "vector", "node")

# --------------------- GENERATE WIND FARM -------------------------------------
println("generating geometry")
wfv.generate_windfarm(rotor_diameter, hub_height, nBlades,
                              turbine_x, turbine_y, turbine_z,
                              yaw,
                              perimeter_points, fdom;
                              NDIVSx=NDIVSx, NDIVSy=NDIVSy, NDIVSz=NDIVSz,
                              save_path=save_path, spl_s=0.01,
                              data_path=data_path, paraview=false);

nothing
