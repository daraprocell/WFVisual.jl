module_path = splitdir(@__FILE__)[1]           

# example 3 - static wind farm with n turbines


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
save_path = "example_3/"         # Save path of this example
gt.create_path(save_path, true)


# --------------------- WIND FARM LAYOUT ---------------------------------------
# turbine_x = [  0.0, 1000.0, -1000.0 ]

# turbine_y = [  0.0, 0.0, 0.0]
# turbine_z = [ 0.,  0.,  0.]

# hub_height = [ 100.,  100.,  100.]
# rotor_diameter = [ 100.,  100.,  100.]
# nBlades = [3, 3, 3]

# nturbs = length(turbine_x)
# hub_height = zeros(nturbs) .+ 100.0
# rotor_diameter = zeros(nturbs) .+ 100.0
# nBlades = zeros(Int8,nturbs) .+ 3

# wind_direction = 228.0

# yaw = [ 0.,  0.,  0.] .- wind_direction

turbine_x = [  909.98349606,  1005.0174903 ,   900.40238835,  -479.57607866,
                 937.92551703,   461.20472344,   123.06165965, -1073.3529325 ,
                -698.76200523, -1083.53094471,  1365.9338773 ,   204.10557045,
                 306.99725222,  -585.35767689,    10.9706543 ,   543.54927959,
                 -20.50212678, -1313.07762466, -1327.00852816,   137.06316521,
                1378.99410072,  -645.66250771,  -535.7876181 ,  -100.66538419,
               -1331.29278587]

turbine_y = [  500.59889721,  -974.65186327,    76.51586507,  1315.29789208,
                1039.37304649, -1321.85187113,  -300.57817461,  -765.82650713,
               -1213.14956771,   886.54973742,  -306.96866239,   194.51410928,
                 735.69027357,    81.95252777,  -892.14222889,  -930.41798152,
                 970.9117144 ,   398.75277565,    -3.52236011,  1393.27798468,
                 241.61586052,  -791.92492366,   919.26560593, -1396.37638078,
                -433.19849727]
turbine_z = [ 0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,
                0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.]

# hub_height = [ 100.,  100.,  100.,  100.,  100.,  100.,  100.,  100.,  100.,
#                 100.,  100.,  100.,  100.,  100.,  100.,  100.,  100.,  100.,
#                 100.,  100.,  100.,  100.,  100.,  100.,  100.]

# rotor_diameter = [ 100.,  100.,  100.,  100.,  100.,  100.,  100.,  100.,  100.,
#                 100.,  100.,  100.,  100.,  100.,  100.,  100.,  100.,  100.,
#                 100.,  100.,  100.,  100.,  100.,  100.,  100.]

# nBlades = [3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
#                3, 3]
nturbs = length(turbine_x)
hub_height = zeros(nturbs) .+ 100.0
rotor_diameter = zeros(nturbs) .+ 100.0
nBlades = zeros(Int8,nturbs) .+ 3

wind_direction = 228.0

yaw = [ 0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,
        0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.] .- wind_direction




# --------------------- PERIMETER AND FLUID DOMAIN -----------------------------
NDIVSx = 20              # Cells in the parametric x-direction
NDIVSy = 20            # Cells in the parametric y-direction
NDIVSz = 10              # Cells in the geometric z-direction

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
invVinf_Oaxis = Vinf_Oaxis'             # Inverse transformation
refH = mean(hub_height)                 # (m) reference height for shear layer


function wake(X)

  k = 0.0325 # we can let this vary
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
  windspeed *= (X[3]/refH)^alpha

  # Return velocity vector in global coordinate system
  return Vinf_Oaxis*[ windspeed*(1-loss), 0, 0 ]
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
