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
turbine_x = [ 909.98349606,  1005.0174903 ,   900.40238835,  -479.57607866,
937.92551703,   461.20472344,   123.06165965, -1073.3529325 ,
-698.76200523, -1083.53094471]

turbine_y = [  500.59889721,  -974.65186327,    76.51586507,  1315.29789208,
1039.37304649, -1321.85187113,  -300.57817461,  -765.82650713,
-1213.14956771,   886.54973742]

turbine_z = [ 0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0., 0. ]

hub_height = [ 100.,  100.,  100., 100.,  100.,  100., 100.,  100.,  100., 100. ]

rotor_diameter = [ 100.,  100.,  100., 100.,  100.,  100., 100.,  100.,  100., 100. ]

nBlades = [3, 3, 3, 3, 3, 3, 3, 3, 3, 3]

wind_direction = 270.0

yaw = [ 0.,  0.,  0., 0., 0., 0., 0., 0., 0., 0. ] .- wind_direction #direction against wind. We want it yawing back and forth


# --------------------- PERIMETER AND FLUID DOMAIN -----------------------------
NDIVSx = 50              # Cells in the parametric x-direction
NDIVSy = 50            # Cells in the parametric y-direction
NDIVSz = 50               # Cells in the geometric z-direction

# Dummy perimeter
Rper = 1500.0
perimeter_points = Rper.*[ [cos(a), sin(a), 0] for a in range(0, 2*pi, 179)]

# Dummy wake function
wake(X) = 1.0*[cos(wind_direction*pi/180), sin(wind_direction*pi/180), 0]



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

rotation_angle_start = [0,59,92, 73, 45, 66, 50, 43, 21, 88]
rot_add = 2


# --------------------- GENERATE WIND FARM -------------------------------------
x = range(0,100,length=100)
for i=1:100
  rotation_angle = rotation_angle_start .+ rot_add*i

  turbine_x[1] = x[i]
  wfv.generate_windfarm(rotor_diameter, hub_height, nBlades,
                                turbine_x, turbine_y, turbine_z,
                                yaw,
                                perimeter_points, fdom;
                                NDIVSx=NDIVSx, NDIVSy=NDIVSy, NDIVSz=NDIVSz,
                                save_path=save_path, spl_s=0.01,
                                data_path=data_path, paraview=false,time_step=i,rotation_angle=rotation_angle);
end

nothing
