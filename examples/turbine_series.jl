module_path = splitdir(@__FILE__)[1]            # Path to this file


# NOTE: In Julia v0.7, JLD throws the error "Cannot `convert` an object of type
# getfield(JLD, Symbol("##JLD.AssociativeWrapper" unless, imported outside of
# WFVisual through the following line. I'm not sure, but I think it's because
# JLD must be imported before GeometricTools to be able to read Grid JLDs.
# Also, notice that `generate_windfarm(...)` takes an abnormal long time to run
# in v0.7.

import WFVisualr                
wfv=WFVisual

# Load GeometricTools: https://github.com/byuflowlab/GeometricTools.jl
import GeometricTools
gt=GeometricTools

D = 200.0
H = 150.0
nblades = 3

blade = "NREL5MW"
hub = "hub"
tower = "tower1"

module_path = splitdir(@__FILE__)[1]
data_path = joinpath(module_path, "../datav07/")

save_path = "turbine/"         # Save path of this example
gt.create_path(save_path, true)

# 199 time steps
N = 200
rotation_angle = range(0,720,length=N)
for i=1:N
    wind_turbine = wfv.generate_windturbine(D/2, H, blade, hub, tower;
                                nblades=nblades, data_path=data_path,
                                save_path=save_path,file_name="windturbine",
                                paraview=false, rot=rotation_angle[i], pitch=0,
                                time_step=i)
end

nothing
