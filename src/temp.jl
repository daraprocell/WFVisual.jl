



# function generate_windturbine(Rtip::Float64, h::Float64, blade_name::String,
#     hub_name::String, tower_name::String;
#     nblades::Int64=3, data_path::String=def_data_path,
#     save_path="/Users/dprocell/WF/files", file_name="windturbine",
#     paraview=true,
#     rot=nothing, random_rot::Bool=true, pitch=70)

H = 150.0
D = 200.0
nblades = 3

blade = "NREL5MW"
hub = "hub"
tower = "tower1"

module_path = splitdir(@__FILE__)[1]
data_path = joinpath(module_path, "../datav07/")

save_path = "test"
gt.create_path(save_path, true)

generate_windturbine(D/2, H, blade, hub, tower;
                    nblades=nblades, data_path=data_path,
                    save_path=save_path, file_name="windturbine",
                    paraview=false)

#end
