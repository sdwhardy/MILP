#Includes all dependancies
include("eqp/eqp_data.jl")
include("cost/cst_data.jl")
include("cost/cst_functions.jl")
include("eqp/eqp_functions.jl")
include("cost/cst_structure.jl")
include("eqp/eqp_structure.jl")
include("wind/wnd_data.jl")
include("wind/wnd_structure.jl")
include("wind/wnd_functions.jl")
include("eens/eens_functions.jl")
include("post_process/pp_functions.jl")
#A main function is used to encasulate code as outa function is global scope and should be avoided
function main()
    l=25
    S_min=50
    S_max=4000
    kv=138.0
    wp=wndD_prof()
    o2o=false
    #returns [b,m,Smax,Smin] of mx+b where x is power
    #Smin and Smax are same as initial argument if they were within range else max or min of range is returned
    ab=cstF_linearize_cbl(l,S_min,S_max,kv,wp,o2o)
    println(ab)
end

main()
