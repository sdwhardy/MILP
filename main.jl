#Includes all dependancies
#All data structures
include("cost/cst_structure.jl")
include("wind/wnd_structure.jl")
include("eqp/eqp_structure.jl")#must be included after cst and wind strucs

#input data
include("eqp/eqp_data.jl")
include("cost/cst_data.jl")
include("wind/wnd_data.jl")
#functions
include("cost/cst_functions.jl")#cost
include("eqp/eqp_functions.jl")#equipment
include("wind/wnd_functions.jl")#wind profile
include("eens/eens_functions.jl")#EENS calc
include("post_process/pp_functions.jl")#Post processing

function cable_cost(l,S,kv,wp,o2o)
    cb=cstF_cbl_ttl(l,S,kv,wp,o2o)
    print("Cable: ")
    println(cb)
end
function xfmr_cost(S,wp)
    xfm=cstF_oss_ttl(S,wp)
    print("Xfm: ")
    println(xfm)
end
function xfmr_cbl_cost(l,S,kv,wp,o2o)
    arc=cstF_cblWT_ttl(l,S,kv,wp,o2o)
    print("arc cost: ")
    println(arc)
    #println(S)
end
#linearization
function lcbl_cost(l,kv,wp,o2o)
    S_min=50
    S_max=4000
    #returns [b,m,Smax,Smin] of mx+b where x is power
    #Smin and Smax are same as initial argument if they were within range else max or min of range is returned
    ab=cstF_linearize_cbl(l,S_min,S_max,kv,wp,o2o)
    print("linear cable model: ")
    println(ab)
end
#A main function is used to encasulate code as outa function is global scope and should be avoided
function main()
    l=200
    S=500
    kv=138.0
    wp=wndD_prof()
    o2o=false
    #cable_cost(l,S,kv,wp,o2o)
    #xfmr_cost(S,wp)
    #for S=500:1:1500
        xfmr_cbl_cost(l,S,kv,wp,o2o)
    #end
#cstF_linearize_cbl
    #lcbl_cost(l,kv,wp,o2o)
end
main()
