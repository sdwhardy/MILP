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
#A main function is used to encasulate code as outa function is global scope and should be avoided
function main()
    l=250
    S=200
    kv=138.0
    wp=wndD_prof()
    oss2oss=false
    cb=cbl()
    cb=cstF_cbl_ttl(l,S,kv,wp,oss2oss)
    println(cb)
end

main()
