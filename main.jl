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
    l=250
    S_min=50
    S_max=2000
    S=S_max
    kv=138.0
    csvfile=ppf_cblFileName(S_min,S_max,l,kv)
    wp=wndD_prof()
    oss2oss=false
    cb=cbl()
    cb=cstF_cbl_ttl(l,S,kv,wp,oss2oss)
    ppf_wrt2CblFile(csvfile,cb,S)
    close(csvfile)
end

main()
