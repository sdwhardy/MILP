include("eqp/eqp_data.jl")
include("cost/cst_data.jl")
include("cost/cst_functions.jl")
include("eqp/eqp_functions.jl")
include("eqp/eqp_structure.jl")
include("cost/cst_structure.jl")
include("wind/wnd_data.jl")
include("wind/wnd_structure.jl")
include("wind/wnd_functions.jl")
include("eens/eens_functions.jl")
function main()
    l=250
    S=2000
    kv=138.0
    wp=wndD_prof()
    oss2oss=false
    cst=results()
    cst=cstF_cbl_ttl(l,S,kv,wp,oss2oss)
    println(cst)
end

main()
