#Includes all dependancies
#All data structures
include("cost/cst_structure.jl")#cost
include("wind/wnd_structure.jl")#wind
include("eqp/eqp_structure.jl")#must be included after cst and wind strucs
include("layout/lo_structure.jl")#layout

#input data
include("eqp/eqp_data.jl")#equipment
include("cost/cst_data.jl")#cost
include("wind/wnd_data.jl")#wind
include("layout/lo_data.jl")#layout
#functions
include("cost/cst_functions.jl")#cost
include("eqp/eqp_functions.jl")#equipment
include("wind/wnd_functions.jl")#wind profile
include("eens/eens_functions.jl")#EENS calc
include("post_process/pp_functions.jl")#Post processing
include("layout/lo_functions.jl")#layout

using Distributions
using StatsPlots
using SpecialFunctions
using Polynomials
##################################################################
#################### Cost of cable with no transformers #################
function cbl_cost(l,S,kv,wp,o2o)
    cb=cstF_cbl_ttl(l,S,kv,wp,o2o)
    print("Cable: ")
    println(cb.results.ttl)
end

#################### Cost of transformer with no cables #################
function xfmr_cost(S,wp,o2o)
    xfm=cstF_xfm_ttl(S,wp,o2o)
    print("Xfm: ")
    println(xfm.results.ttl)
end

#################### Cost of transformer and cables #################
function xfmr_cbl_cost(l,S,kv,wp,o2o)
    arc=cstF_cblWT_ttl(l,S,kv,wp,o2o)
    eqpD_xfoXR(kv,arc.xfm)
    arc.cable.ohm=((arc.cable.ohm*l)/arc.cable.num)/eqpD_pu()[4]
    arc.cable.xl=((arc.cable.xl*l)/arc.cable.num)/eqpD_pu()[4]
    arc.cable.yc=arc.cable.yc*l*arc.cable.num*eqpD_pu()[4]
    print("Xfm/Cbl: ")
    println(arc.costs.ttl)
    println()
    println(arc.cable)
end
###################### linearization of cost ########################
function lcbl_cost(l,kv,wp,o2o)
    S_min=50
    S_max=4000
    #returns [b,m,Smax,Smin] of mx+b where x is power
    #Smin and Smax are same as initial argument if they were within range else max or min of range is returned
    ab=cstF_linearize_cbl(l,S_min,S_max,kv,wp,o2o)
    print("linear cable model: ")
    println(ab)
end
##################################################################
#A main function is used to encasulate code as outa function is global scope and should be avoided
function main()
    l=1
    S=200
    kv=33

    mn=13
    a=11.08
    k=2.32
    wp=wndF_wndPrf(mn,a,k)
    #wp=wndF_wndPrf()

    #o2o=false#PCC transformer/ compensation 50-50 offshore-onshore
    o2o=true#OSS transformer/ compensation all offshore
    #cbl_cost(l,S,kv,wp,o2o)
    xfmr_cost(S,wp,o2o)

    #for S=500:1:1500
    xfmr_cbl_cost(l,S,kv,wp,o2o)
    #end
    plotly()
    plot(wp.pu,wp.ce)

    #lcbl_cost(l,kv,wp,o2o)#cstF_linearize_cbl
end
########################################################################################
#main()
lof_layoutOcn()
