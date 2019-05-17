using Distributions
using StatsPlots
using SpecialFunctions
using Polynomials
#MILP
using PowerModels
#using Mosek
using JuMP, Gurobi
#using Cbc
#using Ipopt
#using Juniper
#using Pavito
#using MAT

#Includes all dependancies
#All data structures
include("cost/cst_structure.jl")#cost
include("wind/wnd_structure.jl")#wind
include("eqp/eqp_structure.jl")#must be included after cst and wind strucs
include("layout/lo_structure.jl")#layout
include("TNEP/tpp_structure.jl")#layout

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
include("TNEP/tnep_milp.jl")#milp
include("TNEP/tnep_post_process.jl")#Post processing milp


#main()
########################################################################################
############################# MILP #####################################################
########################################################################################
function mn_slvMilp(map)
    idd,fmap,raw,nt=tnep_milp(map)
    mystring=tpp_prnt2Scrn(raw,nt)
    objmn=trunc(Int,ceil(raw["objective"]))
    return idd,objmn,mystring,fmap
end
function mn_buildMap(cnts)
    map=lof_layoutOcn(cnts)
    return map
end
function mn_setUpTnep()
    solmin=Array{Int64,1}()
    objmin=10000
    cntrls=control()
    cntrls.xrad=true
    cntrls.neib1=true
    cntrls.neib3=true
    cntrls.xradPcc=false
    cntrls.xradHlf=false
    cntrls.spcfy=false
    cntrls.xXrad=[1,2,3,4,5,6,7]
    cntrls.xXneib1=[1,2,3,4,5,6,7]
    cntrls.xXneib3=[1,2,3,4,5,6,7]
    mp=mn_buildMap(cntrls)
    tpp_main2mfile(mp,solmin,objmin)
    return mp.asBuilt
end



function mn_tnep_full()
    ab=mn_setUpTnep()
    mv("results/owpp_tnep_map.mat", "results/owpp_tnep_map.m", force=true)
    network_data_all = PowerModels.parse_file("results/owpp_tnep_map.m")
    #network_data_all["bus"]
    #network_data_all["ne_branch"]

    #length(network_data_all["bus"])
    oplssngs=Array{oplossing,1}()
    cmbs=[[1],[2,3],[4,5,6,7]]
        for vec in cmbs
            println(vec)
        #x=2
            #solve initialization problem
            solmin=Array{Int64,1}()
            objmin=10000
            cntrls=control()
            cntrls.xrad=true
            cntrls.neib1=false
            cntrls.neib3=false
            cntrls.xradPcc=false
            cntrls.xradHlf=false
            cntrls.spcfy=false
            cntrls.xXrad=deepcopy(vec)
            cntrls.xXneib1=deepcopy(vec)
            cntrls.xXneib3=deepcopy(vec)
            oplssnga=oplossing()
            map=mn_buildMap(cntrls)
            map.asBuilt=ab
            tpp_main2PartMfile(map,network_data_all,solmin,objmin)
            oplssnga.eyeDs,oplssnga.objBst,oplssnga.sumry,oplssnga.layout=mn_slvMilp(map)

            #solve initialization problem
            cntrls=control()
            cntrls.xrad=true
            cntrls.neib1=true
            cntrls.neib3=true
            cntrls.xradPcc=false
            cntrls.xradHlf=false
            cntrls.spcfy=false
            cntrls.xXrad=deepcopy(vec)
            cntrls.xXneib1=deepcopy(vec)
            cntrls.xXneib3=deepcopy(vec)
            println(vec)
            oplssng=oplossing()
            mp=mn_buildMap(cntrls)
            mp.asBuilt=ab
            tpp_main2PartMfile(mp,network_data_all,oplssnga.eyeDs,oplssnga.objBst)
            oplssng.eyeDs,oplssng.objBst,oplssng.sumry,oplssng.layout=mn_slvMilp(mp)
            push!(oplssngs,oplssng)
        end
        return oplssngs
        #return 1
end
#ab=mn_setUpTnep()
ops=mn_tnep_full()
println(ops[1].layout.pccs)
ppf_printOcn(ops[1].layout)#print ocean












########################################################################################
############################ Cost Functions ############################################
########################################################################################
#################### Cost of cable with no transformers #################
function cbl_cost(l,S,kv,wp,o2o)
    arc=cstF_cbl_ttl(l,S,kv,wp,o2o)
    arc.ohm=((arc.ohm*l)/arc.num)/eqpD_pu()[4]
    arc.xl=((arc.xl*l)/arc.num)/eqpD_pu()[4]
    arc.yc=arc.yc*l*arc.num*eqpD_pu()[4]
    print("Cable: ")
    println(arc)
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
    println(arc)
    println()
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
    l=100
    S=1000
    kv=132

    mn=13
    a=11.08
    k=2.32
    #wp=wndF_wndPrf(mn,a,k)
    wp=wndF_wndPrf()

    #o2o=false#PCC transformer/ compensation 50-50 offshore-onshore
    o2o=false#OSS transformer/ compensation all offshore
    #cbl_cost(l,S,kv,wp,o2o)
    #xfmr_cost(S,wp,o2o)

    #for S=500:1:1500
    xfmr_cbl_cost(l,S,kv,wp,o2o)
    #end
    #plotly()
    #plot(wp.pu,wp.ce)

    #lcbl_cost(l,kv,wp,o2o)#cstF_linearize_cbl
end
########################################################################################
########################################################################################
########################################################################################
