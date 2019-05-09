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
function mn_buildMap(cnts,sol,obm)
    map=lof_layoutOcn(cnts)
    tpp_main2mfile(map,sol,obm)
    ppf_printOcn(map)#print ocean
    return map
end
function mn_slvMilp(map)
    idd,fmap,raw,nt=tnep_milp(map)
    tpp_prnt2Scrn(raw,nt)
    ppf_printOcn(fmap)
    objmn=trunc(Int,ceil(raw["objective"]))
    return idd,objmn
end
function mn_stUpTnep()
    solmin=Array{Int64,1}()
    objmin=10000
    cntrls=control()
    cntrls.xrad=true
    cntrls.neib1=false
    cntrls.neib3=false
    cntrls.xradPcc=false
    cntrls.xradHlf=false
    cntrls.spcfy=false
    cntrls.xXrad=[3,5,7]
    cntrls.xXneib1=[1,7]
    cntrls.xXneib3=[1,7]
    mp=mn_buildMap(cntrls,solmin,objmin)
    ppf_printOcn(mp)#print ocean
    #=solmin=[]
    objmin=0
    solmin,objmin=mn_slvMilp(mp)
    cntrls.neib1=true
    mp=mn_buildMap(cntrls,solmin,objmin)
    solmin=[]
    objmin=0
    solmin,objmin=mn_slvMilp(mp)
    cntrls.xradPcc=true
    mp=mn_buildMap(cntrls,solmin,objmin)
    solmin=[]
    objmin=0
    solmin,objmin=mn_slvMilp(mp)
    cntrls.xXrad=2
    cntrls.xXneib1=2
    cntrls.xXneib3=2
    mp=mn_buildMap(cntrls,solmin,objmin)
    solmin=[]
    objmin=0
    solmin,objmin=mn_slvMilp(mp)
    ppf_printOcn(mp)#print ocean=#
end
mn_stUpTnep()


network_data = PowerModels.parse_file("results/owpp_tnep_map.m")
network_data["ne_branch"]["3935"]["br_status"]
network_data["ne_branch"]["3935"]["br_status"]=1.0
network_data["ne_branch"]["3935"]["br_status"]
fileName="results/owpp_tnep_nw.mat"
mfile = open(fileName,"w")
print(mfile,network_data)
occursin("1422450110", "1341422450110")
occursin("122450110", "1341422450110")













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
    l=45
    S=1517
    kv=220

    #mn=13
    ka=Array{Tuple,1}()
    push!(ka,(2.32,11.08))#[(k,a)]
    wp=wndF_wndPrf(ka,turb())
    #wp=wndF_wndPrf()

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
