###############################################################################
####################### Problem size adjustments ##############################
###############################################################################
########################### control booleans ##################################
###############################################################################
#place oss in radius on route to neighbouring oss
function lod_xrad10()
    return false
end
###############################################################################
#place oss in middle between neighbouring oss
function lod_oss1neib10()
    return true
end
###############################################################################
#place 2nd and 3rd oss between neighbouring oss
function lod_oss3neibs10()
    return false
end
###############################################################################
#place oss in near radius of oss towards pcc
function lod_ossXradPcc10()
    return false
end
###############################################################################
#place oss half way (adjustable with  lod_frcNum()) to pcc from owpp
function lod_ossXradPccHalf10()
    return false
end
###############################################################################
# set if close osss should be combined
function lod_ossSpcfy()
    return false
end
###############################################################################
########################## distance adjustments ###############################
###############################################################################
#fraction of distance to place oss on lod_frcNum() of furthest owpp (mid point compensation)
function lod_rdFrc()
    return 0.5
end
###############################################################################
#furthest # of owpp to add lod_rdFrc() point oss (midpoint compensation)
function lod_frcNum()
    return 4
end
###############################################################################
#distance from owpp to place nearest oss on path to neighbouring owpp
function lod_rad()
    return 1
end
###############################################################################
#the number of neighbours to place closest oss between (those on radius)
function lod_xXneibs()
    return 1
end
###############################################################################
#the number of closest neighbours to place 1 oss at centre point
function lod_x1neibs()
    return 1
end
###############################################################################
#place 2 and 3 in middle of x number of neighbouring owpp
function lod_x3neibs()
    return 1
end
###############################################################################
#set maximum distance to connect the gens to pccs with MV cable
#values are set compared to 220kV grid
function lod_mxMv2PccKm(cn)
    if cn.kv == 33.0
        km=10
    elseif cn.kv == 66.0
        km=25
    else
        error("Cable MV does not match option!")
    end
    return km
end
###############################################################################
#set maximum distance to connect the gens to oss with MV cable
#values are set compared to 220kV grid
function lod_mxMvKm(cn)
    if cn.kv == 33.0
        km=3
    elseif cn.kv == 66.0
        km=7
    else
        error("Cable MV does not match option!")
    end
    return km
end
###############################################################################
#set min oss to oss arc length
function lod_mnKm()
    return 2.0
end
################################################################################
#set minimum distance between any neighbouring OSS (only used if lod_spcfy()=true)
function lod_mnDist()
    return 2
end
################################################################################
#sets max distance upstream to connect gen to oss and which wind profiles are included at OSS
function lod_gen2Noss()
    return 2
end
################################################################################
#sets offset of sourounding OSS from center generator
#=function lod_genSpc()
    return 2.5
end=#
################################################################################
#set west buffer on domain
function loD_wbuff()
    buffer=0
    return buffer
end
################################################################################
#set east buffer on domain
function loD_ebuff()
    buffer=0
    return buffer
end
################################################################################
#set south buffer on domain
function loD_sbuff()
    buffer=2.5
    return buffer
end
################################################################################
#set north buffer on domain
function loD_nbuff()
    buffer=0
    return buffer
end
################################################################################
#sets multiple of turbine diameter for oss spacing
#function lod_ossSpc()
#    return 100
#end
###############################################################################
####################### Voltage,power,wind adjustments ########################
###############################################################################
#set collector voltage
function lod_cncsKv()
    return 66.0
end
################################################################################
#set onshore transmission voltage
function lod_pccKv()
    return 220.0
end
################################################################################
#set oss transmission voltage
function lod_ossKv()
    return 220.0
end
################################################################################
function lod_cnceMva()
    return 250.0
end
################################################################################

###############################################################################
####################### OWPP/pcc location #####################################
###############################################################################
function lod_cncesGps()
    c=Array{Tuple,1}()
    wnd=Array{Tuple,1}()
    p=Array{Float64,1}()
    trbs=Array{turb,1}()
    #concessions
    #Norther
    push!(c,(3.015833,51.52806))
    push!(p,250.0)
    push!(wnd,(2.32,11.08))
    trb=turb()
    wndD_TrqCrv(trb)
    push!(trbs,trb)
    #Thornton
    push!(c,((2.97+2.919972)/2,(51.56+51.53997)/2))
    push!(p,250.0)
    push!(wnd,(2.32,11.08))
    trb=turb()
    wndD_TrqCrv(trb)
    push!(trbs,trb)
    #Rentel
    push!(c,(2.939972,51.59))
    push!(p,250.0)
    push!(wnd,(2.32,11.08))
    trb=turb()
    wndD_TrqCrv(trb)
    push!(trbs,trb)
    #Northwind
    push!(c,(2.900972,51.61897))
    push!(p,250.0)
    push!(wnd,(2.32,11.08))
    trb=turb()
    wndD_TrqCrv(trb)
    push!(trbs,trb)
    #Seastar
    push!(c,(2.859972,51.63))
    push!(p,250.0)
    push!(wnd,(2.32,11.08))
    trb=turb()
    wndD_TrqCrv(trb)
    push!(trbs,trb)
    #Nobelwind/Belwind
    push!(c,((2.819972+2.799972)/2,(51.664+51.67)/2))
    push!(p,250.0)
    push!(wnd,(2.32,11.08))
    trb=turb()
    wndD_TrqCrv(trb)
    push!(trbs,trb)
    #Northwester
    push!(c,(2.757,51.68597))
    push!(p,250.0)
    push!(wnd,(2.32,11.08))
    trb=turb()
    wndD_TrqCrv(trb)
    push!(trbs,trb)
    #Mermaid
    push!(c,(2.74,51.71997))
    push!(p,250.0)
    push!(wnd,(2.32,11.08))
    trb=turb()
    wndD_TrqCrv(trb)
    push!(trbs,trb)
    return c,p,wnd,trbs
end
################################################################################
function lod_pccGps()
    pcc=Array{Tuple,1}()
    #PCCs,
    #push!(pcc,(2.939692,51.239737))
    push!(pcc,(3.183611,51.32694))
    base_lg=pcc[1][1]#2.941944
    base_lt=pcc[1][2]#51.24306
    print("reference longitude: ")
    println(base_lg)
    print("reference latitude: ")
    println(base_lt)
    return pcc
end
###############################################################################
####################### Turbine Data ##########################################
###############################################################################
#Default torque curve of the turbine
#used only in case of specified wind profile other than the default and no specific turbine
function wndD_TrqCrv(trb)
trb.mva=8#turbine power
trb.dia=0.1#km turbine diameter
trb.cin=4#turbine cut in speed
trb.cout=25#turbine cut out
#wind speed
trb.wsp=[2.11
2.57
3.06
3.49
4.01
4.49
5.03
5.52
5.97
6.49
6.99
7.48
7.99
8.5
8.96
9.51
10.04
10.48
10.99
11.48
11.98
12.41
12.96
13.53
14.05
14.6
15
15.5
15.95
16.57
17.03
17.51
17.92
18.47
18.95
19.5
20
20.5
21
21.5
22
22.5
23
23.5
24
24.5
25
25.5
26
26.5
27
27.5
28]
#power production
trb.pwr=[0
0
0
3.6
191.2
421.7
716.5
992.9
1267.1
1606.8
1971
2552.9
2928.5
3772.3
4443.9
5418.1
6041.2
6510.9
7098.8
7523.9
8726.7
8065.9
8253.8
8366.1
8359.3
8425.7
8423.8
8436.4
8434.8
8441.5
8443.1
8343.6
8252.3
8241.6
8144.6
8066.1
7953.8
7865.9
7726.7
7523.9
7298.8
6710.9
6041.2
5418.1
4443.9
3772.3
2928.5
2552.9
1971
1606.8
1267.1
992.9
516.5
]
return nothing
end
#######################################################





################################################################################
######################## Removed Functions #####################################
################################################################################
#=function lod_concessionAreas()
    areas=[38.0,19.0,23.0,14.0,18.0,35.0,12.0,16.0]
    return areas
end=#
################################################################################
#place OSS north of generation?
#=function lod_noss()
    return true
end
###############################################################################
#place OSS south of generation?
function lod_soss()
    return true
end
###############################################################################
#place OSS swest of generation?
function lod_woss()
    return true
end
###############################################################################
#place OSS east of generation?
function lod_eoss()
    return true
end
###############################################################################
#place OSS on generation?
function lod_goss()
    return true
end
###############################################################################
#extend OSS to the east of generation?
function lod_eosss()
    return false
end
###############################################################################
#extend OSS to the west of generation?
function lod_wosss()
    return true
end=#
###############################################################################
#sets the power level of each concession
#=function lod_cncesMva(nc)
    mva=Array{Float64,1}()
    for i=1:nc
        push!(mva,lod_cnceMva())
    end
    return mva
end=#
################################################################################
#sets wind of each concession
#=function lod_cncesWnd(nc)
    wnd=Array{Tuple,1}()
    #include special functions for gamma
    #=mn=13
    a=mn/gamma(1.5)=#
    a=11.08
    k=2.32
    for i=1:nc
        push!(wnd,(k,a))
    end
    return wnd
end
################################################################################
#sets wind of each concession
function lod_cncesTrbs(nc)
    trbs=Array{turb,1}()
    for i=1:nc
        trb=turb()
        wndD_TrqCrv(trb)
        push!(trbs,trb)
    end
    return trbs
end=#
################################################################################
