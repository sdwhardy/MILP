################################################################################
function lod_concessionAreas()
    areas=[38.0,19.0,23.0,14.0,18.0,35.0,12.0,16.0]
    return areas
end
################################################################################
function lod_cncesGps()
    c0=[]
    c1=[]
    c2=[]
    c3=[]
    c4=[]
    c5=[]
    c6=[]
    c7=[]
    c=[]
    #concessions
    #Norther
    push!(c0,3.015833)
    push!(c0,51.52806)
    push!(c,c0)
    #Thornton
    push!(c1,(2.97+2.919972)/2)
    push!(c1,(51.56+51.53997)/2)
    push!(c,c1)
    #Rentel
    push!(c2,2.939972)
    push!(c2,51.59)
    push!(c,c2)
    #Northwind
    push!(c3,2.900972)
    push!(c3,51.61897)
    push!(c,c3)
    #Seastar
    push!(c4,2.859972)
    push!(c4,51.63)
    push!(c,c4)
    #Nobelwind/Belwind
    push!(c5,(2.819972+2.799972)/2)
    push!(c5,(51.664+51.67)/2)
    push!(c,c5)
    #Northwester
    push!(c6,2.757)
    push!(c6,51.68597)
    push!(c,c6)
    #Mermaid
    push!(c7,2.74)
    push!(c7,51.71997)
    push!(c,c7)#==#
    return c
end
################################################################################
function lod_pccGps()
    pcc0=Array{Float64,1}()
    pcc1=Array{Float64,1}()
    pcc=[]
    #PCCs
    base_lg=2.941944
    base_lt=51.24306
    push!(pcc0,base_lg)
    push!(pcc0,base_lt)
    push!(pcc,pcc0)#==#
    push!(pcc1,3.183611)
    push!(pcc1,51.32694)
    push!(pcc,pcc1)
    return pcc
end
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
    buffer=0
    return buffer
end
################################################################################
#set north buffer on domain
function loD_nbuff()
    buffer=0
    return buffer
end
################################################################################
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
