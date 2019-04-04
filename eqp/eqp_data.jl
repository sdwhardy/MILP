#=
This file contains the input data for equipment used
=#
########################################################
#sets owpp power factor
function eqpF_pf()
    return 1.0
end
####################transformers########################
function  eqpD_xfo_opt()
    xfos=Array{Float64,1}()
    for i=50:10:1000
        push!(xfos,i)
    end
    return xfos
end
########################################################
#Set maximum of transformers possible in parallel
function eqpD_MAXxfos()
    return 4
end
########################################################
#the efficiency of transformers
function eqpD_xEFF()
    eta=0.994
    return eta
end
########################################################
#failure data for transformers
function eqpD_xfo_fail(x)
    x.fr=0.03#/yr
    x.mttr=6.0#month
    x.mc=2.8#
    return nothing
end
########################################################
########################################################
#set the system AC frequency
function eqpD_freq()
    return 50.0
end
########################################################
#Set maximum of cables possible in parallel
function eqpD_MAXcbls()
    return 12
end

########################################################
#138kV cables
function eqpD_138cbl_opt(cbls,km)
#exchange rate
    p2e=cstD_xchg()
#%kV,cm^2,mohms/km,nF/km,Amps,10^3 euros/km, capacity at km
    a=[132,185,100,165,501,424*p2e]
    b=[132,300,76.1,175,600,504*p2e]
    c=[132,400,60.6,185,677,568*p2e]
    d=[132,500,49.3,192,739,635*p2e]
    e=[132,630,39.5,209,818,685*p2e]
    f=[132,800,32.4,217,888,795*p2e]
    g=[132,1000,27.5,238,949,860*p2e]
    alphbt=[a,b,c,d,e,f,g]
    alphbt=eqpF_cbls_caps(alphbt,km)
    cbls=eqpF_pushArray(cbls,alphbt)
    return cbls
end
########################################################
#220kV cables
function eqpD_220cbl_opt(cbls,km)
#exchange rate
    p2e=cstD_xchg()
#%kV,cm^2,mohms/km,nF/km,Amps,10^3 euros/km, capacity at km
    a=[220,400,60.1,122,665,728*p2e]
    b=[220,500,48.9,136,732,815*p2e]
    c=[220,630,39.1,151,808,850*p2e]
    d=[220,800,31.9,163,879,975*p2e]
    e=[220,1000,27,177,942,1000*p2e]
    alphbt=[a,b,c,d,e]
    alphbt=eqpF_cbls_caps(alphbt,km)
    cbls=eqpF_pushArray(cbls,alphbt)
    return cbls
end
########################################################
#400kV cables
function eqpD_400cbl_opt(cbls,km)
#exchange rate
    p2e=cstD_xchg()
#%kV,cm^2,mohms/km,nF/km,Amps,10^3 euros/km, capacity at km
    a=[400,500,40,117,776,1239*p2e,]
    b=[400,630,36,125,824,1323*p2e]
    c=[400,800,31.4,130,870,1400*p2e]
    d=[400,1000,26.5,140,932,1550*p2e]
    e=[400,1200,22.1,170,986,1700*p2e]
    f=[400,1400,18.9,180,1015,1850*p2e]
    g=[400,1600,16.6,190,1036,2000*p2e]
    h=[400,2000,13.2,200,1078,2150*p2e]
    alphbt=[a,b,c,d,e,f,g,h]
    alphbt=eqpF_cbls_caps(alphbt,km)
    cbls=eqpF_pushArray(cbls,alphbt)
    return cbls

end
########################################################
#150kV hvdc cables
function eqpD_150cbl_opt(cbls)
end
########################################################
#300kV hvdc cables
function eqpD_300cbl_opt(cbls)
end
########################################################
#Sets the limits that cables will be sized as a % of OWPP capacity
function eqpD_eqp_lims()
    return [0.9,1.5]
end
########################################################
#failure data for cables
function eqpD_cbl_fail(cbl)
    cbl.fr=0.04#/yr/100km
    cbl.mttr=2.0#/yr/100km
    cbl.mc=0.56
    return nothing
end
