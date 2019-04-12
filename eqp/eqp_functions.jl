#=
This file contains functions that manipulate the equipment data
=#

########################################################
function eqpF_puChgBs(Sn,z)
    z=z*100*10^6/Sn
    return z
end
########################################################
#loads values into the end of an array.
function eqpF_cbls_caps(cbls,km)
    for i in cbls
        push!(i,eqpF_km_cap(km,i[1],i[4],i[5]))
    end
    return cbls
end
########################################################
#Calculates the new hvac cable capacity after 50-50 compensation at distance km.
function eqpF_km_cap(l,v,q,a)
#get system frequency
    f=eqpD_freq()
#Calculates the square of new hvac cable capacity after 50-50 compensation at distance km.
    mva=(sqrt(3)*v*10^3*a/10^6)^2-((0.5*((v*10^3)^2*2*pi*f*l*q*10^-9))/10^6)^2
#takes square root if positive returns zero if negative
    if mva>=0
        mva=sqrt(mva)
    else
        mva=0.0
    end
 return mva
end
########################################################
#loads a value into the last position of each array in an araay of arrays
function eqpF_pushArray(eqp,array)
    for i in array
        push!(eqp,i)
    end
    return eqp
end
########################################################
#Logic function that gets the cable data of appropriae voltage level
function eqpF_cbl_opt(kv,cbls,km)
    if kv==33.0
        opt=eqpD_33cbl_opt(cbls,km)
    elseif kv==66.0
        opt=eqpD_66cbl_opt(cbls,km)
    elseif kv==132.0
        opt=eqpD_132cbl_opt(cbls,km)
    elseif kv==220.0
        opt=eqpD_220cbl_opt(cbls,km)
    elseif kv==400.0
        opt=eqpD_400cbl_opt(cbls,km)
    elseif kv==150.0
        opt=eqpD_150cbl_opt(cbls,km)
    elseif kv==300.0
        opt=eqpD_300cbl_opt(cbls,km)
    else
        error("No matching cables for specified KV.")
    end
    return opt
end
########################################################
#Fills in the physical data of a cable into the cable structure
function eqpF_cbl_struct(cb,km,num)
    cbl_data=cbl()
    cbl_data.volt=cb[1]
    cbl_data.size=cb[2]
    cbl_data.ohm=cb[3]*10^-3
    cbl_data.farrad=cb[4]*10^-9
    cbl_data.amp=cb[5]
    cbl_data.cost=cb[6]*10^-3
    cbl_data.length=km
    cbl_data.henry=cb[7]*10^-3
    cbl_data.xl=eqpD_xl(cbl_data.henry)
    cbl_data.yc=eqpD_yc(cbl_data.farrad)
    cbl_data.mva=cb[8]
    cbl_data.num=num
#Set failure data
    eqpD_cbl_fail(cbl_data)
#scale for return
    #cbl_data.ohm=cbl_data.ohm*10^-3
    #cbl_data.farrad=cbl_data.farrad*10^-9
    #cbl_data.cost=cbl_data.cost*10^-3
    return cbl_data
end
########################################################
#return cable inductive reactance
function eqpD_xl(l)
    xl=2*pi*eqpD_freq()*l
    return xl
end
########################################################
#return cable capacitive reactance
function eqpD_yc(c)
    yc=2*pi*eqpD_freq()*c
    return yc
end
########################################################
#Selects sets of cables that satisfy ampacity requirements given by limits
function eqpF_cbl_sel(cbls,S,l)
    cbls_2use=[]
#Get limits and max cables possible in parallel - specified in eqp_data.jl
    lims=eqpD_eqp_lims()
    parCmax=eqpD_MAXcbls(cbls[1][1])
    for i in cbls
        for j=1:parCmax
            if ((j*i[8])>lims[1]*S && (j*i[8])<lims[2]*S)
                push!(cbls_2use,eqpF_cbl_struct(i,l,j))
            end
        end
    end
    return cbls_2use
end
########################################################################################################################################################################
#built chosen sizes into transformer structured array
function eqpF_xfo_struct(s,num)
    xfm=xfo()
    xfm.mva=s
    xfm.num=num
    xfm.eta=eqpD_xEFF()
#Set failure data
    eqpD_xfo_fail(xfm)
    return xfm
end
###########################################################################
#Selects sets of transformers that satisfy power requirements given limits
function eqpF_xfo_sel(xfos,S)
    xfms_2use=Array{xfo,1}()
#Get limits and max cables possible in parallel - specified in eqp_data.jl
    lims=eqpD_eqp_lims()
    parXmax=eqpD_MAXxfos()
    for i in xfos
        for j=1:parXmax
            if ((j*i)>lims[1]*S && (j*i)<lims[2]*S)
                push!(xfms_2use,eqpF_xfo_struct(i,j))
            end
        end
    end
    return xfms_2use
end
###########################################################################


#################################################################################################################
################################# Removed #######################################################################
#################################################################################################################
#cable suceptance
#=function eqpF_cblB(cb)
    print("C: ")
    println(2*pi*50*cb.length*cb.farrad)

    print("X: ")
    println(cb.henry)
    b=(-1*cb.henry)/abs(cb.henry^2+(cb.ohm*cb.length)^2)
    return b
end=#
##################################################################################################################
