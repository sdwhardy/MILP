function eqpF_cbls_caps(cbls,km)
    for i in cbls
        push!(i,eqpF_km_cap(km,i[1],i[4],i[5]))
    end
    return cbls
end
########################################################
function eqpF_km_cap(l,v,q,a)
    f=50
    mva=(sqrt(3)*v*10^3*a/10^6)^2-((0.5*((v*10^3)^2*2*pi*f*l*q*10^-9))/10^6)^2
    if mva>=0
        mva=sqrt(mva)
    else
        mva=0.0
    end
 return mva
end
########################################################
function eqpF_pushArray(eqp,array)
    #for i=1:length(array)
    for i in array
        push!(eqp,i)
    end
    return eqp
end
########################################################
function eqpF_cbl_opt(kv,cbls,km)
    if kv==138.0
        opt=eqpD_138cbl_opt(cbls,km)
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
function eqpF_cbl_struct(cb,km,num)
    cbl_data=cbl()
    cbl_data.volt=cb[1]
    cbl_data.size=cb[2]
    cbl_data.ohm=cb[3]
    cbl_data.farrad=cb[4]
    cbl_data.amp=cb[5]
    cbl_data.cost=cb[6]
    cbl_data.length=km
    cbl_data.mva=cb[7]
    cbl_data.num=num
    #Set failure data
    eqpD_cbl_fail(cbl_data)
    #scale for return
    cbl_data.ohm=cbl_data.ohm*10^-3
    cbl_data.farrad=cbl_data.farrad*10^-9
    cbl_data.cost=cbl_data.cost*10^-3
    return cbl_data
end
########################################################
function eqpF_cbl_sel(cbls,S,l)
    cbls_2use=[]
    lims=eqpD_cbl_lims()
    #number of cables in parallel
    for i in cbls
        for j=1:10
            if ((j*i[7])>lims[1]*S && (j*i[7])<lims[2]*S)
                push!(cbls_2use,eqpF_cbl_struct(i,l,j))
            end
        end
    end
    return cbls_2use
end
########################################################
