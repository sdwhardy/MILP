#=
This file contains the functions that calculate
the cost of the OWPP and manipulate cost data.
It is devided into subsections based on type of
equipment being costed.
=#
#######################################################################################################################################
############# Transformers #######################
function cstF_oss_ttl(S,wp)
#get the cost factors
    ks=cstD_cfs()
#create transformer object
    xfm=xfo()
    xfm.results.ttl=Inf
#get all available xformer sizes
    xfos_all=eqpD_xfo_opt()
#selsct combinations to calculate
    xfos_2use=eqpF_xfo_sel(xfos_all,S)
    for x in xfos_2use
#capex
        x.results.cpx=cstF_oss_cpx(x,ks)
#cost of losses
        x.results.tlc=cstF_xfo_tlc(x,S,ks,wp)
#corrective maintenance
        x.results.cm=cstF_eqp_cm(x,ks)
#eens calculation
        x.results.eens=eensF_eqpEENS(x,S,ks,wp)
#totals the cable cost
        x.results.ttl=cstF_Xttl(x.results)
#store lowest cost option
        if x.results.ttl<xfm.results.ttl
            xfm=deepcopy(x)
        end
    end
    return xfm
end
#############################################
function cstF_oss_cpx(x,k)
#OSS platform CAPEX Calculation
    A=(1+k.dc*(x.num-2))
    B=(k.f_ct+k.p_ct)
    cpx=k.FC_ac+A*B*x.num*x.mva
    return cpx
end
#############################################
function cstF_xfo_tlc(xfo,S,ks,wp)
#OSS tlc Calculation
    pf=eqpF_pf()
    tlc=S*pf*(1-xfo.eta)*ks.T_op*ks.E_op*wp.delta
    return tlc
end
#=############# CABLES #######################
Functions related to the cost of cables are in this section
#############################################
the main calculation flow of the cost of an HVAC cable is performed.
l is the length of the cable, S is the power to be transmitted,
kv is the voltage rating of the cable, wp is an object that describes
the wind profile, os is a binary variable that which is true if a cable
connects 2 OSS and false if from OSS to PCC=#
#######################################################################################################################################
#picks points along the range to calculate
function cstF_points(mn,mx)
    points=Array{Float64,1}()
#higher number increases the points chosen but increases computation time (max tried is 6 which is 4x slower than 3)
    acrcy=3
#selets points equally allong the range
    div=(mx-mn)/acrcy
#sets mini range around points selected
    offset=div/acrcy
#sets amount of increase between points
    unit=2*offset/acrcy
#loads all points into an array
    for i=1:(acrcy-1)
        push!(points,(mn+i*div)-offset)
        for j=1:acrcy
            push!(points,((mn+i*div)-offset)+j*unit)
        end
    end
    return points
end
#############################################
function avrg(a)
    ttl=sum(a)
    av=ttl/length(a)
    return av
end
#############################################
#This function eliminates points were no cable is really suitable so cost is unrealistically high
function cstF_rmvCblOutLrs(cbls,pnts)
    lngth=length(cbls)
    cap_rt=[]
#Calculate ratio of cable capacity to owpp capacity
    for i=1:lngth
        push!(cap_rt,(cbls[i].num*cbls[i].mva)/pnts[i])
    end
#Checks if selcted points average ratio is acceptable, gives warning if not
    avg=avrg(cap_rt)
    if avg<0.93 !! avg>0.97
        println("Cable capacity ratio is outside normal limits!")
        println(avg)
    end
#Eliminate outliers and keep only the cables that are "well" sized
    while length(cbls)>trunc(Int, lngth/2)
        mx=findmax(cap_rt)
        deleteat!(pnts,mx[2])
        deleteat!(cbls,mx[2])
        deleteat!(cap_rt,mx[2])
        mn=findmin(cap_rt)
        deleteat!(pnts,mn[2])
        deleteat!(cbls,mn[2])
        deleteat!(cap_rt,mn[2])
    end
#Checks if selcted points average ratio is acceptable, gives warning if not
    avg=avrg(cap_rt)
    if avg<0.95 || avg>0.99
        println("Cable capacity ratio  after removal of outsiders is outside normal limits!")
        println(avg)
    end
end
#############################################
#Calculates a linear regression model for a kv cable length l
function cstF_linearize_cbl(l,S_min,S_max,kv,wp,o2o)
#checks if limits are reasonable for cable type
    S_min,S_max=cstF_chkCblLms(l,S_min,S_max,kv,wp,o2o)
    ttls=Array{Float64,1}()
    cbls=[]
#Selects appropriate points along the range to calculate
    pnts=cstF_points(S_min,S_max)
    for pnt in pnts
        cb=cbl()
#Finds optimal cable for point
        cb=cstF_cbl_ttl(l,pnt,kv,wp,o2o)
        push!(cbls,cb)
    end
#removes cables at extreme under or over loading
    cstF_rmvCblOutLrs(cbls,pnts)
#seperates results into array of total costs only
    for cbl in cbls
        push!(ttls,cbl.results.ttl)
    end
#fits linear model
    alph_beta=reverse([pnts ones(length(pnts))]\ttls)
#Adds limits to array and returns
    push!(alph_beta,S_min)
    push!(alph_beta,S_max)
    return alph_beta
end
#############################################
#Adjusts demanded limits of cable if outside acceptable linear range
#%%%%%%%%%%%%%%%%%% Future Upgrade Investigate further how to ensure linear range is maintained %%%%%%%%%%%%%%%
function cstF_chkCblLms(l,S_min,S_max,kv,wp,o2o)
    cbls_all=[]
    cbls_all=eqpF_cbl_opt(kv,cbls_all,l)
#Checks upper range
    if S_max>eqpD_MAXcbls()*cbls_all[length(cbls_all)][7]
        S_max=eqpD_MAXcbls()*cbls_all[length(cbls_all)][7]
    end
#Checks lower range
    if 1.5*S_min<cbls_all[1][7]
        S_min=cbls_all[1][7]/1.5
    end
    return S_min,S_max
end
##############################################
function cstF_cbl_ttl(l,S,kv,wp,os)
    cbls_all=[]
    cbls_2use=[]
#create 1 object of type cbl_costs
    cb=cbl()
#Initialize to very high total for comparison
    cb.results.ttl=Inf
#get the cost factors
    ks=cstD_cfs()
#returns all base data available for kv cables
    cbls_all=eqpF_cbl_opt(kv,cbls_all,l)
#Selects 1 to 10 of the cables in parallel appropriate for reuired capacity
    cbls_2use=eqpF_cbl_sel(cbls_all,S,l)
    for i in cbls_2use
#capex
        i.results.cbc=cstF_cbl_cpx(i)
#cost of compensastion
        i.results.qc=cstF_ACcbl_q(i,os,ks)
#cost of losses
        i.results.rlc=cstF_ACcbl_rlc(i,S,ks,wp)
#corrective maintenance
        i.results.cm=cstF_eqp_cm(i,ks)
#eens calculation
        i.results.eens=eensF_eqpEENS(i,S,ks,wp)
#totals the cable cost
        i.results.ttl=cstF_CBLttl(i.results)
#store lowest cost option
        if i.results.ttl<cb.results.ttl
            cb=deepcopy(i)
        end
    end
#return optimal cable costs
    return cb
end
#############################################
#CAPEX of cable
function cstF_cbl_cpx(cbl)
    cbc=cbl.length*cbl.num*cbl.cost
    return cbc
end
#############################################
#Ac compensation cost
function cstF_ACcbl_q(cbl,os,ks)
    cstD_QC(os,ks)
#div sets compensation to 50-50 split
    f=eqpD_freq()
    div=0.5
    A=cbl.farrad*cbl.length*cbl.num
    Q=2*pi*f*cbl.volt^2*A
    Q_oss=Q*div
    Q_pcc=Q*(1-div)
    qc=ks.Qc_oss*Q_oss+ks.Qc_pcc*Q_pcc
    return qc
end
#############################################
#Ac cable loss cost
function cstF_ACcbl_rlc(cbl,s,ks,wp)
    eta=eqpD_xEFF()
#eta is the efficiencu of a feeder transformer, s power to transmit
    A=s*eta
#cable current
    B=cbl.volt*sqrt(3)
    I=A/B
#cable resistance
    R=(cbl.length*cbl.ohm)/(cbl.num)
#I^2R losses times cost factores
#delta is related to wind profile, T_op lifetime hours and E_op cost of energy
    rlc=I^2*R*ks.T_op*ks.E_op*wp.delta
    return rlc
end
#############################################
#generic functions used in calculations for more than 1 specified piece of equipment follow
############## EQUIPMENT ####################
#corrective maintenance of equipment
#eq is the equipment, k are the cost factors
function cstF_eqp_cm(eq,k)
    A=(eq.num*eq.mc)
    B=(1/eq.fr)
#changes time base of mttr of equipment: mnth-hrs/yr-hrs
    C=(eq.mttr*30.417*24.0)/8760.0
    cm=k.cf*(A/(B+C))
    return cm
end
#############################################
#sums all cable costs and returns the total
function cstF_CBLttl(res)
    ttl=res.rlc+res.qc+res.cbc+res.cm+res.eens
    return ttl
end
#############################################
#sums all transformer costs and returns the total
function cstF_Xttl(res)
    ttl=res.cpx+res.tlc+res.cm+res.eens
    return ttl
end
#############################################
