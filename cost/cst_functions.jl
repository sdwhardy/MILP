#=
This file contains the functions that calculate
the cost of the OWPP and manipulate cost data.
It is devided into subsections based on type of
equipment being costed.
=#
#######################################################################################################################################
################################################# Transformers ########################################################################
#######################################################################################################################################
#OSS/PCC transformer cost
function cstF_xfm_ttl(S,wp,o2o)
    ks=cstD_cfs()#get the cost factors
    xfm=xfo()#create transformer object
    xfm.results.ttl=Inf
    xfos_all=eqpD_xfo_opt()#get all available xformer sizes
    xfos_2use=eqpF_xfo_sel(xfos_all,S)#select combinations to calculate

    for x in xfos_2use
        if o2o==true
            x.results.cpx=cstF_oss_cpx(x,ks)#capex oss
        else
            x.results.cpx=cstF_pcc_cpx(x)#capex pcc
        end

        x.results.tlc=cstF_xfo_tlc(x,S,ks,wp)#cost of losses
        x.results.cm=cstF_eqp_cm(x,ks)#corrective maintenance
        x.results.eens=eensF_eqpEENS(x,S,ks,wp)#eens calculation
        x.results.ttl=cstF_Xttl(x.results)#totals the cable cost
        #store lowest cost option
        if x.results.ttl<xfm.results.ttl
            xfm=deepcopy(x)
        end
    end
    return xfm
end
#############################################
#OSS platform CAPEX Calculation
function cstF_oss_cpx(x,k)
    A=(1+k.dc*(x.num-2))
    B=(k.f_ct+k.p_ct)
    cpx=k.FC_ac+A*B*x.num*x.mva
    return cpx
end
#############################################
#PCC platform CAPEX Calculation
function cstF_pcc_cpx(x)
    cpx=0.02621*(x.mva*x.num)^0.7513
    return cpx
end

#############################################
#xfm Losses Calculation
function cstF_xfo_tlc(xfo,S,ks,wp)
    pf=eqpD_pf()
    tlc=S*pf*(1-xfo.eta)*ks.T_op*ks.E_op*wp.delta
    return tlc
end
#############################################
#sums all transformer costs and returns the total
function cstF_Xttl(res)
    ttl=res.cpx+res.tlc+res.cm+res.eens
    return ttl
end
#############################################
#######################################################################################################################################
########################################## CABLES #####################################################################################
#######################################################################################################################################
#=Functions related to the cost of cables are in this section
the main calculation flow of the cost of an HVAC cable is performed.
l is the length of the cable, S is the power to be transmitted,
kv is the voltage rating of the cable, wp is an object that describes
the wind profile, os is a binary variable that which is true if a cable
connects 2 OSS and false if from OSS to PCC=#
##############################################
#Cost of optimal cable under l,S,kv,wp,os without considering an OSS
function cstF_cbl_ttl(l,S,kv,wp,os)
    cbls_all=[]
    cbls_2use=[]
    cb=cbl()#create 1 object of type cbl_costs
    cb.results.ttl=Inf#Initialize to very high total for comparison
    ks=cstD_cfs()#get the cost factors
    cbls_all=eqpF_cbl_opt(kv,cbls_all,l)#returns all base data available for kv cables
    cbls_2use=eqpF_cbl_sel(cbls_all,S,l)#Selects 1 to ...(adjust in data) of the cables in parallel appropriate for reuired capacity

    for i in cbls_2use
        i.results.cbc=cstF_cbl_cpx(i)#capex
        i.results.qc=cstF_ACcbl_q(i,os,ks)#cost of compensastion
        i.results.rlc=cstF_ACcbl_rlc(i,S,ks,wp)#cost of losses
        i.results.cm=cstF_eqp_cm(i,ks)#corrective maintenance
        i.results.eens=eensF_eqpEENS(i,S,ks,wp)#eens calculation
        i.results.ttl=cstF_CBLttl(i.results)#totals the cable cost
    #store lowest cost option
        if i.results.ttl<cb.results.ttl
            cb=deepcopy(i)
        end
    end
#return optimal cable object
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
#sums all cable costs and returns the total
function cstF_CBLttl(res)
    ttl=res.rlc+res.qc+res.cbc+res.cm+res.eens
    return ttl
end
##########################################################################################
############## WIND FARM #################################################################
##########################################################################################
#Cost of optimal cable under l,S,kv,wp,os AND considering an OSS
function cstF_cblWT_ttl(l,S,kv,wp,os)
    cbls_all=[]
    cbls_2use=[]
    cb=cbl()#create 1 object of type cbl_costs
    xfm=xfo()
    xfm=cstF_xfm_ttl(S,wp,os)#get transformer
    cb.results.ttl=Inf#Initialize cable to very high total for comparison
    ks=cstD_cfs()#get the cost factors
    cbls_all=eqpF_cbl_opt(kv,cbls_all,l)#returns all base data available for kv cables
    cbls_2use=eqpF_cbl_sel(cbls_all,S,l)#Selects 1 to 10 of the cables in parallel appropriate for reuired capacity

# loop through all cables
    for i in cbls_2use
        i.results.cbc=cstF_cbl_cpx(i)#capex
        i.results.qc=cstF_ACcbl_q(i,os,ks)#cost of compensastion
        i.results.rlc=cstF_ACcbl_rlc(i,S,ks,wp)#cost of losses
        i.results.cm=cstF_eqp_cm(i,ks)#corrective maintenance
        i.results.eens=eensF_EENS(xfm,i,S,ks,wp)#eens calculation
        i.results.ttl=cstF_CBLttl(i.results)#totals the cable cost
        #store lowest cost option
        if i.results.ttl<cb.results.ttl
            cb=deepcopy(i)
        end
    end

#Store values as an owpp object and return
    arc=owpp()
    arc.mva=S
    arc.km=l
    arc.wp=wp
    arc.cable=cb
    arc.xfm=xfm
    cstF_owpp_results(arc)#summarizes final values
    return arc#return optimal cable and transformer
end
#############################################
#summarizes owpp resulst concisely
function cstF_owpp_results(owp)
    owp.costs.cpx=owp.xfm.results.cpx+owp.cable.results.cbc+owp.cable.results.qc
    owp.costs.loss=owp.xfm.results.tlc+owp.cable.results.rlc
    owp.costs.opex=owp.cable.results.eens+owp.xfm.results.cm+owp.cable.results.cm#don't add xfm eens already included
    owp.costs.ttl=owp.costs.cpx+owp.costs.opex+owp.costs.loss
    owp.cable.results.eens=owp.cable.results.eens-owp.xfm.results.eens#adjusts cable eens to validate sum
end
#############################################
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
#########################################################################################
################################## Linearization ########################################
#########################################################################################
#picks points along the range to calculate
function cstF_points(mn,mx)
    points=Array{Float64,1}()
    acrcy=3#higher number increases the points chosen but increases computation time (max tried is 6 which is 4x slower than 3)
    div=(mx-mn)/acrcy#selets points equally allong the range
    offset=div/acrcy#sets mini range around points selected

    unit=2*offset/acrcy#sets amount of increase between points
#loads points into an array
    for i=1:(acrcy-1)
        push!(points,(mn+i*div)-offset)
        for j=1:acrcy
            push!(points,((mn+i*div)-offset)+j*unit)
        end
    end

    return points
end
#############################################
#finds the average of array a
function cstF_avrg(a)
    ttl=sum(a)
    av=ttl/length(a)
    return av
end
#############################################
#eliminates points were no cable is really suitable so cost is unrealistically high
function cstF_rmvCblOutLrs(cbls,pnts)
    lngth=length(cbls)
    cap_rt=[]

#Calculate ratio of cable capacity to owpp capacity
    for i=1:lngth
        push!(cap_rt,(cbls[i].num*cbls[i].mva)/pnts[i])
    end

#Checks if selcted points average ratio is acceptable, gives warning if not
    avg=cstF_avrg(cap_rt)
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
    avg=cstF_avrg(cap_rt)
    if avg<0.95 || avg>0.99
        println("Cable capacity ratio  after removal of outsiders is outside normal limits!")
        println(avg)
    end
end
#############################################
#Calculates a linear regression model for a kv cable length l
function cstF_linearize_cbl(l,S_min,S_max,kv,wp,o2o)
    S_min,S_max=cstF_chkCblLms(l,S_min,S_max,kv,wp,o2o)#checks if limits are reasonable for cable type
    ttls=Array{Float64,1}()
    cbls=[]
    pnts=cstF_points(S_min,S_max)#Selects appropriate points along the range to calculate

#loops through all desired points
    for pnt in pnts
        cb=cbl()
        cb=cstF_cbl_ttl(l,pnt,kv,wp,o2o)#Finds optimal cable for point
        push!(cbls,cb)
    end

    cstF_rmvCblOutLrs(cbls,pnts)#removes cables at extreme under or over loading

#seperates results into array of total costs only
    for cbl in cbls
        push!(ttls,cbl.results.ttl)
    end

    alpha_beta=reverse([pnts ones(length(pnts))]\ttls)#fits linear model

#Adds limits to array and return
    push!(alpha_beta,S_min)
    push!(alpha_beta,S_max)
    return alpha_beta
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
############################################################################################






############################################################################################
###################################### Removed Functions ###################################
############################################################################################
#############################################
#cost of single transformer considering cable
#=function cstF_ossWC_ttl(l,S,kv,wp,o2o)
#get the cost factors
    ks=cstD_cfs()
#create transformer object
    xfm=xfo()
    xfm_orig=xfo()
#Finds optimal cable for point
    cb=cbl()
    xfm,cb=cstF_cblWT_ttl(l,S,kv,wp,o2o)
    xfm_orig=deepcopy(xfm)
#get all available xformer sizes
    xfos_all=eqpD_xfo_opt()
#selsct combinations to calculate
    s_max=min(cb.mva*cb.num,S)
    xfos_2use=eqpF_xfo_sel(xfos_all,s_max)
    for x in xfos_2use
#capex
        x.results.cpx=cstF_oss_cpx(x,ks)
#cost of losses
        x.results.tlc=cstF_xfo_tlc(x,s_max,ks,wp)
#corrective maintenance
        x.results.cm=cstF_eqp_cm(x,ks)
#eens calculation
        #x.results.eens=eensF_eqpEENS(x,S,ks,wp)
        x.results.eens=eensF_EENS(x,cb,S,ks,wp)
#totals the cable cost
        x.results.ttl=cstF_Xttl(x.results)
#store lowest cost option
        if x.results.ttl<xfm.results.ttl
            xfm=deepcopy(x)
        end
    end
    if xfm_orig.results.ttl != xfm.results.ttl
        println("XFRM Changed!!!!")
        print("original: ")
        print(xfm_orig.num)
        print(" - ")
        println(xfm_orig.mva)
        print("new: ")
        print(xfm.num)
        print(" - ")
        println(xfm.mva)
        println("XFRM Changed!!!!")
    end
    return xfm,cb
end=#
#################################################
#=cost of single offhore platform without considering any cable
function cstF_pcc_ttl(S,wp)
    ks=cstD_cfs()#get the cost factors
    xfm=xfo()#create transformer object
    xfm.results.ttl=Inf
    xfos_all=eqpD_xfo_opt()#get all available xformer sizes
    xfos_2use=eqpF_xfo_sel(xfos_all,S)#select combinations to calculate

    for x in xfos_2use
        x.results.cpx=cstF_pcc_cpx(x)#capex
        x.results.tlc=cstF_xfo_tlc(x,S,ks,wp)#cost of losses
        x.results.cm=cstF_eqp_cm(x,ks)#corrective maintenance
        x.results.eens=eensF_eqpEENS(x,S,ks,wp)#eens calculation
        x.results.ttl=cstF_Xttl(x.results)#totals the cable cost
        #store lowest cost option
        if x.results.ttl<xfm.results.ttl
            xfm=deepcopy(x)
        end
    end

    return xfm
end
#############################################
#cost of single offhore platform without considering any cable
function cstF_oss_ttl(S,wp)
    ks=cstD_cfs()#get the cost factors
    xfm=xfo()#create transformer object
    xfm.results.ttl=Inf
    xfos_all=eqpD_xfo_opt()#get all available xformer sizes
    xfos_2use=eqpF_xfo_sel(xfos_all,S)#select combinations to calculate

    for x in xfos_2use
        x.results.cpx=cstF_oss_cpx(x,ks)#capex
        x.results.tlc=cstF_xfo_tlc(x,S,ks,wp)#cost of losses
        x.results.cm=cstF_eqp_cm(x,ks)#corrective maintenance
        x.results.eens=eensF_eqpEENS(x,S,ks,wp)#eens calculation
        x.results.ttl=cstF_Xttl(x.results)#totals the cable cost
        #store lowest cost option
        if x.results.ttl<xfm.results.ttl
            xfm=deepcopy(x)
        end
    end

    return xfm
end=#
