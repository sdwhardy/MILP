#=
This file contains the functions that calculate
the cost of the OWPP and manipulate cost data.
It is devided into subsections based on type of
equipment being costed.
=#
#=############# CABLES #######################
Functions related to the cost of cables are in this section
#############################################
the main calculation flow of the cost of an HVAC cable is performed.
l is the length of the cable, S is the power to be transmitted,
kv is the voltage rating of the cable, wp is an object that describes
the wind profile, os is a binary variable that which is true if a cable
connects 2 OSS and false if from OSS to PCC=#
function cstF_cbl_ttl(l,S,kv,wp,os)
    cbls_all=[]
    cbls_2use=[]
#create 2 objects of type results
    cst_tmp=results()
    cst=results()
#Initialize to very high total for comparison
    cst.ttl=Inf
#create an object of type ks
    ks=cstD_cfs()
#returns all base data available for kv cables
    cbls_all=eqpF_cbl_opt(kv,cbls_all,l)
#Selects 1 to 10 of the cables in parallel appropriate for reuired capacity
    cbls_2use=eqpF_cbl_sel(cbls_all,S,l)
    for i in cbls_2use
#capex
        cst_tmp.cbc=cstF_cbl_cpx(i)
#cost of compensastion
        cst_tmp.qc=cstF_ACcbl_q(i,os,ks)
#cost of losses
        cst_tmp.rlc=cstF_ACcbl_rlc(i,S,ks,wp)
#corrective maintenance
        cst_tmp.cm=cstF_eqp_cm(i,ks)
#eens calculation
        cst_tmp.eens=eensF_eqpEENS(i,S,ks,wp)
#totals the cable cost
        cst_tmp.ttl=cstF_ttl(cst_tmp)
#store lowest cost option
        if cst_tmp.ttl<cst.ttl
            cst=deepcopy(cst_tmp)
        end
    end
#return optimal cable costs
    return cst
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
    div=0.5
    A=cbl.farrad*cbl.length*cbl.num
    Q=2*pi*50*cbl.volt^2*A
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
#sums all costs and returns the total
function cstF_ttl(res)
    ttl=res.oppc+res.opc+res.tlc_pcc+res.tlc_oss+res.rlc+res.qc+res.cbc+res.cm+res.eens
    return ttl
end
#############################################
