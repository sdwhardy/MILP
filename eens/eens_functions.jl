#=
This file contains all functions used in the calculation of EENS
=#
###############################################################
#EENS of identical equipment in parallel
function eensF_eqpEENS(eqp, S, ks, wp)
#make capacity probability table
    cpt_tbl=eensF_eqpCPT(eqp,S)
#Create eens array
    eens_all=[]
    eens=0.0
#loop through rows of cpt
    for i=1:length(cpt_tbl[:,1])
#find PU curtailment ratio
        ratio_curt=cpt_tbl[i,1]/S
#find closest PU of wind power series to PU curtail ratio
        diff=wp.pu.-ratio_curt
        i_min=argmin(sqrt.((diff[:]).^2))
#check if curt ratio is at or above full power and set ce=curtailed energy to 0
        if i_min == length(diff) && diff[i_min]<0
            ce=0
#check if curt ratio is at or below zero power and set ce to max
        elseif i_min == 1 && diff[i_min]>0
            ce=wp.ce[1]
#if curt ratio is a mid point interpolate ce
        elseif i_min < length(diff) && diff[i_min]<0
            ce=eensF_intPole(ratio_curt,wp.pu[i_min],wp.pu[i_min+1],wp.ce[i_min],wp.ce[i_min+1])
        elseif i_min > 1 && diff[i_min]>0
            ce=eensF_intPole(ratio_curt,wp.pu[i_min-1],wp.pu[i_min],wp.ce[i_min-1],wp.ce[i_min])
#if exact match occurs
        else
            ce=wp.ce[i_min]
        end
#multiply PU curtailed energy with max power and availability, then store
        push!(eens_all, ce*S*cpt_tbl[i,2])
    end
#sum all eens and multiply by cost factors
    eens=sum(eens_all)*ks.life*ks.E_op
    return eens
end
###############################################################
#linearly interpolates 2 points of graph
function eensF_intPole(true_x,min_x,max_x,min_y,max_y)
    slope=(max_y-min_y)/(max_x-min_x)
    b=min_y-slope*min_x
    true_y=slope*true_x+b
    return true_y
end
###############################################################
#creates a blank capacity probability table
function eensF_blank_TBL(rows,clms)
    XFM_CBL=zeros(rows,clms)
#=create all combinations ie
111
110
101
011
100
010
001
000
=#
    round=1
    k=1
    multi=1
    while round<=clms
      while k<rows
        while k<=(multi*2^(clms-round))
          XFM_CBL[k,round]=1
          k=k+1
        end
        multi=multi+2
        k=k+2^(clms-round)
      end
      round=round+1
      k=1
      multi=1
    end
    return XFM_CBL
end
###############################################################
#The calculation of equipment level capacity probability table
function eensF_eqpCPT(eqp,S)
#Calculate failure rate for entire length if cable
    if typeof(eqp)==typeof(cbl())
        eqp.fr=(eqp.fr/100.0)*eqp.length
    end
#Calculate Availability of eqiupment
    A_eqp=1.0/(1.0+eqp.fr*(eqp.mttr*30.0*24.0/8760.0))
#Create combinatorial matrix of 0s and 1s
    clms=trunc(Int,eqp.num)
    rows=trunc(Int, 2.0^clms)
    empty_tbl=eensF_blank_TBL(rows,clms)
#Create blank power and availability tables
    PWR_tbl=zeros(Float64,rows,1)
    AVL_tbl=ones(Float64,rows,1)
#Set powers and availabilities by looping through the CPT
    for k=1:clms
        for j=1:rows
#if 1 the equipment is functional and the power is added to total
#the availability is multiplied
            if trunc(Int,empty_tbl[j,k])==1
                AVL_tbl[j]=AVL_tbl[j]*A_eqp
                PWR_tbl[j]=min(S,PWR_tbl[j]+eqp.mva)
#if 0 the equipment is broken and no power is transmitted
            else
                AVL_tbl[j]=AVL_tbl[j]*(1-A_eqp)
            end
          end
      end
#all unique power levels are extracted
tbl_c1=unique(PWR_tbl)
tbl_c2=zeros(Float64,length(tbl_c1),1)
for k=1:length(tbl_c1)
  for j=1:length(AVL_tbl)
#Availabilities are summed for common power levels
      if PWR_tbl[j]==tbl_c1[k]
          tbl_c2[k]=tbl_c2[k]+AVL_tbl[j]
      end
  end
end
#Checks if probability sums to 1 else erro ris thrown
    if sum(tbl_c2) > 1.00001 || sum(tbl_c2) < 0.99999
        error("probability does not sum to 1")
    elseif maximum(tbl_c1) > S && minimum(tbl_c1) > 1
        error("power is not correct")
    else
      return [tbl_c1 tbl_c2]
    end
end
##########################################################
