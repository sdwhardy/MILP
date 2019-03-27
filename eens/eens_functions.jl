###############################################################
function eensF_eqpEENS(eqp, S, ks, wp)
#make capacity probability table
    cpt_tbl=eensF_eqpCPT(eqp,S)
#Create eens array
    eens_all=[]
    eens=0.0
#find curtailment ratio
    for i=1:length(cpt_tbl[:,1])
        ratio_curt=cpt_tbl[i,1]/S
#find closest wind data to curtail ratio
        diff=wp.pu.-ratio_curt
        i_min=argmin(sqrt.((diff[:]).^2))
        if i_min == length(diff) && diff[i_min]<0
            ce=0
        elseif i_min == 1 && diff[i_min]>0
            ce=wp.ce[1]
        elseif i_min < length(diff) && diff[i_min]<0
            ce=eensF_intPole(ratio_curt,wp.pu[i_min],wp.pu[i_min+1],wp.ce[i_min],wp.ce[i_min+1])
        elseif i_min > 1 && diff[i_min]>0
            ce=eensF_intPole(ratio_curt,wp.pu[i_min-1],wp.pu[i_min],wp.ce[i_min-1],wp.ce[i_min])
        else
            ce=wp.ce[i_min]
        end
        push!(eens_all, ce*S*0.95*cpt_tbl[i,2])

    end
    eens=sum(eens_all)*ks.life*ks.E_op
    return eens
end
###############################################################
function eensF_intPole(true_x,min_x,max_x,min_y,max_y)
    slope=(max_y-min_y)/(max_x-min_x)
    b=min_y-slope*min_x
    true_y=slope*true_x+b
    return true_y
end
###############################################################
function eensF_blank_TBL(rows,clms)
    XFM_CBL=zeros(rows,clms)
#create all combinations
    round=1;
    k=1;
    multi=1;
    while round<=clms
      while k<rows
        while k<=(multi*2^(clms-round))
          XFM_CBL[k,round]=1;
          k=k+1;
        end
        multi=multi+2;
        k=k+2^(clms-round);
      end
      round=round+1;
      k=1;
      multi=1;
    end
    return XFM_CBL
end
###############################################################
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
#Set powers and availabilities
    for k=1:clms
        for j=1:rows
            if trunc(Int,empty_tbl[j,k])==1
                AVL_tbl[j]=AVL_tbl[j]*A_eqp
                PWR_tbl[j]=min(S,PWR_tbl[j]+eqp.mva)
            else
                AVL_tbl[j]=AVL_tbl[j]*(1-A_eqp)
            end
          end
      end
tbl_c1=unique(PWR_tbl)
tbl_c2=zeros(Float64,length(tbl_c1),1)
for k=1:length(tbl_c1)
  for j=1:length(AVL_tbl)
      if PWR_tbl[j]==tbl_c1[k]
          tbl_c2[k]=tbl_c2[k]+AVL_tbl[j]
      end
  end
end
    if sum(tbl_c2) > 1.00001 || sum(tbl_c2) < 0.99999
        error("probability does not sum to 1")
    elseif maximum(tbl_c1) > S && minimum(tbl_c1) > 1
        error("power is not correct")
    else
        #display([tbl_c1 tbl_c2])
      return [tbl_c1 tbl_c2]
    end
end
##########################################################
