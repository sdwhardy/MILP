#######################################################
#returns default profile if no arguments are given
function wndF_wndPrf()
    return wndD_prof()
end
#################################################################################
#
function wndF_wndPrf(mn, a, k)
    wnd=wind()#creates wind object
    trb=turb()#creates turbine object
    #cut_off=trb.cout#change cuttoff speed here
    aque=100
    dis=MixtureModel(Weibull,[(k,a)])#Generates a weibull distribution of shape factors a,k
    spds=[0:1/aque:trb.cout]#builds array of wind speeds within turbine range
    probs=sum(pdf.(dis, spds))/aque#maps the speeds to the wiebull distribution
    wndD_TrqCrv(trb)#calculates the output power of the turbine based on wind distribution
    pwr_prb=zeros(length(spds[]),2)
    #pud=findmax(wnd.pwr)[1]#finds the turbines max output level
#loops through wind speeds and matches probability to power level
    for i=1:length(spds[])
        for j = 1:length(trb.wsp)
            if spds[1][i]<=trb.wsp[1] || spds[1][i]>trb.wsp[length(trb.wsp)]
                pwr_prb[i,1]=0.0
                pwr_prb[i,2]=probs[i]
                break
            elseif spds[1][i]==trb.wsp[length(trb.wsp)]
                pwr_prb[i,1]=trb.pwr[length(trb.wsp)]
                pwr_prb[i,2]=probs[i]
                break
            elseif trb.wsp[j]<=spds[1][i] && spds[][i]<=trb.wsp[j+1]
                pwr_prb[i,1]=eensF_intPole(spds[1][i],trb.wsp[j],trb.wsp[j+1],trb.pwr[j],trb.pwr[j+1])
                pwr_prb[i,2]=probs[i]
                break
            end
        end
    end
    pud=findmax(pwr_prb[:,1])[1]#finds the turbines max output level
    pwr_prb=wndD_unic(pwr_prb,pud)#reduces array to unic matches
    graph=wndD_puVShrs(pwr_prb)#puts results in terms of PU and and hrs per year
    wndD_conENG(graph,wnd,pud)#constraint energy calc
    wndD_LDlss(pwr_prb, wnd)#calculates loss factor
    return wnd
end
####################################################################
#eliminates duplicate information from wind profile array
function wndD_unic(pwr_prb,max)
    tbl_c1=[0:0.01:1]#creates array of 0 to 1
    tbl_c2=zeros(Float64,length(tbl_c1[]),1)#creates array of zeroes to match
    tbl_c3=pwr_prb[:,1]./max#Puts power in PU
#calculates energy above capacity that is curtailed
    for i=1:length(tbl_c1[])
        for j=1:length(tbl_c3)
            if (tbl_c3[j]==0.0)
                 tbl_c2[i]=tbl_c2[i]+pwr_prb[j,2]
            elseif (i != length(tbl_c1[]) && tbl_c3[j]>=tbl_c1[][i] && tbl_c3[j]<=tbl_c1[][i+1])
                    tbl_c2[i]=tbl_c2[i]+pwr_prb[j,2]
            else
            end
        end
     end
  return [tbl_c1[] tbl_c2]
end
################################################################################
#calculates load loss factor
function wndD_LDlss(div, wind)
  half_hours=div[:,2]/sum(div[:,2])*(365*48)#half hourly probabilities
  pow_hrs=[div[:,1] half_hours]
  pow_hrs3=pow_hrs[:,1].*pow_hrs[:,2]
  pwr_hrs=[pow_hrs[:,1] pow_hrs[:,2] pow_hrs3]#power-probability-half hourly array
  peak=findmax(pwr_hrs[:,1])[1]#fid the value of max power
  llf=0.0
  lf=0.0
 #load and loss load factor calculation
  for k=1:length(pwr_hrs[:,2])
    for j=1:pwr_hrs[k,2]
      llf=llf+(pwr_hrs[k,1]^2)
      lf=lf+(pwr_hrs[k,1])
    end
  end
  wind.delta=((llf)/(365*48))#saves load loss factor
  wind.lf=((lf)/(365*48))#saves loss factor
end
################################################################################
#Caculates the PU out hours/yr curve
function wndD_puVShrs(div)
  hours=div[:,2]/sum(div[:,2])*(365*24)#hourly probabilities
  pow_hrs=[div[:,1] hours]#power at set hours
  puVShrs=zeros(length(div[:,2]),2)

#create PU vs Hrs array
  for i=1:length(pow_hrs[:,2])
      puVShrs[i,2]=div[i,1]
      for j=1:length(pow_hrs[:,2])
          if pow_hrs[j,1]>=pow_hrs[i,1]
              puVShrs[i,1]=puVShrs[i,1]+pow_hrs[j,2]
          end
      end
  end
  return puVShrs
end
#################################################################################
#calc for constrained energy
function wndD_conENG(graph,wind,max)
#create sized arrays
    B=zeros(length(graph[:,2]),2)
    conENG=zeros(length(graph[:,2]),2)
    p_div=polyfit(graph[:,1],graph[:,2],3)#make a polynomial approximation
    area=zeros(length(graph[:,2]),1)
    integral=polyint(p_div)#set up the integral
    area=(polyval(integral,graph[:,1])-(graph[:,1].*graph[:,2]))#take integral to find area under curve
    x_axis=reverse(graph[:,2],2)
    y_axis=reverse(area[:,1],2)#reverse x and y axis
    B=[x_axis y_axis]
    conENG=sortslices(B,dims=1)#sort by x axis
    wind.ce=conENG[:,2]
    wind.pu=conENG[:,1]#store pu and constraind energy in wind object
    return nothing
end
#################################################################################
