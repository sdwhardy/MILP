#this file contains functions that calculate the layout of the area
###############################################################################
############################ Main #############################################
###############################################################################
function lof_layoutOcn(cnt)
    ocean=eez()#build the eez in the ocean
    lof_shoreConnect(ocean.pccs)#add the gps of points of common coupling
    println("PCCs positioned...")
    ocean.cnces=lof_layoutZone(length(ocean.pccs))#add the region with concessions
    println("OWPPs positioned...")
    base=lof_bseCrd(ocean)#find base coordinates
    lof_dists(ocean,base)#set layout in terms of base
    println("GPS coordinates projected onto cartesian plane...")
    os=lof_rotateOcn(ocean)#apply rotation to align n-s with y
    lof_slideOcn(ocean,os)#slide to align western most point with x=0
    println("Axis transformed...")
    lof_ewbnd(ocean)#find the boundary of region containnig oss
    #=println(ocean.bnd.wbnd.lims)
    println(ocean.bnd.wbnd.lmodel)
    println(ocean.bnd.ebnd.lims)
    println(ocean.bnd.ebnd.lmodel)=#
    println("boundary drawn...")
    lof_osss(ocean,cnt)#add all osss within boundary
    println("OSSs positioned...")
    lof_GoArcs(ocean)#add all gen to oss arcs within boundary
    println("OWPP to OSS arcs complete...")
    lof_GpArcs(ocean)#add all gen to pcc arcs within boundary
    println("OWPP to PCC arcs complete...")
    lof_OpArcs(ocean)#add all oss to pcc arcs within boundary
    println("OSS to PCC arcs complete...")
    lof_OoArcs(ocean)#add all oss to oss arcs within boundary
    println("OSS to OSS arcs complete...")
    #println(length(ocean.gOarcs)+length(ocean.oOarcs)+length(ocean.oParcs))
    #println(length(ocean.cnces)+length(ocean.osss)+length(ocean.pccs))
    #ppf_printOcn(ocean)#print ocean
    return ocean
end
################################################################################

################################################################################
########################### General purpose ####################################
################################################################################
#returns the hypotenuse distance between 2 points
function lof_pnt2pnt_dist(pnt1,pnt2)
    hyp=sqrt((pnt2.x-pnt1.x)^2+(pnt2.y-pnt1.y)^2)
    if hyp < 1
        hyp=1
        println("Length less than 1km set to 1km.")
    end
    return hyp
end
################################################################################
#Change degrees to radians
function lof_d2r(deg)
    return deg*pi/180
end
################################################################################
#calculates length of 1 deg og longitude at given lattitude
function lof_lg1deg(lat)
    return cos(lof_d2r(lat))*111
end
################################################################################
#changes angle to a length
function lof_deg2lgth(d,dPl)
    return d*dPl
end
################################################################################

################################################################################
########################### Region Layout ######################################
################################################################################
#rotates axis to align n-s with y
function lof_rotation(x,y,theta)
    co_od=[x y]
    rotated=co_od*[cos(theta) -1*sin(theta);sin(theta) cos(theta)]
    return rotated
end
################################################################################
#shifts western most to be at x=0
function lof_shift(x,os)
    return x+os
end
################################################################################
#base layout of the zone consisting of all concessions
function lof_layoutZone(num)
    km=1#scale
    gpss,mvas,wnds,trbs=lod_cncesGps()#get gps coords/power/winds/turbines of each concession
    #mvas=lod_cncesMva(length(gpss))#get power of each concession
    #wnds=lod_cncesWnd(length(gpss))#get wind profiles for each concession
    #trbs=lod_cncesTrbs(length(gpss))#get turbine info for each concession
    #areas=lod_concessionAreas()
    zone=Array{cnce,1}()
    #loop through all concessions
    for i=1:length(gpss)
        concession=cnce()
        concession.gps.lng=gpss[i][1]#set longitude
        concession.gps.lat=gpss[i][2]#set latittude
        concession.mva=mvas[i]#set concession power
        concession.wnd=wnds[i]#set wind profiles
        concession.trb=trbs[i]#set turbine type
        concession.kv=lod_cncsKv()#set collector kv
        concession.num=num+i
        #concession.area=areas[i]
        push!(zone,deepcopy(concession))
    end
    return zone
end
################################################################################
#sets the gps coords that are the reference coords
function lof_bseCrd(ocean)
    base=gps()
    base.lat=ocean.pccs[length(ocean.pccs)].gps.lat#base lat
    base.lng=ocean.pccs[length(ocean.pccs)].gps.lng#base long
    #println(base)
    return base
end
################################################################################
#loops through all coords to get relative lengths
function lof_dist(location,base,lnthLT,km)
    for i=1:length(location)
        location[i].coord.x=lof_deg2lgth(location[i].gps.lng-base.lng,lof_lg1deg(location[i].gps.lat+base.lng)*km)
        #location[i].coord.y=abs(lof_deg2lgth(location[i].gps.lat-base.lat,lnthLT))
        location[i].coord.y=lof_deg2lgth(location[i].gps.lat-base.lat,lnthLT)
    end
end
################################################################################
#calculates lengths based on latitude
#as lattitude changes number of km should be updated
function lof_dists(ocean,base)
    km=1
    lnth_1deg_LT=111*km
    lof_dist(ocean.cnces,base,lnth_1deg_LT,km)
    lof_dist(ocean.pccs,base,lnth_1deg_LT,km)
end
################################################################################
#loops through to apply appropriate rotations for coords
function lof_rot(location,theta,os)
    for i=1:length(location)
        xy=lof_rotation(location[i].coord.x,location[i].coord.y,theta)
        location[i].coord.x=xy[1]
        location[i].coord.y=xy[2]
        if location[i].coord.x<os
            os=location[i].coord.x
        end
    end
    return os
end
################################################################################
#rotates the entire region
function lof_rotateOcn(ocean)
    theta=atan((ocean.pccs[length(ocean.pccs)].coord.x-ocean.cnces[length(ocean.cnces)].coord.x)/(ocean.cnces[length(ocean.cnces)].coord.y-ocean.pccs[length(ocean.pccs)].coord.y))
    os=0.0
    os=lof_rot(ocean.cnces,theta,os)
    os=lof_rot(ocean.pccs,theta,os)
    return os
end
################################################################################
#translates the entire region
function lof_slideOcn(ocean,os)
    for cnce in ocean.cnces
        cnce.coord.x=lof_shift(cnce.coord.x,abs(os))
        cnce.id="1"*string(trunc(Int,10*cnce.coord.x))*string(trunc(Int,10*cnce.coord.y))
    end
    for pcc in ocean.pccs
        pcc.coord.x=lof_shift(pcc.coord.x,abs(os))
        pcc.id="0"*string(trunc(Int,10*pcc.coord.x))*string(trunc(Int,10*pcc.coord.y))
    end
end
###############################################################################

###############################################################################
############################ border ###########################################
###############################################################################
#find the boundary on the western side
function lof_wBnd(bd,f,pnts)

    thetas=Array{Float64,1}()
    while (bd[length(bd)].x != f.x || bd[length(bd)].y != f.y)#loop until furthest point is used
        for i in pnts
            ang=atan((i.y-bd[length(bd)].y),(i.x-bd[length(bd)].x))
            #if ang>=0
                push!(thetas,ang)#calc angles between current point and remainder
            #end
        end
        if length(thetas) != 0

            pnt=findmax(thetas)[2]#select th eminimum angle

            push!(bd,pnts[pnt])#chose this point as next boundary point
            if pnts[pnt].x != f.x && pnts[pnt].y != f.y
                deleteat!(pnts,pnt)#remove the used point
            end
            thetas=[]
        end
    end
    return bd,pnts
end
################################################################################
#find the boundary on the eastern side
function lof_eBnd(bd,f,pnts)
    thetas=Array{Float64,1}()
    while (bd[length(bd)].x != f.x || bd[length(bd)].y != f.y)#loop until furthest point is used
        for i in pnts
            ang=atan((i.y-bd[length(bd)].y),(i.x-bd[length(bd)].x))
            if ang>=0
                push!(thetas,ang)#calc angles between current point and remainder
            else
                push!(thetas,2*pi)#ensure negative angles are not selected
            end
        end
        if length(thetas) != 0
            #println(thetas)
            pnt=findmin(thetas)[2]#select the minimum angle
            #println(pnts[pnt])
            push!(bd,pnts[pnt])#chose this point as next boundary point
            if pnts[pnt].x != f.x || pnts[pnt].y != f.y
                deleteat!(pnts,pnt)#remove the used point
            end
            thetas=[]
        end
    end
    return bd
end
###############################################################################
#adds specified buffer to the west
function lof_addWBuff(bnd)
    wbuff=loD_wbuff()
    #sbuff=loD_sbuff()
    nbuff=loD_nbuff()
    for i in bnd
        i.x=i.x-wbuff
    end
    bnd[length(bnd)].y=bnd[length(bnd)].y+nbuff
    #bnd[1].y=bnd[1].y-sbuff
    return bnd
end
################################################################################
#adds specified buffer to the east
function lof_addEBuff(bnd)
    ebuff=loD_ebuff()
    #sbuff=loD_sbuff()
    nbuff=loD_nbuff()
    for i in bnd
        i.x=i.x+ebuff
    end
    bnd[length(bnd)].y=bnd[length(bnd)].y+nbuff
    #bnd[1].y=bnd[1].y-sbuff
    return bnd
end
################################################################################
#makes the linear model given boundary points
function lof_lnrBnd(xbnd)
    for i=1:length(xbnd.lims)-1
        push!(xbnd.lmodel,line())
        if xbnd.lims[i].y == xbnd.lims[i+1].y
            #xbnd.lims[i+1].x=xbnd.lims[i+1].x
            alpha_beta=[xbnd.lims[i].y,0]
        else
            alpha_beta=reverse([[xbnd.lims[i].y,xbnd.lims[i+1].y] ones(2)]\[xbnd.lims[i].x,xbnd.lims[i+1].x])#fits linear model
        end
        xbnd.lmodel[i].b=alpha_beta[1]
        xbnd.lmodel[i].m=alpha_beta[2]
    end
    return nothing
end
################################################################################
#finds the full boundary points around the concession and makes linear boundary model
function lof_ewbnd(ocn)
    #dummy arrays
    all=Array{xy,1}()
    all_oss=Array{xy,1}()
    Wbnd=Array{xy,1}()
    Ebnd=Array{xy,1}()

    strt=ocn.pccs[length(ocn.pccs)].coord#define 0 point
    fnsh=ocn.cnces[length(ocn.cnces)].coord#define furthest concession
    push!(Wbnd,deepcopy(strt))
    push!(Ebnd,deepcopy(strt))
    push!(all,deepcopy(ocn.pccs[1].coord))
    for i in ocn.cnces[1:length( ocn.cnces)]
        push!(all,deepcopy(i.coord))

    end

    #sets east and west boundary

    Wbnd,all=lof_wBnd(Wbnd,fnsh,all)#finds western boundary
    Ebnd=lof_eBnd(Ebnd,fnsh,deepcopy(all))#finds eastern boundary
    Wbnd=lof_addWBuff(Wbnd)#add buffer to western boundary
    Ebnd=lof_addEBuff(Ebnd)#add buffer to eastern boundary
    ocn.bnd.wbnd.lims=Wbnd
    ocn.bnd.ebnd.lims=reverse(Ebnd,1)

    #set calculated n-s boundary points
    push!(ocn.bnd.nbnd.lims,ocn.bnd.wbnd.lims[length(ocn.bnd.wbnd.lims)])
    push!(ocn.bnd.nbnd.lims,ocn.bnd.ebnd.lims[1])
    push!(ocn.bnd.sbnd.lims,ocn.bnd.ebnd.lims[length(ocn.bnd.ebnd.lims)])
    push!(ocn.bnd.sbnd.lims,ocn.bnd.wbnd.lims[1])

    #calculates linear model for boundary
    lof_lnrBnd(ocn.bnd.ebnd)
    lof_lnrBnd(ocn.bnd.wbnd)
    lof_lnrBnd(ocn.bnd.sbnd)
    lof_lnrBnd(ocn.bnd.nbnd)

    return nothing
end
###############################################################################

##############################################################################
############################ laying PCC nodes ################################
##############################################################################
#Places the pccs
function lof_shoreConnect(location)
    gpss=lod_pccGps()
    for i=1:length(gpss)
        shore=pcc()
        shore.gps.lng=gpss[i][1]
        shore.gps.lat=gpss[i][2]
        shore.kv=lod_pccKv()
        shore.num=i
        push!(location,deepcopy(shore))
    end
end
##############################################################################

###############################################################################
########################## Laying OSS nodes ###################################
###############################################################################
#finds the western boundary at y=y
function lof_wLim(y,bd)
    wlx=0
    for i=1:length(bd.lims)-1
        if (bd.lims[i].y <= y && y <= bd.lims[i+1].y)
            wlx=bd.lmodel[i].m*y+bd.lmodel[i].b
        end
    end
    return wlx
end
###############################################################################
#finds the eastern boundary at y=y
function lof_eLim(y,bd)
    elx=0
    for i=1:length(bd.lims)-1
        if (bd.lims[i+1].y <= y && y <= bd.lims[i].y)
            elx=bd.lmodel[i].m*y+bd.lmodel[i].b
        end
    end
    return elx
end
function lof_mrgOss(xy1,xy2)
    xy1.x=(xy1.x+xy2.x)/2
    xy1.y=(xy1.y+xy2.y)/2
end
###############################################################################
function lof_renumOss(nm,j,osss)
    for i in osss[j:length(osss)]
        i.num=deepcopy(nm)
        nm=nm+1
    end
end
###############################################################################
function lof_ossSprcfy(osss)
    mnDist=lod_mnDist()
    for i=1:length(osss)-1
        j=i+1
        while j<=length(osss)
            if lof_pnt2pnt_dist(osss[i].coord,osss[j].coord) <= mnDist
                lof_mrgOss(osss[i].coord,osss[j].coord)
                osss[j].id="2"*string(trunc(Int,10*osss[j].coord.x))*string(trunc(Int,10*osss[j].coord.y))
                #nm=osss[j].num
                deleteat!(osss,j)#remove the used point
                #println(length(osss))
                #lof_renumOss(nm,j,osss)
            end
            j=j+1
        end
    end
end
###############################################################################
function lof_wndPfOss(osss,ocn)
    avePcc=lof_avePnt(ocn.pccs)
    nos=lod_gen2Noss()
    for i in osss
        for j in ocn.cnces
            for k in ocn.pccs
                if (lof_pnt2pnt_dist(i.coord,k.coord) <= lof_pnt2pnt_dist(j.coord,k.coord)+nos)
                    push!(i.wnds,j.wnd)
                    push!(i.mvas,j.mva)
                    @goto wnd_stored
                else
                end
            end
            @label wnd_stored
        end
    end
end
###############################################################################
function lof_osss(ocn,cnt)
    num=length(ocn.pccs)+length(ocn.cnces)+1
    cns=reverse(ocn.cnces,1)
    osss=Array{oss,1}()


    for i=1:length(cns)
        if cnt.xrad==true
            num=lof_ossXradius(i,num,cns,osss,cnt)#lay oss version 3
        end
        if cnt.neib1==true
            num=lof_oss1neibs(i,num,cns,osss,cnt)#lay oss version 3
        end
        if cnt.neib3==true
            num=lof_oss3neibs(i,num,cns,osss,cnt)
        end
        if cnt.xradPcc==true
            num=lof_ossXradPcc(i,num,cns,ocn.pccs,osss,lod_rad(),false)#lay oss version 3
        end
        if (i<=lod_frcNum() && cnt.xradHlf==true)
            num=lof_ossXradPcc(i,num,cns,ocn.pccs,osss,lod_rdFrc(),true)
        end
    end
    nl=1
    ol=2
    if lod_ossSpcfy()==true
        while nl != ol
            ol=deepcopy(length(osss))
            lof_ossSprcfy(osss)
            nl=deepcopy(length(osss))
        end
    end
    osss=lof_ossOrder(osss,ocn)
    lof_wndPfOss(osss,ocn)
    ocn.osss=deepcopy(osss)
end
###############################################################################
function lof_ossOrder(osss,ocn)
    lnths=Array{Float64,1}()
    ordrdOsss=Array{oss,1}()
    num=length(ocn.pccs)+length(ocn.cnces)
    for oss in osss
        pcc_close=lof_xClosestPcc(oss,ocn.pccs)
        push!(lnths,lof_pnt2pnt_dist(oss.coord,pcc_close.coord))
    end
    for lps=1:1:length(lnths)
        push!(ordrdOsss,osss[findmax(lnths)[2]])
        ordrdOsss[length(ordrdOsss)].num=deepcopy(num+lps)
        lnths[findmax(lnths)[2]]=0
    end
    osss = ordrdOsss
    return osss
end
###############################################################################
########################## Oss Layout schemes #################################
###############################################################################
#find closest PCC
function lof_xClosestPcc(i,pc)
    x_close=pcc()
    lnths=Array{Float64,1}()
    for j=1:length(pc)
        push!(lnths,lof_pnt2pnt_dist(i.coord,pc[j].coord))
    end
    mn=findmin(lnths)
    x_close=pc[mn[2]]
    return x_close
end
###############################################################################
#find x closest owpp neighbours
#=function lof_xClosestCns(i,cns,x)
    x_close=Array{cnce,1}()
    lnths=Array{Float64,1}()
    for j=i+1:1:length(cns)
        push!(lnths,lof_pnt2pnt_dist(cns[i].coord,cns[j].coord))
    end
    if x>length(lnths)
        x=length(lnths)
    end
    while length(x_close)<x && length(lnths) != 0
        mn=findmin(lnths)
        push!(x_close,cns[mn[2]+i])
        lnths[mn[2]]=Inf
    end
    return x_close
end=#
###############################################################################
#find x closest owpp neighbours
function lof_xClosestCns(i,cns,x)
    x_close=Array{cnce,1}()
    lnths=Array{Float64,1}()
    for j=i+1:1:length(cns)
        push!(lnths,lof_pnt2pnt_dist(cns[i].coord,cns[j].coord))
    end
    while (length(lnths) != 0 && length(x) != 0 && findmax(x)[1]>length(lnths))
        deleteat!(x,findmax(x)[2])
    end
    lp=1
    while length(x_close)<length(x) && length(lnths) != 0
        mn=findmin(lnths)
        if lp in x
            push!(x_close,cns[mn[2]+i])
        end
        lnths[mn[2]]=Inf
        lp=lp+1
    end
    return x_close
end
############################### Start #########################################
######################## sets OSS between owpp and pccs  ######################
###############################################################################
function lof_ossXradPcc(i,num,cns,pccs,osss,rad,frac)
    x=lof_xClosestPcc(cns[i],pccs)
    if frac == true
        rad=rad*lof_pnt2pnt_dist(cns[i].coord,x.coord)
    end
    osub=oss()
    alpha_beta=reverse([[cns[i].coord.y,x.coord.y] ones(2)]\[cns[i].coord.x,x.coord.x])#fits linear model
    osub.coord=lof_atXPcc(alpha_beta,cns[i],x,rad)
    osub.id="2"*string(trunc(Int,10*osub.coord.x))*string(trunc(Int,10*osub.coord.y))
    push!(osss,deepcopy(osub))
    num=num+1
    return num
end
###############################################################################
#Sorts the relative position of OWPP to pcc
function lof_atXPcc(mb,p1,pcc,rad)
    xy1=xy()
    if (p1.coord.x == pcc.coord.x && p1.coord.y > pcc.coord.y)
        xy1.y=p1.coord.y-rad
        xy1.x=p1.coord.x
    elseif (p1.coord.x == pcc.coord.x && p1.coord.y < pcc.coord.y)
        xy1.y=p1.coord.y+rad
        xy1.x=p1.coord.x
    elseif (p1.coord.y == pcc.coord.y && p1.coord.x > pcc.coord.x)
        xy1.x=p1.coord.x-rad
        xy1.y=p1.coord.y
    elseif (p1.coord.y == pcc.coord.y && p1.coord.x < pcc.coord.x)
        xy1.x=p1.coord.x+rad
        xy1.y=p1.coord.y
    else
         lof_solvIntersPcc(xy1,p1,rad,mb)
    end
    return xy1
end
###############################################################################
function lof_solvIntersPcc(xy1,p1,os,mb)
    xy1.y=p1.coord.y-os
    xy1.x=xy1.y*mb[2]+mb[1]
    while lof_pnt2pnt_dist(p1.coord,xy1)<os
        xy1.y=xy1.y-0.1
        xy1.x=xy1.y*mb[2]+mb[1]
    end
end
############################### Start #########################################
###############sets OSS X OSS on the line to neighbouring owpp  ###############
###############################################################################
function lof_mdOss(c1,c2)
    xy0=xy()
    xy0.x=(c1.coord.x+c2.coord.x)/2
    xy0.y=(c1.coord.y+c2.coord.y)/2
    return xy0
end
###############################################################################
function lof_oss1neibs(i,num,cns,osss,cnt)
    x=cnt.xXneib1
    x_close=lof_xClosestCns(i,cns,x)
    for x in x_close
        osub=oss()
        osub.coord=lof_mdOss(x,cns[i])
        osub.id="2"*string(trunc(Int,10*osub.coord.x))*string(trunc(Int,10*osub.coord.y))
        push!(osss,deepcopy(osub))
        num=num+1
    end
    return num
end
###############################################################################
function lof_oss3neibs(i,num,cns,osss,cnt)
    x=cnt.xXneib3
    x_close=lof_xClosestCns(i,cns,x)
    xy0=cnce()
    for x in x_close
        osub1=oss()
        osub2=oss()
        xy0.coord=lof_mdOss(x,cns[i])
        osub1.coord=lof_mdOss(xy0,cns[i])
        osub2.coord=lof_mdOss(x,xy0)
        osub1.id="2"*string(trunc(Int,10*osub1.coord.x))*string(trunc(Int,10*osub1.coord.y))
        num=num+1
        osub2.id="2"*string(trunc(Int,10*osub2.coord.x))*string(trunc(Int,10*osub2.coord.y))
        push!(osss,deepcopy(osub1))
        push!(osss,deepcopy(osub2))
        num=num+1
    end
    return num
end
############################### Start #########################################
###sets OSS on the radius x cicle on the line to neighbouring owpp  ###########
###############################################################################
#finds the point on the circle of radius r around the owpp that intersects the line to a neighbouring OWPP
function lof_solvIntersect(xy1,xy2,p1,p2,os,mb)
    xy1.y=p1.coord.y+os
    xy1.x=xy1.y*mb[2]+mb[1]
    while lof_pnt2pnt_dist(p1.coord,xy1)>2
        xy1.y=xy1.y-0.1
        xy1.x=xy1.y*mb[2]+mb[1]
    end

    xy2.y=p2.coord.y-os
    xy2.x=xy2.y*mb[2]+mb[1]
    while lof_pnt2pnt_dist(p2.coord,xy2)>2
        xy2.y=xy2.y+0.1
        xy2.x=xy2.y*mb[2]+mb[1]
    end
end
###############################################################################
#Sorts the relative position of one OWPP to another
function lof_atX(mb,p1,p2)
    os=lod_rad()
    xy1=xy()
    xy2=xy()
    if (p1.coord.x == p2.coord.x && p1.coord.y > p2.coord.y)
        xy1.y=p1.coord.y-os
        xy2.y=p2.coord.y+os
    elseif (p1.coord.x == p2.coord.x && p1.coord.y < p2.coord.y)
        xy1.y=p1.coord.y+os
        xy2.y=p2.coord.y-os
    elseif (p1.coord.y == p2.coord.y && p1.coord.x > p2.coord.x)
        xy1.x=p1.coord.x-os
        xy2.x=p2.coord.x+os
    elseif (p1.coord.y == p2.coord.y && p1.coord.x < p2.coord.x)
        xy1.x=p1.coord.x+os
        xy2.x=p2.coord.x-os
    elseif (p1.coord.x < p2.coord.x && p1.coord.y < p2.coord.y)
         lof_solvIntersect(xy1,xy2,p1,p2,os,mb)
    elseif (p1.coord.x < p2.coord.x && p1.coord.y > p2.coord.y)
        lof_solvIntersect(xy2,xy1,p2,p1,os,mb)
    elseif (p1.coord.x > p2.coord.x && p1.coord.y < p2.coord.y)
        lof_solvIntersect(xy1,xy2,p1,p2,os,mb)
    elseif (p1.coord.x > p2.coord.x && p1.coord.y > p2.coord.y)
        lof_solvIntersect(xy2,xy1,p2,p1,os,mb)
    else
        println("Caution: No OSS Xkm radius placement matched!")
    end
    return xy1,xy2
end
###############################################################################
#main control func to lay OSS on radius around the OWPP
function lof_ossXradius(i,num,cns,osss,cnt)
    x=cnt.xXrad
    x_close=lof_xClosestCns(i,cns,x)
    for x in x_close
        osub1=oss()
        osub2=oss()
        alpha_beta=reverse([[cns[i].coord.y,x.coord.y] ones(2)]\[cns[i].coord.x,x.coord.x])#fits linear model
        osub1.coord,osub2.coord=lof_atX(alpha_beta,cns[i],x)
        osub1.id="2"*string(trunc(Int,10*osub1.coord.x))*string(trunc(Int,10*osub1.coord.y))
        num=num+1
        osub2.id="2"*string(trunc(Int,10*osub2.coord.x))*string(trunc(Int,10*osub2.coord.y))
        push!(osss,deepcopy(osub1))
        push!(osss,deepcopy(osub2))
        num=num+1
    end
    return num
end
############################### Start #########################################
####### sets x number of OSS on the line between neighbouring owpp  ###########
###############################################################################
###############################################################################
#Layout scheme 1
function lof_ossLine(xy,num,rg,osss)
    spc=lod_genSpc()
    mxy=rg.cnces[length(rg.cnces)].coord.y+loD_nbuff()
    mny=rg.cnces[1].coord.y-loD_sbuff()
    elx=lof_eLim(xy.y,rg.bnd.ebnd)
    wlx=lof_wLim(xy.y,rg.bnd.wbnd)

    #add oss to the north of generation
    if (xy.y+spc <= mxy && lod_noss()==true)
        x=max(xy.x, lof_wLim(xy.y+spc,rg.bnd.wbnd))
        x=min(x, lof_eLim(xy.y+spc,rg.bnd.ebnd))
        osub=oss()
        osub.coord.x=deepcopy(x)
        osub.coord.y=deepcopy(xy.y+spc)
        osub.num=deepcopy(num)
        push!(osss,osub)
        num=num+1
    else
    end


    #add oss to the east of generation
    if (xy.x+spc <= elx && lod_eoss()==true)
        osub=oss()
        osub.coord.x=deepcopy(xy.x+spc)
        osub.coord.y=deepcopy(xy.y)
        osub.num=deepcopy(num)
        push!(osss,osub)
        num=num+1
    else
    end
    #add oss to the east of generation
    if (xy.x != elx && lod_eosss()==true)
        i=1
        while (xy.x+i*spc <= elx)
            x=min(xy.x+(i+1)*spc, elx)
            osub=oss()
            osub.coord.x=deepcopy(x)
            osub.coord.y=deepcopy(xy.y)
            osub.num=deepcopy(num)
            push!(osss,osub)
            num=num+1
            i=i+1
        end
    else
    end

    #Oss at gen location
    #=if (lod_goss()==true)
        osub=oss()
        osub.coord.x=deepcopy(xy.x)
        osub.coord.y=deepcopy(xy.y)
        osub.num=deepcopy(num)
        push!(osss,osub)
        num=num+1
    else
    end=#

    #add oss to the south of generation
    if (xy.y-spc >= mny && lod_soss()==true)
        osub=oss()
        x=max(xy.x, lof_wLim(xy.y-spc,rg.bnd.wbnd))
        x=min(x, lof_eLim(xy.y-spc,rg.bnd.ebnd))
        osub.coord.x=deepcopy(x)
        osub.coord.y=deepcopy(xy.y-spc)
        osub.num=deepcopy(num)
        push!(osss,osub)
        num=num+1
    end
    #add oss to the west of generation
    if (xy.x-spc >= wlx && lod_woss()==true)
        osub=oss()
        osub.coord.x=deepcopy(xy.x-spc)
        osub.coord.y=deepcopy(xy.y)
        osub.num=deepcopy(num)
        push!(osss,osub)
        num=num+1
    else
    end
    #add oss to the west of generation
    if (xy.x != wlx && lod_wosss()==true)
        i=1
        while (xy.x-i*spc >= wlx)
            x=max(xy.x-(i+1)*spc, wlx)
            osub=oss()
            osub.coord.x=deepcopy(x)
            osub.coord.y=deepcopy(xy.y)
            osub.num=deepcopy(num)
            push!(osss,osub)
            num=num+1
            i=i+1
        end
    else
    end
    return num
end
###############################################################################
#Layout scheme 2
function lof_ossLine2(xy,num,rg,osss)
    spc=lod_genSpc()
    mxy=rg.cnces[length(rg.cnces)].coord.y+loD_nbuff()
    mny=rg.cnces[1].coord.y-loD_sbuff()
    elx=lof_eLim(xy.y,rg.bnd.ebnd)
    wlx=lof_wLim(xy.y,rg.bnd.wbnd)



    #Oss at gen location
    if (lod_goss()==true)
        osub=oss()
        osub.coord.x=deepcopy(xy.x)
        osub.coord.y=deepcopy(xy.y)
        osub.num=deepcopy(num)
        push!(osss,osub)
        num=num+1
    else
    end
    for i in rg.cnces
        osub=oss()
        osub.coord=lof_avePnt([xy,i.coord])
        osub.num=deepcopy(num)
        push!(osss,osub)
        num=num+1
    end
    return num
end
###############################################################################
###############################################################################
#Layout scheme 3
function lof_ossLine3(i,num,rg,cns,osss)

    xy=cns[i].coord
    elx=lof_eLim(xy.y,rg.bnd.ebnd)
    wlx=lof_wLim(xy.y,rg.bnd.wbnd)



    #Oss at gen location
    if (lod_goss()==true)
        osub=oss()
        osub.coord.x=deepcopy(xy.x)
        osub.coord.y=deepcopy(xy.y)
        osub.num=deepcopy(num)
        push!(osss,osub)
        num=num+1
    else
    end
    if i<length(cns)
        osub=oss()
        osub.coord=lof_avePnt([xy,cns[i+1].coord])
        osub.num=deepcopy(num)
        push!(osss,osub)
        num=num+1
    end
    return num
end
###############################################################################
###############################################################################


###########################################################################################################
########################################### arcs ##########################################################
###########################################################################################################
function lof_buldGoArc(tl,hd,km)
    #push!(hd.wnds,tl.wnd)
    #push!(hd.mvas,tl.mva)
    a=gOarc()
    a.head=hd
    a.tail=tl
    #a.mva=tl.mva
    #a.kv=tl.kv
    a.lngth=deepcopy(km)
    #a.wnd=tl.wnd
    return a
end
###############################################################################
function lof_buldOpArc(tl,hd,km)
    a=oParc()
    a.head=hd
    a.tail=tl
    #a.mva=tl.mva
    #a.kv=tl.kv
    a.lngth=deepcopy(km)
    #a.wnd=tl.wnd
    return a
end
###############################################################################
function lof_buldGpArc(tl,hd,km)
    a=gParc()
    a.head=hd
    a.tail=tl
    #a.mva=tl.mva
    #a.kv=tl.kv
    a.lngth=deepcopy(km)
    #a.wnd=tl.wnd
    return a
end
###############################################################################
function lof_buldOoArc(tl,hd,km)
    a=oOarc()
    a.head=hd
    a.tail=tl
    #a.mva=tl.mva
    #a.kv=tl.kv
    a.lngth=deepcopy(km)
    #a.wnd=tl.wnd
    return a
end
###############################################################################
#OSS to OSS connection
function lof_OpArcs(ocn)
    for i in ocn.osss
        for j in ocn.pccs
            km=lof_pnt2pnt_dist(i.coord,j.coord)
            push!(ocn.oParcs,lof_buldOpArc(i,j,km))
        end
    end
end
###############################################################################
#average point of PCCs
function lof_avePnt(vec)
    XY=xy()
    XY.x=0
    XY.y=0
    tot=0
    for i in vec
        if typeof(i) != typeof(xy())
            XY.x=XY.x+i.coord.x
            XY.y=XY.y+i.coord.y
        else
            XY.x=XY.x+i.x
            XY.y=XY.y+i.y
        end
        tot=tot+1
    end
    XY.x=XY.x/tot
    XY.y=XY.y/tot
    return XY
end
###############################################################################
#generator to OSS connection
function lof_GoArcs(ocn)
    avePcc=lof_avePnt(ocn.pccs)
    nos=lod_gen2Noss()
    for i in ocn.cnces
        mxKm=lod_mxMvKm(i)
        for j in ocn.osss
            km=lof_pnt2pnt_dist(i.coord,j.coord)
            if (km <= mxKm && lof_pnt2pnt_dist(i.coord,avePcc)+nos >= lof_pnt2pnt_dist(j.coord,avePcc))
            #if km <= mxKm
                push!(ocn.gOarcs,lof_buldGoArc(i,j,km))
            else
            end
        end
    end
end
###############################################################################
#OSS to OSS connection
function lof_OoArcs(ocn)
    for i=1:length(ocn.osss)
        mnKm=lod_mnKm()
        for j=(i+1):length(ocn.osss)
            km=lof_pnt2pnt_dist(ocn.osss[i].coord,ocn.osss[j].coord)
            if mnKm <= km
                push!(ocn.oOarcs,lof_buldOoArc(ocn.osss[i],ocn.osss[j],km))
            else
            end
        end
    end
end
##############################################################################
#generator to Pcc connection
function lof_GpArcs(ocn)
    for i in ocn.cnces
        mxKm=lod_mxMv2PccKm(i)
        for j in ocn.pccs
            km=lof_pnt2pnt_dist(i.coord,j.coord)
            if km <= mxKm
                push!(ocn.gParcs,lof_buldGpArc(i,j,km))
            else
            end
        end
    end
end
###############################################################################



















################################################################################
############################## Removed Functions ###############################
###############################################################################

###############################################################################
############################ border removed #######################################
###############################################################################
#finds the western boundary at y=y
#=function lof_wLim(y,bd)
    wlx=0
    for i=1:length(bd.lims)-1
        if (bd.lims[i].y <= y && y <= bd.lims[i+1].y)
            wlx=bd.lmodel[i].m*y+bd.lmodel[i].b
        end
    end
    return wlx
end
###############################################################################
#finds the eastern boundary at y=y
function lof_eLim(y,bd)
    elx=0
    for i=1:length(bd.lims)-1
        if (bd.lims[i+1].y <= y && y <= bd.lims[i].y)
            elx=bd.lmodel[i].m*y+bd.lmodel[i].b
        end
    end
    return elx
end
#################################################################################
#=finds the full boundary points around the concession and makes linear boundary model
function lof_ewbnd(ocn)
    #dummy arrays
    all=Array{xy,1}()
    all_oss=Array{xy,1}()
    Wbnd=Array{xy,1}()
    Ebnd=Array{xy,1}()

    strt=ocn.cnces[1].coord#define 0 point
    fnsh=ocn.cnces[length(ocn.cnces)].coord#define furthest concession
    push!(Wbnd,deepcopy(strt))
    push!(Ebnd,deepcopy(strt))
    for i in ocn.cnces[2:length( ocn.cnces)]
        push!(all,deepcopy(i.coord))
    end

    #sets east and west boundary
    Wbnd,all=lof_wBnd(Wbnd,fnsh,all)#finds western boundary
    Ebnd=lof_eBnd(Ebnd,fnsh,deepcopy(all))#finds eastern boundary
    Wbnd=lof_addWBuff(Wbnd)#add buffer to western boundary
    Ebnd=lof_addEBuff(Ebnd)#add buffer to eastern boundary
    ocn.bnd.wbnd.lims=Wbnd
    ocn.bnd.ebnd.lims=reverse(Ebnd,1)

    #set calculated n-s boundary points
    push!(ocn.bnd.nbnd.lims,ocn.bnd.wbnd.lims[length(ocn.bnd.wbnd.lims)])
    push!(ocn.bnd.nbnd.lims,ocn.bnd.ebnd.lims[1])
    push!(ocn.bnd.sbnd.lims,ocn.bnd.ebnd.lims[length(ocn.bnd.ebnd.lims)])
    push!(ocn.bnd.sbnd.lims,ocn.bnd.wbnd.lims[1])

    #calculates linear model for boundary
    lof_lnrBnd(ocn.bnd.ebnd)
    lof_lnrBnd(ocn.bnd.wbnd)
    lof_lnrBnd(ocn.bnd.sbnd)
    lof_lnrBnd(ocn.bnd.nbnd)

    return nothing
end=#
###############################################################################
#=function lof_cnsPeri(cns)
    hgts=[]
    wdts=[]
    sps=[]
    push!(hgts, cns[2].coord.y-cns[1].coord.y)
    push!(wdts, (cns[1].area/hgts[1])/2)
    push!(sps, cns[2].coord.y-cns[1].coord.y)
    for i=1:length(cns)-2
        push!(sps, cns[i+1].coord.y-cns[i].coord.y)
        push!(hgts, (cns[i+1].coord.y-cns[i].coord.y+cns[i+2].coord.y-cns[i+1].coord.y)/2)
        push!(wdts, (cns[i+1].area/hgts[i+1])/2)
    end
    push!(sps, cns[length(cns)].coord.y-cns[length(cns)-1].coord.y)
    push!(sps, cns[length(cns)].coord.y-cns[length(cns)-1].coord.y)
    push!(hgts, cns[length(cns)].coord.y-cns[length(cns)-1].coord.y)
    push!(wdts, (cns[length(cns)].area/hgts[length(cns)])/2)
    maxl=cns[1].coord.x-wdts[1]
    maxr=cns[1].coord.x+wdts[1]
    for i=2:length(cns)
        println(cns[i].coord.x)
        if maxl>cns[i].coord.x-wdts[i]
            maxl=cns[i].coord.x-wdts[i]
        end
        if maxr<cns[i].coord.x+wdts[i]
            maxr=cns[i].coord.x+wdts[i]
        end
    end
    for i=1:length(cns)
        cns[i].coord.nw.x=cns[i].coord.x-wdts[i]
        cns[i].coord.ne.x=cns[i].coord.x+wdts[i]
        cns[i].coord.sw.x=cns[i].coord.x-wdts[i]
        cns[i].coord.se.x=cns[i].coord.x+wdts[i]
        cns[i].coord.nw.y=cns[i].coord.y+sps[i+1]/2
        cns[i].coord.ne.y=cns[i].coord.y+sps[i+1]/2
        cns[i].coord.sw.y=cns[i].coord.y-sps[i]/2
        cns[i].coord.se.y=cns[i].coord.y-sps[i]/2
        println(cns[i].coord)
    end
end=#
###############################################################################
#=function lof_sbnd(ocn)
    cns=Array{Float64,1}()
    for i in ocn.cnces
        push!(cns,i.coord.y)
    end
    tmp=findmin(cns)[2]
    tmp1=findmax(cns)[2]
    cns[tmp]=cns[tmp1]+10
    cls0=deepcopy(ocn.cnces[tmp].coord)
    cls0.y=cls0.y-buffer
    tmp=findmin(cns)[2]
    cls1=deepcopy(ocn.cnces[tmp].coord)
    cls1.y=cls1.y-buffer
    cns=[]
    push!(cns,cls0)
    push!(cns,cls1)
    return cns
end
##########################################################################################
#find the full boundary around the concessions
function lof_outerbnd(ocn)
    all=Array{xy,1}()
    all_pcc=Array{xy,1}()
    all_oss=Array{xy,1}()
    Wbnd=Array{xy,1}()
    Ebnd=Array{xy,1}()
    strt=ocn.pccs[length(ocn.pccs)].coord#define 0 point
    fnsh=ocn.cnces[length(ocn.cnces)].coord#define furthest concession
    push!(Wbnd,strt)
    push!(Ebnd,strt)
#combine all points into single array for processing
    for i=1:length(ocn.pccs)-1
        push!(all_pcc,deepcopy(ocn.pccs[i].coord))
        push!(all,all_pcc[length(all_pcc)])
    end

    for i in ocn.cnces
        #push!(all_oss,deepcopy(i.coord))
        push!(all,deepcopy(i.coord))
    end

    Wbnd,all=lof_wBnd(Wbnd,fnsh,all)#finds western boundary
    println(Wbnd)
    Ebnd,all=lof_eBnd(Ebnd,fnsh,deepcopy(all))#finds eastern boundary

    buffer=loD_ewbuff()
    Wbnd=lof_addBuff(Wbnd,-1*buffer,all_pcc)#add buffer to western boundary

    Ebnd=lof_addBuff(Ebnd,buffer,all_pcc)#add buffer to eastern boundary

    Ebnd=reverse(Ebnd,1)
    for i in Ebnd
        push!(Wbnd,i)
    end
    return Wbnd
end=#
#################################################################################


###############################################################################
########################### oss placement old #################################
###############################################################################
#=function lof_farthest(pccs,rg)
    fthst=xy()
    l1=0
    l2=0

    for i in pccs
        l2=l2+lof_pnt2pnt_dist(i.coord,rg.bnd.nbnd.lims[2])
        l1=l1+lof_pnt2pnt_dist(i.coord,rg.bnd.nbnd.lims[1])
    end
    if l1>=l2
        spc=rg.cnces[1].trb.dia*lod_ossSpc()
        fthst=rg.bnd.nbnd.lims[1]
    else
        fthst=rg.bnd.nbnd.lims[2]
        spc=-1*(rg.cnces[1].trb.dia)*lod_ossSpc()
    end
    return fthst, spc
end=#
###############################################################################
function lof_ossLine(y,elx,wlx,spc,num,osss)
    if spc>0
        x=wlx
    elseif spc<0
        x=elx
    else
        error("OSS spacing is 0!")
    end
    while (x <= elx && wlx <= x)
        osub=oss()
        osub.coord.x=deepcopy(x)
        osub.coord.y=deepcopy(y)
        osub.num=num
        push!(osss,osub)
        x=x+spc
        num=num+1
    end
    return num
end
###############################################################################
function lof_osss(ocn)
    strtPnt,spc=lof_farthest(ocn.pccs,ocn)#selects either nw or ne to start and sets the spcing
    ns=strtPnt.y
    num=length(ocn.pccs)+length(ocn.cnces)+1
    while ns>ocn.bnd.sbnd.lims[1].y
        elx=lof_eLim(ns,ocn.bnd.ebnd)
        wlx=lof_wLim(ns,ocn.bnd.wbnd)
        num=lof_ossLine(ns,elx,wlx,spc,num,ocn.osss)#lay a line of oss
        ns=ns-abs(spc)
    end
end=#
###############################################################################















###############################################################################
#############################$# Repeats #######################################
###############################################################################
#=function lof_farthest(pccs,rg)
    fthst=xy()
    l1=0
    l2=0

    for i in pccs
        l2=l2+lof_pnt2pnt_dist(i.coord,rg.bnd.nbnd.lims[2])
        l1=l1+lof_pnt2pnt_dist(i.coord,rg.bnd.nbnd.lims[1])
    end
    if l1>=l2
        spc=rg.cnces[1].trb.dia*lod_ossSpc()
        fthst=rg.bnd.nbnd.lims[1]
    else
        fthst=rg.bnd.nbnd.lims[2]
        spc=-1*(rg.cnces[1].trb.dia)*lod_ossSpc()
    end
    return fthst, spc
end=#
###############################################################################
#finds the western boundary at y=y
#=function lof_wLim(y,bd)
    wlx=0
    for i=1:length(bd.lims)-1
        if (bd.lims[i].y <= y && y <= bd.lims[i+1].y)
            wlx=bd.lmodel[i].m*y+bd.lmodel[i].b
        end
    end
    return wlx
end
###############################################################################
#finds the eastern boundary at y=y
function lof_eLim(y,bd)
    elx=0
    for i=1:length(bd.lims)-1
        if (bd.lims[i+1].y <= y && y <= bd.lims[i].y)
            elx=bd.lmodel[i].m*y+bd.lmodel[i].b
        end
    end
    return elx
end
###############################################################################
function lof_ossLine(y,elx,wlx,spc,num,osss)
    if spc>0
        x=wlx
    elseif spc<0
        x=elx
    else
        error("OSS spacing is 0!")
    end
    while (x <= elx && wlx <= x)
        osub=oss()
        osub.coord.x=deepcopy(x)
        osub.coord.y=deepcopy(y)
        osub.num=num
        push!(osss,osub)
        x=x+spc
        num=num+1
    end
    return num
end
###############################################################################
function lof_osss(ocn)
    strtPnt,spc=lof_farthest(ocn.pccs,ocn)#selects either nw or ne to start and sets the spcing
    ns=strtPnt.y
    num=length(ocn.pccs)+length(ocn.cnces)+1
    while ns>ocn.bnd.sbnd.lims[1].y
        elx=lof_eLim(ns,ocn.bnd.ebnd)
        wlx=lof_wLim(ns,ocn.bnd.wbnd)
        num=lof_ossLine(ns,elx,wlx,spc,num,ocn.osss)#lay a line of oss
        ns=ns-abs(spc)
    end
end=#
###############################################################################
