#this file contains functions that calculate the layout of the area
################################################################################
#returns the hypotenuse
function lof_pnt2pnt_dist(pnt1,pnt2)
    hyp=sqrt((pnt2.x-pnt1.x)^2+(pnt2.y-pnt1.y)^2)
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
#changes degrees to a length
function lof_deg2lgth(d,dPl)
    return d*dPl
end
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
    gpss=lod_cncesGps()#get gps coords of each concession
    mvas=lod_cncesMva(length(gpss))#get power of each concession
    wnds=lod_cncesWnd(length(gpss))#get wind profiles for each concession
    trbs=lod_cncesTrbs(length(gpss))#get turbine info for each concession
    #areas=lod_concessionAreas()
    zone=region()
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
        push!(zone.cnces,deepcopy(concession))
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
#loops through all coords to get relative lengthd
function lof_dist(location,base,lnthLT,km)
    for i=1:length(location)
        location[i].coord.cnt.x=lof_deg2lgth(location[i].gps.lng-base.lng,lof_lg1deg(location[i].gps.lat+base.lng)*km)
        location[i].coord.cnt.y=abs(lof_deg2lgth(location[i].gps.lat-base.lat,lnthLT))
    end
end
################################################################################
#calculates lengths based on latitude
#as lattitude changes number of km should be updated
function lof_dists(ocean,base)
    km=1
    lnth_1deg_LT=111*km
    lof_dist(ocean.reg.cnces,base,lnth_1deg_LT,km)
    lof_dist(ocean.pccs,base,lnth_1deg_LT,km)
end
################################################################################
#loops through to apply appropriate rotations for coords
function lof_rot(location,theta,os)
    for i=1:length(location)
        xy=lof_rotation(location[i].coord.cnt.x,location[i].coord.cnt.y,theta)
        location[i].coord.cnt.x=xy[1]
        location[i].coord.cnt.y=xy[2]
        if location[i].coord.cnt.x<os
            os=location[i].coord.cnt.x
        end
    end
    return os
end
################################################################################
#rotates the entire region
function lof_rotateOcn(ocean)
    theta=atan((ocean.pccs[length(ocean.pccs)].coord.cnt.x-ocean.reg.cnces[length(ocean.reg.cnces)].coord.cnt.x)/(ocean.reg.cnces[length(ocean.reg.cnces)].coord.cnt.y-ocean.pccs[length(ocean.pccs)].coord.cnt.y))
    os=0.0
    os=lof_rot(ocean.reg.cnces,theta,os)
    os=lof_rot(ocean.pccs,theta,os)
    return os
end
################################################################################
#translates the entire region
function lof_slideOcn(ocean,os)
    for i=1:length(ocean.reg.cnces)
        ocean.reg.cnces[i].coord.cnt.x=lof_shift(ocean.reg.cnces[i].coord.cnt.x,abs(os))
    end
    for i=1:length(ocean.pccs)
        ocean.pccs[i].coord.cnt.x=lof_shift(ocean.pccs[i].coord.cnt.x,abs(os))
    end
end
################################################################################
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
################################################################################
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
            pnt=findmin(thetas)[2]#select th eminimum angle
            #println(pnts[pnt])
            push!(bd,pnts[pnt])#chose this point as next boundary point
            if pnts[pnt].x != f.x && pnts[pnt].y != f.y
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
    sbuff=loD_sbuff()
    nbuff=loD_nbuff()
    for i in bnd
        i.x=i.x-wbuff
    end
    bnd[length(bnd)].y=bnd[length(bnd)].y+nbuff
    bnd[1].y=bnd[1].y-sbuff
    return bnd
end
################################################################################
#adds specified buffer to the east
function lof_addEBuff(bnd)
    ebuff=loD_ebuff()
    sbuff=loD_sbuff()
    nbuff=loD_nbuff()
    for i in bnd
        i.x=i.x+ebuff
    end
    bnd[length(bnd)].y=bnd[length(bnd)].y+nbuff
    bnd[1].y=bnd[1].y-sbuff
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

    strt=ocn.pccs[2].coord.cnt#define 0 point
    fnsh=ocn.reg.cnces[length(ocn.reg.cnces)].coord.cnt#define furthest concession
    push!(Wbnd,deepcopy(strt))
    push!(Ebnd,deepcopy(strt))
    push!(all,deepcopy(ocn.pccs[1].coord.cnt))
    for i in ocn.reg.cnces[1:length( ocn.reg.cnces)]
        push!(all,deepcopy(i.coord.cnt))
    end

    #sets east and west boundary
    Wbnd,all=lof_wBnd(Wbnd,fnsh,all)#finds western boundary
    Ebnd=lof_eBnd(Ebnd,fnsh,deepcopy(all))#finds eastern boundary
    Wbnd=lof_addWBuff(Wbnd)#add buffer to western boundary
    Ebnd=lof_addEBuff(Ebnd)#add buffer to eastern boundary
    ocn.reg.bnd.wbnd.lims=Wbnd
    ocn.reg.bnd.ebnd.lims=reverse(Ebnd,1)

    #set calculated n-s boundary points
    push!(ocn.reg.bnd.nbnd.lims,ocn.reg.bnd.wbnd.lims[length(ocn.reg.bnd.wbnd.lims)])
    push!(ocn.reg.bnd.nbnd.lims,ocn.reg.bnd.ebnd.lims[1])
    push!(ocn.reg.bnd.sbnd.lims,ocn.reg.bnd.ebnd.lims[length(ocn.reg.bnd.ebnd.lims)])
    push!(ocn.reg.bnd.sbnd.lims,ocn.reg.bnd.wbnd.lims[1])

    #calculates linear model for boundary
    lof_lnrBnd(ocn.reg.bnd.ebnd)
    lof_lnrBnd(ocn.reg.bnd.wbnd)
    lof_lnrBnd(ocn.reg.bnd.sbnd)
    lof_lnrBnd(ocn.reg.bnd.nbnd)

    return nothing
end
###############################################################################
function lof_farthest(pccs,rg)
    fthst=xy()
    l1=0
    l2=0

    for i in pccs
        l2=l2+lof_pnt2pnt_dist(i.coord.cnt,rg.bnd.nbnd.lims[2])
        l1=l1+lof_pnt2pnt_dist(i.coord.cnt,rg.bnd.nbnd.lims[1])
    end
    if l1>=l2
        spc=rg.cnces[1].trb.dia*lod_ossSpc()
        fthst=rg.bnd.nbnd.lims[1]
    else
        fthst=rg.bnd.nbnd.lims[2]
        spc=-1*(rg.cnces[1].trb.dia)*lod_ossSpc()
    end
    return fthst, spc
end
##############################################################################################################################################################
############################################################################### arcs #########################################################################
##############################################################################################################################################################
function lof_mxMvKm(cn)
    if cn.kv == 33.0
        km=10
    elseif cn.kv==66.0
        km=20
    else
        error("Cable MV does not match option!")
    end
    return km
end
function lof_mnKm()
    return 5.0
end
###############################################################################
function lof_buldGoArc(tl,hd,km)
    push!(hd.wnds,tl.wnd)
    push!(hd.mvas,tl.mva)
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
    for i in ocn.reg.osss
        for j in ocn.pccs
            km=lof_pnt2pnt_dist(i.coord.cnt,j.coord.cnt)
            push!(ocn.reg.oParcs,lof_buldOpArc(i,j,km))
        end
    end
end
###############################################################################
#generator to OSS connection
function lof_GoArcs(ocn)
    for i in ocn.reg.cnces
        mxKm=lof_mxMvKm(i)
        for j in ocn.reg.osss
            km=lof_pnt2pnt_dist(i.coord.cnt,j.coord.cnt)
            if km <= mxKm
                push!(ocn.reg.gOarcs,lof_buldGoArc(i,j,km))
            else
            end
        end
    end
end
###############################################################################
#OSS to OSS connection
function lof_OoArcs(ocn)
    for i=1:length(ocn.reg.osss)
        mnKm=lof_mnKm()
        for j=(i+1):length(ocn.reg.osss)
            km=lof_pnt2pnt_dist(ocn.reg.osss[i].coord.cnt,ocn.reg.osss[j].coord.cnt)
            if mnKm <= km
                push!(ocn.reg.oOarcs,lof_buldOoArc(ocn.reg.osss[i],ocn.reg.osss[j],km))
            else
            end
        end
    end
end
###############################################################################
############################ border ###########################################
###############################################################################
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
        osub.coord.cnt.x=deepcopy(x)
        osub.coord.cnt.y=deepcopy(y)
        osub.num=num
        push!(osss,osub)
        x=x+spc
        num=num+1
    end
    return num
end
###############################################################################
function lof_osss(ocn)
    strtPnt,spc=lof_farthest(ocn.pccs,ocn.reg)#selects either nw or ne to start and sets the spcing
    ns=strtPnt.y
    num=length(ocn.pccs)+length(ocn.reg.cnces)+1
    while ns>ocn.reg.bnd.sbnd.lims[1].y
        elx=lof_eLim(ns,ocn.reg.bnd.ebnd)
        wlx=lof_wLim(ns,ocn.reg.bnd.wbnd)
        num=lof_ossLine(ns,elx,wlx,spc,num,ocn.reg.osss)#lay a line of oss
        ns=ns-abs(spc)
    end
end=#
###############################################################################
#finds the western boundary at y=y
function lof_wLim2(y,bd)
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
function lof_eLim2(y,bd)
    elx=0
    for i=1:length(bd.lims)-1
        if (bd.lims[i+1].y <= y && y <= bd.lims[i].y)
            elx=bd.lmodel[i].m*y+bd.lmodel[i].b
        end
    end
    return elx
end
###############################################################################
function lof_ossLine2(xy,num,rg,osss)
    spc=lod_genSpc()
    mxy=rg.cnces[length(rg.cnces)].coord.cnt.y+loD_nbuff()
    mny=rg.cnces[1].coord.cnt.y-loD_sbuff()
    elx=lof_eLim2(xy.y,rg.bnd.ebnd)
    wlx=lof_wLim2(xy.y,rg.bnd.wbnd)

    #add oss to the north of generation
    if xy.y+spc <= mxy
        x=max(xy.x, lof_wLim2(xy.y+spc,rg.bnd.wbnd))
        x=min(xy.x, lof_eLim2(xy.y+spc,rg.bnd.ebnd))
        print(xy.x)
        print(" - ")
        println(xy.y)
        println(xy.y+spc)
        println(lof_wLim2(xy.y+spc,rg.bnd.wbnd))
        println(x)
        osub=oss()
        osub.coord.cnt.x=deepcopy(x)
        osub.coord.cnt.y=deepcopy(xy.y+spc)
        osub.num=num
        push!(osss,osub)
        num=num+1
    else
    end


    #add oss to the east of generation
    if xy.x+spc <= elx
        osub=oss()
        osub.coord.cnt.x=deepcopy(xy.x+spc)
        osub.coord.cnt.y=deepcopy(xy.y)
        osub.num=num
        push!(osss,osub)
        num=num+1
    else
    end
    #add oss to the east of generation
    if xy.x != elx
        x=min(xy.x+2*spc, elx)
        osub=oss()
        osub.coord.cnt.x=deepcopy(x)
        osub.coord.cnt.y=deepcopy(xy.y)
        osub.num=num
        push!(osss,osub)
        num=num+1
    else
    end

    #Oss at gen location
    osub=oss()
    osub.coord.cnt.x=deepcopy(xy.x)
    osub.coord.cnt.y=deepcopy(xy.y)
    osub.num=num
    push!(osss,osub)
    num=num+1

    #add oss to the south of generation

    if xy.y-spc >= mny
        osub=oss()
        x=max(xy.x, lof_wLim2(xy.y-spc,rg.bnd.wbnd))
        x=min(xy.x, lof_eLim2(xy.y-spc,rg.bnd.ebnd))
        osub.coord.cnt.x=deepcopy(x)
        osub.coord.cnt.y=deepcopy(xy.y-spc)
        osub.num=num
        push!(osss,osub)
        num=num+1
    end
    #add oss to the west of generation
    if xy.x-spc >= wlx
        osub=oss()
        osub.coord.cnt.x=deepcopy(xy.x-spc)
        osub.coord.cnt.y=deepcopy(xy.y)
        osub.num=num
        push!(osss,osub)
        num=num+1
    else
    end
    #add oss to the west of generation
    if xy.x != wlx
        x=max(xy.x-2*spc, wlx)
        osub=oss()
        osub.coord.cnt.x=deepcopy(x)
        osub.coord.cnt.y=deepcopy(xy.y)
        osub.num=num
        push!(osss,osub)
        num=num+1
    else
    end
    return num
end
###############################################################################
function lof_osss2(ocn)

    num=length(ocn.pccs)+length(ocn.reg.cnces)+1
    cns=reverse(ocn.reg.cnces,1)
    osss=Array{oss,1}()

    for i in cns
        gen=i.coord.cnt
        num=lof_ossLine2(gen,num,ocn.reg,osss)#lay a line of oss
    end
    ocn.reg.osss=deepcopy(osss)
end
###############################################################################
############################ Main #############################################
###############################################################################
function lof_layoutOcn()
    ocean=eez()#build the eez in the ocean
    lof_shoreConnect(ocean.pccs)#add the gps of points of common coupling
    ocean.reg=lof_layoutZone(length(ocean.pccs))#add the region with concessions
    base=lof_bseCrd(ocean)#find base coordinates
    lof_dists(ocean,base)#set layout in terms of base
    os=lof_rotateOcn(ocean)#apply rotation to align n-s with y
    lof_slideOcn(ocean,os)#slide to align western most point with x=0
    lof_ewbnd(ocean)#find the boundary of region containnig oss
    lof_osss2(ocean)#add all osss within boundary
    #lof_GoArcs(ocean)#add all gen to oss arcs within boundary
    #lof_OpArcs(ocean)#add all oss to pcc arcs within boundary
    #lof_OoArcs(ocean)#add all oss to oss arcs within boundary
    #println(length(ocean.reg.gOarcs)+length(ocean.reg.oOarcs)+length(ocean.reg.oParcs))
    ppf_printOcn(ocean)#print ocean
end
################################################################################
















################################################################################
############################## Removed Functions ###############################
################################################################################
#=function lof_cnsPeri(cns)
    hgts=[]
    wdts=[]
    sps=[]
    push!(hgts, cns[2].coord.cnt.y-cns[1].coord.cnt.y)
    push!(wdts, (cns[1].area/hgts[1])/2)
    push!(sps, cns[2].coord.cnt.y-cns[1].coord.cnt.y)
    for i=1:length(cns)-2
        push!(sps, cns[i+1].coord.cnt.y-cns[i].coord.cnt.y)
        push!(hgts, (cns[i+1].coord.cnt.y-cns[i].coord.cnt.y+cns[i+2].coord.cnt.y-cns[i+1].coord.cnt.y)/2)
        push!(wdts, (cns[i+1].area/hgts[i+1])/2)
    end
    push!(sps, cns[length(cns)].coord.cnt.y-cns[length(cns)-1].coord.cnt.y)
    push!(sps, cns[length(cns)].coord.cnt.y-cns[length(cns)-1].coord.cnt.y)
    push!(hgts, cns[length(cns)].coord.cnt.y-cns[length(cns)-1].coord.cnt.y)
    push!(wdts, (cns[length(cns)].area/hgts[length(cns)])/2)
    maxl=cns[1].coord.cnt.x-wdts[1]
    maxr=cns[1].coord.cnt.x+wdts[1]
    for i=2:length(cns)
        println(cns[i].coord.cnt.x)
        if maxl>cns[i].coord.cnt.x-wdts[i]
            maxl=cns[i].coord.cnt.x-wdts[i]
        end
        if maxr<cns[i].coord.cnt.x+wdts[i]
            maxr=cns[i].coord.cnt.x+wdts[i]
        end
    end
    for i=1:length(cns)
        cns[i].coord.nw.x=cns[i].coord.cnt.x-wdts[i]
        cns[i].coord.ne.x=cns[i].coord.cnt.x+wdts[i]
        cns[i].coord.sw.x=cns[i].coord.cnt.x-wdts[i]
        cns[i].coord.se.x=cns[i].coord.cnt.x+wdts[i]
        cns[i].coord.nw.y=cns[i].coord.cnt.y+sps[i+1]/2
        cns[i].coord.ne.y=cns[i].coord.cnt.y+sps[i+1]/2
        cns[i].coord.sw.y=cns[i].coord.cnt.y-sps[i]/2
        cns[i].coord.se.y=cns[i].coord.cnt.y-sps[i]/2
        println(cns[i].coord)
    end
end=#
###############################################################################
#=function lof_sbnd(ocn)
    cns=Array{Float64,1}()
    for i in ocn.reg.cnces
        push!(cns,i.coord.cnt.y)
    end
    tmp=findmin(cns)[2]
    tmp1=findmax(cns)[2]
    cns[tmp]=cns[tmp1]+10
    cls0=deepcopy(ocn.reg.cnces[tmp].coord.cnt)
    cls0.y=cls0.y-buffer
    tmp=findmin(cns)[2]
    cls1=deepcopy(ocn.reg.cnces[tmp].coord.cnt)
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
    strt=ocn.pccs[length(ocn.pccs)].coord.cnt#define 0 point
    fnsh=ocn.reg.cnces[length(ocn.reg.cnces)].coord.cnt#define furthest concession
    push!(Wbnd,strt)
    push!(Ebnd,strt)
#combine all points into single array for processing
    for i=1:length(ocn.pccs)-1
        push!(all_pcc,deepcopy(ocn.pccs[i].coord.cnt))
        push!(all,all_pcc[length(all_pcc)])
    end

    for i in ocn.reg.cnces
        #push!(all_oss,deepcopy(i.coord.cnt))
        push!(all,deepcopy(i.coord.cnt))
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
#=finds the full boundary points around the concession and makes linear boundary model
function lof_ewbnd(ocn)
    #dummy arrays
    all=Array{xy,1}()
    all_oss=Array{xy,1}()
    Wbnd=Array{xy,1}()
    Ebnd=Array{xy,1}()

    strt=ocn.reg.cnces[1].coord.cnt#define 0 point
    fnsh=ocn.reg.cnces[length(ocn.reg.cnces)].coord.cnt#define furthest concession
    push!(Wbnd,deepcopy(strt))
    push!(Ebnd,deepcopy(strt))
    for i in ocn.reg.cnces[2:length( ocn.reg.cnces)]
        push!(all,deepcopy(i.coord.cnt))
    end

    #sets east and west boundary
    Wbnd,all=lof_wBnd(Wbnd,fnsh,all)#finds western boundary
    Ebnd=lof_eBnd(Ebnd,fnsh,deepcopy(all))#finds eastern boundary
    Wbnd=lof_addWBuff(Wbnd)#add buffer to western boundary
    Ebnd=lof_addEBuff(Ebnd)#add buffer to eastern boundary
    ocn.reg.bnd.wbnd.lims=Wbnd
    ocn.reg.bnd.ebnd.lims=reverse(Ebnd,1)

    #set calculated n-s boundary points
    push!(ocn.reg.bnd.nbnd.lims,ocn.reg.bnd.wbnd.lims[length(ocn.reg.bnd.wbnd.lims)])
    push!(ocn.reg.bnd.nbnd.lims,ocn.reg.bnd.ebnd.lims[1])
    push!(ocn.reg.bnd.sbnd.lims,ocn.reg.bnd.ebnd.lims[length(ocn.reg.bnd.ebnd.lims)])
    push!(ocn.reg.bnd.sbnd.lims,ocn.reg.bnd.wbnd.lims[1])

    #calculates linear model for boundary
    lof_lnrBnd(ocn.reg.bnd.ebnd)
    lof_lnrBnd(ocn.reg.bnd.wbnd)
    lof_lnrBnd(ocn.reg.bnd.sbnd)
    lof_lnrBnd(ocn.reg.bnd.nbnd)

    return nothing
end=#
#################################################################################
