#this file contains functions that calculate the layout of the area
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
function lof_layoutZone()
    km=1#scale
    gpss=lod_cncesGps()
    areas=lod_concessionAreas()
    zone=region()
    for i=1:length(gpss)
        concession=cnce()
        concession.gps.lng=gpss[i][1]
        concession.gps.lat=gpss[i][2]
        concession.area=areas[i]
        push!(zone.cnces,deepcopy(concession))
    end
    return zone
end
################################################################################
#sets the gps coords that are the reference coords
function lof_bseCrd(ocean)
    base=gps()
    base.lat=ocean.pccs[length(ocean.pccs)].gps.lat
    base.lng=ocean.pccs[length(ocean.pccs)].gps.lng
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
end
###############################################################################
function lof_layoutOcn()
    ocean=eez()
    ocean.reg=lof_layoutZone()
    lof_shoreConnect(ocean.pccs)
    base=lof_bseCrd(ocean)
    lof_dists(ocean,base)
    os=lof_rotateOcn(ocean)
    lof_slideOcn(ocean,os)
    lof_ewbnd(ocean)
    ppf_printOcn(ocean)
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
################################################################################
#find the full boundary around the concessions
#=function lof_ewbnd(ocn)
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
    #println(Wbnd)
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
end=#
