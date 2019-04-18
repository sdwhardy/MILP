###################################################################
mutable struct xy
      x::Float64
      y::Float64
end
xy()=xy(69.69,69.69)
###################################################################
mutable struct line
      m::Float64
      b::Float64
end
line()=line(69.69,69.69)
###################################################################
mutable struct gps
      lat::Float64
      lng::Float64
end
gps()=gps(69.69,69.69)
###################################################################
mutable struct coord
      #nw::xy
      #ne::xy
      #sw::xy
      #se::xy
      cnt::xy
end
#coord()=coord(xy(),xy(),xy(),xy(),xy())
coord()=coord(xy())
###################################################################
mutable struct turb
      wsp::Array{Float64}
      pwr::Array{Float64}
      dia::Float64
      mva::Float64
      cin::Float64
      cout::Float64
end
turb()=turb([],[],69.69,69.69,69.69,69.69)
###################################################################
mutable struct cnce
      gps::gps
      #name::String
      #area::Float64
      coord::coord
      mva::Float64
      wnd::wind
      trb::turb
      kv::Float64
      num::Int64
end
#cnce()=cnce(gps(),"colruyt",69.69,coord())
cnce()=cnce(gps(),coord(),69.69,wind(),turb(),69.69,69)
###################################################################
mutable struct oss
      coord::coord
      mvas::Array{Float64}
      wnds::Array{wind}
      num::Int64
end
oss()=oss(coord(),[],[],69)
###################################################################
mutable struct pcc
      gps::gps
      coord::coord
      kv::Float64
      num::Int64
end
pcc()=pcc(gps(),coord(),69.69,69)
##############################################################################################################################
##################################################### arcs ###################################################################
##############################################################################################################################
mutable struct oOarc
      head::oss
      tail::oss
      #mva::Float64
      #kv::Float64
      cost::Float64
      xr::Float64
      xl::Float64
      yb::Float64
      lngth::Float64
      #wnd::wind
end
oOarc()=oOarc(oss(),oss(),69.69,69.69,69.69,69.69,69.69)
##############################################################################################################################
mutable struct oParc
      head::pcc
      tail::oss
      #mva::Float64
      #kv::Float64
      cost::Float64
      xr::Float64
      xl::Float64
      yb::Float64
      lngth::Float64
      #wnd::wind
end
oParc()=oParc(pcc(),oss(),69.69,69.69,69.69,69.69,69.69)
###################################################################
mutable struct gOarc
      head::oss
      tail::cnce
      #mva::Float64
      #kv::Float64
      cost::Float64
      xr::Float64
      xl::Float64
      yb::Float64
      lngth::Float64
      #wnd::wind
end
gOarc()=gOarc(oss(),cnce(),69.69,69.69,69.69,69.69,69.69)
###################################################################
mutable struct lnbnd
      lims::Array{xy}
      lmodel::Array{line}
end
lnbnd()=lnbnd([],[])
###################################################################
mutable struct domain
      ebnd::lnbnd
      wbnd::lnbnd
      sbnd::lnbnd
      nbnd::lnbnd
end
domain()=domain(lnbnd(),lnbnd(),lnbnd(),lnbnd())
###################################################################
mutable struct region
      #gps::gps
      #coord::coord
      #area::Float64
      #ncons::Float64
      osss::Array{oss}
      cnces::Array{cnce}
      gOarcs::Array{gOarc}
      oOarcs::Array{oOarc}
      oParcs::Array{oParc}
      #gParcs::Array{gParc}
      #oOarcs::Array{oOarc}
      #oParcs::Array{oParc}
      bnd::domain
      #sth::Array{xy}
end
#region()=region(gps(),coord(),69.69,69.69,[],domain(),[])
region()=region([],[],[],[],[],domain())
###################################################################
mutable struct eez
      #gps::gps
      reg::region
      #area::Float64
      pccs::Array{pcc}
      #coord::coord
end
#eez()=eez(gps(),region(),69.69,[],coord())
eez()=eez(region(),[])
###################################################################
