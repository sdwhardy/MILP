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
      nw::xy
      ne::xy
      sw::xy
      se::xy
      cnt::xy
end
coord()=coord(xy(),xy(),xy(),xy(),xy())
###################################################################
mutable struct cnce
      gps::gps
      name::String
      area::Float64
      coord::coord
end
cnce()=cnce(gps(),"colruyt",69.69,coord())
###################################################################
mutable struct pcc
      gps::gps
      coord::coord
end
pcc()=pcc(gps(),coord())
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
      fbnd::lnbnd
end
domain()=domain(lnbnd(),lnbnd(),lnbnd(),lnbnd(),lnbnd())
###################################################################
mutable struct region
      gps::gps
      coord::coord
      area::Float64
      ncons::Float64
      cnces::Array{cnce}
      bnd::domain
      sth::Array{xy}
end
region()=region(gps(),coord(),69.69,69.69,[],domain(),[])
###################################################################
mutable struct eez
      gps::gps
      reg::region
      area::Float64
      pccs::Array{pcc}
      coord::coord
end
eez()=eez(gps(),region(),69.69,[],coord())
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
