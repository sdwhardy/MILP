###################################################################
mutable struct xy
      x::Float64
      y::Float64
end
xy()=xy(69.69,69.69)
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
mutable struct region
      gps::gps
      coord::coord
      area::Float64
      ncons::Float64
      cnces::Array{cnce}
      bnd::Array{xy}
end
region()=region(gps(),coord(),69.69,69.69,[],[])
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
