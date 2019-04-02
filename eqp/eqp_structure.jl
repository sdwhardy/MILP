#data structures used for equipment are specified in this file
##################################################################
mutable struct xfo_costs
   cpx::Float64
   tlc::Float64
   cm::Float64
   eens::Float64
   ttl::Float64
end
xfo_costs()=xfo_costs(0.0,0.0,0.0,0.0,0.0)
###################################################################
mutable struct xfo
   mva::Float64
   num::Float64
   eta::Float64
   fr::Float64
   mttr::Float64
   mc::Float64
   results::xfo_costs
end
xfo()=xfo(0.0,0.0,0.0,0.0,0.0,0.0,xfo_costs())
######################################################################################################################################
######################################################################################################################################
#the structure of costs for a cable
mutable struct cbl_costs
   qc::Float64
   cbc::Float64
   rlc::Float64
   cm::Float64
   eens::Float64
   ttl::Float64
end
cbl_costs()=cbl_costs(0.0,0.0,0.0,0.0,0.0,0.0)
###################################################
#the structure used for a cable
mutable struct cbl
   mva::Float64
   length::Float64
   size::Float64
   amp::Float64
   volt::Float64
   ohm::Float64
   farrad::Float64
   cost::Float64
   num::Float64
   fr::Float64
   mttr::Float64
   mc::Float64
   results::cbl_costs
end
cbl()=cbl(0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,cbl_costs())
