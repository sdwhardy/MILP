#structures used for equipment are specified in this file
###################################################################
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
