###################################################################
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
end
cbl()=cbl(0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0)
