#=
This file defines the structure of any objects associated with costs
=###################################################
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
#an object that contains cost components and totals calculated
mutable struct results
     oppc::Float64
     opc::Float64
     tlc_pcc::Float64
     tlc_oss::Float64
     qc::Float64
     cbc::Float64
     rlc::Float64
     cm::Float64
     eens::Float64
     ttl::Float64
end
results()=results(0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0)
###################################################
#an object that contains all cost factors used in the calculations
mutable struct cstS_ks
   FC_ac::Float64
   FC_dc::Float64
   dc::Float64
   f_ct::Float64
   p_ct::Float64
   c_ct::Float64
   Qc_oss::Float64
   Qc_pcc::Float64
   life::Float64
   T_op::Float64
   E_op::Float64
   cf::Float64
end
cstS_ks()=cstS_ks(0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0)
###################################################################
