mutable struct control
      xrad::Bool
      neib1::Bool
      neib3::Bool
      xradPcc::Bool
      xradHlf::Bool
      spcfy::Bool
      xXrad::Array{Int64,1}
      xXneib1::Array{Int64,1}
      xXneib3::Array{Int64,1}
end
control()=control(false,false,false,false,false,false,[],[],[])
###############################################################################
mutable struct oplossing
      eyeDs::Array{Int64,1}
      objBst::Float64
end
oplossing()=oplossing([],69.69)
###############################################################################
