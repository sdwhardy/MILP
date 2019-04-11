#this file describes the structure of wind objects
###################################################
mutable struct wind
      wsp::Array{Float64}
      pwr::Array{Float64}
      cp::Array{Float64}
      pu::Array{Float64}
      ce::Array{Float64}
      delta::Float64
      lf::Float64
      co::Float64
end
wind()=wind([],[],[],[],[],0.0,0.0,0.0)
###################################################
