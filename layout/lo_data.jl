################################################################################
function lod_concessionAreas()
    areas=[38.0,19.0,23.0,14.0,18.0,35.0,12.0,16.0]
    return areas
end
################################################################################
function lod_cncesGps()
    c0=[]
    c1=[]
    c2=[]
    c3=[]
    c4=[]
    c5=[]
    c6=[]
    c7=[]
    c=[]
    #concessions
    #Norther
    push!(c0,3.015833)
    push!(c0,51.52806)
    push!(c,c0)
    #Thornton
    push!(c1,(2.97+2.919972)/2)
    push!(c1,(51.56+51.53997)/2)
    push!(c,c1)
    #Rentel
    push!(c2,2.939972)
    push!(c2,51.59)
    push!(c,c2)
    #Northwind
    push!(c3,2.900972)
    push!(c3,51.61897)
    push!(c,c3)
    #Seastar
    push!(c4,2.859972)
    push!(c4,51.63)
    push!(c,c4)
    #Nobelwind/Belwind
    push!(c5,(2.819972+2.799972)/2)
    push!(c5,(51.664+51.67)/2)
    push!(c,c5)
    #Northwester
    push!(c6,2.757)
    push!(c6,51.68597)
    push!(c,c6)
    #Mermaid
    push!(c7,2.74)
    push!(c7,51.71997)
    push!(c,c7)#==#
    return c
end
################################################################################
function lod_pccGps()
    pcc0=Array{Float64,1}()
    pcc1=Array{Float64,1}()
    pcc=[]
    #PCCs
    base_lg=2.941944
    base_lt=51.24306
    push!(pcc0,base_lg)
    push!(pcc0,base_lt)
    push!(pcc,pcc0)
    push!(pcc1,3.183611)
    push!(pcc1,51.32694)
    push!(pcc,pcc1)
    return pcc
end
################################################################################
