#########################################################
function ppf_cblFileName(S_min,S_max,l,kv)
	fileName="results/ACCBL"*string(trunc(Int,kv))*"kv"*string(trunc(Int,S_min))*"to"*string(trunc(Int,S_max))*"mw"*string(trunc(Int,l))*"km.csv"
	csvfile = open(fileName,"w")
	ppf_cblTtle2f(csvfile)
    return csvfile
end
#########################################################
function ppf_cblTtle2f(outFile)
	strwidth = 13
	# Print headers
	print(outFile, rpad("#1-CAP [MW]", strwidth), "\t")
	print(outFile, rpad("2-DIST [KM]", strwidth), "\t")
	print(outFile, rpad("3-TOTAL [ME]", strwidth), "\t")
	print(outFile, rpad("4-CAPEX [ME]", strwidth), "\t")
	print(outFile, rpad("5-LOSSES [ME]", strwidth), "\t")
	print(outFile, rpad("6-Q [ME]", strwidth), "\t")
	print(outFile, rpad("7-CM [ME]", strwidth), "\t")
	print(outFile, rpad("8-EENS [ME]", strwidth), "\t")
	print(outFile, rpad("9-CBL1 [#]", strwidth), "\t")
	print(outFile, rpad("10-CBL1 [kV]", strwidth), "\t")
	print(outFile, rpad("11-CBL1 [MM]", strwidth), "\t")
	print(outFile, rpad("12-CBL1 [MW]", strwidth), "\t")
end
#########################################################
function  ppf_wrt2CblFile(outFile,cbl,S)
	strwidth = 13
	pad=4
	println(outFile)
	print(outFile, rpad(round(S, sigdigits=strwidth-pad), strwidth), "\t")
	print(outFile, rpad(round(cbl.length, sigdigits=strwidth-pad), strwidth), "\t")
	print(outFile, rpad(round(cbl.results.ttl, sigdigits=strwidth-pad), strwidth), "\t")
	print(outFile, rpad(round(cbl.results.cbc, sigdigits=strwidth-pad), strwidth), "\t")
	print(outFile, rpad(round(cbl.results.rlc, sigdigits=strwidth-pad), strwidth), "\t")
	print(outFile, rpad(round(cbl.results.qc, sigdigits=strwidth-pad), strwidth), "\t")
	print(outFile, rpad(round(cbl.results.cm, sigdigits=strwidth-pad), strwidth), "\t")
	print(outFile, rpad(round(cbl.results.eens, sigdigits=strwidth-pad), strwidth), "\t")
	print(outFile, rpad(round(cbl.num, sigdigits=strwidth-pad), strwidth), "\t")
	print(outFile, rpad(round(cbl.volt, sigdigits=strwidth-pad), strwidth), "\t")
	print(outFile, rpad(round(cbl.size, sigdigits=strwidth-pad), strwidth), "\t")
	print(outFile, rpad(round(cbl.mva, sigdigits=strwidth-pad), strwidth), "\t")
end
#########################################################
function ppf_printOcn(ocean)
	x=Array{Float64,1}()
	y=Array{Float64,1}()
	for i in ocean.pccs
		push!(x,i.coord.x)
		push!(y,i.coord.y)
		#println(i.num)
	end
	for i in ocean.cnces
		push!(x,i.coord.x)
		push!(y,i.coord.y)
		#println(i.num)
	end
	xb=Array{Float64,1}()
	yb=Array{Float64,1}()
	#println(ocean.bnd)
	for i in ocean.bnd.wbnd.lims
		push!(xb,i.x)
		push!(yb,i.y)
	end
	#println(ocean.bnd)
	for i in ocean.bnd.nbnd.lims
		push!(xb,i.x)
		push!(yb,i.y)
	end
	for i in ocean.bnd.ebnd.lims
		push!(xb,i.x)
		push!(yb,i.y)
	end
	for i in ocean.bnd.sbnd.lims
		push!(xb,i.x)
		push!(yb,i.y)
	end
	xoss=Array{Float64,1}()
	yoss=Array{Float64,1}()
	for i in ocean.osss
		push!(xoss,i.coord.x)
		push!(yoss,i.coord.y)
		#println(i.coord)
	end
	#println(lof_pnt2pnt_dist(ocean.osss[2].coord,ocean.osss[3].coord))
	os=1
	xlimax=trunc(Int,findmax(x)[1])+os
	ylimax=trunc(Int,findmax(y)[1])+os
	xlimin=trunc(Int,findmin(x)[1])-os
	ylimin=trunc(Int,findmin(y)[1])-os
	plotly()

	p=plot(xoss,yoss,seriestype=:scatter,xticks = ylimin:5:ylimax,xlims=(ylimin,ylimax),yticks = ylimin:5:ylimax)
	plot!(p,x,y,seriestype=:scatter,xticks = ylimin:5:ylimax,xlims=(ylimin,ylimax),yticks = ylimin:5:ylimax)
	plot!(p,xb,yb,xticks = ylimin:5:ylimax,xlims=(ylimin,ylimax),yticks = ylimin:5:ylimax)

	xd=Array{Float64,1}()
	yd=Array{Float64,1}()
	for i in ocean.gOarcs
		push!(xd,i.tail.coord.x)
		push!(xd,i.head.coord.x)
		push!(yd,i.tail.coord.y)
		push!(yd,i.head.coord.y)
		plot!(p,xd,yd,xticks = ylimin:5:ylimax,xlims=(ylimin,ylimax),yticks = ylimin:5:ylimax,label="")
		xd=[]
		yd=[]
	end
	for i in ocean.oParcs
		push!(xd,i.tail.coord.x)
		push!(xd,i.head.coord.x)
		push!(yd,i.tail.coord.y)
		push!(yd,i.head.coord.y)
		plot!(p,xd,yd,xticks = ylimin:5:ylimax,xlims=(ylimin,ylimax),yticks = ylimin:5:ylimax,label="")
		xd=[]
		yd=[]
	end
	for i in ocean.oOarcs
		push!(xd,i.tail.coord.x)
		push!(xd,i.head.coord.x)
		push!(yd,i.tail.coord.y)
		push!(yd,i.head.coord.y)
		plot!(p,xd,yd,xticks = ylimin:5:ylimax,xlims=(ylimin,ylimax),yticks = ylimin:5:ylimax,label="")
		xd=[]
		yd=[]
	end#==#
	p
	end
################################################################################
