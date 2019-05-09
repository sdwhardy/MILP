##################################################################################################################
######################## printing to screen ######################################################################
##################################################################################################################
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
	xpcc=Array{Float64,1}()
	ypcc=Array{Float64,1}()
	pcc=Array{Tuple,1}()
	for i in ocean.pccs
		push!(xpcc,i.coord.x)
		push!(ypcc,i.coord.y)
		tx=text(string(i.num),12,:black,:right)
		push!(pcc,(i.coord.x,i.coord.y,tx))
		#println(i.num)
	end
	xgen=Array{Float64,1}()
	ygen=Array{Float64,1}()
	gen=Array{Tuple,1}()
	for i in ocean.cnces
		push!(xgen,i.coord.x)
		push!(ygen,i.coord.y)
		tx=text(string(i.num),12,:red,:right)
		push!(gen,(i.coord.x,i.coord.y,tx))
		#println(lof_pnt2pnt_dist(ocean.pccs[2].coord,i.coord))
		#println(lof_pnt2pnt_dist(ocean.pccs[1].coord,i.coord))
	end


	#=xb=Array{Float64,1}()
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
	end=#
	#anoss=Array{Int64,1}()
	xoss=Array{Float64,1}()
	yoss=Array{Float64,1}()
	oss=Array{Tuple,1}()
	for i in ocean.osss
		push!(xoss,i.coord.x)
		push!(yoss,i.coord.y)
		tx=text(string(i.num),12,:blue,:left)
		push!(oss,(i.coord.x,i.coord.y,tx))
		#push!(anoss,(i.num))
		#println(i.coord)
	end
	#println(lof_pnt2pnt_dist(ocean.osss[2].coord,ocean.osss[3].coord))
	os=1
	xlimax=trunc(Int,findmax(findmax([xoss,xpcc,xgen])[1])[1])+os
	ylimax=trunc(Int,findmax(findmax([yoss,ypcc,ygen])[1])[1])+os
	xlimin=trunc(Int,findmin(findmin([xoss,xpcc,xgen])[1])[1])-os
	ylimin=trunc(Int,findmin(findmin([yoss,ypcc,ygen])[1])[1])-os
	plotly()



	p=plot(xoss,yoss,annotations=oss,color = :blue,seriestype=:scatter,xticks = ylimin:1:ylimax,xlims=(ylimin,ylimax),yticks = ylimin:5:ylimax,label="")
	plot!(p,xpcc,ypcc,annotations=pcc,color = :black,seriestype=:scatter,xticks = ylimin:1:ylimax,xlims=(ylimin,ylimax),yticks = ylimin:5:ylimax,label="")
	#annotate!(oss)
	plot!(p,xgen,ygen,annotations=gen,color = :red,seriestype=:scatter,xticks = ylimin:1:ylimax,xlims=(ylimin,ylimax),yticks = ylimin:5:ylimax,label="")
	#p=plot(xoss,yoss,seriestype=:scatter)
	#plot!(p,x,y,seriestype=:scatter)

	#plot!(p,xb,yb,xticks = ylimin:5:ylimax,xlims=(ylimin,ylimax),yticks = ylimin:5:ylimax)

	xd=Array{Float64,1}()
	yd=Array{Float64,1}()
	op=Array{Tuple,1}()
	#=for i in ocean.oParcs
		push!(xd,i.tail.coord.x)
		push!(xd,i.head.coord.x)
		push!(yd,i.tail.coord.y)
		push!(yd,i.head.coord.y)
		plot!(p,xd,yd,color = :red,xticks = ylimin:1:ylimax,xlims=(ylimin,ylimax),yticks = ylimin:5:ylimax,label="")
		xd=[]
		yd=[]
	end
	for i in ocean.oOarcs
		push!(xd,i.tail.coord.x)
		push!(xd,i.head.coord.x)
		push!(yd,i.tail.coord.y)
		push!(yd,i.head.coord.y)
		plot!(p,xd,yd,color = :black,xticks = ylimin:1:ylimax,xlims=(ylimin,ylimax),yticks = ylimin:5:ylimax,label="")
		xd=[]
		yd=[]
	end
	for i in ocean.gOarcs
		push!(xd,i.tail.coord.x)
		push!(xd,i.head.coord.x)
		push!(yd,i.tail.coord.y)
		push!(yd,i.head.coord.y)
		plot!(p,xd,yd,color = :blue,xticks = ylimin:1:ylimax,xlims=(ylimin,ylimax),yticks = ylimin:5:ylimax,label="")
		xd=[]
		yd=[]
	end
	for i in ocean.gParcs
		push!(xd,i.tail.coord.x)
		push!(xd,i.head.coord.x)
		push!(yd,i.tail.coord.y)
		push!(yd,i.head.coord.y)
		plot!(p,xd,yd,color = :blue,xticks = ylimin:1:ylimax,xlims=(ylimin,ylimax),yticks = ylimin:5:ylimax,label="")
		xd=[]
		yd=[]
	end=#
	p
	end
################################################################################
