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
		push!(x,i.coord.cnt.x)
		push!(y,i.coord.cnt.y)
	end
	for i in ocean.reg.cnces
		push!(x,i.coord.cnt.x)
		push!(y,i.coord.cnt.y)
	end
	xb=Array{Float64,1}()
	yb=Array{Float64,1}()
	#println(ocean.reg.bnd)
	for i in ocean.reg.bnd
		push!(xb,i.x)
		push!(yb,i.y)
	end
	
	#=xc=Array{Float64,1}()
	yc=Array{Float64,1}()
	for i in ocean.reg.sth
		push!(xc,i.x)
		push!(yc,i.y)
	end=#
    #x=[ocean.pccs[1].coord.cnt.x,ocean.pccs[2].coord.cnt.x,ocean.reg.cnces[1].coord.cnt.x,ocean.reg.cnces[2].coord.cnt.x]#,ocean.reg.cnces[3].coord.cnt.x,ocean.reg.cnces[4].coord.cnt.x,ocean.reg.cnces[5].coord.cnt.x,ocean.reg.cnces[6].coord.cnt.x,ocean.reg.cnces[7].coord.cnt.x,ocean.reg.cnces[8].coord.cnt.x]
    #y=[ocean.pccs[1].coord.cnt.y,ocean.pccs[2].coord.cnt.y,ocean.reg.cnces[1].coord.cnt.y,ocean.reg.cnces[2].coord.cnt.y]#,ocean.reg.cnces[3].coord.cnt.y,ocean.reg.cnces[4].coord.cnt.y,ocean.reg.cnces[5].coord.cnt.y,ocean.reg.cnces[6].coord.cnt.y,ocean.reg.cnces[7].coord.cnt.y,ocean.reg.cnces[8].coord.cnt.y]
	#println(findmax(y)[1])
	#println(findmax(y)[2])
	os=1
	xlimax=trunc(Int,findmax(x)[1])+os
	ylimax=trunc(Int,findmax(y)[1])+os
	xlimin=trunc(Int,findmin(x)[1])-os
	ylimin=trunc(Int,findmin(y)[1])-os
	plotly()
	plot(x,y,seriestype=:scatter,xticks = xlimin:5:xlimax,xlims=(xlimin,xlimax),yticks = ylimin:5:ylimax)
	plot!(xb,yb,xticks = xlimin:5:xlimax,xlims=(xlimin,xlimax),yticks = ylimin:5:ylimax)
	#plot!(xc,yc,xticks = xlimin:5:xlimax,xlims=(xlimin,xlimax),yticks = ylimin:5:ylimax)

    #=i=1
    x1=[ocean.reg.cnces[i].coord.nw.x,ocean.reg.cnces[i].coord.ne.x,ocean.reg.cnces[i].coord.se.x,ocean.reg.cnces[i].coord.sw.x,ocean.reg.cnces[i].coord.nw.x]
    y1=[ocean.reg.cnces[i].coord.nw.y,ocean.reg.cnces[i].coord.ne.y,ocean.reg.cnces[i].coord.se.y,ocean.reg.cnces[i].coord.sw.y,ocean.reg.cnces[i].coord.nw.y]
    i=2
    x2=[ocean.reg.cnces[i].coord.nw.x,ocean.reg.cnces[i].coord.ne.x,ocean.reg.cnces[i].coord.se.x,ocean.reg.cnces[i].coord.sw.x,ocean.reg.cnces[i].coord.nw.x]
    y2=[ocean.reg.cnces[i].coord.nw.y,ocean.reg.cnces[i].coord.ne.y,ocean.reg.cnces[i].coord.se.y,ocean.reg.cnces[i].coord.sw.y,ocean.reg.cnces[i].coord.nw.y]
	i=3
    x3=[ocean.reg.cnces[i].coord.nw.x,ocean.reg.cnces[i].coord.ne.x,ocean.reg.cnces[i].coord.se.x,ocean.reg.cnces[i].coord.sw.x,ocean.reg.cnces[i].coord.nw.x]
    y3=[ocean.reg.cnces[i].coord.nw.y,ocean.reg.cnces[i].coord.ne.y,ocean.reg.cnces[i].coord.se.y,ocean.reg.cnces[i].coord.sw.y,ocean.reg.cnces[i].coord.nw.y]
    i=4
    x4=[ocean.reg.cnces[i].coord.nw.x,ocean.reg.cnces[i].coord.ne.x,ocean.reg.cnces[i].coord.se.x,ocean.reg.cnces[i].coord.sw.x,ocean.reg.cnces[i].coord.nw.x]
    y4=[ocean.reg.cnces[i].coord.nw.y,ocean.reg.cnces[i].coord.ne.y,ocean.reg.cnces[i].coord.se.y,ocean.reg.cnces[i].coord.sw.y,ocean.reg.cnces[i].coord.nw.y]
    i=5
    x5=[ocean.reg.cnces[i].coord.nw.x,ocean.reg.cnces[i].coord.ne.x,ocean.reg.cnces[i].coord.se.x,ocean.reg.cnces[i].coord.sw.x,ocean.reg.cnces[i].coord.nw.x]
    y5=[ocean.reg.cnces[i].coord.nw.y,ocean.reg.cnces[i].coord.ne.y,ocean.reg.cnces[i].coord.se.y,ocean.reg.cnces[i].coord.sw.y,ocean.reg.cnces[i].coord.nw.y]
    i=6
    x6=[ocean.reg.cnces[i].coord.nw.x,ocean.reg.cnces[i].coord.ne.x,ocean.reg.cnces[i].coord.se.x,ocean.reg.cnces[i].coord.sw.x,ocean.reg.cnces[i].coord.nw.x]
    y6=[ocean.reg.cnces[i].coord.nw.y,ocean.reg.cnces[i].coord.ne.y,ocean.reg.cnces[i].coord.se.y,ocean.reg.cnces[i].coord.sw.y,ocean.reg.cnces[i].coord.nw.y]
    i=7
    x7=[ocean.reg.cnces[i].coord.nw.x,ocean.reg.cnces[i].coord.ne.x,ocean.reg.cnces[i].coord.se.x,ocean.reg.cnces[i].coord.sw.x,ocean.reg.cnces[i].coord.nw.x]
    y7=[ocean.reg.cnces[i].coord.nw.y,ocean.reg.cnces[i].coord.ne.y,ocean.reg.cnces[i].coord.se.y,ocean.reg.cnces[i].coord.sw.y,ocean.reg.cnces[i].coord.nw.y]
    i=8
    x8=[ocean.reg.cnces[i].coord.nw.x,ocean.reg.cnces[i].coord.ne.x,ocean.reg.cnces[i].coord.se.x,ocean.reg.cnces[i].coord.sw.x,ocean.reg.cnces[i].coord.nw.x]
    y8=[ocean.reg.cnces[i].coord.nw.y,ocean.reg.cnces[i].coord.ne.y,ocean.reg.cnces[i].coord.se.y,ocean.reg.cnces[i].coord.sw.y,ocean.reg.cnces[i].coord.nw.y]=#
    #plot!(x1,y1,xticks = 0:5:xlim,xlims=(0,xlim),yticks = 0:5:ylim,ylims=(0,ylim))
    #plot!(x2,y2,xticks = 0:5:xlim,xlims=(0,xlim),yticks = 0:5:ylim,ylims=(0,ylim))

    #=plot!(x3,y3,xticks = 0:5:20,xlims=(0,20),yticks = 0:5:55)
    plot!(x4,y4,xticks = 0:5:20,xlims=(0,20),yticks = 0:5:55)
    plot!(x5,y5,xticks = 0:5:20,xlims=(0,20),yticks = 0:5:55)
    plot!(x6,y6,xticks = 0:5:20,xlims=(0,20),yticks = 0:5:55)
    plot!(x7,y7,xticks = 0:5:20,xlims=(0,20),yticks = 0:5:55)
    plot!(x8,y8,xticks = 0:5:20,xlims=(0,20),yticks = 0:5:55)=#
end
################################################################################
