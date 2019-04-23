##################################################################################################################
################### parse to .m file ################################################################################
##################################################################################################################
function ppf_main2mfile(mp)
	mf=ppf_dotMfile(mp.cnces[1].mva)
	ppf_prntBuss(mf,mp)
	ppf_prntGens(mf,mp)
	ppf_prntBrns(mf,mp)
	ppf_prntNeBrns(mf,mp)
	close(mf)
end
########################################################
function ppf_prntBuss(mf,mp)
	println(mf, "%bus data")
	println(mf, "%	bus_i	type	Pd	Qd	Gs	Bs	area	Vm	Va	baseKV	zone	Vmax	Vmin")
	println(mf, "mpc.bus = [")
	ppf_busNode(mf,mp.pccs)
	ppf_busNode(mf,mp.cnces)
	ppf_busNode(mf,mp.osss)
	Binf=cnce()
	Binf.num=mp.osss[length(mp.osss)].num+1
	Binf.mva=length(mp.cnces)*lod_cnceMva()
	ppf_busInf(mf,Binf)
	println(mf, "];")
	println(mf, "")
end
########################################################
function ppf_prntGens(mf,mp)
	println(mf, "%generator data")
	println(mf, "%	bus	Pg	Qg	Qmax	Qmin	Vg	mbase	status	Pmax	Pmin")
	println(mf, "mpc.gen = [")
	ppf_genNode(mf,mp.cnces)
	println(mf, "];")
	println(mf, "")
	println(mf, "%generator cost data")
	println(mf, "mpc.gencost = [];")
	println(mf, "")
end
########################################################
function ppf_genNode(mf,nds)
	for n in nds
		print(mf,n.num,"\t")
		print(mf,n.mva,"\t")
		print(mf,0.0,"\t")
		print(mf,0.0,"\t")
		print(mf,0.0,"\t")
		print(mf,1.1,"\t")
		print(mf,n.mva,"\t")
		print(mf,1.0,"\t")
		print(mf,n.mva,"\t")
		print(mf,0.0,"\t")
		println(mf,";")
	end
end
########################################################
function ppf_busInf(mf,n)
	print(mf,n.num,"\t")
	print(mf,1.0,"\t")
	print(mf,n.mva,"\t")
	print(mf,0.0,"\t")
	print(mf,0.0,"\t")
	print(mf,1.0,"\t")
	print(mf,1.0,"\t")
	print(mf,1.0,"\t")
	print(mf,1.05,"\t")
	print(mf,lod_pccKv(),"\t")
	print(mf,1.0,"\t")
	print(mf,1.1,"\t")
	print(mf,0.9)
	println(mf,";")
end
########################################################
function ppf_busNode(mf,nds)
	if typeof(nds[1]) == typeof(cnce())
		tp=3.0
	elseif typeof(nds[1]) == typeof(pcc())
		tp=1.0
	elseif typeof(nds[1]) == typeof(oss())
		tp=1.0
	else
		println("No type match for node!")
	end
	for n in nds
		print(mf,n.num,"\t")
		print(mf,tp,"\t")
		print(mf,0.0,"\t")
		print(mf,0.0,"\t")
		print(mf,0.0,"\t")
		print(mf,1.0,"\t")
		print(mf,1.0,"\t")
		print(mf,1.0,"\t")
		print(mf,1.05,"\t")
		print(mf,lod_pccKv(),"\t")
		print(mf,1.0,"\t")
		print(mf,1.1,"\t")
		print(mf,0.9)
		println(mf,";")
		if typeof(nds[1]) == typeof(cnce())
			tp=2.0
		else
		end
	end
end
########################################################
function ppf_prntBrn(mf,a,cb)
	print(mf,a.tail.num,"\t")
	print(mf,a.head.num,"\t")
	print(mf,cb.ohm,"\t")
	print(mf,cb.xl,"\t")
	print(mf,cb.yc,"\t")
	print(mf,cb.mva*cb.num,"\t")
	print(mf,cb.mva*cb.num,"\t")
	print(mf,cb.mva*cb.num,"\t")
	print(mf,0.0,"\t")
	print(mf,0.0,"\t")
	print(mf,1.0,"\t")
	print(mf,-30.0,"\t")
	print(mf,30.0,"\t")
end
########################################################
function ppf_ooBrn(mf,as)

	for a in as
		wp=wndF_wndPrf([a.head.wnds[1]])
		cb=cstF_cbl_ttl(a.lngth,a.head.mvas[1]/2,lod_ossKv(),wp,false)
		cb.ohm=((cb.ohm*a.lngth)/cb.num)/eqpD_pu()[4]
		cb.xl=((cb.xl*a.lngth)/cb.num)/eqpD_pu()[4]
		cb.yc=cb.yc*a.lngth*cb.num*eqpD_pu()[4]
		ppf_prntBrn(mf,a,cb)
		print(mf,cb.results.ttl+2*cstD_cfs().FC_ac)
		println(mf,";")

		#=cb=cstF_cbl_ttl(a.lngth,a.head.mvas[1]/3,lod_ossKv(),wp,false)
		cb.ohm=((cb.ohm*a.lngth)/cb.num)/eqpD_pu()[4]
		cb.xl=((cb.xl*a.lngth)/cb.num)/eqpD_pu()[4]
		cb.yc=cb.yc*a.lngth*cb.num*eqpD_pu()[4]
		ppf_prntBrn(mf,a,cb)
		print(mf,cb.results.ttl+2*cstD_cfs().FC_ac)
		println(mf,";")=#

		mva=0.0
		ka=Array{Tuple,1}()
		for j=1:length(a.head.mvas)
			mva=mva+a.head.mvas[j]
			push!(ka,a.head.wnds[j])
			wp=wndF_wndPrf(ka)
			cb=cstF_cbl_ttl(a.lngth,mva,lod_ossKv(),wp,false)
		    cb.ohm=((cb.ohm*a.lngth)/cb.num)/eqpD_pu()[4]
		    cb.xl=((cb.xl*a.lngth)/cb.num)/eqpD_pu()[4]
		    cb.yc=cb.yc*a.lngth*cb.num*eqpD_pu()[4]
			ppf_prntBrn(mf,a,cb)
			print(mf,cb.results.ttl+2*cstD_cfs().FC_ac)
			println(mf,";")
		end
		ka=[]
	end
end
########################################################
function ppf_prntBrns(mf,mp)
	println(mf, "%branch data")
	println(mf, "%	fbus	tbus	r	x	b	rateA	rateB	rateC	ratio	angle	status	angmin	angmax")
	println(mf, "mpc.branch = [")
	println(mf, "];")
	println(mf, "")
end
########################################################
function ppf_prntNeBrns(mf,mp)
	println(mf, "%candidate branch data")
	println(mf, "%column_names%	f_bus	t_bus	br_r	br_x	br_b	rate_a	rate_b	rate_c	tap	shift	br_status	angmin	angmax	construction_cost")
	println(mf, "mpc.ne_branch = [")
	#ppf_goBrn(mf,mp.gOarcs)
	ppf_ooBrn(mf,mp.oOarcs)
	#ppf_opBrn(mf,mp.oParcs)
	#ppf_gpBrn(mf,mp.gParcs)
	println(mf, "];")
end
########################################################
function ppf_dotMfile(S_pu)
	fileName="results/owpp_tnep_map.m"
	mfile = open(fileName,"w")
	ppf_cblTtle2m(mfile,S_pu)
    return mfile
end
########################################################
function ppf_addGens()

end
#########################################################
function ppf_cblTtle2m(mf,S_pu)
	# Print headers
	println(mf, "%TNEP optimizor input file test")
	println(mf, "")
	println(mf, "function mpc = owpp_tnep_map")
	println(mf, "mpc.version = '2';")
	print(mf, "mpc.baseMVA = ", "\t")
	print(mf, S_pu)
	println(mf, ";")
	println(mf, "")
end
#########################################################

##################################################################################################################
################### printing ################################################################################
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
		#println(lof_pnt2pnt_dist(ocean.pccs[2].coord,i.coord))
		#println(lof_pnt2pnt_dist(ocean.pccs[1].coord,i.coord))
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
	#p=plot(xoss,yoss,seriestype=:scatter)
	#plot!(p,x,y,seriestype=:scatter)

	#plot!(p,xb,yb,xticks = ylimin:5:ylimax,xlims=(ylimin,ylimax),yticks = ylimin:5:ylimax)

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
	#=for i in ocean.oParcs
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
	end
	for i in ocean.gParcs
		push!(xd,i.tail.coord.x)
		push!(xd,i.head.coord.x)
		push!(yd,i.tail.coord.y)
		push!(yd,i.head.coord.y)
		plot!(p,xd,yd,xticks = ylimin:5:ylimax,xlims=(ylimin,ylimax),yticks = ylimin:5:ylimax,label="")
		xd=[]
		yd=[]
	end=#
	p
	end
################################################################################
