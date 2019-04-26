################################################################################################################
################################## Print to optimization .m file ###############################################
################################################################################################################
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
	println(mf, "%bus_i	type	Pd	Qd	Gs	Bs	area	Vm	Va	baseKV	zone	Vmax	Vmin")
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
	println(mf, "%bus	Pg	Qg	Qmax	Qmin	Vg	mbase	status	Pmax	Pmin")
	println(mf, "mpc.gen = [")
	ppf_genNode(mf,mp.cnces)
	println(mf, "];")
	println(mf, "")
	println(mf, "%generator cost data")
	println(mf, "mpc.gencost = [")
	for i=1:length(mp.cnces)
		println(mf, "2	 0.0	 0.0	 0	   0	   0	   0;")
	end
	println(mf, "];")
	println(mf, "")
end
########################################################
function ppf_genNode(mf,nds)
	for n in nds
		print(mf,n.num,"\t")
		print(mf,trunc(Int,n.mva),"\t")
		print(mf,0.0,"\t")
		print(mf,0.0,"\t")
		print(mf,0.0,"\t")
		print(mf,1.0,"\t")
		print(mf,trunc(Int,n.mva),"\t")
		print(mf,1.0,"\t")
		print(mf,trunc(Int,n.mva),"\t")
		print(mf,0.0)
		println(mf,";")
	end
end
########################################################
function ppf_busInf(mf,n)
	print(mf,n.num,"\t")
	print(mf,1.0,"\t")
	print(mf,trunc(Int,n.mva),"\t")
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

################################################################################################################
################################ Arcs ##########################################################################
################################################################################################################
function ppf_prntNeBrns(mf,mp)
	println(mf, "%candidate branch data")
	println(mf, "%column_names%	f_bus	t_bus	br_r	br_x	br_b	rate_a	rate_b	rate_c	tap	shift	br_status	angmin	angmax	construction_cost")
	println(mf, "mpc.ne_branch = [")
	ppf_gpBrn(mf,mp.gParcs,mp.asBuilt)
	ppf_goBrn(mf,mp.gOarcs,mp.asBuilt)
	ppf_ooBrn(mf,mp.oOarcs,mp.asBuilt)
	ppf_opBrn(mf,mp.oParcs,mp.asBuilt)
	println(mf, "];")
end
################################################################################
#converts transformer and cable values to PU then calls .m printing
function ppf_xcXrlPu(a,arc,mf)
	eqpf_puImped(arc.cable,a.lngth)
	xxr=(arc.xfm.ohm/arc.xfm.num)/eqpD_pu()[4]
	xxl=(arc.xfm.xl/arc.xfm.num)/eqpD_pu()[4]
	ppf_prntBrn(mf,a,arc.cable,xxr,xxl)
end
################################################################################
function ppf_asBuilt(a,abs,mva,arc,cst,km)
	ab=lnks()
	ab.hd=deepcopy(a.head)
	ab.tl=deepcopy(a.tail)
	#arc.mva=deepcopy(mva)
	#arc.km=deepcopy(km)
	#arc.costs.ttl=deepcopy(cst)
	ab.elec=deepcopy(arc)
	ab.elec.costs.ttl=cst
	ab.elec.km=deepcopy(km)
	ab.elec.mva=deepcopy(mva)
	push!(abs,ab)
end
################################################################################
function ppf_gpBrn(mf,as,abs)
	mvas=Array{Float64,1}()
	for a in as
		#push!(mvas,a.tail.mva/3)
		push!(mvas,a.tail.mva/2)
		push!(mvas,a.tail.mva)
		wp=wndF_wndPrf([a.tail.wnd],a.tail.trb)

		for mva in mvas
			arc=cstF_cblWT_ttl(a.lngth,mva,a.tail.kv,wp,false)
		    eqpD_xfoXR(a.tail.kv,arc.xfm)
			if arc.cable.num != 0
				ppf_xcXrlPu(a,arc,mf)
				cst=arc.costs.ttl+2*cstD_cfs().FC_ac
				println(cst)
				ppf_asBuilt(a,abs,mva,arc,cst,a.lngth)
				print(mf,cst)
				println(mf,";")
			else
				println("No suitable "*string(a.tail.kv)*"Kv, "*string(mva)*"MVA cable for "*string(trunc(Int,a.lngth))*"Km. -removing")
			end
		end
		mvas=[]
	end
end
################################################################################
function ppf_goBrn(mf,as,abs)
	mvas=Array{Float64,1}()
	#println("length of cable: ")
	#println(length(as))
	for a in as
		wp=wndF_wndPrf([a.tail.wnd],a.tail.trb)
		#push!(mvas,a.tail.mva/3)
		push!(mvas,a.tail.mva/2)
		push!(mvas,a.tail.mva)
		wp=wndF_wndPrf([a.tail.wnd],a.tail.trb)

		for mva in mvas
			arc=cstF_cblWT_ttl(a.lngth,mva,a.tail.kv,wp,true)
		    eqpD_xfoXR(a.tail.kv,arc.xfm)
			if arc.cable.num != 0
				ppf_xcXrlPu(a,arc,mf)
				cst=arc.costs.ttl+(-1*cstD_cfs().FC_ac)
				ppf_asBuilt(a,abs,mva,arc,cst,a.lngth)
				print(mf,cst)
				println(mf,";")
			else
				println("No suitable "*string(a.tail.kv)*"Kv, "*string(mva)*"MVA cable for "*string(trunc(Int,a.lngth))*"Km. -removing")
			end
		end
		mvas=[]
	end
end
################################################################################
function ppf_ooBrn(mf,as,abs)
#TODO: turb() should be updated to carry concession turbines
	mvas=Array{Float64,1}()
	arc=owpp()
	for a in as
		wp=wndF_wndPrf([a.tail.wnds[1]],turb())
		#push!(mvas,a.tail.mvas[1]/3))
		push!(mvas,a.tail.mvas[1]/2)
		for mva in mvas
			cb=cstF_cbl_ttl(a.lngth,mva,lod_ossKv(),wp,true)
			if cb.num != 0
				eqpf_puImped(cb,a.lngth)
				ppf_prntBrn(mf,a,cb,0.0,0.0)
				arc.cable=cb
				cst=cb.results.ttl+2*cstD_cfs().FC_ac
				ppf_asBuilt(a,abs,mva,arc,cst,a.lngth)
				print(mf,cst)
				println(mf,";")
			else
				println("No suitable "*string(a.tail.kv)*"Kv, "*string(mva)*"MVA cable for "*string(trunc(Int,a.lngth))*"Km. -removing")
			end
		end
		mvas=[]
		mva=0.0
		ka=Array{Tuple,1}()
		for j=1:length(a.tail.mvas)
			mva=mva+a.tail.mvas[j]
			push!(ka,a.tail.wnds[j])
			wp=wndF_wndPrf(ka,turb())
			#println(a.lngth)
			cb=cstF_cbl_ttl(a.lngth,mva,lod_ossKv(),wp,true)
			if cb.num != 0
			    eqpf_puImped(cb,a.lngth)
				ppf_prntBrn(mf,a,cb,0.0,0.0)
				arc.cable=cb
				cst=cb.results.ttl+2*cstD_cfs().FC_ac
				ppf_asBuilt(a,abs,mva,arc,cst,a.lngth)
				print(mf,cst)
				println(mf,";")
			else
				println("No suitable "*string(lod_ossKv())*"Kv, "*string(mva)*"MVA cable for "*string(trunc(Int,a.lngth))*"Km. -removing")
			end
		end
		ka=[]
	end
end
################################################################################
function ppf_opBrn(mf,as,abs)
	#TODO: turb() should be updated to carry concession turbines
		mvas=Array{Float64,1}()
		for a in as
			#println(a)
			wp=wndF_wndPrf([a.tail.wnds[1]],turb())
			#push!(mvas,a.tail.mvas[1]/3))
			push!(mvas,a.tail.mvas[1]/2)
			for mva in mvas
				arc=cstF_cblWT_ttl(a.lngth,mva,lod_ossKv(),wp,false)
				eqpD_xfoXR(lod_ossKv(),arc.xfm)
				if arc.cable.num != 0
					ppf_xcXrlPu(a,arc,mf)
					cst=arc.costs.ttl+cstD_cfs().FC_ac
					ppf_asBuilt(a,abs,mva,arc,cst,a.lngth)
					print(mf,cst)
					println(mf,";")
				else
					println("No suitable "*string(lod_ossKv())*"Kv, "*string(mva)*"MVA cable for "*string(trunc(Int,a.lngth))*"Km. -removing")
				end
			end
			mvas=[]
			mva=0.0
			ka=Array{Tuple,1}()
			for j=1:length(a.tail.mvas)
				mva=mva+a.tail.mvas[j]
				push!(ka,a.tail.wnds[j])
				wp=wndF_wndPrf(ka,turb())
				#println(a.lngth)
				arc=cstF_cblWT_ttl(a.lngth,mva,lod_ossKv(),wp,false)
				eqpD_xfoXR(lod_ossKv(),arc.xfm)
				if arc.cable.num != 0
				    ppf_xcXrlPu(a,arc,mf)
					cst=arc.costs.ttl+cstD_cfs().FC_ac
					ppf_asBuilt(a,abs,mva,arc,cst,a.lngth)
					print(mf,cst)
					println(mf,";")
				else
					println("No suitable "*string(lod_ossKv())*"Kv, "*string(mva)*"MVA cable for "*string(trunc(Int,a.lngth))*"Km. -removing")
				end
			end
			ka=[]
		end
	end
########################################################
function ppf_prntBrn(mf,a,cb,xr,xxl)
	print(mf,a.tail.num,"\t")
	print(mf,a.head.num,"\t")
	print(mf,cb.ohm+xr,"\t")
	print(mf,cb.xl+xxl,"\t")
	print(mf,cb.yc,"\t")
	print(mf,trunc(Int,cb.mva*cb.num),"\t")
	print(mf,trunc(Int,cb.mva*cb.num),"\t")
	print(mf,trunc(Int,cb.mva*cb.num),"\t")
	print(mf,0.0,"\t")
	print(mf,0.0,"\t")
	print(mf,1.0,"\t")
	print(mf,-30.0,"\t")
	print(mf,30.0,"\t")
end
########################################################
function ppf_prntBrnOnShre(mf,fn,tn,p)
	print(mf,fn,"\t")
	print(mf,tn,"\t")
	print(mf,0.0093,"\t")
	print(mf,0.0222,"\t")
	print(mf,0.2217,"\t")
	print(mf,p,"\t")
	print(mf,p,"\t")
	print(mf,p,"\t")
	print(mf,0.0,"\t")
	print(mf,0.0,"\t")
	print(mf,1.0,"\t")
	print(mf,-30.0,"\t")
	print(mf,30.0)
	println(mf,";")
end
########################################################
function ppf_prntBrns(mf,mp)
	println(mf, "%branch data")
	println(mf, "%fbus	tbus	r	x	b	rateA	rateB	rateC	ratio	angle	status	angmin	angmax")
	println(mf, "mpc.branch = [")
	mva=0.0
	tn=length(mp.pccs)+length(mp.cnces)+length(mp.osss)+1
	for i in mp.cnces
		mva=mva+i.mva
	end
	for i in mp.pccs
		ppf_prntBrnOnShre(mf,i.num,tn,mva)
	end
	println(mf, "];")
	println(mf, "")
end
################################################################################################################
################################ Nodes #########################################################################
################################################################################################################
function ppf_dotMfile(S_pu)
	fileName="results/owpp_tnep_map.mat"
	mfile = open(fileName,"w")
	ppf_cblTtle2m(mfile,S_pu)
    return mfile
end
########################################################
function ppf_cblTtle2m(mf,S_pu)
	# Print headers
	println(mf, "%TNEP optimizor input file test")
	println(mf, "")
	println(mf, "function mpc = owpp_tnep_map")
	println(mf, "mpc.version = '2';")
	print(mf, "mpc.baseMVA = ")
	print(mf, S_pu)
	println(mf, ";")
	println(mf, "")
end
#########################################################

##################################################################################################################
######################## printing ################################################################################
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
	anoss=Array{Int64,1}()
	xoss=Array{Float64,1}()
	yoss=Array{Float64,1}()
	for i in ocean.osss
		push!(xoss,i.coord.x)
		push!(yoss,i.coord.y)
		push!(anoss,i.num)
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
	end
	for i in ocean.gParcs
		push!(xd,i.tail.coord.x)
		push!(xd,i.head.coord.x)
		push!(yd,i.tail.coord.y)
		push!(yd,i.head.coord.y)
		plot!(p,xd,yd,xticks = ylimin:5:ylimax,xlims=(ylimin,ylimax),yticks = ylimin:5:ylimax,label="")
		xd=[]
		yd=[]
	end
	p
	end
################################################################################
