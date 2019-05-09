################################################################################
######################## reconstruction of optimal solution ####################
################################################################################
function tpp_1build(eyed,solId)
	yeyne=0.0
	for i=1:1:length(solId)
		if string(eyed) == string(solId[i])
			yeyne=1.0
			deleteat!(solId,i)
			@goto end_1build
		end
	end
	@label end_1build
	return yeyne
end
################################################################################
#compares solution to original data file to reconstruct OSS and cables
function tpp_reCnstrctSol(res,nt,ocn)
    s=tpp_xtrctBrch(res,nt)
    #fileName="results/owpp_tnep_map.m"
	#mfile = open(fileName)
	idd, syAsBuilt=tpp_xtrctOwpp(s,ocn.asBuilt)
	#close(mfile)
	return idd, syAsBuilt
end
################################################################################
function tpp_xtrctOwpp(s,ab)
	idd=Array{Int64,1}()
	tmpSys=eez()
	for i in s
		for j in ab
			if (j.tl.num==i["f_bus"] && j.hd.num==i["t_bus"] && j.elec.costs.ttl==i["construction_cost"])
				push!(idd,i["brnID"])
				if typeof(j.tl) == typeof(oss())
					push!(tmpSys.osss,j.tl)
					if typeof(j.hd) == typeof(oss())
						push!(tmpSys.osss,j.hd)
						a=oOarc()
						a.tail=j.tl
						a.head=j.hd
						push!(tmpSys.oOarcs,a)
					elseif typeof(j.hd) == typeof(pcc())
						push!(tmpSys.pccs,j.hd)
						a=oParc()
						a.tail=j.tl
						a.head=j.hd
						push!(tmpSys.oParcs,a)
					else
					end
				elseif typeof(j.tl) == typeof(cnce())
					push!(tmpSys.cnces,j.tl)
					if typeof(j.hd) == typeof(oss())
						push!(tmpSys.osss,j.hd)
						a=gOarc()
						a.tail=j.tl
						a.head=j.hd
						push!(tmpSys.gOarcs,a)
					elseif typeof(j.hd) == typeof(pcc())
						push!(tmpSys.pccs,j.hd)
						a=gParc()
						a.tail=j.tl
						a.head=j.hd
						push!(tmpSys.gParcs,a)
					else
					end
				else
				end
			end
		end
	end
	return idd, tmpSys
end
################################################################################
function tpp_xtrctBrch(r,n)
        sol=Array{Dict,1}()
		for (b, branch) in r["solution"]["ne_branch"]
		        if isapprox(branch["built"], 1.0, atol = 0.01) == 1
                        push!(sol,n["ne_branch"]["$b"])
                end
        end
		#println(sol)
	return sol
end
################################################################################
function tpp_prnt2Scrn(r,nt)
	for (b, branch) in r["solution"]["ne_branch"]
	        if isapprox(branch["built"], 1.0, atol = 0.01) == 1
				print(b)
				print(" - ")
				print(nt["ne_branch"]["$b"]["rate_a"]*lod_cnceMva())
				print("MVA: ")
				print(" Nodes:")
	            print(nt["ne_branch"]["$b"]["f_bus"])
	            print(" - ")
	            print(nt["ne_branch"]["$b"]["t_bus"])
	            print(" Cost:")
	            print(nt["ne_branch"]["$b"]["construction_cost"])
				print(" Built:")
				println(branch["built"])
	        end
	end
end
################################################################################
################################################################################################################
################################## opens and titles .m file ####################################################
################################################################################################################
#main function to create the .m file
function tpp_main2mfile(mp,s,ob)
	mf=tpp_dotMfile()#open the .m file
	tpp_cblTtle2m(mf,mp.cnces[1].mva)#print top function data
	tpp_prntBuss(mf,mp)#prints the bus data
	tpp_prntGens(mf,mp)#prints all generator (OWPP) data
	tpp_prntBrns(mf,mp,ob)#prints any pre-existing branches (onshore connections)
	tpp_prntNeBrns(mf,mp,s)#prints all candiadate branch data
	close(mf)
end
########################################################
#opens the .m file
function tpp_dotMfile()
	fileName="results/owpp_tnep_map.mat"
	mfile = open(fileName,"w")
    return mfile
end
########################################################
# file description and function headers
function tpp_cblTtle2m(mf,S_pu)
	println(mf, "%TNEP optimizor input file test")
	println(mf, "")
	println(mf, "function mpc = owpp_tnep_map")
	println(mf, "mpc.version = '2';")
	print(mf, "mpc.baseMVA = ")
	print(mf, S_pu)
	println(mf, ";")
	println(mf, "")
end
################################################################################

################################################################################################################
################################## Prints Node data ############################################################
################################################################################################################
#prints bus titles and data for pccs, concessions, oss and onshore star point
function tpp_prntBuss(mf,mp)
	println(mf, "%bus data")
	println(mf, "%bus_i	type	Pd	Qd	Gs	Bs	area	Vm	Va	baseKV	zone	Vmax	Vmin id_bus")
	println(mf, "mpc.bus = [")
	tpp_busNode(mf,mp.pccs)#prints all pccs to the bus
	tpp_busNode(mf,mp.cnces)#prints all concessions to the bus
	tpp_busNode(mf,mp.osss)#prints all oss to the bus
	Binf=cnce()#creates an onshore connection point to connect all PCCS
	Binf.num=mp.osss[length(mp.osss)].num+1#gets the next available bus number
	Binf.mva=length(mp.cnces)*lod_cnceMva()#adds a load to the bus equal to all generation
	tpp_busInf(mf,Binf)#builds the infinite onshore bus and load
	println(mf, "];")
	println(mf, "")
end
########################################################
#builds infinite bus onshore data
function tpp_busInf(mf,n)
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
	print(mf,0.9,"\t")
	print(mf,"3"*string(n.num))
	println(mf,";")
end
########################################################
#general node printing
function tpp_busNode(mf,nds)
	if typeof(nds[1]) == typeof(cnce())
		tp=3.0#designates the first generator the reference node
	elseif (typeof(nds[1]) == typeof(pcc()) || typeof(nds[1]) == typeof(oss()))
		tp=1.0#designates the pccs pq nodes
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
		print(mf,0.9,"\t")
		print(mf,n.id)
		println(mf,";")
		if typeof(nds[1]) == typeof(cnce())
			tp=2.0#designates all other generators pv nodes
		else
		end
	end
end
########################################################
#prints all the generators (OWPP) titles/data
function tpp_prntGens(mf,mp)
	println(mf, "%generator data")
	println(mf, "%bus	Pg	Qg	Qmax	Qmin	Vg	mbase	status	Pmax	Pmin")
	println(mf, "mpc.gen = [")
	tpp_genNode(mf,mp.cnces)#the concessions are the only sources of generation
	println(mf, "];")
	println(mf, "")
	println(mf, "%generator cost data")
	println(mf, "mpc.gencost = [")
	for i=1:length(mp.cnces)
		println(mf, "2	 0.0	 0.0	 0	   0	   0	   0;")#Adds 0 cost function for each generator
	end
	println(mf, "];")
	println(mf, "")
end
########################################################
#prints the generator data
function tpp_genNode(mf,nds)
	for n in nds
		print(mf,n.num,"\t")#adds bus ID
		print(mf,trunc(Int,n.mva),"\t")#rounds power to integer
		print(mf,0.0,"\t")#Qg
		print(mf,0.0,"\t")#Qmax
		print(mf,0.0,"\t")#Qmin
		print(mf,1.0,"\t")#Vg
		print(mf,trunc(Int,n.mva),"\t")#mbase
		print(mf,1.0,"\t")#status
		print(mf,trunc(Int,n.mva),"\t")#Pmax
		print(mf,0.0)#Pmin
		println(mf,";")
	end
end
########################################################

################################################################################################################
################################ Arcs ##########################################################################
################################################################################################################
function tpp_prntNeBrns(mf,mp,s)
	println(mf, "%candidate branch data")
	println(mf, "%column_names%	f_bus	t_bus	brnID	br_r	br_x	br_b	rate_a	rate_b	rate_c	tap	shift	br_status	angmin	angmax	construction_cost	branch_tnep_start")
	println(mf, "mpc.ne_branch = [")
	tpp_gpBrn(mf,mp.gParcs,mp.asBuilt,s)
	tpp_goBrn(mf,mp.gOarcs,mp.asBuilt,s)
	tpp_ooBrn(mf,mp.oOarcs,mp.asBuilt,s)
	tpp_opBrn(mf,mp.oParcs,mp.asBuilt,s)
	println(mf, "];")
end
################################################################################
#converts transformer and cable values to PU then calls .m printing
function tpp_xcXrlPu(a,arc,mf)
	eqpf_puImped(arc.cable,a.lngth)
	xxr=(arc.xfm.ohm/arc.xfm.num)/eqpD_pu()[4]
	xxl=(arc.xfm.xl/arc.xfm.num)/eqpD_pu()[4]
	tpp_prntBrn(mf,a,arc.cable,xxr,xxl)
end
################################################################################
function tpp_asBuilt(a,abs,mva,arc,cst,km)
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
function tpp_gpBrn(mf,as,abs,s)
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
				tpp_xcXrlPu(a,arc,mf)
				cst=arc.costs.ttl+2*cstD_cfs().FC_ac
				#println(cst)
				tpp_asBuilt(a,abs,mva,arc,cst,a.lngth)
				print(mf,cst,"\t")
				strt=tpp_1build(string(trunc(Int,arc.cable.mva*arc.cable.num))*string(a.tail.id)*string(a.head.id),s)
				print(mf,Float64(strt))
				println(mf,";")
			else
				println("No suitable "*string(a.tail.kv)*"Kv, "*string(mva)*"MVA cable for "*string(trunc(Int,a.lngth))*"Km. -removing")
			end
		end
		mvas=[]
	end
end
################################################################################
function tpp_goBrn(mf,as,abs,s)
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
				tpp_xcXrlPu(a,arc,mf)
				cst=arc.costs.ttl+(-1*cstD_cfs().FC_ac)
				tpp_asBuilt(a,abs,mva,arc,cst,a.lngth)
				print(mf,cst,"\t")
				strt=tpp_1build(string(trunc(Int,arc.cable.mva*arc.cable.num))*string(a.tail.id)*string(a.head.id),s)
				print(mf,Float64(strt))
				println(mf,";")
			else
				println("No suitable "*string(a.tail.kv)*"Kv, "*string(mva)*"MVA cable for "*string(trunc(Int,a.lngth))*"Km. -removing")
			end
		end
		mvas=[]
	end
end
################################################################################
function tpp_ooBrn(mf,as,abs,s)
#TODO: turb() should be updated to carry concession turbines
	mvas=Array{Float64,1}()
	arc=owpp()
	for a in as
		wp=wndF_wndPrf([a.tail.wnds[1]],turb())
		#push!(mvas,a.tail.mvas[1]/3))
		#push!(mvas,a.tail.mvas[1]/2)
		for mva in mvas
			cb=cstF_cbl_ttl(a.lngth,mva,lod_ossKv(),wp,true)
			if cb.num != 0
				eqpf_puImped(cb,a.lngth)
				tpp_prntBrn(mf,a,cb,0.0,0.0)
				arc.cable=cb
				cst=cb.results.ttl+2*cstD_cfs().FC_ac
				tpp_asBuilt(a,abs,mva,arc,cst,a.lngth)
				print(mf,cst,"\t")
				strt=tpp_1build(string(trunc(Int,arc.cable.mva*arc.cable.num))*string(a.tail.id)*string(a.head.id),s)
				print(mf,Float64(strt))
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
				tpp_prntBrn(mf,a,cb,0.0,0.0)
				arc.cable=cb
				cst=cb.results.ttl+2*cstD_cfs().FC_ac
				tpp_asBuilt(a,abs,mva,arc,cst,a.lngth)
				print(mf,cst,"\t")
				strt=tpp_1build(string(trunc(Int,cb.mva*cb.num))*string(a.tail.id)*string(a.head.id),s)
				print(mf,Float64(strt))
				println(mf,";")
			else
				println("No suitable "*string(lod_ossKv())*"Kv, "*string(mva)*"MVA cable for "*string(trunc(Int,a.lngth))*"Km. -removing")
			end
		end
		ka=[]
	end
end
################################################################################
function tpp_opBrn(mf,as,abs,s)
	#TODO: turb() should be updated to carry concession turbines
		mvas=Array{Float64,1}()
		for a in as
			#println(a)
			wp=wndF_wndPrf([a.tail.wnds[1]],turb())
			#push!(mvas,a.tail.mvas[1]/3))
			#push!(mvas,a.tail.mvas[1]/2)
			for mva in mvas
				arc=cstF_cblWT_ttl(a.lngth,mva,lod_ossKv(),wp,false)
				eqpD_xfoXR(lod_ossKv(),arc.xfm)
				if arc.cable.num != 0
					tpp_xcXrlPu(a,arc,mf)
					cst=arc.costs.ttl+cstD_cfs().FC_ac
					tpp_asBuilt(a,abs,mva,arc,cst,a.lngth)
					print(mf,cst,"\t")
					strt=tpp_1build(string(trunc(Int,arc.cable.mva*arc.cable.num))*string(a.tail.id)*string(a.head.id),s)
					print(mf,Float64(strt))
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
				    tpp_xcXrlPu(a,arc,mf)
					cst=arc.costs.ttl+cstD_cfs().FC_ac
					tpp_asBuilt(a,abs,mva,arc,cst,a.lngth)
					print(mf,cst,"\t")
					strt=tpp_1build(string(trunc(Int,arc.cable.mva*arc.cable.num))*string(a.tail.id)*string(a.head.id),s)
					print(mf,Float64(strt))
					println(mf,";")
				else
					println("No suitable "*string(lod_ossKv())*"Kv, "*string(mva)*"MVA cable for "*string(trunc(Int,a.lngth))*"Km. -removing")
				end
			end
			ka=[]
		end
	end
########################################################
function tpp_prntBrn(mf,a,cb,xr,xxl)
	eyed=string(trunc(Int,cb.mva*cb.num))*string(a.tail.id)*string(a.head.id)
	print(mf,a.tail.num,"\t")
	print(mf,a.head.num,"\t")
	print(mf,eyed,"\t")
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
function tpp_prntBrnOnShre(mf,fn,tn,p,ob)
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
	print(mf,30.0,"\t")
	print(mf,Float64(ob),"\t")
	println(mf,";")
end
########################################################
function tpp_prntBrns(mf,mp,ob)
	println(mf, "%branch data")
	println(mf, "%fbus	tbus	r	x	b	rateA	rateB	rateC	ratio	angle	status	angmin	angmax objmax")
	println(mf, "mpc.branch = [")
	mva=0.0
	tn=length(mp.pccs)+length(mp.cnces)+length(mp.osss)+1
	for i in mp.cnces
		mva=mva+i.mva
	end
	for i in mp.pccs
		tpp_prntBrnOnShre(mf,i.num,tn,mva,ob)
	end
	println(mf, "];")
	println(mf, "")
end
########################################################
