################################################################################
function tpp_reCnstrctSol(res,nt,ocn)
        s=tpp_xtrctBrch(res,nt)
        fileName="results/owpp_tnep_map.mat"
	mfile = open(fileName)
	syAsBuilt=tpp_xtrctOwpp(s,ocn.asBuilt)
	close(mfile)
	return syAsBuilt
end
################################################################################
function tpp_xtrctOwpp(s,ab)
	tmpSys=eez()
	for i in s
		for j in ab
			if (j.tl.num==i["f_bus"] && j.hd.num==i["t_bus"] && j.elec.costs.ttl==i["construction_cost"])
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
	return tmpSys
end
################################################################################
function tpp_xtrctBrch(r,n)
        sol=Array{Dict,1}()
        for i =1:length(r["solution"]["ne_branch"])
                if r["solution"]["ne_branch"][string(i)]["built"] == 1.0
                        push!(sol,n["ne_branch"][string(i)])
                end
        end
		println(sol)
	return sol
end
################################################################################
function tpp_prnt2Scrn(r,nt)
	for i =1:length(r["solution"]["ne_branch"])
	        if r["solution"]["ne_branch"][string(i)]["built"] == 1.0
	            print(nt["ne_branch"][string(i)]["f_bus"])
	            print(" - ")
	            print(nt["ne_branch"][string(i)]["t_bus"])
	            print(" Cost:")
	            println(nt["ne_branch"][string(i)]["construction_cost"])
	        end
	end
end
################################################################################
