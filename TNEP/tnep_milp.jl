function tnep_milp(ocn)
    mv("results/owpp_tnep_map.mat", "results/owpp_tnep_map.m", force=true)

    mip_solver = MosekSolver()
    mip2_solver=GurobiSolver(Presolve=1)
    mip3_solver=CbcSolver()
    mip_solver_drives = true
    log_level = 3

    result = run_tnep("results/owpp_tnep_map.m", DCPPowerModel, mip2_solver)

    #mip_solver = MosekSolver()
    #result = run_tnep("results/owpp_tnep_map.m", DCPPowerModel, mip_solver)
    network_data = PowerModels.parse_file("results/owpp_tnep_map.m")
    idd,fmap=tpp_reCnstrctSol(result,network_data,ocn)
    return idd, fmap, result, network_data
end
