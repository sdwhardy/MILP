

function tnep_milp(ocn)
    cp("results/owpp_tnep_map.mat", "results/owpp_tnep_map.m", force=true)
    mip_solver = MosekSolver()
    mip_solver_drives = true
    log_level = 3
    result = run_tnep("results/owpp_tnep_map.m", DCPPowerModel, mip_solver)
    network_data = PowerModels.parse_file("results/owpp_tnep_map.m")
    fmap=tpp_reCnstrctSol(result,network_data,ocn)
    return fmap, result, network_data
end
