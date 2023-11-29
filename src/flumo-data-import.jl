function load_datasets(path_online, path_offline, experiments = experiments_total)
    dfs_online = []
    dfs_offline = []
    for i in experiments
        push!(dfs_online, DataFrame(XLSX.readtable(path_online, "run_$i", "A:D")))
        push!(dfs_offline, DataFrame(XLSX.readtable(path_offline, "run_$i", "A:Q")))
    end
    return (online=dfs_online, offline=dfs_offline)
end

function get_LDH_data()
    path_online = "data/LDH/231002_Data_Pulse_ges_online.xlsx"
    path_offline = "data/LDH/231002_Data_Pulse_ges_offline.xlsx"
    dfs_online, dfs_offline = load_datasets(path_online, path_offline, 1:42)
    return (online=dfs_online, offline=dfs_offline)
end