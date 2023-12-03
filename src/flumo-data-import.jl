function load_datasets(path_online, path_offline, experiments = experiments_total)
    dfs_online = []
    dfs_offline = []
    for i in experiments
        push!(dfs_online, DataFrame(XLSX.readtable(path_online, "run_$i", "A:D")))
        push!(dfs_offline, DataFrame(XLSX.readtable(path_offline, "run_$i", "A:Q")))
    end
    return (online=dfs_online, offline=dfs_offline)
end

function load_datasets_Galox(path_online, path_offline, experiments = experiments_total)
    dfs_online = []
    df_offline = DataFrame(XLSX.readtable(path_offline, "All"))
    for i in experiments
        push!(dfs_online, DataFrame(XLSX.readtable(path_online, "Run$i")))
    end
    return (online=dfs_online, offline=df_offline)
end

function get_LDH_data()
    path_online = "data/LDH/231002_Data_Pulse_ges_online.xlsx"
    path_offline = "data/LDH/231002_Data_Pulse_ges_offline.xlsx"
    dfs_online, dfs_offline = load_datasets(path_online, path_offline, 1:42)
    return (online=dfs_online, offline=dfs_offline)
end

function get_GalOx_data()
    path_online = "data/GalOx/231017_Data_Combined_GalOx.xlsx"
    path_offline = "data/GalOx/231019_Galox_Off.xlsx"
    dfs_online, dfs_offline = load_datasets_Galox(path_online, path_offline, 1:7)
    return (online=dfs_online, offline=dfs_offline)
end

function get_HRP_data()
    path_online = "data/231117_Data_Combined_HRP.xlsx"
    #path_offline = "data/GalOx/231019_Galox_Off.xlsx"
    dfs_online = []
    #dfs_offline = []
    for i in 1:4
        push!(dfs_online, DataFrame(XLSX.readtable(path_online, "run$i", "B:D")))
    end
    return dfs_online
end
