function simulate_LDH_experiment(online, offline)
    mD_pulses = vcat(offline.c_GuHCl[1], diff(offline.c_GuHCl[1:end-1])).*offline.V[1:end-1]./1000
    V_pulses = vcat(offline.V[1], diff(offline.V[1:end-1]))./1000
    pulse_times = offline.time[1:end-1]./60
    end_time = online[end,1]./60
    if "c_p" in names(offline)
        mP_pulses = vcat(offline.c_p[1], diff(offline.c_p[1:end-1])).*offline.V[1:end-1]./1000
    else
        mP_pulses = [offline.c_P[1].*offline.V[1]./1000]
        V_pulses = [V_pulses[1]]
        pulse_times = [pulse_times[1]]
    end

    sys = FLUMO_MECH(pulse_times[2:end], mP_pulses[2:end], mD_pulses[2:end], zeros(length(pulse_times)-1), V_pulses[2:end])
    p = (sys.a_n => p_LDH[:a_n], sys.b_n => p_LDH[:b_n], 
        sys.a_a => p_LDH[:a_a], sys.b_a => p_LDH[:b_a], sys.a_ic => 0.0,
        sys.a_nc => 0.0, sys.a_cn => 0.0, sys.p1 => p_LDH[:p_u1][1], sys.p2 => p_LDH[:p_u1][2],
        sys.p3 => p_LDH[:p_u2][1], sys.p4 => p_LDH[:p_u2][2], sys.p5 => p_LDH[:p_u2][3]
    )
    u0 = [sys.D => mD_pulses[1], sys.I => mP_pulses[1], sys.A => 0.0, 
        sys.C => 0.0, sys.IC => 0.0, sys.NC => 0.0,
        sys.N => 0.0, sys.V => V_pulses[1].±0.001
    ]
    tspan = (0.,end_time)
    oprob = ODEProblem(sys, u0, tspan, p)
    osol  = solve(oprob, Tsit5())
    ts = range(tspan..., length=300)

    p = plot(ts, osol(ts, idxs=sys.cI).u, label = "cI(t)", 
    xlabel="Time (h)", ylabel="Concentration (mol/L)")
    p = plot!(ts, osol(ts, idxs=sys.cP+2*sys.cA).u, label = "cP_tot(t)")
    #scatter!(offline.time./60, offline[!,"c_N "], label = "N(t) data")
    if "c_p" in names(offline)
        scatter!(offline.time./60, offline.c_p, label = "P(t) data")
    else
        scatter!(offline.time./60, offline.c_P, label = "P(t) data")
    end
    p2 = plot(ts, osol(ts, idxs=sys.F).u, label = "model")
    p2 = plot!(online[!,1]./60, online[!,3], label = "measured") #online[!,2].-online[1,2]
    return p, p2
end


function simulate_LDH_soft_sensor(online, offline)
    mD_pulses = vcat(offline.c_GuHCl[1], diff(offline.c_GuHCl[1:end-1])).*offline.V[1:end-1]./1000
    V_pulses = vcat(offline.V[1], diff(offline.V[1:end-1]))./1000
    pulse_times = offline.time[1:end-1]./60
    end_time = online[end,1]./60
    if "c_p" in names(offline)
        mP_pulses = vcat(offline.c_p[1], diff(offline.c_p[1:end-1])).*offline.V[1:end-1]./1000
    else
        mP_pulses = [offline.c_P[1].*offline.V[1]./1000]
        V_pulses = [V_pulses[1]]
        pulse_times = [pulse_times[1]]
    end
    tx = Float64.(online[!,1]./60)
    F = Float64.(online.Integral)
    dAEWdt = Float64.(online.diff_AEW)
    function F_fun(t)
        xint = LinearInterpolation(F, tx, extrapolate=true) #DataInterpolations.linear_interpolation(tx, F, extrapolation_bc=Line())
        return xint(t)
    end
    function dAEWdt_fun(t)
        xint = LinearInterpolation(dAEWdt, tx, extrapolate=true)#DataInterpolations.linear_interpolation(tx, dAEWdt, extrapolation_bc=Line())
        return xint(t)
    end

    @named sys = FLUMO_SOFTSENSOR(F_fun, dAEWdt_fun, pulse_times[2:end], mP_pulses[2:end], V_pulses[2:end])
    sys_simp = structural_simplify(sys)
    p = (sys_simp.f_IAEW => -0.1±0.01, sys_simp.F0 => F[1], 
        sys_simp.P0 => mP_pulses[1]/V_pulses[1], 
    )
    u0 = [sys_simp.I => mP_pulses[1]/V_pulses[1]±0.001, sys_simp.V => V_pulses[1]±0.001]

    # sys = FLUMO_MECH(pulse_times[2:end], mP_pulses[2:end], mD_pulses[2:end], zeros(length(pulse_times)-1), V_pulses[2:end])
    # p = (sys.a_n => p_LDH[:a_n], sys.b_n => p_LDH[:b_n], 
    #     sys.a_a => p_LDH[:a_a], sys.b_a => p_LDH[:b_a], sys.a_ic => 0.0,
    #     sys.a_nc => 0.0, sys.a_cn => 0.0, sys.p1 => p_LDH[:p_u1][1], sys.p2 => p_LDH[:p_u1][2],
    #     sys.p3 => p_LDH[:p_u2][1], sys.p4 => p_LDH[:p_u2][2], sys.p5 => p_LDH[:p_u2][3]
    # )
    # u0 = [sys.D => mD_pulses[1], sys.I => mP_pulses[1], sys.A => 0.0, 
    #     sys.C => 0.0, sys.IC => 0.0, sys.NC => 0.0,
    #     sys.N => 0.0, sys.V => V_pulses[1].±0.001
    # ]

    tspan = (0.,end_time)
    oprob = ODEProblem(sys_simp, u0, tspan, p)
    osol  = solve(oprob, Tsit5())
    ts = range(tspan..., length=300)

    #p = plot(ts, osol(ts, idxs=sys_simp.cI).u, label = "cI(t)", 
    # xlabel="Time (h)", ylabel="Concentration (mol/L)")
    # p = plot!(ts, osol(ts, idxs=sys.cP+2*sys.cA).u, label = "cP_tot(t)")
    # #scatter!(offline.time./60, offline[!,"c_N "], label = "N(t) data")
    # if "c_p" in names(offline)
    #     scatter!(offline.time./60, offline.c_p, label = "P(t) data")
    # else
    #     scatter!(offline.time./60, offline.c_P, label = "P(t) data")
    # end
    # p2 = plot(ts, osol(ts, idxs=sys.F).u, label = "model")
    # p2 = plot!(online[!,1]./60, online[!,3], label = "measured") #online[!,2].-online[1,2]
    p = plot(ts, osol(ts, idxs=sys_simp.dAEWdt).u, label = "dAEWdt(t)")
    p2 = plot(ts, osol(ts, idxs=sys_simp.I).u, label = "I(t)")
    p2 = plot!(ts, osol(ts, idxs=sys_simp.N).u, label = "N(t)")
    p2 = plot!(ts, osol(ts, idxs=sys_simp.A).u, label = "A(t)")
    p2 = plot!(ts, osol(ts, idxs=sys_simp.N+sys_simp.I).u, label = "SOL(t)")
    p2 = plot!(ts, osol(ts, idxs=sys_simp.N+sys_simp.A).u, label = "REACT(t)")
    if "c_p" in names(offline)
        scatter!(offline.time./60, offline.c_p, label = "P(t) data")
    else
        scatter!(offline.time./60, offline.c_P, label = "P(t) data")
    end
    if "c_sol" in names(offline)
        scatter!(offline.time./60, offline.c_sol, label="c_sol")
    end
    if "c_P_sol" in names(offline)
        scatter!(offline.time./60, offline.c_P_sol, label="c_sol")
    end
    if "c_N " in names(offline)
        scatter!(offline.time./60, offline[!,"c_N "], label="c_N")
    end
    if "N" in names(offline)
        scatter!(offline.time./60, offline[!,"N"], label="c_N")
    end
    if "A" in names(offline)
        scatter!(offline.time./60, offline[!,"A"], label="c_A")
    end
    if "c_A" in names(offline)
        scatter!(offline.time./60, offline[!,"c_A"], label="c_A")
    end
    p3 = plot(ts, osol(ts, idxs=sys_simp.F).u, label = "F(t)")
    p4 = plot(online[!,1]./60, online[!,2], label = "measured")
    pt = plot(p, p3, p2,p4, layout=(2,2), xlabel="Time (h)", size = (1000,1000))
    display(pt)
    return pt
end

function simulate_GalOx_experiment(online, offline)
    mP = offline.cP_theo
    mD = offline.GuHCl
    V = 1.0
    t_copper = offline.t_cop./60
    c_copper = 10.0 # wie viel annehmen?
    end_time = online[end,1]./60

    sys = FLUMO_MECH([t_copper], [0.0], [0.0], [c_copper], [0.0]) # add cofactor pulse
    p = (sys.a_n => p_GalOx[:a_n], sys.b_n => p_GalOx[:b_n], 
        sys.a_a => p_GalOx[:a_a], sys.b_a => p_GalOx[:b_a], sys.a_ic => p_GalOx[:a_ic],
        sys.a_nc => p_GalOx[:a_nc], sys.a_cn => p_GalOx[:a_cn], sys.p1 => p_GalOx[:p_u1][1], 
        sys.p2 => p_GalOx[:p_u1][2], sys.p3 => p_GalOx[:p_u2][1], sys.p4 => p_GalOx[:p_u2][2], 
        sys.p5 => p_GalOx[:p_u2][3]
    )
    u0 = [sys.D => mD, sys.I => mP, sys.A => 0.0, 
        sys.C => 0.0, sys.IC => 0.0, sys.NC => 0.0,
        sys.N => 0.0, sys.V => V.±0.001
    ]
    tspan = (0.,end_time)
    oprob = ODEProblem(sys, u0, tspan, p)
    osol  = solve(oprob, Tsit5())
    ts = range(tspan..., length=300)

    p = plot(ts, osol(ts, idxs=sys.cN).u, label = "cN(t)", 
    title = "Testset MTK Catalyst Connections", xlabel="Time (h)", ylabel="Concentration (mol/L)")
    #p = plot!(ts, osol(ts, idxs=sys.cI).u, label = "cI(t)")
    #p = plot!(ts, osol(ts, idxs=sys.cA).u, label = "cA(t)")
    p = plot!(ts, osol(ts, idxs=sys.cP+2*sys.cA).u, label = "cP(t)")
    p = plot!(ts, osol(ts, idxs=sys.cNC).u, label = "cNC(t)")
    #scatter!(offline.time./60, offline[!,"c_N "], label = "N(t) data")
    scatter!([0., end_time], [mP, offline.cP_sol], label = "P(t) data")
    p2 = plot(ts, osol(ts, idxs=sys.F).u, label = "model")
    p2 = plot!(online[!,1]./60, online[!,end-1], label = "measured") #online[!,2].-online[1,2]
    return p, p2
    #p = plot(ts, osol(ts, idxs=sys.k_n).u, label = "k_n", title = "Testset MTK Catalyst Connections", 
    #xlabel="Time (h)", ylabel="Concentration (mol/L)")
end