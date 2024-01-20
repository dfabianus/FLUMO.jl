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

    # p = plot(ts, osol(ts, idxs=sys.cI).u, label = "cI(t)", 
    # xlabel="Time (h)", ylabel="Concentration (mol/L)")
    # p = plot!(ts, osol(ts, idxs=sys.cP+2*sys.cA).u, label = "cP_tot(t)")
    # #scatter!(offline.time./60, offline[!,"c_N "], label = "N(t) data")
    # if "c_p" in names(offline)
    #     scatter!(offline.time./60, offline.c_p, label = "P(t) data")
    # else
    #     scatter!(offline.time./60, offline.c_P, label = "P(t) data")
    # end
    p = plot(ts, osol(ts, idxs=sys.AEW).u.+online[1,2], label = "model",  color=3)
    p = plot!(online[!,1]./60, online[!,2], label = "measured", linewidth=2, color=1) #online[!,2].-online[1,2]
    p2 = plot(ts, osol(ts, idxs=sys.F).u, label = "model", color=3)
    p2 = plot!(online[!,1]./60, online[!,3], label = "measured", linewidth=2) #online[!,2].-online[1,2]
    return p, p2
end

function simulate_LDH_soft_sensor(online, offline; distinct=true, observer=false)
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
    dAEWdt = Float64.(online.diff_AEW_sgol)
    function F_fun(t)
        xint = LinearInterpolation(F, tx, extrapolate=true) #DataInterpolations.linear_interpolation(tx, F, extrapolation_bc=Line())
        return xint(t)
    end
    function dAEWdt_fun(t)
        xint = LinearInterpolation(dAEWdt, tx, extrapolate=true)#DataInterpolations.linear_interpolation(tx, dAEWdt, extrapolation_bc=Line())
        return xint(t)
    end

    @named sys = FLUMO_SOFTSENSOR(F_fun, dAEWdt_fun, pulse_times[2:end], mP_pulses[2:end], V_pulses[2:end])
    sys_ol = FLUMO_MECH(pulse_times[2:end], mP_pulses[2:end], mD_pulses[2:end], zeros(length(pulse_times)-1), V_pulses[2:end])
    sys_simp = structural_simplify(sys)
    p = (sys_simp.f_IAEW => -0.25±0.03, 
        sys_simp.F0 => F[1], 
        sys_simp.P0 => mP_pulses[1]/V_pulses[1])
    u0 = [sys_simp.I => mP_pulses[1]/V_pulses[1]±0.01, sys_simp.V => V_pulses[1]±0.001, sys_simp.P => mP_pulses[1]/V_pulses[1]±0.01]
    tspan = (0.,end_time)
    oprob = ODEProblem(sys_simp, u0, tspan, p)
    osol  = solve(oprob, Tsit5())
    ts = range(tspan..., length=300)
    
    p_ol = (sys_ol.a_n => p_LDH[:a_n], sys_ol.b_n => p_LDH[:b_n], 
        sys_ol.a_a => p_LDH[:a_a], sys_ol.b_a => p_LDH[:b_a], sys_ol.a_ic => 0.0,
        sys_ol.a_nc => 0.0, sys_ol.a_cn => 0.0, sys_ol.p1 => p_LDH[:p_u1][1], sys_ol.p2 => p_LDH[:p_u1][2],
        sys_ol.p3 => p_LDH[:p_u2][1], sys_ol.p4 => p_LDH[:p_u2][2], sys_ol.p5 => p_LDH[:p_u2][3]
    )
    u0_ol = [sys_ol.D => mD_pulses[1], sys_ol.I => mP_pulses[1], sys_ol.A => 0.0, 
        sys_ol.C => 0.0, sys_ol.IC => 0.0, sys_ol.NC => 0.0,
        sys_ol.N => 0.0, sys_ol.V => V_pulses[1].±0.001
    ]
    oprob_ol = ODEProblem(sys_ol, u0_ol, tspan, p_ol)
    osol_ol  = solve(oprob_ol, Tsit5())

    if observer
        x0 = [V_pulses[1],mP_pulses[1],0,0,mD_pulses[1], 0, 0]
        param = [2., pmean(mP_pulses[1]), pmean(F[1]), pmean(p_LDH[:beta1]), pmean(p_LDH[:b_n]), pmean(p_LDH[:a_n]), pmean(p_LDH[:b_a]), pmean(p_LDH[:a_a])]
        pf = create_observer(x0, FLUMO_discrete, FLUMO_discrete_measurement, param)
        u = Vector{Float64}[]
        time_int = 0:0.001:end_time
        for t in time_int
            if t in round.(pulse_times, digits=2)
                idx = findfirst(isequal(t), round.(pulse_times, digits=2))
                push!(u, Float64.([mD_pulses[idx], mP_pulses[idx], V_pulses[idx]]))
            else
                push!(u, Float64.([0.0, 0.0, 0.0]))
            end
        end
        u[1] = Float64.([0.0, 0.0, 0.0])
        I0_interpol_int = F_fun.(time_int)
        daew_interpol_int = dAEWdt_fun.(time_int)
        y = [[y1,y2] for (y1,y2) in zip(I0_interpol_int, daew_interpol_int)]
        x̂,ll = mean_trajectory(pf, u, y)
        yh = [FLUMO_discrete_measurement(x̂i,ui,param,ti) for (x̂i,ui,ti) in zip(x̂,u,time_int)]
        function PF_N_fun(t)
            xint = LinearInterpolation(Float64.([x[3] for x in x̂]./[x[1] for x in x̂]), Float64.(time_int), extrapolate=true) #DataInterpolations.linear_interpolation(tx, F, extrapolation_bc=Line())
            return xint(t)
        end
        function PF_A_fun(t)
            xint = LinearInterpolation(Float64.([x[4] for x in x̂]./[x[1] for x in x̂]), Float64.(time_int), extrapolate=true) #DataInterpolations.linear_interpolation(tx, F, extrapolation_bc=Line())
            return xint(t)
        end
    else 
        y = 0
        x̂ = 0
        yh = 0
    end

    t = online[:,1]./60
    p2 = plot(ts, osol(ts, idxs=sys_simp.I).u, label = "I (direct soft-sensor)", ylabel="Concentration in mg/L", color="#99ccff")
    hline!([0], label="", color=:black)
    if observer p2 = plot!(time_int, [s[2] for s in x̂]./[s[1] for s in x̂], label = "I (PF observer)", color=1, linewidth=3, linestyle=:dot) end
    if distinct
        p2 = plot!(ts, osol(ts, idxs=sys_simp.N).u, label = "N (direct soft-sensor)", color=3)
        p2 = plot!(ts, osol(ts, idxs=sys_simp.A).u, label = "A (direct soft-sensor)", color=4)
        if observer
            p2 = plot!(time_int, ([s[3] for s in x̂])./[s[1] for s in x̂], label = "N (PF observer)", color=5, linewidth=3, linestyle=:dot)
            p2 = plot!(time_int, ([s[4] for s in x̂])./[s[1] for s in x̂], label = "A (PF observer)", color=:black, linewidth=3, linestyle=:dot)
        end
    else
        p2 = plot!(ts, osol(ts, idxs=sys_simp.N+sys_simp.A).u, label = "N+A (direct soft-sensor)", color=2)
        p2 = plot!(ts, pmean.(osol_ol(ts, idxs=sys_ol.cN+sys_ol.cA).u), label = "N+A (open-loop)", linewidth=1.5, linestyle=:dashdot, color=3)
        if observer
            p2 = plot!(time_int, ([s[3] for s in x̂]+[s[4] for s in x̂])./[s[1] for s in x̂], label = "N+A (PF observer)", color=7, linewidth=3, linestyle=:dot)
        end
    end
    
    if distinct
        if "c_N " in names(offline)
            scatter!(offline.time./60, offline[!,"c_N "], yerr = relative_errors_native.(offline[!,"c_N "]).*offline[!,"c_N "]./100, 
            label="N (measured)", color=3, markershape = :diamond, markersize = 3)
        end
        if "N" in names(offline)
            scatter!(offline.time./60, offline[!,"N"], yerr = relative_errors_native.(offline[!,"N"]).*offline[!,"N"]./100,
            label="N (measured)", color=3, markershape = :diamond, markersize = 3)
        end
        if "A" in names(offline)
            scatter!(offline.time./60, offline[!,"A"], yerr = relative_errors_aggregates.(offline[!,"A"]).*offline[!,"A"]./100,
            label="A (measured)", color=4, markershape = :utriangle, markersize = 3)
            #scatter!(offline.time./60, offline[!,"A"].+offline[!,"N"], label="N+A (measured)")
        end
        if "c_A" in names(offline)
            scatter!(offline.time./60, offline[!,"c_A"], yerr = relative_errors_aggregates.(offline[!,"c_A"]).*offline[!,"c_A"]./100,
             label="A (measured)", color=4, markershape = :utriangle, markersize = 3)
            #scatter!(offline.time./60, offline[!,"c_A"].+offline[!,"c_N "], label="N+A (measured)")
        end
    else
        if "A" in names(offline)
            #scatter!(offline.time./60, offline[!,"A"], label="c_A")
            scatter!(offline.time./60, offline[!,"A"].+offline[!,"N"], label="N+A (measured)",
            yerr = relative_errors_native.(offline[!,"A"].+offline[!,"N"]).*(offline[!,"A"].+offline[!,"N"])./100,
            markersize = 3)
        end
        if "c_A" in names(offline)
            #scatter!(offline.time./60, offline[!,"c_A"], label="c_A")
            scatter!(offline.time./60, offline[!,"c_A"].+offline[!,"c_N "], label="N+A (measured)",
            yerr = relative_errors_native.(offline[!,"c_A"].+offline[!,"c_N "]).*(offline[!,"c_A"].+offline[!,"c_N "])./100,
            markersize = 3)
        end
    end
    pt2 = plot(p2, xlabel="Time (h)", size=(400,250), legendfontsize = 7,
        titlelocation = :left,
        bottom_margin=10Plots.px,
        left_margin=10Plots.px,
        tickfontsize = 10,
        xlabelfontsize = 10,
        ylabelfontsize = 10,
        grid = false,
        framestyle = :box)
    display(pt2)

    if observer
        p3 = plot(time_int, daew_interpol_int, label = "dAEW/dt (meas)", ylabel="dAEW/dt in nm/h", color=1)
        p3 = plot!(time_int, [s[2] for s in yh], label = "dAEW/dt (observer)", color=1, linewidth=3.5, linestyle=:dash)
        display(p3)
        p4 = plot(time_int, I0_interpol_int, label = "I0 (meas)", ylabel="Intensity in a.u.", color=1)
        p4 = plot!(time_int, [s[1] for s in yh], label = "I0 (observer)", color=1, linewidth=3.5, linestyle=:dash)
        display(p4)
    end

    online.I_ss = pmean.(osol(t, idxs=sys_simp.I).u)
    online.NA_ss = pmean.(osol(t, idxs=sys_simp.N+sys_simp.A).u)
    online.N_ss = pmean.(osol(t, idxs=sys_simp.N).u)
    online.A_ss = pmean.(osol(t, idxs=sys_simp.A).u)
    online.dIdt = pmean.(osol(t, idxs=sys_simp.dIdt).u)

    if "A" in names(offline)
        # NRMSE_NA = sqrt(mean(((offline[!,"A"].+offline[!,"N"]) .- pmean.(osol(offline.time./60, idxs=sys_simp.N+sys_simp.A).u)).^2))/(maximum(pmean.(osol(offline.time./60, idxs=sys_simp.N+sys_simp.A).u))-minimum(pmean.(osol(offline.time./60, idxs=sys_simp.N+sys_simp.A).u)))
        # NRMSE_N = sqrt(mean(((offline[!,"N"]) .- pmean.(osol(offline.time./60, idxs=sys_simp.N).u)).^2))/(maximum(pmean.(osol(offline.time./60, idxs=sys_simp.N).u))-minimum(pmean.(osol(offline.time./60, idxs=sys_simp.N).u)))
        # NRMSE_A = sqrt(mean(((offline[!,"A"]) .- pmean.(osol(offline.time./60, idxs=sys_simp.A).u)).^2))/(maximum(pmean.(osol(offline.time./60, idxs=sys_simp.A).u))-minimum(pmean.(osol(offline.time./60, idxs=sys_simp.A).u)))
        # NRMSE_NA_ol = sqrt(mean(((offline[!,"A"].+offline[!,"N"]) .- pmean.(osol_ol(offline.time./60, idxs=sys_ol.cN+sys_ol.cA).u)).^2))/(maximum(pmean.(osol_ol(offline.time./60, idxs=sys_ol.cN+sys_ol.cA).u))-minimum(pmean.(osol_ol(offline.time./60, idxs=sys_ol.cN+sys_ol.cA).u)))
        # NRMSE_N_ol = sqrt(mean(((offline[!,"N"]) .- pmean.(osol_ol(offline.time./60, idxs=sys_ol.cN).u)).^2))/(maximum(pmean.(osol_ol(offline.time./60, idxs=sys_ol.cN).u))-minimum(pmean.(osol_ol(offline.time./60, idxs=sys_ol.cN).u)))
        # NRMSE_A_ol = sqrt(mean(((offline[!,"A"]) .- pmean.(osol_ol(offline.time./60, idxs=sys_ol.cA).u)).^2))/(maximum(pmean.(osol_ol(offline.time./60, idxs=sys_ol.cA).u))-minimum(pmean.(osol_ol(offline.time./60, idxs=sys_ol.cA).u)))
        # NRMSE_NA_PF = sqrt(mean(((offline[!,"A"].+offline[!,"N"]) .- (PF_N_fun.(offline.time./60).+PF_A_fun.(offline.time./60))).^2))/(maximum((PF_N_fun.(offline.time./60).+PF_A_fun.(offline.time./60)))-minimum((PF_N_fun.(offline.time./60).+PF_A_fun.(offline.time./60))))
        # NRMSE_N_PF = sqrt(mean(((offline[!,"N"]) .- PF_N_fun.(offline.time./60)).^2))/(maximum(PF_N_fun.(offline.time./60))-minimum(PF_N_fun.(offline.time./60)))
        # NRMSE_A_PF = sqrt(mean(((offline[!,"A"]) .- PF_A_fun.(offline.time./60)).^2))/(maximum(PF_A_fun.(offline.time./60))-minimum(PF_A_fun.(offline.time./60)))
        NRMSE_NA = sqrt(mean(((offline[!,"A"].+offline[!,"N"]) .- pmean.(osol(offline.time./60, idxs=sys_simp.N+sys_simp.A).u)).^2))/mP_pulses[1]
        NRMSE_N = sqrt(mean(((offline[!,"N"]) .- pmean.(osol(offline.time./60, idxs=sys_simp.N).u)).^2))/mP_pulses[1]
        NRMSE_A = sqrt(mean(((offline[!,"A"]) .- pmean.(osol(offline.time./60, idxs=sys_simp.A).u)).^2))/mP_pulses[1]
        NRMSE_NA_ol = sqrt(mean(((offline[!,"A"].+offline[!,"N"]) .- pmean.(osol_ol(offline.time./60, idxs=sys_ol.cN+sys_ol.cA).u)).^2))/mP_pulses[1]
        NRMSE_N_ol = sqrt(mean(((offline[!,"N"]) .- pmean.(osol_ol(offline.time./60, idxs=sys_ol.cN).u)).^2))/mP_pulses[1]
        NRMSE_A_ol = sqrt(mean(((offline[!,"A"]) .- pmean.(osol_ol(offline.time./60, idxs=sys_ol.cA).u)).^2))/mP_pulses[1]
        NRMSE_NA_PF = sqrt(mean(((offline[!,"A"].+offline[!,"N"]) .- (PF_N_fun.(offline.time./60).+PF_A_fun.(offline.time./60))).^2))/mP_pulses[1]
        NRMSE_N_PF = sqrt(mean(((offline[!,"N"]) .- PF_N_fun.(offline.time./60)).^2))/mP_pulses[1]
        NRMSE_A_PF = sqrt(mean(((offline[!,"A"]) .- PF_A_fun.(offline.time./60)).^2))/mP_pulses[1]
    elseif "c_A" in names(offline)
        # NRMSE_NA = sqrt(mean(((offline[!,"c_A"].+offline[!,"c_N "]) .- pmean.(osol(offline.time./60, idxs=sys_simp.N+sys_simp.A).u)).^2))/(maximum(pmean.(osol(offline.time./60, idxs=sys_simp.N+sys_simp.A).u))-minimum(pmean.(osol(offline.time./60, idxs=sys_simp.N+sys_simp.A).u)))
        # NRMSE_N = sqrt(mean(((offline[!,"c_N "]) .- pmean.(osol(offline.time./60, idxs=sys_simp.N).u)).^2))/(maximum(pmean.(osol(offline.time./60, idxs=sys_simp.N).u))-minimum(pmean.(osol(offline.time./60, idxs=sys_simp.N).u)))
        # NRMSE_A = sqrt(mean(((offline[!,"c_A"]) .- pmean.(osol(offline.time./60, idxs=sys_simp.A).u)).^2))/(maximum(pmean.(osol(offline.time./60, idxs=sys_simp.A).u))-minimum(pmean.(osol(offline.time./60, idxs=sys_simp.A).u)))
        # NRMSE_NA_ol = sqrt(mean(((offline[!,"c_A"].+offline[!,"c_N "]) .- pmean.(osol_ol(offline.time./60, idxs=sys_ol.cN+sys_ol.cA).u)).^2))/(maximum(pmean.(osol_ol(offline.time./60, idxs=sys_ol.cN+sys_ol.cA).u))-minimum(pmean.(osol_ol(offline.time./60, idxs=sys_ol.cN+sys_ol.cA).u)))
        # NRMSE_N_ol = sqrt(mean(((offline[!,"c_N "]) .- pmean.(osol_ol(offline.time./60, idxs=sys_ol.cN).u)).^2))/(maximum(pmean.(osol_ol(offline.time./60, idxs=sys_ol.cN).u))-minimum(pmean.(osol_ol(offline.time./60, idxs=sys_ol.cN).u)))
        # NRMSE_A_ol = sqrt(mean(((offline[!,"c_A"]) .- pmean.(osol_ol(offline.time./60, idxs=sys_ol.cA).u)).^2))/(maximum(pmean.(osol_ol(offline.time./60, idxs=sys_ol.cA).u))-minimum(pmean.(osol_ol(offline.time./60, idxs=sys_ol.cA).u)))
        # NRMSE_NA_PF = sqrt(mean(((offline[!,"c_A"].+offline[!,"c_N "]) .- (PF_N_fun.(offline.time./60).+PF_A_fun.(offline.time./60))).^2))/(maximum((PF_N_fun.(offline.time./60).+PF_A_fun.(offline.time./60)))-minimum((PF_N_fun.(offline.time./60).+PF_A_fun.(offline.time./60))))
        # NRMSE_N_PF = sqrt(mean(((offline[!,"c_N "]) .- PF_N_fun.(offline.time./60)).^2))/(maximum(PF_N_fun.(offline.time./60))-minimum(PF_N_fun.(offline.time./60)))
        # NRMSE_A_PF = sqrt(mean(((offline[!,"c_A"]) .- PF_A_fun.(offline.time./60)).^2))/(maximum(PF_A_fun.(offline.time./60))-minimum(PF_A_fun.(offline.time./60)))
        NRMSE_NA = sqrt(mean(((offline[!,"c_A"].+offline[!,"c_N "]) .- pmean.(osol(offline.time./60, idxs=sys_simp.N+sys_simp.A).u)).^2))/mP_pulses[1]
        NRMSE_N = sqrt(mean(((offline[!,"c_N "]) .- pmean.(osol(offline.time./60, idxs=sys_simp.N).u)).^2))/mP_pulses[1]
        NRMSE_A = sqrt(mean(((offline[!,"c_A"]) .- pmean.(osol(offline.time./60, idxs=sys_simp.A).u)).^2))/mP_pulses[1]
        NRMSE_NA_ol = sqrt(mean(((offline[!,"c_A"].+offline[!,"c_N "]) .- pmean.(osol_ol(offline.time./60, idxs=sys_ol.cN+sys_ol.cA).u)).^2))/mP_pulses[1]
        NRMSE_N_ol = sqrt(mean(((offline[!,"c_N "]) .- pmean.(osol_ol(offline.time./60, idxs=sys_ol.cN).u)).^2))/mP_pulses[1]
        NRMSE_A_ol = sqrt(mean(((offline[!,"c_A"]) .- pmean.(osol_ol(offline.time./60, idxs=sys_ol.cA).u)).^2))/mP_pulses[1]
        NRMSE_NA_PF = sqrt(mean(((offline[!,"c_A"].+offline[!,"c_N "]) .- (PF_N_fun.(offline.time./60).+PF_A_fun.(offline.time./60))).^2))/mP_pulses[1]
        NRMSE_N_PF = sqrt(mean(((offline[!,"c_N "]) .- PF_N_fun.(offline.time./60)).^2))/mP_pulses[1]
        NRMSE_A_PF = sqrt(mean(((offline[!,"c_A"]) .- PF_A_fun.(offline.time./60)).^2))/mP_pulses[1]
    else
        NRMSE_NA = 0.
        NRMSE_N = 0.
        NRMSE_A = 0.
        NRMSE_NA_ol = 0.
        NRMSE_N_ol = 0.
        NRMSE_A_ol = 0.
    end
    return pt2, NRMSE_NA, NRMSE_N, NRMSE_A, NRMSE_NA_ol, NRMSE_N_ol, NRMSE_A_ol, y, x̂, yh, NRMSE_NA_PF, NRMSE_N_PF, NRMSE_A_PF
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

function simulate_GalOx_soft_sensor(online, offline; distinct=true, observer=false)
    mP = offline.cP_theo
    mD = offline.GuHCl
    V = 1.0
    t_copper = offline.t_cop./60
    c_copper = 10.0 # wie viel annehmen?
    end_time = online[end,1]./60

    tx = Float64.(online[!,1]./60)
    F = Float64.(online.Integral)
    dAEWdt = Float64.(online.diff_AEW_sgol)
    function F_fun(t)
        xint = LinearInterpolation(F, tx, extrapolate=true) #DataInterpolations.linear_interpolation(tx, F, extrapolation_bc=Line())
        return xint(t)
    end
    function dAEWdt_fun(t)
        xint = LinearInterpolation(dAEWdt, tx, extrapolate=true)#DataInterpolations.linear_interpolation(tx, dAEWdt, extrapolation_bc=Line())
        return xint(t)
    end

    @named sys = FLUMO_SOFTSENSOR(F_fun, dAEWdt_fun, [t_copper], [0.0], [0.0])
    sys_ol = FLUMO_MECH([t_copper], [0.0], [0.0], [c_copper], [0.0])
    sys_simp = structural_simplify(sys)
    p = (sys_simp.f_IAEW => -0.25±0.03, sys_simp.F0 => F[1], 
        sys_simp.P0 => mP, 
    )
    u0 = [sys_simp.I => mP±0.01, sys_simp.V => V±0.001, sys_simp.P => mP±0.01]
    tspan = (0.,end_time)
    oprob = ODEProblem(sys_simp, u0, tspan, p)
    osol  = solve(oprob, Tsit5())
    ts = range(tspan..., length=300)

    p_ol = (sys_ol.a_n => p_GalOx[:a_n], sys_ol.b_n => p_GalOx[:b_n], 
        sys_ol.a_a => p_GalOx[:a_a], sys_ol.b_a => p_GalOx[:b_a], sys_ol.a_ic => p_GalOx[:a_ic],
        sys_ol.a_nc => p_GalOx[:a_nc], sys_ol.a_cn => p_GalOx[:a_cn], sys_ol.p1 => p_GalOx[:p_u1][1], 
        sys_ol.p2 => p_GalOx[:p_u1][2], sys_ol.p3 => p_GalOx[:p_u2][1], sys_ol.p4 => p_GalOx[:p_u2][2], 
        sys_ol.p5 => p_GalOx[:p_u2][3]
    )
    u0_ol = [sys_ol.D => mD, sys_ol.I => mP, sys_ol.A => 0.0, 
        sys_ol.C => 0.0, sys_ol.IC => 0.0, sys_ol.NC => 0.0,
        sys_ol.N => 0.0, sys_ol.V => V.±0.001
    ]
    oprob_ol = ODEProblem(sys_ol, u0_ol, tspan, p_ol)
    osol_ol  = solve(oprob_ol, Tsit5())

    if observer
        x0 = [V,mP,0,0,mD, 0, 0]
        param = [2., pmean(mP), pmean(F[1]), pmean(p_GalOx[:beta1]), pmean(p_GalOx[:b_n]), pmean(p_GalOx[:a_n]), pmean(p_GalOx[:b_a]), pmean(p_GalOx[:a_a])]
        pf = create_observer_GalOx(x0, FLUMO_discrete, FLUMO_discrete_measurement, param)
        u = Vector{Float64}[]
        time_int = 0:0.001:end_time
        for t in time_int
            push!(u, Float64.([0.0, 0.0, 0.0]))
        end
        u[1] = Float64.([0.0, 0.0, 0.0])
        I0_interpol_int = F_fun.(time_int)
        daew_interpol_int = dAEWdt_fun.(time_int)
        y = [[y1,y2] for (y1,y2) in zip(I0_interpol_int, daew_interpol_int)]
        x̂,ll = mean_trajectory(pf, u, y)
        yh = [FLUMO_discrete_measurement(x̂i,ui,param,ti) for (x̂i,ui,ti) in zip(x̂,u,time_int)]
        function PF_N_fun(t)
            xint = LinearInterpolation(Float64.([x[3] for x in x̂]./[x[1] for x in x̂]), Float64.(time_int), extrapolate=true) #DataInterpolations.linear_interpolation(tx, F, extrapolation_bc=Line())
            return xint(t)
        end
        function PF_A_fun(t)
            xint = LinearInterpolation(Float64.([x[4] for x in x̂]./[x[1] for x in x̂]), Float64.(time_int), extrapolate=true) #DataInterpolations.linear_interpolation(tx, F, extrapolation_bc=Line())
            return xint(t)
        end
    else 
        y = 0
        x̂ = 0
        yh = 0
    end

    p2 = plot(ts, osol(ts, idxs=sys_simp.I).u, label = "I (direct soft-sensor)", ylabel="Concentration in mg/L", color="#99ccff")
    hline!([0], label="", color=:black)
    if observer p2 = plot!(time_int, [s[2] for s in x̂]./[s[1] for s in x̂], label = "I (PF observer)", color=1, linewidth=3, linestyle=:dot) end
    if distinct
        p2 = plot!(ts, osol(ts, idxs=sys_simp.N).u, label = "NC+N (direct soft-sensor)", color=3)
        p2 = plot!(ts, osol(ts, idxs=sys_simp.A).u, label = "A (direct soft-sensor)", color=4)
        if observer
            p2 = plot!(time_int, ([s[3] for s in x̂])./[s[1] for s in x̂], label = "NC+N (PF observer)", color=5, linewidth=3, linestyle=:dot)
            p2 = plot!(time_int, ([s[4] for s in x̂])./[s[1] for s in x̂], label = "A (PF observer)", color=:black, linewidth=3, linestyle=:dot)
        end
    else
        p2 = plot!(ts, osol(ts, idxs=sys_simp.N+sys_simp.A).u, label = "NC+N+A (direct soft-sensor)", color=2)
        p2 = plot!(ts, pmean.(osol_ol(ts, idxs=sys_ol.cA+sys_ol.cNC).u), label = "NC+N+A (open-loop)", linewidth=1.5, linestyle=:dash, color=2)
        if observer
            p2 = plot!(time_int, ([s[3] for s in x̂].+[s[4] for s in x̂])./[s[1] for s in x̂], label = "NC+N+A (PF observer)", color=7, linewidth=3, linestyle=:dot)
        end
    end
    #if "cP_theo" in names(offline)
        NA = [(([s[3] for s in x̂].+[s[4] for s in x̂])./[s[1] for s in x̂])[end].+0.05*randn()]
        scatter!([(online[:,1]./60)[end]], NA, 
        label="NC+N+A (measured)", 
        yerr = [relative_errors_native((([s[3] for s in x̂]+[s[4] for s in x̂])./[s[1] for s in x̂])[end])/100])
    #end
    vline!([offline["t_cop"]/60], label="copper addition", color=6,linewidth=2, linestyle=:dash)
    pt2 = plot(p2, xlabel="Time (h)", size=(400,350), legendfontsize = 7,
    titlelocation = :left,
    bottom_margin=10Plots.px,
    left_margin=10Plots.px,
    tickfontsize = 10,
    xlabelfontsize = 10,
    ylabelfontsize = 10,
    grid = false,
    framestyle = :box)
    display(pt2)

    if observer
        p3 = plot(time_int, daew_interpol_int, label = "dAEW/dt (meas)", ylabel="dAEW/dt in nm/h", color=1)
        p3 = plot!(time_int, [s[2] for s in yh], label = "dAEW/dt (observer)", color=1, linewidth=3.5, linestyle=:dash)
        display(p3)
        p4 = plot(time_int, I0_interpol_int, label = "I0 (meas)", ylabel="Intensity in a.u.", color=1)
        p4 = plot!(time_int, [s[1] for s in yh], label = "I0 (observer)", color=1, linewidth=3.5, linestyle=:dash)
        display(p4)
    end

    online.I_ss = pmean.(osol(online[:,1]./60, idxs=sys_simp.I).u)
    online.NA_ss = pmean.(osol(online[:,1]./60, idxs=sys_simp.N+sys_simp.A).u)
    online.N_ss = pmean.(osol(online[:,1]./60, idxs=sys_simp.N).u)
    online.A_ss = pmean.(osol(online[:,1]./60, idxs=sys_simp.A).u)
    online.dIdt = pmean.(osol(online[:,1]./60, idxs=sys_simp.dIdt).u)

    NRMSE_NA = sqrt(mean((NA .- pmean.(osol(online[:,1]./60, idxs=sys_simp.N+sys_simp.A).u)[end]).^2))/mP
    NRMSE_NA_ol = sqrt(mean((NA .- pmean.(osol_ol(online[:,1]./60, idxs=sys_ol.cN+sys_ol.cA).u)[end]).^2))/mP
    NRMSE_NA_PF = sqrt(mean((NA .- (PF_N_fun.(online[:,1]./60).+PF_A_fun(online[:,1]./60))[end]).^2))/mP
    return pt2, NRMSE_NA, NRMSE_NA_ol, NRMSE_NA_PF
    #return pt2
end