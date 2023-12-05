function diff_AEW_LDH!(online)
    t = deepcopy(online[:,1]./60)
    #F = deepcopy(online.Integral)
    AEW = deepcopy(Float64.(online[:,2]))
    dAEWdt = diff(AEW)./diff(t)
    idxInf = findall(x->x>5, dAEWdt)
    for i in idxInf
        AEW[i+1:end] = AEW[i+1:end].-diff(AEW)[i]
    end
    #AEW = hampel(AEW, 1, 3)

    dAEWdt = diff(AEW)./diff(t)
    dAEWdt_smooth = savitzky_golay(AEW, 15, 5, deriv=1).y[2:end]./diff(t)
    # idxInf = findall(x->x==-Inf, dAEWdt)
    # for i in idxInf
    #     dAEWdt[i] = (dAEWdt[i+1] + dAEWdt[i-1])/2
    #     dAEWdt_smooth[i] = (dAEWdt_smooth[i+1] + dAEWdt_smooth[i-1])/2
    # end
    #dAEWdt = hampel(dAEWdt, 2, 2)
    #dAEWdt_smooth = hampel(dAEWdt_smooth, 2, 2)
    dAEWdt_smooth_loess = predict(loess(t[2:end], dAEWdt, span=0.09), t[2:end])
    dAEWdt = [dAEWdt[1]; dAEWdt]
    dAEWdt_smooth = [dAEWdt_smooth[1]; dAEWdt_smooth]
    dAEWdt_smooth_loess = [dAEWdt_smooth_loess[1]; dAEWdt_smooth_loess]
    
    replace!(dAEWdt_smooth_loess, -Inf=>0, NaN=>0)
    replace!(dAEWdt_smooth, -Inf=>0, NaN=>0)
    replace!(dAEWdt, -Inf=>0, NaN=>0)

    online.diff_AEW_raw = [maximum([0, -d]) for d in dAEWdt]
    online.diff_AEW_sgol = [maximum([0, -d]) for d in dAEWdt_smooth]
    online.diff_AEW_loess = [maximum([0, -d]) for d in dAEWdt_smooth_loess]
    online.AEW_adapted = AEW
    online.dAEW_adapted = AEW .- AEW[1]
end