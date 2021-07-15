"""
    $(SIGNATURES)
    quality function constructor for least square optimization
"""
function make_qfunc(model, x, y)
    function qfunc(p)
        mask = (y .!= 0)
        model_res = model(x, p)
        sum(((model_res[mask] .- y[mask]).^2) ./ y[mask])
    end
    qfunc
end


"""
    $(SIGNATURES)
    performs pre fit
    # Arguments
    - `chargedist => ChargeDist`: charge distribution on which pre fit is
                                  performed
"""
function pre_fit(chargedist::ChargeDist)
    x = chargedist.x
    y = chargedist.y
    qfunc_ped = make_qfunc(gauss, x, y)
    mxval_ped, mxidx_ped = findmax(y)
    p0_ped = [x[mxidx_ped], 0.01, mxval_ped]
    fit_ped = optimize(qfunc_ped, p0_ped, NewtonTrustRegion())
    popt_ped = Optim.minimizer(fit_ped)
    mask = (x .> (popt_ped[1] + 3 * popt_ped[2]))
    x_spe = x[mask]
    y_spe = y[mask]
    mxval_spe, mxidx_spe = findmax(y_spe)
    mask = (x_spe .< x_spe[mxidx_spe] * 2)
    x_spe = x_spe[mask]
    y_spe = y_spe[mask]
    qfunc_spe = make_qfunc(gauss, x_spe, y_spe)
    p0_spe = [x_spe[mxidx_spe], 0.1, mxval_spe]
    fit_spe = optimize(qfunc_spe, p0_spe, NewtonTrustRegion())
    popt_spe = Optim.minimizer(fit_spe)
    PreFit(popt_ped..., popt_spe...)
end

"""
    $(SIGNATURES)
    performs pmt response fit
    # Arguments
    - `charges => ChargeDist`: charge distribution on which fit is performed
    - `prefit_results => PreFitResults`: results from pre fit used as starting
                                         values for fit
"""
function pmtresp_fit(chargedist::ChargeDist, prefit::PreFit; mod=:default)
    p0 = [prefit.μₚ,
          prefit.σₚ,
          prefit.μₛ,
          prefit.σₛ,
          prefit.Aₛ / prefit.Aₚ,
          prefit.Aₚ + prefit.Aₛ]
    if mod == :uap
        push!(p0, prefit.μₛ / 5)
        push!(p0, prefit.σₛ / 5)
        push!(p0, prefit.Aₛ / 5)
        func = pmtresp_uap
        ResultStruct = PMTRespUapFit
    else
        func = pmtresp
        ResultStruct = PMTRespFit
    end
    qfunc = make_qfunc(func, chargedist.x, chargedist.y)
    fit = optimize(qfunc, p0, Newton())
    popt = Optim.minimizer(fit)
    
    ResultStruct(popt...)
end