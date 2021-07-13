
@recipe function f(prefit::PreFit; labels=("pedestal", "spe"))
    xs=-.2:0.001:1
    xguide := "charge [A.U.]"
    yguide := "counts"
    yscale := :log10
    ylims := (1, 1e5)
    @series begin
        linewidth := 2
        label := labels[1]
        xs, gauss(xs, [prefit.μₚ, prefit.σₚ, prefit.Aₚ])
    end
    @series begin
        linewidth := 2
        label := labels[2]
        xs, gauss(xs, [prefit.μₛ, prefit.σₛ, prefit.Aₛ])
    end
end

@recipe function f(pmtresp::PMTRespFit)
    xs=-.2:0.001:1
    xguide := "charge [A.U.]"
    yguide := "counts"
    yscale := :log10
    ylims := (1, 1e5)
    @series begin
        linewidth := 2
        xs, pmtresp(xs, [pmtresp.μₚ, pmtresp.σₚ,
                         pmtresp.μₛ, pmtresp.σₛ, pmtresp.nₚₑ,
                         pmtresp.A])
    end
end

@recipe function f(pmtresp::PMTRespUapFit)
    xs=-.2:0.001:1
    xguide := "charge [A.U.]"
    yguide := "counts"
    yscale := :log10
    ylims := (1, 1e5)
    @series begin
        linewidth := 2
        xs, pmtresp_uap(xs, [pmtresp.μₚ, pmtresp.σₚ,
                         pmtresp.μₛ, pmtresp.σₛ, pmtresp.nₚₑ,
                         pmtresp.A, pmtresp.μᵤₐₚ, pmtresp.σᵤₐₚ,
                         pmtresp.Aᵤₐₚ])
    end
end