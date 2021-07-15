
@recipe function f(prefit::PreFit, xs; labels=("pedestal", "spe"))
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

@recipe function f(pmtrespfit::PMTRespFit, xs; show_single_gaussians=false, n_gaussians=5)
    xguide := "charge [A.U.]"
    yguide := "counts"
    yscale := :log10
    ylims := (1, 1e5)
    @series begin
        linewidth := 2
        xs, pmtresp(xs, [pmtrespfit.μₚ, pmtrespfit.σₚ,
                         pmtrespfit.μₛ, pmtrespfit.σₛ, pmtrespfit.nₚₑ,
                         pmtrespfit.A])
    end
    if show_single_gaussians
        poisson = Poisson(abs(pmtrespfit.nₚₑ))
        for i in 1:n_gaussians
            @series begin
                linewidth := 2
                linestyle := :dash
                seriescolor := :black
                xs, pmtrespfit.A .* gauss(xs,
                                        [i * pmtrespfit.μₛ,
                                         sqrt(i) * pmtrespfit.σₛ,
                                         pdf(poisson, i)])
            end
        end
    end
end

@recipe function f(pmtrespfit::PMTRespUapFit, xs; show_single_gaussians=false)
    xguide := "charge [A.U.]"
    yguide := "counts"
    yscale := :log10
    ylims := (1, 1e5)
    @series begin
        linewidth := 2
        xs, pmtresp_uap(xs, [pmtrespfit.μₚ, pmtrespfit.σₚ,
                         pmtrespfit.μₛ, pmtrespfit.σₛ, pmtrespfit.nₚₑ,
                         pmtrespfit.A, pmtrespfit.μᵤₐₚ, pmtrespfit.σᵤₐₚ,
                         pmtrespfit.Aᵤₐₚ])
    end
    if show_single_gaussians
        println("fdsagf")
        poisson = Poisson(abs(pmtrespfit.nₚₑ))
        for i in 1:n_gaussians
            @series begin
                linewidth := 2
                linestyle := :dash
                seriescolor := :black
                xs, pmtrespfit.A .* gauss(xs,
                                        [i * pmtrespfit.μₛ,
                                         sqrt(i) * pmtrespfit.σₛ,
                                         pdf(poisson, i)])
            end
        end
    end
end