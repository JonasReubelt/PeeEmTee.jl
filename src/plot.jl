
@recipe function f(prefit::PreFit)
    xs=-1:0.001:1
    seriestype := :scatter
    xguide := "charge [A.U.]"
    yguide := "counts"
    y = gauss(xs, [prefit.μₚ, prefit.σₚ, prefit.Aₚ])
    @series begin
        collect(xs), y
    end
end
