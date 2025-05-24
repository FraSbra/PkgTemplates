module figure1

export do_figure1

include("Custom_fct.jl")

using DataFrames, DataFramesMeta, StatFiles, Plots, Statistics, StatsBase, CategoricalArrays, Distributions, .Custom


function do_figure1(m)



function mea_ci(x; alpha = 0.05)
    data = collect(skipmissing(x))          # <- materialise the iterator
    n = length(data)
    if n == 0
        return (mean = missing, lb = missing, ub = missing)
    end
    μ  = mean(data)
    se = std(data) / sqrt(n)
    z  = quantile(Normal(), 1 - alpha/2)
    (mean = μ,
     lb   = μ - z*se,
     ub   = μ + z*se)
end

#----------------------------------------------------------------------------
#Replication of Figure 1
#----------------------------------------------------------------------------

#----------------------------------------------------------------------------
##Lines from 27 to 82 2_cohort
#----------------------------------------------------------------------------

#Parameters
alow = 40
ahigh = 43
maxobs = ahigh-alow +1    # age range for lifetime income
minobs = 3     		# minimum observation requirement
cohortlength = 9 - minobs	#max cohort length to stay within new sample periods (2009 onward) 
ng=5 				# number of income groups

df1 = DataFrame(load(m*"/data/data_clean.dta"));

# control inflation (2012 price)
df_cpi = DataFrame(load(m*"/data/cpi.dta"));
df_merged = leftjoin(df1, df_cpi, on=:year);
df = dropmissing(df_merged, Symbol.(setdiff(names(df_cpi), ["year"])));

sort!(df, [:hhid, :year])

#completed fertility
compfert = combine(groupby(df, :hhid), :n_child_head => maximum => :comp_fertility)
df = leftjoin(df, compfert, on = :hhid)

# Setting age restriction to construct lifetime income measure
df = filter(row -> !ismissing(row.age_head_f) && alow <= row.age_head_f <= ahigh, df)

# drop if hhid has two wives
birthyear_stats = combine(groupby(df, :hhid), 
    :birth_year_f => mean => :birth_year_f_mean)
df = leftjoin(df, birthyear_stats, on = :hhid)
df.birth_year_f_dev = abs.(df.birth_year_f .- df.birth_year_f_mean)
df = filter(:birth_year_f_dev => x -> x < 2, df)

# maximum & minimum age of wife in the data
age_stats = combine(groupby(df, :hhid), 
    :age_head_f => minimum => :minage, 
    :age_head_f => maximum => :maxage,
    nrow => :numobs)
df = leftjoin(df, age_stats, on = :hhid)

# keep only "married with spouse"
df.couple = ifelse.(.!ismissing.(df.age_head_m) .&& .!ismissing.(df.age_head_f), 1, 0)

df = filter(row -> row.couple == 1 &&
              (row.married_f isa Missing || 
			  (row.married_f != 4 && row.married_f != 5)),
            df)


# including samples observed any x periods 
df = filter(:numobs => x -> x >= minobs, df)

lifeinc_stats = combine(groupby(df, :hhid), :inc => mean => :lifeinc)
df = leftjoin(df, lifeinc_stats, on = :hhid)

# Define cohort age band (based on 2017 being most recent year)
cohort1max = 2017 - alow - minobs + 1
cohort1min = cohort1max - cohortlength + 1
cohort0max = 2008 - alow - minobs + 1
cohort0min = cohort0max - cohortlength + 1

# generate birth cohort (automatically select)
df.birth_cohort = assign_cohort.(
	df.birth_year_f,
	Ref(cohort0min),
	Ref(cohort0max),
	Ref(cohort1min),
	Ref(cohort1max)
)
	
# Remove repetition of households
df = filter(:birth_cohort => x -> !ismissing(x) && x in (0, 1), df)
df = combine(groupby(df, :hhid), first) #cohort=0: Old cohort, cohort=1 Young cohort

# Create childless variable
df.childless = Int.(df.comp_fertility .== 0)


#----------------------------------------------------------------------------
# Lines 112 to 168
#----------------------------------------------------------------------------
results = DataFrame(
    birth_year = Int[],
    group = Int[],
    totfert_inc_b = Float64[],
    totfert_inc_b_ub = Float64[],
    totfert_inc_b_lb = Float64[],
    ln_totfert_inc_b = Float64[],
    cless_inc_b = Float64[],
    cless_inc_b_ub = Float64[],
    cless_inc_b_lb = Float64[],
    inc_b = Float64[],
    ln_inc_b = Float64[],
    n_total = Int[],
    n_single = Int[],
)


	row = 1

for b in (0, 1)
    df_b = filter(:birth_cohort => ==(b), df)

    # Create income quantiles (line 127)
    df_b.inc_both = cut(df_b.lifeinc, quantile(skipmissing(df_b.lifeinc), 0:1/ng:1), 	labels=1:ng, extend=true)

    for g in 1:ng
        df_bg = filter(row -> !ismissing(row.inc_both) && row.inc_both == g, df_b)

        fert_stats = mea_ci(df_bg.comp_fertility)
        childless_stats = mea_ci(df_bg.childless)
        inc_stats = mea_ci(df_bg.lifeinc)

        total_n = nrow(df_bg)
        n_single = count(r -> r.couple != 1, eachrow(df_bg))

        push!(results, (
            b,
            g,
            fert_stats.mean,
            fert_stats.ub,
            fert_stats.lb,
            log(fert_stats.mean),
            childless_stats.mean,
            childless_stats.ub,
            childless_stats.lb,
            inc_stats.mean,
            log(inc_stats.mean),
            total_n,
            n_single
        ))
    end
end


#----------------------------------------------------------------------------
# FIGURE 1!
#----------------------------------------------------------------------------
# Filter data for birth_year == 1
fig1a_data = filter(:birth_year => ==(1), results)

# Create the plot

p1 = plot(
    fig1a_data.ln_inc_b,
    fig1a_data.totfert_inc_b,
    seriestype = :scatter,
    markershape = :diamond,
    markersize = 6,
    label = "",
    color = :blue,
    xlabel = "Log family income in 2012 KRW",
    ylabel = "Number of children",
    xtickfont = font(12),
    ytickfont = font(12),
    guidefont = font(14),
    ylims = (1.7, 2.15),
    yticks = 1.7:0.1:2.1,
    xlims = (7.5, 9.5),
	size = (500, 500),
)

# Add line connecting the points
plot!(
    fig1a_data.ln_inc_b,
    fig1a_data.totfert_inc_b,
    seriestype = :line,
    linewidth = 1,
    color = :blue,
    label = ""
)

p1
# Save
savefig(p1, m*"/output/Figure_1_A.pdf")

#----------------------------------------------------------------------------
# FIGURE 2 (Wrong)
#----------------------------------------------------------------------------


# Create new column with percentage
results.cless_inc_b_pct = results.cless_inc_b .* 100

# Filter for birth_year == 1
fig1b_data = filter(:birth_year => ==(1), results)

# Create the plot
p2 = plot(
    fig1b_data.ln_inc_b,
    fig1b_data.cless_inc_b_pct,
    seriestype = :scatter,
    markershape = :diamond,
    markersize = 6,
    line = (:solid, 0.8),
    label = "",
    color = :red,
    xlabel = "Log family income in 2012 KRW",
    ylabel = "Childlessness rate (%)",
    xtickfont = font(12),
    ytickfont = font(12),
    guidefont = font(14),
	ylims = (0,5.5),
	xlims = (7.5, 9.5),
	size = (500, 500),
)


plot!(
    fig1b_data.ln_inc_b,
    fig1b_data.cless_inc_b_pct,
    seriestype = :line,
    linewidth = 1,
    color = :red,
    label = ""
)


# Save 
savefig(p2, m*"/output/Figure_1_B.pdf")

end

end