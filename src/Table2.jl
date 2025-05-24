
module table2

    export do_table2

    include("Custom_fct.jl")

    #using Pkg

    #Pkg.add("CategoricalArrays")


    using DataFrames, PrettyTables, LaTeXStrings, StatFiles, FreqTables, Statistics, StatsBase, CategoricalArrays, Distributions, CSV, GLM, .Custom

    function do_table2(m)

    function cycle_stats(df, stage; quintile_col=:inc_both, cohort=1)
        colnames = [:cost, :ln_cost, :inc, :ln_inc, :share,
                    :nhh, :birth_year, :group, :cndtl_cost,
                    :zero_share, :expdt, :share_expdt]
    
        out = DataFrame()                       # will grow to 5 × 12
    
        for g in 1:5
            sub = filter(r ->
                    coalesce(r[quintile_col] == g, false) &&
                    coalesce(r.birth_cohort  == cohort, false),
                    df)
    
            cost = mean(skipmissing(sub[!, Symbol("cost_sch$(stage)")] ))
            inc  = mean(skipmissing(sub.inc))
    
            ln_cost = (!ismissing(cost) && cost > 0) ? log(cost) : missing
            ln_inc  = (!ismissing(inc)  && inc  > 0) ? log(inc)  : missing
            share   = (!ismissing(cost) && !ismissing(inc)) ? cost / (inc/12) * 100 : missing
            nhh     = nrow(sub)
    
            # conditional mean (cost > 0)
            pos = filter(r -> begin
                         v = r[Symbol("cost_sch$(stage)")]
                         !ismissing(v) && v > 0
                     end, sub)
            cndtl_cost = mean(skipmissing(pos[!, Symbol("cost_sch$(stage)")] ))
            zero_share = nhh == 0 ? missing : (nhh - nrow(pos)) / nhh
    
            expdt = mean(skipmissing(sub.expenditure_month))
            share_expdt = (!ismissing(cost) && !ismissing(expdt)) ?
                           cost / expdt * 100 : missing
    
            # --- push as **NamedTuple** so columns match by name, not position ----
            push!(out, (;  cost,
                          ln_cost,
                          inc    = inc / 12,        # monthly
                          ln_inc,
                          share,
                          nhh,
                          birth_year = cohort,
                          group      = g,
                          cndtl_cost,
                          zero_share,
                          expdt,
                          share_expdt))
        end
    
        rename!(out, colnames)                  # give Stata-style names
        return out
    end
    
    # Parameters
    alow, ahigh       = 40, 43
    minobs            = 3
    cohortlength      = 9 - minobs        
    ng                = 5                

    # Load panel data and CPI, merge & deflate
    df_raw = DataFrame(load(joinpath(m, "data/data_clean.dta")))

    df_cpi = DataFrame(load(joinpath(m, "data/cpi.dta")))

    df = leftjoin(df_raw, df_cpi, on = :year)
    df = dropmissing(df, Symbol.(setdiff(names(df_cpi), ["year"])))  

    df.inc               .= df.inc               .* 100 ./ df.cpi
    df.expenditure_month .= df.expenditure_month .* 100 ./ df.cpi
    df.educost_pub_month .= df.educost_pub_month .* 100 ./ df.cpi
    df.educost_priv_month.= df.educost_priv_month.* 100 ./ df.cpi

    # Define birth-year cohorts
    cohort1max = 2017 - alow - minobs + 1        
    cohort1min = cohort1max - cohortlength + 1   
    cohort0max = 2008 - alow - minobs + 1        
    cohort0min = cohort0max - cohortlength + 1   

    df = filter(:birth_year_f => x -> !ismissing(x), df)     # keep only cohorts 0 & 1

    df.birth_cohort = assign_cohort.(
        df.birth_year_f,
        Ref(cohort0min), Ref(cohort0max),
        Ref(cohort1min), Ref(cohort1max)
    )
    df = filter(:birth_cohort => x -> !ismissing(x) && x in (0,1), df)     # keep only cohorts 0 & 1


    # Create income-quintile variables within each calendar year
    ng = 5                                    # number of quantile groups
    years = 4:20                              # survey years 4–20 (2004–2020)

    n = nrow(df)
    df.inc_both  = Vector{Union{Missing,Int}}(missing, n)
    df.inc_both0 = Vector{Union{Missing,Int}}(missing, n)
    df.inc_both1 = Vector{Union{Missing,Int}}(missing, n)

    for y in years
        idx = findall(df.year .== y)
        df.inc_both[idx] = xtile_ado(df.inc[idx], 5)
    end

    for y in years, b in (0,1)
        idx = findall(df.year .== y .&& df.birth_cohort .== b)
        q = xtile_ado(df.inc[idx], 5)
        if b == 0
            df.inc_both0[idx] = q
        else
            df.inc_both1[idx] = q
        end
    end

    ft = freqtable(df, :year, :inc_both)
    #println(ft)

    println("\nCohort 0:")
    println(freqtable(df[df.birth_cohort .== 0, :], :year, :inc_both0))

    println("\nCohort 1:")
    println(freqtable(df[df.birth_cohort .== 1, :], :year, :inc_both1))



    #println(describe(df)) #I am not able to get the same results as in stata!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    # Drop households where oldest child >24 or none children
    df = filter(row -> (ismissing(row.max_age_child) || row.max_age_child ≤ 24) &&
                   row.n_child_head > 0, df)

    # Per-child education costs
    df.educost_priv_cap = df.educost_priv_month ./ df.n_child_head
    df.educost_pub_cap  = df.educost_pub_month  ./ df.n_child_head
    df.educost_tot_cap  = df.educost_priv_cap .+ df.educost_pub_cap
    replace!(df.educost_priv_cap, missing=>0); replace!(df.educost_pub_cap, missing=>0)

    # Build result matrix for one-child households 
    results = DataFrame(age=Int[], educost_age=Float64[],
                    income_age=Float64[], share_age=Float64[],
                    nhh_age=Int[])

    for a in 0:24
        base   = filter(r -> r.n_child_head == 1 &&
                           r.max_age_child == a &&
                           r.birth_cohort  == 1, df)

        # --- ❶ mean education cost (missing costs automatically skipped) ----------
        mean_edu = mean(skipmissing(base.educost_tot_cap))

        # --- ❷ restrict to rows where cost is present -----------------------------
        with_cost = filter(r -> !ismissing(r.educost_tot_cap), base)

        mean_inc = mean(skipmissing(with_cost.inc)) / 12       # monthly
        n_hh     = nrow(with_cost)

        share = (!ismissing(mean_edu) && !ismissing(mean_inc)) ?
             mean_edu / mean_inc * 100 : missing

        push!(results, (a, mean_edu, mean_inc, share, n_hh))
    end
    # println(results) <- correspond to stata's result!

    # -------------------------------------------------------------------	
    # Fraction of total life-time education spending per child in income
    # -------------------------------------------------------------------
    total_edu   = sum(skipmissing(results.educost_age))
    total_inc   = sum(skipmissing(results.income_age))
    lifetime_pct = total_edu / total_inc * 100
	
    outp = DataFrame(lifetime_pct = lifetime_pct)   
    CSV.write(m*"/output/Lifetime_share.csv", outp) #Correct!

    # -------------------------------------------------------------------

    # -------------------------------------------------------------------

    # -------------------------------------------------------------------
	
    # merge education cost data for each children
    wide = DataFrame(load(joinpath(m, "data/data_educost_detail_wide.dta")))
    df = leftjoin(df, wide, on = [:hhid, :year])          
    # rows existing only in the wide file automatically discarded
    for i in 1:5
        c = Symbol("edu_totcost$(i)")
        df[!, c] .= df[!, c] .* 100 ./ df.cpi
    end

    # Remove info for child < 24 y.o.
    for i in 1:5
        costcol = Symbol("edu_totcost$(i)")
        agecol  = Symbol("chage$(i)")
        # set to missing if age > 24
        over24  = .!ismissing.(df[!, agecol]) .&& (df[!, agecol] .> 24)
        df[over24, costcol] .= missing
        df[over24, agecol]  .= missing
    end

    # Recode each child's age to school cycle 1-4 (pre=1,elem=2,middle=3,high=4)
    for i in 1:5
        df[!, Symbol("school$(i)")] = age_to_cycle.(df[!, Symbol("chage$(i)")])
    end

    # costs for children before going to a college
    for sch in 1:4
        cost_vec = zeros(Float64, nrow(df))
        n_vec    = zeros(Int,     nrow(df))

    for i in 1:5
        sc = Symbol("school$(i)")          # stage (1‥4 or missing)
        cc = Symbol("edu_totcost$(i)")     # cost for that child

        belongs = df[!, sc] .== sch        # Bool or Missing
        # Stata treats Missing as FALSE when it adds the indicator
        belongs_bool = coalesce.(belongs, false)

        # add cost only when child belongs to this stage
        cost_vec .+= ifelse.(belongs_bool,
                             coalesce.(df[!, cc], 0.0), 0.0)

        # ---- denominator ---------------------------------------------------
        # Every row contributes one observation; it contributes +1 if this
        # particular child is in the stage, +0 otherwise – exactly what
        #   rowtotal(ind*) , missing
        # does in Stata.
        n_vec .+= belongs_bool
    end

    # if a household has no child in the stage, Stata stores "."  (missing)
    cost_ave = ifelse.(n_vec .== 0, missing, cost_vec ./ n_vec);

    df[!, Symbol("cost_sch$(sch)")] = cost_ave;
    df[!, Symbol("n_sch$(sch)")]    = n_vec;
    end

    println(describe(df))

    # Build age specific results
    result1 = cycle_stats(df, 1)   # pre-school
    result2 = cycle_stats(df, 2)   # elementary
    result3 = cycle_stats(df, 3)   # middle
    result4 = cycle_stats(df, 4)   # high


    # Weighted average across all ages
    weights    = [7, 6, 3, 3] ./ 19                    
    stage_dfs  = [result1, result2, result3, result4]

    #convert DataFrame to Matrix	
    weighted_mat = reduce(+,  map((df, w) -> Matrix(df) .* w, stage_dfs, weights)) 
    result_avg = DataFrame(weighted_mat, names(result1))   # keep original col-names

    # recompute derived columns that depend on logs / ratios
    for r in 1:5
        result_avg.ln_cost[r] = log(result_avg.cost[r])
        result_avg.ln_inc[r]  = log(result_avg.inc[r])
        result_avg.share[r]   = result_avg.cost[r] / (result_avg.inc[r]) * 100
    end


    #------------------------


    elasticity_avg = coef(lm(@formula(ln_cost ~ ln_inc), result_avg))[2]
    println("Income elasticity of investment (weighted): ",
        round(elasticity_avg, digits = 3))
    ou = DataFrame(elasticity_avg = elasticity_avg)   
    CSV.write(m*"/output/Income_elasticity.csv", ou)


    #---------------------------
    table2 = DataFrame(
        pre     = result1.share,
        elem    = result2.share,
        middle  = result3.share,
        high    = result4.share,
        average = result_avg.share,
    )



    # Helper function to format values
fmt(x; d=1) = string(round(x, digits=d))

# Build table data (assuming `table2` already exists and contains correct `share` columns)
data = Matrix{Any}(undef, 5, 6)
for r in 1:5
    data[r,1] = "$(r)$(r==1 ? "st" : r==2 ? "nd" : r==3 ? "rd" : "th")"
    data[r,2] = fmt(table2.pre[r])
    data[r,3] = fmt(table2.elem[r])
    data[r,4] = fmt(table2.middle[r])
    data[r,5] = fmt(table2.high[r])
    data[r,6] = fmt(table2.average[r])
end

# Define header with subheaders
header = (
    ["Income quintile", "", "", "", "", ""],
    ["", "Preschool", "Elementary school", "Middle school", "High school", "Weighted average"]
)

# Horizontal lines
hlines = [1]

latex_path = joinpath(m, "output/Table2.tex")

# Write LaTeX file
open(latex_path, "w") do io
    println(io, "\\begin{center}")
    println(io, "\\textbf{Table 2—Percentage of Income Spent on Private Education}")
    println(io, "\\end{center}")
    println(io, "")

    pretty_table(io,
        data;
        header = header,
        backend = Val(:latex),
        tf = tf_latex_booktabs,
        alignment = :r,
        hlines = hlines
    )
end

    end

end
