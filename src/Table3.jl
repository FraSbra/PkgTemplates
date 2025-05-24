

module table3

export do_table3

using Pkg
using StatsModels
using DataFrames, StatFiles, StatsBase, Statistics, CSV, FixedEffectModels, CategoricalArrays, ShiftedArrays
using Econometrics, RegressionTables, WildBootTests, Distributions
using StatsModels: coefnames, @formula 
using GLM, Random, Statistics
using Random: MersenneTwister
using PrettyTables
using LaTeXStrings

function do_table3(m)

#Pkg.add("StatsModels")
"""
Pkg.add("Econometrics")
Pkg.add("ShiftedArrays")
Pkg.add("DataFrames")
Pkg.add("StatFiles")
Pkg.add("StatsBase")
Pkg.add("Statistics")

Pkg.add("CSV")
Pkg.add("FixedEffectModels")
Pkg.add("CategoricalArrays")
Pkg.add("PrettyTables")
Pkg.add("RegressionTables")
Pkg.add("WildBootTests")
Pkg.add("Distributions")
Pkg.add("PrettyTables")     
Pkg.add("LaTeXStrings")
"""

# Load packages




#funcspace
# Assignment function
function assign_curfews!(df)
    for row in eachrow(df)
        p, y = row.province, row.year
        if p == 1 && y >= 11  # Seoul
            row.curfew_mid10 = 1; row.curfew_mid11 = 1; row.curfew_mid12 = 1
        elseif p == 2 && y >= 11  # Busan
            row.curfew_mid10 = 1; row.curfew_mid11 = 1; row.curfew_mid12 = 1
        elseif p == 3 && y >= 14  # Daegu: first two
            row.curfew_mid10 = 1; row.curfew_mid11 = 1
        end
        if p == 3 && y >= 3      # Daegu: third
            row.curfew_mid12 = 1
        end
        if p == 4 && y >= 12     # Daejeon
            row.curfew_mid11 = 1; row.curfew_mid12 = 1
        end
        if p == 5 && y >= 15     # Incheon
            row.curfew_mid10 = 1; row.curfew_mid11 = 1
        end
        if p == 5 && y >= 11
            row.curfew_mid12 = 1
        end
        if p == 6 && y >= 11     # Gwangju
            row.curfew_mid10 = 1; row.curfew_mid11 = 1; row.curfew_mid12 = 1
        end
        if p == 7 && y >= 11     # Ulsan
            row.curfew_mid12 = 1
        end
        if p == 8 && y >= 14     # Gyeonggi
            row.curfew_mid10 = 1
        end
        if p == 8 && y >= 11
            row.curfew_mid11 = 1; row.curfew_mid12 = 1
        end
        if p == 9 && y >= 15     # Gangwon
            row.curfew_mid11 = 1
        end
        if p == 9 && y >= 4
            row.curfew_mid12 = 1
        end
        if p == 10 && y >= 10    # Chungbuk
            row.curfew_mid11 = 1; row.curfew_mid12 = 1
        end
        if p == 11 && y >= 15    # Chungnam
            row.curfew_mid11 = 1
        end
        if p == 11 && y >= 11
            row.curfew_mid12 = 1
        end
        if p == 12 && y >= 15    # Jeonbuk
            row.curfew_mid10 = 1
        end
        if p == 12 && y >= 12
            row.curfew_mid11 = 1; row.curfew_mid12 = 1
        end
        if p == 13 && y >= 14    # Jeonnam
            row.curfew_mid10 = 1; row.curfew_mid11 = 1
        end
        if p == 13 && y >= 10
            row.curfew_mid12 = 1
        end
        if p == 14 && (y >= 3 && y <= 10)  # Gyeongbuk
            row.curfew_mid10 = 1
        end
        if p == 14 && y >= 3
            row.curfew_mid11 = 1; row.curfew_mid12 = 1
        end
        if p == 15 && y >= 16    # Gyeongnam
            row.curfew_mid11 = 1
        end
        if p == 15 && y >= 11
            row.curfew_mid12 = 1
        end
        if p == 16 && y >= 15    # Jeju
            row.curfew_mid11 = 1
        end
        if p == 16 && y >= 10
            row.curfew_mid12 = 1
        end
    end
    return df
end

function cluster_boot_se(fitfun, formula::FormulaTerm, df::DataFrame, term::Symbol;
                         cluster_var::Symbol = :province,
                         B::Int            = 999,
                         rng::Integer      = 2023)

    # 1) Initial fit to build the name→coef map
    m0    = fitfun(formula, df)
    names = coefnames(m0)
    vals  = coef(m0)
    coefmap0 = Dict(names[i] => vals[i] for i in eachindex(names))

    key = string(term)
    @assert haskey(coefmap0, key) "Term $term not found in $(names)"

    # 2) Set up cluster bootstrap
    clusters = unique(df[!, cluster_var])
    nclus    = length(clusters)
    draws    = Vector{Float64}(undef, B)
    rng      = MersenneTwister(rng)

    for b in 1:B
        # sample whole provinces
        samp = rand(rng, clusters, nclus)
        idx  = reduce(vcat, [findall(df[!, cluster_var] .== c) for c in samp])
        boot = df[idx, :]

        # re‐fit and rebuild the map
        mb     = fitfun(formula, boot)
        nmb    = coefnames(mb)
        vmb    = coef(mb)
        cmap_b = Dict(nmb[i] => vmb[i] for i in eachindex(nmb))

        # extract the same key
        draws[b] = cmap_b[key]
    end

    return std(draws; corrected=true)
end

function collect_boot!(triples, terms::Vector{Symbol};
                       cluster_var::Symbol = :province,
                       B::Int             = 999)

    out = Dict{Int, Dict{Symbol,Float64}}()
    for (i, (fitfun, f, df)) in enumerate(triples)
        sdmap = Dict{Symbol,Float64}()
        for t in terms
            sdmap[t] = cluster_boot_se(fitfun, f, df, t;
                                       cluster_var = cluster_var,
                                       B           = B)
        end
        out[i] = sdmap
    end
    return out
end

function coef_se_cell(c, s)
    c_rounded = round(c, digits=3)
    s_rounded = round(s, digits=3)
    return "\$$c_rounded\$ \\\\ \\small(\$$s_rounded\$)"
end

#----------------------------------------
# Set paths
#----------------------------------------
data_path= "/Users/fra/VS_CCA/Replication/PkgTemplates/data/"
output_path="/Users/fra/VS_CCA/Replication/PkgTemplates/output/"

#----------------------------------------
# Load and merge data
#----------------------------------------
# Master dataset
df = DataFrame(load(data_path*"data_clean.dta"))
# CPI for deflation
df_cpi = DataFrame(load(data_path*"cpi.dta"))
# Children
df_school = DataFrame(load(data_path*"children_school.dta"))

# Drop families with no children under head
filter!(:n_child_head => x -> x != 0, df)

# Merge school-age indicator (1:1 by hhid, year)
df = leftjoin(df, df_school, on = [:hhid, :year])
# Coerce missing child flags to zero and select middle-school
df.child_middle = coalesce.(df.child_middle, 0)
df.sel_school = df.child_middle .== 1

# carry forward last year’s values *within* each household
df = transform(
    groupby(df, :hhid, sort=true),
    :inc               => (x -> lead(x))               => :new_inc,
    :expenditure_month => (x -> lead(x))               => :new_expenditure_month,
    :educost_pub_month => (x -> lead(x))               => :new_educost_pub_month,
    :educost_priv_month=> (x -> lead(x))               => :new_educost_priv_month,
)

# overwrite the originals:
df.inc                .= df.new_inc
df.expenditure_month  .= df.new_expenditure_month
df.educost_pub_month  .= df.new_educost_pub_month
df.educost_priv_month .= df.new_educost_priv_month

# then drop the helpers
select!(df, Not([:new_inc, :new_expenditure_month, 
                 :new_educost_pub_month, :new_educost_priv_month]))

# Deflate prices
df = leftjoin(df, df_cpi, on = :year)
filter!(:cpi => x -> !ismissing(x), df)
df.inc .*= 100 ./ df.cpi
df.expenditure_month .*= 100 ./ df.cpi
df.educost_pub_month .*= 100 ./ df.cpi
df.educost_priv_month .*= 100 ./ df.cpi

# Compute per-child costs
df.educost_priv_cap     = df.educost_priv_month ./ df.n_child_head
df.educost_priv_cap_con = df.educost_priv_cap ./ df.expenditure_month

# Regional averages & logs
df.region_priv_edu_mid_top1 = df.micro_priv_edu_top1
df.region_priv_edu_mid_top2 = df.micro_priv_edu_top2
df.region_priv_edu_mid_top4 = df.micro_priv_edu_top4

for k in (1,2,4)
    df[!, Symbol("ln_region_priv_edu_mid_top$(k)")] = log.(df[!, Symbol("region_priv_edu_mid_top$(k)")])
end


# Curfew dummies
df.curfew_mid10 = zeros(Int, nrow(df))
df.curfew_mid11 = zeros(Int, nrow(df))
df.curfew_mid12 = zeros(Int, nrow(df))

assign_curfews!(df)

# Boolean mask for “both parents’ education ≤5 (with missing→99)”
mask = (coalesce.(df.edu_head_m, 99) .<= 5) .& (coalesce.(df.edu_head_f, 99) .<= 5)
# turn mask into 0/1
df.loweduc = Int.(mask)

# Average parents' age
df.avg_age_parents = map(
  (f, m) -> begin
    if ismissing(f) && ismissing(m)
      missing
    elseif ismissing(f)
      m
    elseif ismissing(m)
      f
    else
      (f + m)/2
    end
  end,
  df.age_head_f, df.age_head_m
)

# Recode into seven bins:
n = nrow(df)
df.age_group = Vector{Union{Missing,Int}}(missing, n)
# group 1: 16 ≤ avg ≤ 35
df.age_group[(df.avg_age_parents .>= 16) .& (df.avg_age_parents .<= 35)] .= 1
# group 2: 35.1 ≤ avg ≤ 40
df.age_group[(df.avg_age_parents .>  35) .& (df.avg_age_parents .<= 40)] .= 2
# group 3: 40.1 ≤ avg ≤ 45
df.age_group[(df.avg_age_parents .>  40) .& (df.avg_age_parents .<= 45)] .= 3
# group 4: 45.1 ≤ avg ≤ 50
df.age_group[(df.avg_age_parents .>  45) .& (df.avg_age_parents .<= 50)] .= 4
# group 5: 50.1 ≤ avg ≤ 55
df.age_group[(df.avg_age_parents .>  50) .& (df.avg_age_parents .<= 55)] .= 5
# group 6: 55.1 ≤ avg ≤ 60
df.age_group[(df.avg_age_parents .>  55) .& (df.avg_age_parents .<= 60)] .= 6
# group 7: 60.1 ≤ avg ≤ 95
df.age_group[(df.avg_age_parents .>  60) .& (df.avg_age_parents .<= 95)] .= 7

# Income quartiles by year & age_group
df.group = Vector{Union{Missing, Int}}(missing, nrow(df))

for y in 2:18
    for ag in 1:7
        mask_year = df.year .== y
        mask_age  = coalesce.(df.age_group .== ag, false)
        mask_inc  = .!ismissing.(df.inc)

        mask_all = mask_year .& mask_age .& mask_inc

        idx = findall(mask_all)

        if !isempty(idx)
            vals = df.inc[idx]
            qs   = quantile(vals, [0.25, 0.5, 0.75])
            for (i, v) in enumerate(vals)
                b = findfirst(q -> v <= q, qs)
                df.group[idx[i]] = isnothing(b) ? 4 : b
            end
        end
    end
end

# Low-income
n = nrow(df)
df.lowinc  = zeros(Int, n)
df.lowinc2 = zeros(Int, n)

mask1 = coalesce.(df.group .<= 1, false)
mask2 = coalesce.(df.group .<= 2, false)

df.lowinc[mask1]  .= 1
df.lowinc2[mask2] .= 1


#----------------------------------------
# Prepare regressions
#----------------------------------------
df.province   = categorical(df.province)
df.age_group  = categorical(df.age_group)
df.edu_head_m = categorical(df.edu_head_m)
df.edu_head_f = categorical(df.edu_head_f)

# ------------------------------------------------------------------
# LOW-INCOME 
# ------------------------------------------------------------------
 DEP_VAR  = :educost_priv_cap_con
 LEFT_VAR = :ln_region_priv_edu_mid_top2
 INST     = [:curfew_mid10, :curfew_mid11]
 CONTROLS = [:inc, :edu_head_m, :edu_head_f, :age_group]

sel_period = (df.year .>= 12) .& (df.year .<= 18)
df_lowinc  = df[df.lowinc2 .== 1 .& sel_period .& df.sel_school, :]

formula1     = Term(LEFT_VAR) ~ sum(Term.(INST)) + sum(Term.(CONTROLS)) + fe(:province)
formula2     = Term(LEFT_VAR) ~ sum(Term.(INST)) + sum(Term.(CONTROLS)) + Term.(:year) + fe(:province)
formula3     = Term(LEFT_VAR) ~ sum(Term.(INST)) + sum(Term.(CONTROLS)) + fe(:province) + fe(:year)

ols1         = reg(df_lowinc, formula1, Vcov.cluster(:province))
ols2         = reg(df_lowinc, formula2, Vcov.cluster(:province))
ols3         = reg(df_lowinc, formula3, Vcov.cluster(:province))

iv_formula1  = Term(DEP_VAR) ~ sum(Term.(CONTROLS)) + (Term(LEFT_VAR) ~ sum(Term.(INST))) + fe(:province)
iv_formula2  = Term(DEP_VAR) ~ sum(Term.(CONTROLS)) + (Term(LEFT_VAR) ~ sum(Term.(INST))) + Term.(:year) + fe(:province)
iv_formula3  = Term(DEP_VAR) ~ sum(Term.(CONTROLS)) + (Term(LEFT_VAR) ~ sum(Term.(INST))) + fe(:province) + fe(:year)

iv1          = reg(df_lowinc, iv_formula1, Vcov.cluster(:province))
iv2          = reg(df_lowinc, iv_formula2, Vcov.cluster(:province))
iv3          = reg(df_lowinc, iv_formula3, Vcov.cluster(:province))

# ------------------------------------------------------------------
# LOW-EDUCATION  
# ------------------------------------------------------------------
df_lowedu = df[df.loweduc .== 1 .& sel_period .& df.sel_school, :]

olsE1     = reg(df_lowedu, formula1, Vcov.cluster(:province))
olsE2     = reg(df_lowedu, formula2, Vcov.cluster(:province))
olsE3     = reg(df_lowedu, formula3, Vcov.cluster(:province))

ivE1      = reg(df_lowedu, iv_formula1, Vcov.cluster(:province))
ivE2      = reg(df_lowedu, iv_formula2, Vcov.cluster(:province))
ivE3      = reg(df_lowedu, iv_formula3, Vcov.cluster(:province))

ols_fit = (f,df) -> lm(f, df)
iv_fit  = (f,df) -> reg(df, f, Vcov.cluster(:province))

# ------------------------------------------------------------------- 
# Bootstrap SE
#--------------------------------------------------------------------
y = :ln_region_priv_edu_mid_top2
x_vars1 = [:curfew_mid10, :curfew_mid11, :inc, :edu_head_m, :edu_head_f, :age_group, :province]
x_vars2 = [:curfew_mid10, :curfew_mid11, :inc, :edu_head_m, :edu_head_f, :age_group, :province, :year]
x_vars3 = [:curfew_mid10, :curfew_mid11, :inc, :edu_head_m, :edu_head_f, :age_group, :province, :year]


df_lowinc.z .= 0  
df_lowedu.z .= 0

f1= Term(y) ~ sum(Term.(x_vars1)) 
f2= Term(y) ~ sum(Term.(x_vars2))
f3= Term(y) ~ Term(:z) + sum(Term.(x_vars3))

triples_lowinc_ols = [
  (ols_fit,    f1,    df_lowinc),
  (ols_fit,    f2,    df_lowinc),
  (ols_fit,    f3,    df_lowinc),
]

triples_lowinc_iv  = [
  (iv_fit,     iv_formula1, df_lowinc),
  (iv_fit,     iv_formula2, df_lowinc),
  (iv_fit,     iv_formula3, df_lowinc),
]


triples_lowedu_ols = [
  (ols_fit,    f1,    df_lowedu),
  (ols_fit,    f2,    df_lowedu),
  (ols_fit,    f3,    df_lowedu),
]

triples_lowedu_iv  = [
  (iv_fit,     iv_formula1, df_lowedu),
  (iv_fit,     iv_formula2, df_lowedu),
  (iv_fit,     iv_formula3, df_lowedu),
]

#Perform bootstrap
boot_lowinc    = collect_boot!(triples_lowinc_ols, INST;      cluster_var = :province, B = 999)
boot_lowinc_iv = collect_boot!(triples_lowinc_iv,  [LEFT_VAR]; cluster_var = :province, B = 999)

boot_lowedu    = collect_boot!(triples_lowedu_ols, INST;      cluster_var = :province, B = 999)
boot_lowedu_iv = collect_boot!(triples_lowedu_iv,  [LEFT_VAR]; cluster_var = :province, B = 999)


# -------------------------------------------------------------------
# TABLE 3
#--------------------------------------------------------------------


ols_models = [ols1, ols2, ols3, olsE1, olsE2, olsE3]
iv_models  = [iv1, iv2, iv3, ivE1, ivE2, ivE3]

boot_ols = boot_lowinc
boot_iv  = boot_lowinc_iv

beta_vals = [round(coef(m)[coefnames(m) .== "ln_region_priv_edu_mid_top2"][1], digits=3) for m in iv_models]
beta_se1 = [round(v[:ln_region_priv_edu_mid_top2], digits=3) for (k, v) in boot_lowinc_iv]
beta_se2 = [round(v[:ln_region_priv_edu_mid_top2], digits=3) for (k, v) in boot_lowedu_iv]

delta10_vals = [round(coef(m)[coefnames(m) .== "curfew_mid10"][1], digits=3) for m in ols_models]
delta10_se1 = [round(v[:curfew_mid10], digits=3) for (k, v) in boot_lowinc]
delta10_se2 = [round(v[:curfew_mid10], digits=3) for (k, v) in boot_lowedu]

delta11_vals = [round(coef(m)[coefnames(m) .== "curfew_mid11"][1], digits=3) for m in ols_models]
delta11_se1 = [round(v[:curfew_mid11], digits=3) for (k, v) in boot_lowinc]
delta11_se2 = [round(v[:curfew_mid10], digits=3) for (k, v) in boot_lowedu]


header = (
    ["", "", "", "", "", "", ""],
    ["Samples:", "", "Low income", "", "", "Low education", ""],
    ["", "(1)", "(2)", "(3)", "(4)", "(5)", "(6)"]
)

data = [
    "Panel A. 2nd stage"  ""       ""       ""       ""       ""       "";
    latexstring("\\beta")                  beta_vals[1]    beta_vals[2]    beta_vals[3]    beta_vals[4]    beta_vals[5]    beta_vals[6];
    "SE"                "($(beta_se1[1]))" "($( beta_se1[2]))" "($(beta_se1[3]))" "($(beta_se2[1]))" "($(beta_se2[2]))" "($(beta_se2[3]))";
    ""                   ""       ""       ""       ""       ""       "";
    "Panel B. 1st stage" ""       ""       ""       ""       ""       "";
    latexstring("\\delta_{10}")              delta10_vals[1]   delta10_vals[2]   delta10_vals[3]   delta10_vals[4]   delta10_vals[5]   delta10_vals[6];
    "SE"                "($(delta10_se1[1]))" "($(delta10_se1[2]))" "($(delta10_se2[1]))" "($(delta10_se2[1]))" "($(delta10_se2[2]))" "($(delta10_se2[3]))";
    latexstring("\\delta_{11}")              delta11_vals[1]   delta11_vals[2]   delta11_vals[3]   delta11_vals[4]   delta11_vals[5]   delta11_vals[6];
    "SE"                "($(delta11_se1[1]))" "($(delta11_se1[2]))" "($(delta11_se2[1]))" "($(delta11_se2[1]))" "($(delta11_se2[2]))" "($(delta11_se2[3]))";
    "Trendᵃ"            "No"     "Yes"     "No"     "No"    "Yes"    "No";
    "Year FE"           "No"     "No"     "Yes"     "Yes"    "No"    "Yes";
    "Observations"      ""     14031     ""     ""     12750     "";
]

hl = [
    :header,     
    4,           
    10,           
    13           
]



# Render table
latex_path = output_path*"table3.tex"

open(latex_path, "w") do io
    println(io, "\\begin{center}")
            println(io, "\\textbf{Table 3—Estimation of Private Education Spillovers Using Instrumental Variables}")
            println(io, "\\end{center}")
            println(io, "")
    
    pretty_table(io,
        data;
        header = header,
        backend = Val(:latex),
        title = "Table 3—Estimation of Private Education Spillovers Using Instrumental Variables",
        alignment = :l,
        tf = tf_latex_booktabs,
        hlines = hl,

    )


end


end
end



