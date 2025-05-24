module sbrana
    export run

    include("Figure1.jl")
    include("Table1.jl")
    include("Table2.jl")
    include("Table3.jl")
    include("Table45.jl")

    using .figure1, .table1, .table2, .table3, .ProgMainSS
    using StatsModels, Pkg

    function run()
        main_dir = "/Users/fra/VS_CCA/Replication/PkgTemplates"
        Pkg.add("TestItemRunner")
        Pkg.add("Test")
        Pkg.add("FreqTables")
        Pkg.add("ShiftedArrays")
        Pkg.add("DataFrames")
        Pkg.add("StatFiles")
        Pkg.add("StatsBase")
        Pkg.add("Statistics")
        Pkg.add("StatsModels")
        Pkg.add("CSV")
        Pkg.add("FixedEffectModels")
        Pkg.add("CategoricalArrays")
        Pkg.add("PrettyTables")
        Pkg.add("Econometrics")
        Pkg.add("RegressionTables")
        Pkg.add("WildBootTests")
        Pkg.add("Distributions")
        Pkg.add("LaTeXStrings")
        Pkg.add("GLM")
        Pkg.add("DataFramesMeta")
        Pkg.add("Plots")
        Pkg.add("XLSX")
        Pkg.add("LinearAlgebra")
        Pkg.add("Random")
        Pkg.add("Printf")
        do_figure1(main_dir)
        do_table1(main_dir)
        do_table2(main_dir)
        do_table3(main_dir)
        do_table45(main_dir)
    end
end

using FreqTables, LinearAlgebra, Printf, DataFramesMeta, Plots, DataFrames, StatFiles, StatsBase, CSV
using FixedEffectModels, CategoricalArrays, ShiftedArrays, GLM, Random, Statistics, Econometrics
using RegressionTables, WildBootTests, Distributions, PrettyTables, LaTeXStrings
using StatsModels: coefnames 
using Random: MersenneTwister


