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
        do_figure1(main_dir)
        do_table1(main_dir)
        do_table2(main_dir)
        do_table3(main_dir)
        do_table45(main_dir)
    end
end

