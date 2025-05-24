
module Package
    export inst

    function inst()
        #using Pkg
        Pkg.add("Distributions")
        Pkg.add("XLSX")
        Pkg.add("DataFrames")
        Pkg.add("PrettyTables")
        Pkg.add("CSV")
    end

end