#Test the functions used to build the model (to be completed)

using Pkg
Pkg.add("TestItemRunner")
Pkg.add("Test")

include("/Users/fra/VS_CCA/Replication/PkgTemplates/src/Custom_fct.jl")

using TestItemRunner, Test, .Custom

@testitem "HKprodfn functioning" begin
    # Test the HKprodfn function
    @test isapprox(HKprodfn(3, 4 , 2.3 , 0.5, 2, 0.2), 17.3700, atol=1e-4)
    @test isapprox(HKprodfn(0.8, 2, 0.1 , 2, 0.5, 0.2), 1.8529, atol=1e-4)
    @test isapprox(HKprodfn(1,0.6, 0.3, 4, 0.3, 2), 4.6968, atol=1e-4)
end