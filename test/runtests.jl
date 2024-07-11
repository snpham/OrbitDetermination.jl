using OrbitDetermination
using Test

@testset "OrbitDetermination.jl" begin
    # Write your tests here.

    @test OrbitDetermination.greet_your_package_name() == "Hello OrbitDetermination!"
    @test OrbitDetermination.greet_your_package_name() != "Hello world!"

end
