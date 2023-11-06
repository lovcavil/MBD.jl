include("..\\src\\mathfunction.jl")

# Begin a test set for the add_to_submatrix! function
@testset "add_to_submatrix! Tests" begin
    A = [1 2; 3 4]
    B = [5 6; 7 8]
    m = 1
    n = 1
    
    expected = [1+5 2+6; 3+7 4+8]
    result = add_to_submatrix!(A, B, m, n)
    @test result == expected

    A = [1 2 3; 4 5 6; 7 8 9]
    B = [10 11; 12 13]
    m = 2
    n = 2
    
    expected = [1 2 3; 4 5+10 6+11; 7 8+12 9+13]
    result = add_to_submatrix!(A, B, m, n)
    @test result == expected
end
