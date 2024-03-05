using BenchmarkTools
@benchmark sort(data) setup=(data=rand(10))