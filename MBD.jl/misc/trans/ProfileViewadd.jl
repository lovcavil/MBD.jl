#import Pkg; Pkg.add("ProfileView")
#using Pkg
#Pkg.add("BlockDiagonals")

using BlockDiagonals

bm = BlockDiagonal([rand(2, 3), ones(3, 2)])



using Pkg
Pkg.add("IterativeRefinement")