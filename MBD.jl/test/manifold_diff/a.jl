using ManifoldDiffEq, OrdinaryDiffEq, Manifolds
using GLMakie, LinearAlgebra, Colors

n = 25

θ = [0;(0.5:n-0.5)/n;1]
φ = [(0:2n-2)*2/(2n-1);2]
x = [cospi(φ)*sinpi(θ) for θ in θ, φ in φ]
y = [sinpi(φ)*sinpi(θ) for θ in θ, φ in φ]
z = [cospi(θ) for θ in θ, φ in φ]

function f2(x, y, z)
    Iv = [1.6, 1.0, 2/3]
    p = [x, y, z]
    A = [0 -z y; z 0 -x; -y x 0]
    return A * (p./Iv)
end

tans = f2.(vec(x), vec(y), vec(z))
u = [a[1] for a in tans]
v = [a[2] for a in tans]
w = [a[3] for a in tans]

f = Figure();
Axis3(f[1,1])

arr = GLMakie.arrows!(
           vec(x), vec(y), vec(z), u, v, w;
           arrowsize = 0.02, linecolor = (:gray, 0.7), linewidth = 0.0075, lengthscale = 0.1
)
save("first_example_vector_field.png", f)