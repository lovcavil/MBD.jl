### A Pluto.jl notebook ###
# v0.19.37

using Markdown
using InteractiveUtils

# ╔═╡ 4fd7e7d2-0a31-4be7-b651-da1817262a0e
struct MyStruct
    field::Union{Int, Nothing}
end

# ╔═╡ 2afed058-5ee6-4192-8f4b-2c88b5a3dc09
my_struct_instance = MyStruct(nothing)

# ╔═╡ f1b7ad86-2ee4-420a-8d37-f40400d35d2b
(0.0, 5)


# ╔═╡ 650ee2e6-7fc6-49d9-aa21-097d74c3ef0d
typeof((0.0, 5))

# ╔═╡ Cell order:
# ╠═4fd7e7d2-0a31-4be7-b651-da1817262a0e
# ╠═2afed058-5ee6-4192-8f4b-2c88b5a3dc09
# ╠═f1b7ad86-2ee4-420a-8d37-f40400d35d2b
# ╠═650ee2e6-7fc6-49d9-aa21-097d74c3ef0d
