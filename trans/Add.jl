# 使用循环的方法
function add_to_submatrix_loop!(A, B, m, n)
    r, s = size(B)
    for i in 1:r
        for j in 1:s
            A[m+i, n+j] += B[i, j]
        end
    end
    return A
end

# 使用广播的方法
function add_to_submatrix_broadcast!(A, B, m, n)
    A[m+1:m+size(B, 1), n+1:n+size(B, 2)] .+= B
    return A
end

function add_to_submatrix_view!(A, B, m, n)
    r, s = size(B)
    sub_A = view(A, m+1:m+r, n+1:n+s)
    sub_A .+= B
    return A
end

function add_to_submatrix_inbounds!(A, B, m, n)
    r, s = size(B)
    @inbounds begin
        for i in 1:r
            for j in 1:s
                A[m+i, n+j] += B[i, j]
            end
        end
    end
    return A
end

function add_to_submatrix_inbounds_view!(A, B, m, n)
    r, s = size(B)
    sub_A = view(A, m+1:m+r, n+1:n+s)
    @inbounds for i in 1:r
        for j in 1:s
            sub_A[i, j] += B[i, j]
        end
    end
    return A
end

using BenchmarkTools

# 创建一些测试数据
A = rand(10000, 10000)
B = rand(200, 200)
m = n = 500

# 测试循环方法
println("1")
@btime add_to_submatrix_loop!(copy($A), $B, $m, $n)

# 测试广播方法
println("2")
@btime add_to_submatrix_broadcast!(copy($A), $B, $m, $n)
# 测试广播方法3
println("3")
@btime add_to_submatrix_view!(copy($A), $B, $m, $n)
# 测试广播方法4
println("4")
@btime add_to_submatrix_inbounds!(copy($A), $B, $m, $n)
println("5add_to_submatrix_inbounds_view!")
@btime add_to_submatrix_inbounds_view!(copy($A), $B, $m, $n)