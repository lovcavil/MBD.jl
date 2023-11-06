function compute_nilpotency(A::Array{T}) where T
    m = 0
    B = copy(A)
    while !all(iszero, B)
        B *= A
        m += 1
        if m > size(A, 1)  # 超过矩阵的维度，矩阵不是nilpotent
            return nothing
        end
    end
    return m
end

# 示例矩阵
A = [0 1 0; 0 2 0; 1 0 0]
nilpotency = compute_nilpotency(A)
println("The nilpotency of A is: ", nilpotency)