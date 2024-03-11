function create_updater()
    initial_v = 0 # 初始值
    prev_q = 0    # 上一次的q值

    # 定义并返回闭包，该闭包接受q和v作为输入
    return (q, v) -> begin
        # 如果prev_q <= 0 且当前q > 0，则更新initial_v
        if prev_q <= 0 && q > 0
            initial_v = v
        end
        prev_q = q # 更新prev_q为当前的q值

        # 返回当前有效的initial_v
        return initial_v
    end
end

# 使用闭包创建一个更新器
updater = create_updater()

# 测试不同的q和v值
println(updater(-1, 5))  # q <= 0，不更新initial_v，输出初始的0
println(updater(1, 10))  # prev_q <= 0 且当前q > 0，更新initial_v为10，输出10
println(updater(2, 20))  # q > 0，但initial_v已在上一步更新，输出10
println(updater(-2, 30)) # q <= 0，不更新，输出10
println(updater(3, 40))  # prev_q <= 0 且当前q > 0，更新initial_v为40，输出40