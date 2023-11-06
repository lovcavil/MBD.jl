module MyPackage

# 使用export来导出函数，使它们可以被外部调用
export function1, function2

# 包含子模块文件
include("mathfunction.jl")


# 包内的函数定义
function function1()
    # ...
end

function function2()
    # ...
end

end # module