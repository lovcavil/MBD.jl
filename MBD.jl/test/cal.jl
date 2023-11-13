function subfunction1(n)
    sleep(1)  # Simulating some computation
    return sum([i for i in 1:n])
end

function subfunction2(n)
    sleep(2)  # Simulating more computation
    return prod([i for i in 1:n])
end

function main_function(n)
    a = subfunction1(n)
    b = subfunction2(n)
    return a + b
end

using Profile

@profile main_function(10)  # Profiling the main function
Profile.print(format=:tree)
using ProfileView
ProfileView.view()
