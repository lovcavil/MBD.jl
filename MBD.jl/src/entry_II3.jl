include("mathfunction.jl")
include("AppData.jl")

# AA Spatial Kinematic Analysis
# User Input
qtol = 1e-6     # Tolerance in solving for q
maxiter = 25    # Maximum number of iterations in Newton-Raphson iteration
maxCond = 1e4   # Maximum condition number for Jacobian
h = 0.001       # Time step
tfinal = 2 * π/10
# Application Data
# app=1, Spatial 2-Bar with 2 rotational drivers
# app=2, 4-Bar with 1 rotational driver
# app=3, Spatial Slider-Crank, RotDr Bod1-grnd
# app=4, Fly-Ball Governor
# app=5, Excavator
app = 1
nb, ngc, nh, nhc, nc, nd, SJDT, q0e = AppData(app)    # Data from AppData function
par = Any[nb, ngc, nh, nhc, nd, qtol, app]     # Parameter vector for analysis control

# Make vectors available to all functions
global ux = [1, 0, 0]
global uy = [0, 1, 0]
global uz = [0, 0, 1]
global z3 = zeros(3)


# Data Storage Arrays
# Preallocate a large array with more columns than currently needed.
# For example, 100 columns as a starting point:
Q = zeros(ngc, 100)
Qd = zeros(ngc, 100)
Qdd = zeros(ngc, 100)
filled_columns = 0

Phiq0e = mathfunction.PhiqEval(q0e, SJDT, par)
CondPhiq0e = cond(Phiq0e)

if app==1
    # Define the arrays to store the output data of interest
    y2 = Float64[]
    y2d = Float64[]
    y2dd = Float64[]
    x2 = Float64[]
    x2d = Float64[]
    x2dd = Float64[]
    x2p = Float64[]
    y2p = Float64[]
    z2p = Float64[]
    x2pd = Float64[]
    y2pd = Float64[]
    z2pd = Float64[]
    x2pdd = Float64[]
    y2pdd = Float64[]
    z2pdd = Float64[]
end

# Kinematic Analysis
n = 1    # Time step counter
t = []
tn = 0.0
while tn < tfinal
    global t,tn,n,filled_columns,Q,Qd,Qdd

    # When you need to add a new column:
    if filled_columns == size(Q, 2)
        # If the array is full, double its size.
        Q = hcat(Q, zeros(ngc, size(Q, 2)))
        Qd = hcat(Qd, zeros(ngc, size(Qd, 2)))
        Qdd = hcat(Qdd, zeros(ngc, size(Qdd, 2)))
    end


    filled_columns += 1
    t = push!(t, tn)
    # q-estimate
    if n == 1
        q = q0e    
    else
        q = Q[:, n-1]         
    end
    if n-1 > 1
        q = q + h * Qd[:, n-1] + 0.5 * h^2 * Qdd[:, n-1]
    end
    i = 1    # Position iteration counter
    err = qtol + 1
    errrpt = []
    while err > qtol
        Phi = mathfunction.PhiEval(tn, q, SJDT, par)
        Phiq = mathfunction.PhiqEval(q, SJDT, par)
        if i == 1
            CondPhiq = cond(Phiq)
            if CondPhiq > maxCond   # Check for ill-conditioned constraint Jacobian
                println("Warning: Constraint Jacobian Ill Conditioned")
                break
            end
        end
        delq = -Phiq \ Phi
        q = q + delq
        err = norm(Phi)
        push!(errrpt, err)
        i += 1
        if i > maxiter   # Check for failure to converge in Newton-Raphson
            println("Warning: Newton-Raphson convergence failure")
            break
        end
    end    
    iter = i - 1    # Report number of iterations
    Q[:, n] = q
    # Evaluate qd
    Phiq = mathfunction.PhiqEval(q, SJDT, par)
    P, Pst, Pstt = mathfunction.P5Eval(tn, q, SJDT, par)
    qd = -Phiq \ Pst
    Qd[:, n] = qd
    # Evaluate qdd
    P2 = mathfunction.P2Eval(q, qd, SJDT, par)
    Gam = P2 * qd + Pstt
    qdd = -Phiq \ Gam
    Qdd[:, n] = qdd


    # ... (Rest of the application-specific output calculations)

    # Calculate output data of interest (Enter for each application)
    if app == 1  # 2 Bar with 2 rotational drivers
        global y2, y2d, y2dd, x2, x2d, x2dd, x2p, y2p, z2p, x2pd, y2pd, z2pd, x2pdd, y2pdd, z2pdd 
        push!(y2, q[9])
        push!(y2d, qd[9])
        push!(y2dd, qdd[9])
        push!(x2, q[8])
        push!(x2d, qd[8])
        push!(x2dd, qdd[8])
        r1, p1 = mathfunction.qPart(q, 1)
        r2, p2 = mathfunction.qPart(q, 2)
        r1p = mathfunction.ATran(p1) * ux
        r2p = r1p + mathfunction.ATran(p2) * ux
        r1d, p1d = mathfunction.qPart(qd, 1)
        r2d, p2d = mathfunction.qPart(qd, 2)
        r1pd = mathfunction.BTran(p1, ux) * p1d
        r2pd = r1pd + mathfunction.BTran(p2, ux) * p2d
        r1dd, p1dd = mathfunction.qPart(qdd, 1)
        r2dd, p2dd = mathfunction.qPart(qdd, 2)
        r1pdd = mathfunction.BTran(p1, ux) * p1dd + mathfunction.BTran(p1d, ux) * p1d
        r2pdd = r1pdd + mathfunction.BTran(p2, ux) * p2dd + mathfunction.BTran(p2d, ux) * p2d
        push!(x2p, r2p[1])
        push!(y2p, r2p[2])
        push!(z2p, r2p[3])
        push!(x2pd, r2pd[1])
        push!(y2pd, r2pd[2])
        push!(z2pd, r2pd[3])
        push!(x2pdd, r2pdd[1])
        push!(y2pdd, r2pdd[2])
        push!(z2pdd, r2pdd[3])
    end


    n += 1
    tn += h
end

using Plots

# Assuming t is your time vector and y2, x2, etc. are as defined previously
f=plot(t, y2, label="y2", title="Output Data Plot", xlabel="Time", ylabel="y2", legend=:topright)
plot!(t, x2, label="x2")



display(f)
# To plot position, velocity, and acceleration in subplots
p1 = plot(t, x2p, label="x2 Position")
p2 = plot(t, x2pd, label="x2 Velocity")
p3 = plot(t, x2pdd, label="x2 Acceleration")

# Combine into a single figure with 3 vertically aligned subplots
f2=plot(p1, p2, p3, layout=(3, 1), legend=:topright)
display(f2)

f2_5=plot(x2pd, x2pdd, label="x2pdx2pdd", title="Output Data Plot", xlabel="x2pd", ylabel="x2pdd", legend=:topright)
display(f2_5)
# To visualize the 3D trajectory of a point
f3=scatter3d(x2p, y2p, z2p, label="3D Trajectory", xlabel="x2p", ylabel="y2p", zlabel="z2p")
display(f3)
# To save plots to files
savefig("output_data_plot.png") # Saves the last plot displayed

using CSV, DataFrames

# Assuming x2pdd is a vector. If it's not, you might need to reshape or adjust it.
df = DataFrame(X2pdd = x2pdd)
CSV.write("x2pdd.csv", df)

println("end")
# ... Add more plots for other variables as needed
