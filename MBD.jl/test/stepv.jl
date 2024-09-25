using DifferentialEquations

# Define the system of ODEs
function f!(du, u, p, t)
    du[1] = u[2]
    du[2] = -u[1]
    du[3] = u[1] * u[2]
end

# Initial conditions and time span
u0 = [1.0, 0.0, 0.0]
tspan = (0.0, 10.0)

# Create the problem
prob = ODEProblem(f!, u0, tspan)

# Initialize the integrator
integrator = init(prob, Tsit5())

# Store the endpoint
t_end = prob.tspan[2]

# Loop through the integration steps
while integrator.t < t_end
    step!(integrator)
    
    # Update u[1] if t > 5
    if integrator.t > 5
        integrator.u[1] = 100.0
    end

    println("Time: $(integrator.t), Solution: $(integrator.u)")
end
