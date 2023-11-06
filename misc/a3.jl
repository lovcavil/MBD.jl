using Plots

# Define the function
f(x) = x^2 - 2x + 2

# Generate x values
xvals = -2:0.01:4  # Change this range if you want a different domain for x

# Plot the function
plot(xvals, f.(xvals), label="f(x) = x^2 - 2x + 2", xlabel="x", ylabel="f(x)", title="Plot of f(x) = x^2 - 2x + 2")