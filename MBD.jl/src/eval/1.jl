using Dierckx
using Plots

# Your data points
x = [0.0, 1.0, 2.0, 5.0, 10.0]
y = [1.0, 1.0, 1.0, 1.25, 1.0]

# Create the spline
spl = Spline1D(x, y)

# Generate a finer set of x values for plotting the spline smoothly
x_fine = range(minimum(x), maximum(x), length=300)

# Evaluate the spline at the finer set of x values
y_fine = [spl(x) for x in x_fine]

# Plot the original data points
a=plot(x, y, seriestype=:scatter, label="Data Points")

# Plot the spline curve
plot!(x_fine, y_fine, label="Spline Interpolation")

# Display the plot
display(a)
