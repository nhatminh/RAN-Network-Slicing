using PyPlot
using Distributions

###################
##  Create Data  ##
###################
n = 100
x = linspace(-3, 3, n)
y = linspace(-3,3,n)

xgrid = repmat(x',n,1)
ygrid = repmat(y,1,n)

z = zeros(n,n)

for i in 1:n
    for j in 1:n
        z[i:i,j:j] = pdf(MvNormal(eye(2)),[x[i];y[j]])
    end
end

############
##  Plot  ##
############
fig = figure("pyplot_surfaceplot",figsize=(10,10))
ax = fig[:add_subplot](2,1,1, projection = "3d")
ax[:plot_surface](xgrid, ygrid, z, rstride=2,edgecolors="k", cstride=2, cmap=ColorMap("gray"), alpha=0.8, linewidth=0.25)
xlabel("X")
ylabel("Y")
title("Surface Plot")

subplot(212)
ax = fig[:add_subplot](2,1,2)
cp = ax[:contour](xgrid, ygrid, z, colors="black", linewidth=2.0)
ax[:clabel](cp, inline=1, fontsize=10)
xlabel("X")
ylabel("Y")
title("Contour Plot")
tight_layout()


using PyPlot

######################################
#  Generate an hour of data at 10Hz  #
######################################
x = [DateTime(2013,10,4):Dates.Millisecond(100):DateTime(2013,10,4,1);] # Generate time array
x = map(Float64,x)/1000/60/60/24 # Convert time from milliseconds from day 0 to days from day 0
y = sin(2*pi*collect(0:2*pi/length(x):2*pi-(2*pi/length(x))))
dx = maximum(x) - minimum(x)
dy = maximum(y) - minimum(y)

y2 = 30*(1+sin(2*pi*collect(pi:2*pi/length(x):3*pi-(2*pi/length(x)))))-10
x2 = [minimum(x):dx/20:maximum(x);]
y2 = 10rand(21)-3
x3 = [minimum(x):dx/20:maximum(x);]
y3 = 10rand(21)-3

##########
#  Plot  #
##########
fig = figure("pyplot_annotation",figsize=(10,10)) # Create a figure and save its handle
x_start = 1
y_start = 1
annotate("Starting Point",
          xy=[x_start;y_start],# Arrow tip
          xytext=[0;0], # Text offset from tip
          xycoords="data",# Coordinates in in "data" units
          arrowprops=Dict("facecolor"=>"black")) #
