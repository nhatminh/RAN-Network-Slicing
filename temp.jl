# using Convex
# using ECOS
#
# x = Variable()
# # problem = minimize(max(x,2*x-1, 3*x-3), [x>=0, x<=3])
# # solve!(problem)
# # println(x.value)
# # x = Variable()
# problem = minimize(invpos(3-x), [x>=1, x<=2])
# solve!(problem)
# println(x.value)
#
# function exhaustic_search()
#   for beta =0.1:0.1:1
#     println("a")
#   end
# end
#
# # exhaustic_search()
# a = 1
# println("a $a")

# using DataFrames
# b=[1,2]
# println(size(b)[1])
#
# using Seaborn
#
# x = rand(1000)
# # df = convert(DataFrame, x)
# # rename!(df, [:x1], [:Value])
# ones(Int64,500)
# t = vcat(1*ones(Int64,500),2*ones(Int64,500))
# df1 = DataFrame(Val = x, Type = t)
# df1[[:Val,:Type]]
# # Seaborn.boxplot(x="Type",y="Val",data=df1)
# # Seaborn.boxplot(df1[[:Val,:Type]])
# Seaborn.seaborn[:boxplot](df1[[:Val,:Type]])
# using PyPlot
# using DataFrames
# # x = rand(1000, 4)
# # df = convert(DataFrame, x)
# # #################
# # #  Create Data  #
# # #################
# x = randn(1000) # Values
# nbins = 50 # Number of bins
#
# ##########
# #  Plot  #
# ##########
# fig = figure("pyplot_histogram",figsize=(10,10)) # Not strictly required
#
# h = plt[:hist](x,nbins) # Histogram
#
# grid("on")
# xlabel("X")
# ylabel("Y")
# title("Histogram")

arr = rand(18:20,4)
N_UE_s=[20,19,19,18]
println(sum(N_UE_s[1:3]))
println(sum(N_UE_s))
print(size(arr))
print(size(N_UE_s))

style.use('ggplot')
