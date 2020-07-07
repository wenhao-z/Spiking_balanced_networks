#rand(10)
# using Distributions
# expdist = Exponential()
# print(Exponential(2))
# print(rand(expdist))

# print(zeros(3))
# using Debug
# @debug function test()
#     parts = {}
#     @bp
#     for j=1:3
#         for i=1:3
#             push!(parts,"($i,$j) ")
#         end
#     end
#     @bp
#     println(parts...)
# end
#
# test()

for iter = 1: 5000
  clf
  plot(VArray[iter,:])
  pause
  drawnow
end
