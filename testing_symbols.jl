begin
    using DifferentialEquations
    using Plots
end

v = [-1., +1., 0.1, 0.2]
rhs = (dx,x,p,t) -> dx .= v

syms = [Symbol(:omega, "_" ,i) for i in 1:2]

syms = [Symbol(:omega, "_" ,1), Symbol(:phi, "_" ,1), Symbol(:omega, "_" ,2), Symbol(:phi, "_" ,2)]
syms = [:omega, :phi, :omega,:phi]

rhs_ODE = ODEFunction(rhs, syms=syms)

ic = rand(4)

prob = ODEProblem(rhs_ODE,ic,(0.,5.))

sol = solve(prob)

plot(sol, vars=[:omega])
