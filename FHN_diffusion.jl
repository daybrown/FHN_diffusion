## ODE and PDE discretization.
using BandedMatrices
using DifferentialEquations
using DiffEqOperators
function FHN!(F,u,p)
    x,y = u
    ϵ,a,b = p

    F[1] = x-x^3-y
    F[2] = ϵ*(x-b*y+a)
end

function FHN(u,p)
    x,y = u
    ϵ,a,b = p
    F = similar(u)

    F[1] = x-x^3-y
    F[2] = ϵ*(x-b*y+a)

    return F
end

function FHN_pde(du,u,p,t)
    D = 0
    order = 2
    deriv = 2
    Δx = 1.
    N = 400
    A = BandedMatrix(CenteredDifference(order, deriv, Δx, N))
    A = A[:,2:end-1]

    B = similar(u)
    for i = 1:size(u,1)
        B[i,:] = FHN(u[i,:],p)
    end

    du .= D.*(A*u) + B
end
## Initial conditions.
N = 400
u0 = zeros(N,2)
u0[:,1] .= -0.6
u0[:,2] .= -0.3
u0[Int64(N/2-10):Int64(N/2+10),1] .= 1.
u0
## Solve and plot.
# alg = ORK256()
alg = Trapezoid()
tstep = 0.05
tspan = (0.,400.)
p1 = [0.01,0.1,2.]
prob = ODEProblem(FHN_pde,u0,tspan,p1)
sol = solve(prob,alg,dtmax=tstep);
uxt = zeros(N,length(sol))
for i = 1:N
    uxt[i,:] = getindex.(sol[:],i)
end
fig1 = heatmap(uxt[41:end-40,:],c=:lightrainbow)

p1 = [0.01,0.1,1.5]
prob = ODEProblem(FHN_pde,u0,tspan,p1)
sol = solve(prob,alg,dtmax=tstep);
uxt = zeros(N,length(sol))
for i = 1:N
    uxt[i,:] = getindex.(sol[:],i)
end
fig2 = heatmap(uxt[41:end-40,:],c=:lightrainbow)

p1 = [0.01,0.1,1.]
prob = ODEProblem(FHN_pde,u0,tspan,p1)
sol = solve(prob,alg,dtmax=tstep);
uxt = zeros(N,length(sol))
for i = 1:N
    uxt[i,:] = getindex.(sol[:],i)
end
fig3 = heatmap(uxt[41:end-40,:],c=:lightrainbow)


lo = @layout [a{0.5w} b{0.5w} c{0.5w}]
fig = plot(fig1,fig2,fig3,layout = 3)
