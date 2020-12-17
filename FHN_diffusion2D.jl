## ODE and PDE discretization.
@time begin
using BlockBandedMatrices
using BandedMatrices
using DifferentialEquations
using DiffEqOperators
function FHN_2D(u,p)
    v,w = u
    ϵ,a,b = p
    F = similar(u)

    F[1] = ϵ*(w-a*v)
    F[2] = -v + b*w - w^3

    return F
end

function FHN_pde_2D(du,u,p,t)
    D = Diagonal([0.05,0.00028])
    Δx = 1.
    N = 40
    dv = range(-2,-2,length=N)
    ev = range(1,1,length=N-1)
    inM = SymTridiagonal(dv, ev)
    inM = Array(inM)

    I_N = Matrix(1.0I, N,N)
    outM = BlockTridiagonal(fill(I_N,N-1),fill(inM,N),fill(I_N,N-1))
    outM = BandedMatrix(outM)./Δx^2

    B = similar(u)
    for i = 1:size(u,1)
        B[i,:] = FHN_2D(u[i,:],p)
    end

    du .= outM*u*D + B
end
## Initial conditions.
N = 40
u0 = zeros(N^2,2)
u0 = u0 + randn(N^2,2)
# u0[:,1] .= -0.6
# u0[:,2] .= -0.3
# u0[Int64(N/2-10):Int64(N/2+10),1] .= 1.
# u0
## Solve and plot.
# alg = ORK256()
# alg = Trapezoid()
using LSODA
alg = lsoda()
tstep = 0.1
tspan = (0.,10.)
p1 = [10.,1.,1.]
prob = ODEProblem(FHN_pde_2D,u0,tspan,p1)
sol = solve(prob,alg,dtmax=tstep);

# uxy = sol[end][:,2]
# uxyt =  reshape(uxy,(N,N))
# fig1 = heatmap(uxyt,c=:lightrainbow)



## Animation
anim = @animate for i ∈ 1:length(sol)
    uxy = sol[i][:,2]
    uxyt =  reshape(uxy,(N,N))
    heatmap(uxyt,c=:lightrainbow)
end
gif(anim, "Project_mmi_2\\figures\\anim1_fps15.gif", fps = 15)

# @gif for i ∈ 1:length(sol)
#     uxy = sol[i][:,2]
#     uxyt =  reshape(uxy,(N,N))
#     heatmap(uxyt,c=:lightrainbow)
# end

end
