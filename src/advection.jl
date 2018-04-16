abstract type AbstractAdvection end

struct UpWind2ndOrder <: AbstractAdvection
end

function (a::UpWind2ndOrder)(rhs, ρC, buvec, uvec, Tc, Tbf, ∇Tc, bf2c, bfn, bcv, f2c, fn, cv, fc, ccenter, floops::Face2CellLoop)
    NBF = nbfaces(floops)
    @inbounds for i=1:NBF
        j = bf2c[i][1]
        rhs[j] -= (ρC * Tbf[i] * (buvec[i] ⋅ bfn[i]))*bcv[i]
    end
    NF = nfaces(floops)
    @inbounds for i=1:NF
        j1, j2 = f2c[i]
        n = fn[i]
        v1, v2 = cv[i]
        udotn = uvec[i] ⋅ n
        if udotn >= 0
            Tf = Tc[j1] + (fc[i] - ccenter[i][1])⋅∇Tc[j1]
        else
            Tf = Tc[j2] + (fc[i] - ccenter[i][2])⋅∇Tc[j2]
        end
        qf = ρC*Tf*udotn
        rhs[j1] -= qf*v1
        rhs[j2] += qf*v2 
    end

end

function (a::UpWind2ndOrder)(rhs,p)
    floops = p.mesh.f2cloops
    a(rhs, p.ρC, p.u, p.u, p.Tc, p.Tbf, p.∇Tc, floops.bf2c, floops.bfn, floops.bcv, floops.f2c, floops.fn, floops.cv, floops.fc, floops.ccenter, floops)
end

struct FullUpWind <: AbstractAdvection
end

function (a::FullUpWind)(rhs, ρC, buvec, uvec, Tc, Tbf, bf2c, bfn, bcv, f2c, fn, cv, fc, ccenter, floops::Face2CellLoop)
    NBF = nbfaces(floops)
    @inbounds for i=1:NBF
        j = bf2c[i][1]
        rhs[j] -= (ρC * Tbf[i] * (buvec[i] ⋅ bfn[i]))*bcv[i]
    end
    NF = nfaces(floops)
    @inbounds for i=1:NF
        j1, j2 = f2c[i]
        n = fn[i]
        v1, v2 = cv[i]
        udotn = uvec[i] ⋅ n
        if udotn >= 0
            Tf = Tc[j1]
        else
            Tf = Tc[j2]
        end
        qf = ρC*Tf*udotn
        rhs[j1] -= qf*v1
        rhs[j2] += qf*v2 
    end

end

function (a::FullUpWind)(rhs,p)
    floops = p.mesh.f2cloops
    a(rhs, p.ρC, p.u, p.u, p.Tc, p.Tbf, floops.bf2c, floops.bfn, floops.bcv, floops.f2c, floops.fn, floops.cv, floops.fc, floops.ccenter, floops)
end
