# -*- coding: utf-8 -*-
"""
Created on June 15, 2024

@author: thanhqtran

AK60.jl
deterministic simple direct computation of AK60
based on Heer (2009)'s book [Ch.9] 
"""


using Plots
using NLsolve

# ========================
# == Parameters ==========
# ========================
beta = 0.98
gamma = 2
alpha = 0.3
delta = 0.1
eta = 2
repl = 0.3
T = 40
TR = 20
tau = repl / (2 + repl)
psi = 0.001
r = 0.045
phi = 0.5  # Learning rate for nbar update


# steady state distribution
function update_w(k, n)
    return (1 - alpha) * k^alpha * n^(-alpha)
end

function update_kbar(n, r)
    return n * (alpha / (r + delta))^(1 / (1 - alpha))
end


# ==========================================================
# == Decision Solvers ======================================
# == these functions solve the decision rules for ks and ns=
# ==========================================================

# 1. Going to use for s = 59, 58, ..., 41
function solve_ks_old(ks1, ks2)
    ks = (1 / (1 + r)) * ((beta * (1 + r))^(1 / eta) * ((1 + r) * ks1 + b - ks2 + psi) - (b - ks1 + psi))
    return ks
end

# 2. At s = 40, given ks41 and ns41, solve ks40 and ns40
function fs40!(F, X, ks1, ks2)
    ns1 = 0
    ks, ns = X
    F[1] = (1 - tau) * w * (1 - ns) / gamma - ((1 + r) * ks + (1 - tau) * w * ns - ks1 + psi)
    F[2] = 1 / beta - (((1 + r) * ks1 + (1 - tau) * w * ns1 - ks2 + psi) / ((1 + r) * ks + (1 - tau) * w * ns - ks1 + psi))^(-eta) * (((1 - ns1) / (1 - ns))^(gamma / (1 - eta))) * (1 + r)
end

# 3. At s = 39, 38, ... , 1, given ks[s+1], ns[s+1], ks[s+2], solve ks[s] and ns[s]
function fs_young!(F, X, ks1, ns1, ks2)
    ks, ns = X
    F[1] = (1 - tau) * w * (1 - ns) / gamma - ((1 + r) * ks + (1 - tau) * w * ns - ks1 + psi)
    F[2] = 1 / beta - (((1 + r) * ks1 + (1 - tau) * w * ns1 - ks2 + psi) / ((1 + r) * ks + (1 - tau) * w * ns - ks1 + psi))^(-eta) * (((1 - ns1) / (1 - ns))^(gamma / (1 - eta))) * (1 + r)
end

# =====================================================================
# == Backward Iteration Function ======================================
# == input: given nbar, solve the decision rules and ss age-profile ===
# == output: series of steady state age-profile ks_true and ns_true  ==
# == in the inner loop, first guess k[60] and iterate until k[1] = 0 ==
# =====================================================================
# Target: loop backward iteration until we get k[1] = 0
# change the guess of k[60] and iterate again if k[1] is not 0
max_iter = 30

function backward_iteration(nbar)
    k60_guess = zeros(max_iter + 1)
    k1_res = zeros(max_iter + 1)

    # true value
    ks_true = zeros(61)
    ns_true = zeros(61)

    # set tol value
    i = 1
    tol = 1e-6
    err = 0.1

    while i <= max_iter && abs(err) > tol
        # an empty array to store the results
        ks = zeros(61)
        ns = zeros(61)
        # update k60 guess
        if i == 1
            k60_guess[i] = 0.15
        elseif i == 2
            k60_guess[i] = 0.2
        else
            # update by secant method
            k60_guess[i] = k60_guess[i-1] - (k1_res[i-1] - 0) * (k60_guess[i-1] - k60_guess[i-2]) / (k1_res[i-1] - k1_res[i-2])
        end
        # initiate a guess for k[60]
        ks[60] = k60_guess[i]
        # calculate k[s] and n[s] for s = 59, 58, ..., 1
        for s in 59:-1:1
            if s >= 41
                # at s = 59, 58, ..., 41, given k[s+1] and k[s+2], solve k[s]
                ks[s] = solve_ks_old(ks[s+1], ks[s+2])
                ns[s] = 0
            elseif s == 40
                # at s = 40, given k41 and n41, solve k40 and n40
                result = nlsolve((F, X) -> fs40!(F, X, ks[s+1], ks[s+2]), [ks[s+1], ns[s+1]])
                ks[s] = result.zero[1]
                ns[s] = result.zero[2]
            else
                # at s = 39, 38, ..., 1, given k[s+1] and n[s+1], solve k[s] and n[s]
                result = nlsolve((F, X) -> fs_young!(F, X, ks[s+1], ns[s+1], ks[s+2]), [ks[s+1], ns[s+1]])
                ks[s] = result.zero[1]
                ns[s] = result.zero[2]
            end
        end
        # store the results
        ks_true = ks
        ns_true = ns
        # store k1 and calculate the error
        k1_res[i] = ks[1]
        # update error value
        err = ks[1]
        # increase the iteration if the error is still greater than the tolerance, otherwise break the loop
        if abs(err) > tol
            i += 1
        else
            break
        end
    end

    return ks_true, ns_true
end

# =====================================================================
# == Outer Loop Algorithm =============================================
# == input: a guess of nbar distribution ==============================
# == output: the true nbar with its associated ks_true and ns_true  ===
# =====================================================================

outer_max_iter = 30
outer_tol = 1e-6
outer_err = 1.0
outer_i = 1

# results
K = 0
N = 0
Ny = 0
ks_true = zeros(61)
ns_true = zeros(61)

# Initial guess for nbar
nbar = 0.2
kbar = update_kbar(nbar, r)
w = update_w(kbar, nbar)
b = repl*(1-tau)*w*nbar*3/2

# the loop algorithm
while outer_i <= outer_max_iter && abs(outer_err) > outer_tol
    kbar = update_kbar(nbar, r)
    w = update_w(kbar, nbar)

    # solve the decision rules steady state age-profile
    ks_true, ns_true = backward_iteration(nbar)

    # calculate aggregate capital
    K = sum(ks_true) / (T + TR)
    
    # calculate aggregate labor
    N = sum(ns_true) / (T + TR)
    Ny = sum(ns_true) / T   #workers only
    
    # nbar error
    nbar_error = N - nbar
    #println("nbar error: ", nbar_error)

    # update nbar
    if abs(nbar_error) > outer_tol
        nbar_new = phi * nbar + (1 - phi) * N
        nbar = nbar_new
        outer_err = nbar_error
        outer_i += 1
    else
        break
    end
end

# =====================================================================
# == Calculate other variables ========================================
# == earnings, consumption, welfare ===================================
# =====================================================================

# earnings
function cal_income(s, k, n)
    if s <= 40
        return (1 + r) * k + (1 - tau) * n * ((1 - alpha) * k^alpha * n^(-alpha))
    elseif s <= 60
        return (1 + r) * k + b
    else
        return 0
    end
end

# consumption
function cal_c(s, n, w, k1)
    if s <= 40
        return (1 - tau) * w * (1-n) / gamma - psi
    elseif s <= 60
        return w - k1
    else
        return 0
    end
end

# welfare
function cal_welfare(s, c, n)
    return beta^(s - 1) * (((c + psi) * ((1 - n)^gamma))^(1 - eta) - 1) / (1 - eta)
end

# using ks_true and ns_true to calculate the variables
ws_true = zeros(61)
cs_true = zeros(61)
us_true = zeros(61)
age_true = zeros(61)

for s in 1:60
    age_true[s] = s + 20
    ws_true[s] = cal_income(s, ks_true[s], ns_true[s])
    cs_true[s] = cal_c(s, ns_true[s], ws_true[s], ks_true[s+1])
    us_true[s] = cal_welfare(s, cs_true[s], ns_true[s])
end

ws_true[61] = 0
cs_true[61] = 0
us_true[61] = 0
age_true[61] = 81

# =======================================
# == Print the steady states ============
# =======================================

println("Final nbar: ", nbar)
println("Final K: ", kbar)
println("Final N: ", N)
println("Final Ny: ", Ny)

# ========================
# == Plotting ============
# ========================
# show 4 plots at the same time

# plot of ks

plotk = plot(age_true, ks_true, label="ks", title="Capital Distribution", xlabel="Age", ylabel="Capital", legend=false, color=:blue, lw=2, dpi=300)
plotn = plot(age_true, ns_true, label="ns", title="Labor Distribution", xlabel="Age", ylabel="Labor", legend=false, color=:red, lw=2, dpi=300)
plotw = plot(age_true, ws_true, label="ws", title="Income Distribution", xlabel="Age", ylabel="Wage", legend=false, color=:green, lw=2, dpi=300)
plotc = plot(age_true, cs_true, label="cs", title="Consumption Distribution", xlabel="Age", ylabel="Consumption", legend=false, color=:purple, lw=2, dpi=300)

fig = plot(plotk, plotn, plotw, plotc, layout=(2, 2), legend=false, size=(800, 600))

savefig(fig, "AK60.png")
