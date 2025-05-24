module ProgMainSS

export do_table45

include("Custom_fct.jl")

#using Pkg
#Pkg.add("LaTeXStrings")


using Statistics, LaTeXStrings, PrettyTables, FixedEffectModels, StatsBase, LinearAlgebra, Random, Printf, Distributions, .Custom

    function do_table45(m)
################################################################################
# ─────────────────────────────  Helper Stubs  ──────────────────────────────── #
################################################################################

# Simple linear interpolation on a grid (1‑D, monotone grid)
# Equivalent to MATLAB’s custom lini helper used in the original file.
#   xgrid  : Vector of grid points (increasing)
#   ygrid  : Vector of function values same length as xgrid
#   xquery : Point to interpolate at

################################################################################
# ──────────────────────────────  Main routine  ─────────────────────────────── #
################################################################################

        function Prog_MainSS(option, Ah, psicsd, psil, λ, φ1, φ2add, φ3add, γ, κρ, κσ, α1, θ0, χ, τx, ω, distoff)

            ############################  Housekeeping  ################################
            title = "______________________________ Compute the Steady State Equilibrium ______________________________ "
            println(title)

            nh        = 36          # HK grid for value, decision
            nh_dist   = 500         # HK grid for distribution
            hhigh0    = 15.0        # sanity bound for simulated h

            # Policy: means‑tested pronatal transfers
            inccut    = 100_000.0   # effectively no means test

            hbarhat   = 5.783218901406690
            pct_target = 0.85       # top 15 %

            σ    = 0.5             # not used
            A    = 1.0
            Yss  = 1.0
            α2   = 0.0             # not used
            ψ̅   = 1.0             # psic mean

            # Meta params
            tol      = 1e-4
            toldist  = 1e-6
            tolout   = 1e-4
            damph    = 0.0

            nn   = 4               # max children = 3 (grid entries 0:3)
            nk   = 20              # grid for kappa
            npsi = 20              # grid for psi
            nt   = 50              # grid for leisure‑supply time choice
            EVmat = zeros(nh, npsi, nk, nn)  # expected value function

            ngrid = collect(0:nn-1)   # 0,1,2,3

            # Custom non‑linear time grid (first 90 % dense on [tol,0.7])
            tgrid = zeros(nt)
            nt0   = round(Int, nt*0.9)
            tgrid[1:nt0]   = range(tol,   stop=0.7, length=nt0)
            tgrid[nt0:end] = range(0.7,   stop=1.0, length=nt-nt0+1)

            ############################################################################
            # 1. Stochastic processes & grids                                          #
            ############################################################################

            κgrid, κP = tauchenmethod(0.0, κσ, κρ, nk,    3.0)
            κgrid = exp.(κgrid)

            ψgrid, ψP = tauchenmethod(0.0, psicsd, 0.0, npsi, 3.0)
            ψgrid = ψgrid .* ψ̅
            ψP = ψP[1,:]           # first‑row stationary distribution as in MATLAB

            # Non‑parametric fertility utility grid φ(n)
            φgrid = zeros(nn)
            φgrid[2] = φ1
            φgrid[3] = φ1 + φ2add
            φgrid[4] = φ1 + φ2add + φ3add

            # Human‑capital grids
            hlow  = max(0.5, HKprodfn(Ah*κgrid[1], 0.0, 0.0, θ0, α1, α2))
            hhigh = Yss*hhigh0
            hgrid = gengrid(hlow, hhigh, nh; option=2, degree=2.0)
            hgrid_dist  = collect(range(hlow, hhigh; length=nh_dist))

            # Simulation sample size
            NS = 500_000

            ############################################################################
            # 2. Pre‑allocate large arrays                                             #
            ############################################################################

            Vmat   = zeros(nh, npsi, nk, nn)
            cmat   = similar(Vmat)
            lmat   = similar(Vmat)
            hwmat  = similar(Vmat)   # hours worked
            xmat   = similar(Vmat)   # education spending

            nindmat      = zeros(Int, nh, npsi, nk)
            nmat         = zeros(nh, npsi, nk)
            nindmat_dist = zeros(Int, nh_dist, npsi, nk)

            # Initial guesses
            L      = Yss
            T      = 0.0
            Hmean  = L * 3.0
    H15    = L * 5.5

    mumatpre = fill(1.0/(nh_dist*nk), nh_dist, nk)   # distribution over (h,kappa)

    ############################################################################
    # 3. Outer loop: general‑equilibrium prices                                #
    ############################################################################

    distanceout   = 1e4
    iteration_out = 1
    avgfert       = 0.0
    w = A   # wage = productivity (unit‑normalized)

    THs = 0.0; Ls = 0.0; xmean = 0.0; txmean = 0.0; Tnew = 0.0
    Y  = A*Ls

    while distanceout > tolout && iteration_out ≤ 30
        @info ">>> OUTER ITERATION $iteration_out"
        

        println(@sprintf("H15 =%8.2f   Hmean =%8.2f   hhigh =%8.2f   hlow =%8.2f", H15, Hmean, hhigh, hlow))

        ########################  INNER LOOP  #################################
        distance = 1e4
        iteration_in = 1
        mumatpre0 = copy(mumatpre)

        while distance > toldist && iteration_in ≤ 30
            ################  3.1  VALUE FUNCTION GIVEN N=0  ##################
            in = 1   # index for n=0
            nnow = ngrid[in]
            Threads.@threads for ih in 1:nh
                hpnow = hgrid[ih]
                for ib in 1:npsi, jk in 1:nk
                    ψ = ψgrid[ib]
                    κ = κgrid[jk]

                    lmax = 1.0
                    if w*hpnow*lmax > T
                        nt0 = count(<=(lmax), tgrid)
                        values = fill(-Inf, nt0)
                        for (idx,it) in enumerate(1:nt0)
                            les = tgrid[it]
                            fun = ψ*log(HES(w*hpnow*(lmax-les) + T, 0)) + psil*(les^(1-γ))/(1-γ)
                            isfinite(fun) && (values[idx] = fun)
                        end
                        Nval, idx = findmax(values)
                        lmat[ih,ib,jk,in]  = tgrid[idx]
                        hwmat[ih,ib,jk,in] = lmax - tgrid[idx]
                        Vmat[ih,ib,jk,in]  = Nval
                        cmat[ih,ib,jk,in]  = w*hpnow*hwmat[ih,ib,jk,in] + T
                        xmat[ih,ib,jk,in]  = 0.0
                    else
                        Vmat[ih,ib,jk,in] = -Inf
                    end
                end
            end  # end threads

            ################  3.2  VALUE FUNCTION GIVEN N>0  ##################
            Threads.@threads for ih in 1:nh
                hpnow = hgrid[ih]
                for in in 2:nn
                    nnow   = ngrid[in]
                    φnow   = φgrid[in]
                    for ib in 1:npsi, ik in 1:nk
                        ψ   = ψgrid[ib]
                        κ   = κgrid[ik]

                        htilde = distoff == 1 ? hbarhat : H15
                        lmax   = 1 - λ*nnow
                        nt0    = count(<=(lmax), tgrid)

                        xmin = begin
                            qualutil10(HKprodfn(Ah*κ, hpnow, 0.0, θ0, α1, α2), htilde, χ, σ) > -Inf ? 0.0 :
                            max(((χ*htilde)/(Ah*κ) - θ0)^(1/α1)*hpnow^(-α2/α1), 0.0)
                        end

                        ##################  CASE 1 : interior x  ############
                        v1  = fill(-Inf, nt0)
                        c1  = fill(0.0,   nt0)
                        for (idx,it) in enumerate(1:nt0)
                            les = tgrid[it]
                            Trn_temp = w*hpnow*(lmax-les) < inccut ? Trn(ω, nnow) : 0.0
                            cmax = w*hpnow*(lmax-les) + T + Trn_temp - (1+τx)*xmin*nnow
                            if cmax > 0.0
                                obj(c) = -( ψ*log(HES(c, nnow)) + psil*(les^(1-γ))/(1-γ) +
                                           φnow*qualutil10(HKprodfn(Ah*κ, hpnow,
                                           (w*hpnow*(lmax-les) + T + Trn_temp - c)/((1+τx)*nnow), θ0, α1, α2),
                                           htilde, χ, σ) )
                                c1[idx] = GSS(obj, 0.0, cmax, tol)
                                v1[idx] = -obj(c1[idx])
                            end
                        end
                        Vval1, idx1 = findmax(v1)

                        ##################  CASE 2 : x = xmin  ###############
                        v2 = fill(-Inf, nt0)
                        for (idx,it) in enumerate(1:nt0)
                            les = tgrid[it]
                            Trn_temp = w*hpnow*(lmax-les) < inccut ? Trn(ω, nnow) : 0.0
                            cmax = w*hpnow*(lmax-les) + T + Trn_temp - (1+τx)*xmin*nnow
                            if cmax > 0.0
                                fun = ψ*log(HES(w*hpnow*(lmax-les) + T + Trn_temp - (1+τx)*xmin*nnow, nnow)) +
                                      psil*(les^(1-γ))/(1-γ) +
                                      φnow*qualutil10(HKprodfn(Ah*κ, hpnow, xmin, θ0, α1, α2), htilde, χ, σ)
                                v2[idx] = fun
                            end
                        end
                        Vval2, idx2 = findmax(v2)

                        if Vval1 > Vval2
                            les_star = tgrid[idx1]
                            hw       = lmax - les_star
                            Trn_temp = w*hpnow*hw < inccut ? Trn(ω, nnow) : 0.0
                            x_star   = (w*hpnow*hw + T + Trn_temp - c1[idx1]) / ((1+τx)*nnow)
                            lmat[ih,ib,ik,in]  = les_star
                            hwmat[ih,ib,ik,in] = hw
                            xmat[ih,ib,ik,in]  = x_star
                            cmat[ih,ib,ik,in]  = c1[idx1]
                            Vmat[ih,ib,ik,in]  = Vval1
                        else
                            les_star = tgrid[idx2]
                            hw       = lmax - les_star
                            Trn_temp = w*hpnow*hw < inccut ? Trn(ω, nnow) : 0.0
                            lmat[ih,ib,ik,in]  = les_star
                            hwmat[ih,ib,ik,in] = hw
                            xmat[ih,ib,ik,in]  = xmin
                            cmat[ih,ib,ik,in]  = w*hpnow*hw + T + Trn_temp - (1+τx)*xmin*nnow
                            Vmat[ih,ib,ik,in]  = Vval2
                        end
                    end # loop ψ,κ
                end # loop n
            end # end threads IH

            ################  3.3  Expected value over κ′  ###################
            for ih in 1:nh, ib in 1:npsi, in in 1:nn, ik in 1:nk
                EVmat[ih,ib,ik,in] = dot(κP[ik, :], Vmat[ih,ib,:,in])
            end

            ################  3.4  Fertility choice on coarse grid #############
            for ih in 1:nh, ib in 1:npsi, ik in 1:nk
                _, idx = findmax(EVmat[ih,ib,ik,:])
                nindmat[ih,ib,ik] = idx
                nmat[ih,ib,ik]    = ngrid[idx]
            end

            sum(nmat) == 0 && @warn "Nobody wants kids!"

            ################  3.5  Stationary distribution                    #
            distance0  = 1e4
            mumatpre0_ = copy(mumatpre0)
            iteration_dist = 1
            while distance0 > toldist && iteration_dist ≤ 30
                # Data containers reused each pass
                Tmumat = zeros(nh_dist, nk)

                for ih_d in 1:nh_dist
                    hpnow = hgrid_dist[ih_d]
                    for ik in 1:nk, ib in 1:npsi
                        # fertility index at (h,kappa,psi)
                        valn = [lini(hgrid, EVmat[:,ib,ik,in], hpnow) for in in 1:nn]
                        _, idxn = findmax(valn)
                        nindmat_dist[ih_d, ib, ik] = idxn
                        nnow = ngrid[idxn]
                        for jk in 1:nk
                            κkid = κgrid[jk]
                            xtemp = lini(hgrid, xmat[:,ib,jk,idxn], hpnow)
                            hnext = HKprodfn(Ah*κkid, hpnow, xtemp, θ0, α1, α2)
                            # linear weights on hgrid_dist
                            if hnext ≤ hgrid_dist[1]
                                jh  = 1; wgt = 1.0
                            elseif hnext ≥ hgrid_dist[end]
                                jh  = nh_dist-1; wgt = 0.0
                            else
                                jh  = searchsortedfirst(hgrid_dist, hnext) - 1
                                wgt = (hgrid_dist[jh+1] - hnext) / (hgrid_dist[jh+1] - hgrid_dist[jh])
                            end
                            jhmat = jh
                            whmat = wgt
                            muval = mumatpre0_[ih_d,ik] * ψP[ib]
                            Tmumat[jh,   jk] += κP[ik,jk] * wgt * nnow * muval
                            Tmumat[jh+1, jk] += κP[ik,jk] * (1-wgt) * nnow * muval
                        end
                    end
                end # Threads

                avgfert   = sum(Tmumat)          # expected kids per adult * μ * κ sum
                popgrowth = avgfert / 2          # siblings share? keep from MATLAB
                Tmumat ./= avgfert               # normalise

                distance0 = maximum(abs, Tmumat .- mumatpre0_)
                mumatpre0_ .= Tmumat             # overwrite

                @info @sprintf("    stationary‑dist iter %2d  ||Δμ||=%8.7g", iteration_dist, distance0)
                iteration_dist += 1
            end # stationary distribution loop

            mumatpre0 = copy(mumatpre0_)
            distance  = distoff == 1 || χ==0 ? 0.0 : maximum(abs, mumatpre0 .- mumatpre)
            mumatpre  = copy(mumatpre0)

            ##########  3.6  Update Hmean & H15  ##############################
            muh = vec(sum(mumatpre, dims=2))
            Hmean = dot(hgrid_dist, muh)
            cumh = cumsum(muh)
            topind = findfirst(cumh .>= pct_target)  # first index beyond 85 %
            muhtop = muh[topind:end]
            muhtop[1] -= pct_target - (topind > 1 ? cumh[topind-1] : 0.0)
            hgridtop = hgrid_dist[topind:end]
            H15 = dot(hgridtop, muhtop) / sum(muhtop)

            @info @sprintf("  inner‑loop %2d  ||Δμ||=%8.7g   H15=%6.3f", iteration_in, distance, H15)
            iteration_in += 1
        end # inner loop

        ############################################################################
        # 3.7  Aggregate moments & update prices                                 #
        ############################################################################

        
        for ih_d in 1:nh_dist, ik in 1:nk, ib in 1:npsi
            muval_base = mumatpre[ih_d,ik] * ψP[ib]
            in = nindmat_dist[ih_d,ib,ik]
            jk_loop = 1:nk
            for jk in jk_loop
                muval = muval_base * κP[ik,jk]
                hpnow = hgrid_dist[ih_d]
                hw = lini(hgrid, hwmat[:,ib,jk,in], hpnow)
                xt = lini(hgrid,  xmat[:,ib,jk,in], hpnow)
                nnow = ngrid[in]
                THs += hw * muval
                Ls  += hpnow * hw * muval
                Trn_temp = w*hpnow*hw < inccut ? Trn(ω, nnow) : 0.0
                xmean  += xt * muval
                txmean += nnow * xt * muval
                Tnew   += (τx*xt*nnow - Trn_temp) * muval
            end
        end
        
        Y = A*Ls
        distanceout = max(abs(log(Ls) - log(L)), abs(T - Tnew))
        TY = T/Y

        println(@sprintf("outer loop %2d  |Δ|=%8.5f", iteration_out, distanceout))
        println(@sprintf("  L=%.3f  Ls=%.3f   τx=%.3f  ω=%.3f", L, Ls, τx, ω))
        println(@sprintf("  T=%.3f  Tnew=%.3f  T/Y=%.5f", T, Tnew, TY))
        println(@sprintf("  Fertility=%.3f  TotHours=%.3f  mean(h)=%.3f  mean(x)=%.3f  Y=%.3f", avgfert, THs, Hmean, xmean, Y))
        println()

        L = Ls
        T = damph*T + (1-damph)*Tnew
        iteration_out += 1
    end # outer loop

    ########################################################
    # 4. Simulation (micro sample)                         #
    ########################################################

    # convert mumat to vector for sampling (MATLAB used reshape)
    muvec = vec(mumatpre)
    idxvec = 1:length(muvec)
    rng = MersenneTwister(1)
    distind = sample(rng, idxvec, Weights(muvec), NS)

    hdistsimind = [ rem(idx-1, nh_dist)+1 for idx in distind ]
    kappasimind = [ (idx-1) ÷ nh_dist + 1 for idx in distind ]

    bsimind       = sample(rng, 1:npsi, Weights(ψP), NS)
    bsimindkid    = sample(rng, 1:npsi, Weights(ψP), NS)
    kappasimindkid     = [ sample(rng, 1:nk, Weights(κP[k,:])) for k in kappasimind ]
    kappasimindkidkid = [ sample(rng, 1:nk, Weights(κP[kk,:])) for kk in kappasimindkid ]

    # Allocate simulation arrays
    nsim      = zeros(Int, NS)
    nsimpos   = fill(NaN, NS)
    child0    = zeros(Int, NS)
    hwsim     = zeros(NS)
    consim    = zeros(NS)
    bsim      = zeros(NS)
    xsim      = zeros(NS)
    txsim     = zeros(NS)
    nxsim     = zeros(NS)
    hsim      = zeros(NS, 2)
    incsim    = zeros(NS)
    incsim2   = zeros(NS)
    transim   = zeros(NS)

    Threads.@threads for i in 1:NS
        ih_d = hdistsimind[i]
        ib    = bsimind[i]
        ik    = kappasimind[i]
        jk    = kappasimindkid[i]

        in    = nindmat_dist[ih_d, ib, ik]
        nnow  = ngrid[in]
        nsim[i] = nnow
        if nnow == 0
            child0[i] = 1
        else
            nsimpos[i] = nnow
        end

        hp = hgrid_dist[ih_d]
        hsim[i,1] = hp
        hwsim[i]  = lini(hgrid, hwmat[:,ib,jk,in], hp)
        consim[i] = lini(hgrid, cmat[:,ib,jk,in], hp)
        xsim[i]   = lini(hgrid, xmat[:,ib,jk,in], hp)
        nxsim[i]  = (1+τx)*xsim[i]*nsim[i]
        incsim[i] = w*hp*hwsim[i]
        bsim[i]   = ψgrid[ib]
        txsim[i]  = xsim[i]/incsim[i]
        Trn_temp  = incsim[i] < inccut ? Trn(ω, nnow) : 0.0
        transim[i] = T + Trn_temp
        hsim[i,2] = HKprodfn(Ah*κgrid[jk], hp, xsim[i], θ0, α1, α2)

        # Child’s income (generation+1)
        jb  = bsimindkid[i]
        jjk = kappasimindkidkid[i]
        Vtemp = [ lini(hgrid, EVmat[:,jb,jk2,jn], hsim[i,2]) for (jn,jk2) in zip(1:nn, repeat([jk], nn)) ]
        _, in2 = findmax(Vtemp)
        incsim2[i] = w*hsim[i,2]*lini(hgrid, hwmat[:,jb,jjk,in2], hsim[i,2])
    end

    budgetbalance = (mean(incsim) + mean(transim) - mean(nxsim) - mean(consim))
    avg_educinv_inc = mean(xsim[.!isnan.(nsimpos)]) / mean(incsim[.!isnan.(nsimpos)])
    avg_LS = mean(hwsim)

    # Quintile statistics & elasticities (very condensed vs MATLAB)
    qbreaks = quantile(incsim, [0.2, 0.4, 0.6, 0.8])
    println("qbreaks = ", qbreaks)
    incquin = searchsortedfirst.(Ref(qbreaks), incsim) #.+ 1
    incquin .= min.(incquin, 5)   # ensure 1–5
    #println("incquin = ", incquin)

    AvgI5 = [mean(incsim[incquin .== k]) for k in 1:5]
    AvgF5 = [mean(nsim[incquin .== k])  for k in 1:5]
    AvgX5 = [mean(xsim[incquin .== k .&& nsim .> 0]) for k in 1:5]
    AvgI5a = [mean(incsim[incquin .== k .&& nsim .> 0]) for k in 1:5]
    AvgH5 = [mean(hsim[incquin .== k,1]) for k in 1:5]

    pctN   = [sum(incquin .== k) for k in 1:5]
    println("pctN = ", pctN)
    pctNpos = [sum(incquin .== k .&& nsim .> 0) for k in 1:5]

    Wmat = Diagonal(pctN)
    println("Wmat = ", Wmat)
    #df.Wmat = trace(Diagonal(pctN))
    #log_AvgF5 = log.(AvgF5)
    #log_AvgI5 =log.(AvgI5)
    #log_AvgX5 = log.(AvgX5)

    elast0 = wls1(log.(AvgF5), log.(AvgI5), Wmat)
    #elast0 = reg(df, @formula(log_AvgF5 ~ log_AvgI5, weights = :Wmat))
    println(elast0)
    elastx = wls1(log.(AvgX5), log.(AvgI5), Wmat)
    #elastx = reg(df, @formula(log_AvgX5 ~ log_AvgI5, weights = :Wmat))

    AvgN05 = 1 .- pctNpos ./ pctN
    invinc5 = AvgX5 ./ AvgI5a

    prn = zeros(nn)
    for j in 1:nn
        prn[j] = count(nsim .== (j-1)) / NS
    end

    Giniinc = ginicoeff(incsim)[1][1]

    hpvec = incsim[nsim .> 0]
    hkvec = incsim2[nsim .> 0]
    ige = ols1(log.(hkvec), log.(hpvec))
    #ige = reg(df, @formula(log.(hkvec) ~ log.(hpvec)))

    ############################################################################
# 5-bis.  Investment-spillover IV step (MATLAB “Step 1” → “iv_coeff”)      #
############################################################################
xmeandev      = 0.9                         # –10 % counter-factual change
xmean_p_chg   = 100 * (xmeandev - 1)        # percentage change (−10)
println("xmean_p_chg = ", xmean_p_chg)

# ── 5.1  Identify parents in the top 15 % of human-capital distribution ──
p90val = quantile(hsim[:,1], pct_target)    # 85-th percentile of parental h

topmask  = (hsim[:,1] .>= p90val) .& (nsim .> 0)   # top-15 % *and* has children

xsim_top     = skipmissing(xsim[topmask])           # drop NaNs automatically
kappa_top    = κgrid[kappasimindkid[topmask]]

Hpmean0 = HKprodfn(Ah * mean(kappa_top),
                   H15,
                   mean(xsim_top),
                   θ0, α1, α2)

Hpmean1 = HKprodfn(Ah * mean(kappa_top),
                   H15,
                   xmeandev * mean(xsim_top),
                   θ0, α1, α2)

Hdev = Hpmean1 / Hpmean0                # human-capital “shock” multiplier
med_inc = median(incsim)                # median income (baseline)

# ── 5.2  Re-solve household problem under new h̃ = Hdev·H15 ───────────────
#        (We reuse exactly the same loops used for Vmat but with Hdev)     #

Vmat_iv   = similar(Vmat)   ;  fill!(Vmat_iv, 0.0)
cmat_iv   = similar(Vmat_iv);  lmat_iv = similar(Vmat_iv)
hwmat_iv  = similar(Vmat_iv);  xmat_iv = similar(Vmat_iv)

# --- (a) n = 0 ------------------------------------------------------------
for ih in 1:nh, ib in 1:npsi, ik in 1:nk
    hp   = hgrid[ih]
    ψ    = ψgrid[ib]
    κ    = κgrid[ik]
    lmax = 1.0
    if w*hp*lmax > T
        bestV = -Inf
        bestt = 0.0
        for les in tgrid[tgrid .<= lmax]
            v = ψ*log(HES(w*hp*(lmax-les)+T, 0)) + psil*(les^(1-γ))/(1-γ)
            if v > bestV
                bestV = v
                bestt = les
            end
        end
        lmat_iv[ih,ib,ik,1]  = bestt
        hwmat_iv[ih,ib,ik,1] = lmax - bestt
        Vmat_iv[ih,ib,ik,1]  = bestV
        cmat_iv[ih,ib,ik,1]  = w*hp*hwmat_iv[ih,ib,ik,1] + T
        xmat_iv[ih,ib,ik,1]  = 0.0
    else
        Vmat_iv[ih,ib,ik,1] = -Inf
    end
end

# --- (b) n > 0  -----------------------------------------------------------
for ih in 1:nh, in in 2:nn, ib in 1:npsi, ik in 1:nk
    hp     = hgrid[ih]
    nnow   = ngrid[in]
    φnow   = φgrid[in]
    ψ      = ψgrid[ib]
    κ      = κgrid[ik]

    htilde = distoff == 1 ? Hdev * hbarhat : Hdev * H15
    lmax   = 1.0 - λ*nnow
    xmin   = begin
        qualutil10(HKprodfn(Ah*κ, hp, 0.0, θ0, α1, α2), htilde, χ, σ) > -Inf ? 0.0 :
        max(((χ*htilde)/(Ah*κ) - θ0)^(1/α1)*hp^(-α2/α1), 0.0)
    end

    bestV, bestt, bestc, bestx = -Inf, 0.0, 0.0, xmin
    for les in tgrid[tgrid .<= lmax]
        Trn_tmp = w*hp*(lmax-les) < inccut ? Trn(ω, nnow) : 0.0
        cmax    = w*hp*(lmax-les) + T + Trn_tmp - (1+τx)*xmin*nnow
        if cmax ≤ 0.0; continue; end

        # inner optimum over c (Golden-section)
        # inner optimum over c  (Golden-section search)
        obj_inner(c) =  ψ*log(HES(c,nnow)) +
                psil*(les^(1-γ))/(1-γ) +
                φnow*qualutil10(
                HKprodfn(Ah*κ, hp,
               (w*hp*(lmax-les)+T+Trn_tmp-c)/((1+τx)*nnow),
               θ0, α1, α2),
                htilde, χ, σ)

            c_star = GSS(c -> -obj_inner(c), 0.0, cmax, tol)   # minimise –obj ⇒ maximise obj
            val1   = obj_inner(c_star)                          # value at optimum
     # GSS returns (c*, –V)
        if val1 > bestV
            bestV, bestt, bestc = val1, les, c_star
            bestx = (w*hp*(lmax-les) + T + Trn_tmp - bestc)/((1+τx)*nnow)
        end
        # edge case x = xmin
        val2 = ψ*log(HES(w*hp*(lmax-les)+T+Trn_tmp-(1+τx)*xmin*nnow, nnow)) +
               psil*(les^(1-γ))/(1-γ) +
               φnow*qualutil10(HKprodfn(Ah*κ, hp, xmin, θ0, α1, α2), htilde, χ, σ)
        if val2 > bestV
            bestV, bestt = val2, les
            bestc  = w*hp*(lmax-les)+T+Trn_tmp-(1+τx)*xmin*nnow
            bestx  = xmin
        end
    end
    lmat_iv[ih,ib,ik,in]  = bestt
    hwmat_iv[ih,ib,ik,in] = lmax - bestt
    cmat_iv[ih,ib,ik,in]  = bestc
    xmat_iv[ih,ib,ik,in]  = bestx
    Vmat_iv[ih,ib,ik,in]  = bestV
end

# --- (c) Expectation over κ′ ---------------------------------------------
EVmat_iv = zeros(nh, npsi, nk, nn)
for ih in 1:nh, ib in 1:npsi, ik in 1:nk, in in 1:nn
    EVmat_iv[ih,ib,ik,in] = dot(κP[ik,:], Vmat_iv[ih,ib,:,in])
end

# --- (d) Fertility decision on the fine h-grid ----------------------------
nindmat_iv_dist = zeros(Int, nh_dist, npsi, nk)
for ihd in 1:nh_dist, ib in 1:npsi, ik in 1:nk
    hp = hgrid_dist[ihd]
    vals = [lini(hgrid, EVmat_iv[:,ib,ik,in], hp) for in in 1:nn]
    nindmat_iv_dist[ihd,ib,ik] = argmax(vals)
end

# ── 5.3  Micro simulation with baseline fertility but new choices x,c –─
xsim_iv      = zeros(NS)
consim_iv    = zeros(NS)
incsim_iv    = zeros(NS)
nsim_iv      = zeros(NS)
xratio_base  = zeros(NS)
xratio_iv    = zeros(NS)

for i in 1:NS
    ihd = hdistsimind[i];  ib = bsimind[i]; ik = kappasimind[i];  jk = kappasimindkid[i]
    in  = nindmat_iv_dist[ihd, ib, ik]
    nsim_iv[i] = in - 1                      # fertility under iv spec (n-index starts at 1)
    hp = hgrid_dist[ihd]

    consim_iv[i] = lini(hgrid, cmat_iv[:,ib,jk,in], hp)
    xsim_iv[i]   = lini(hgrid, xmat_iv[:,ib,jk,in], hp)
    #println("xsim_iv = ", xsim_iv[i], "  consim_iv = ", consim_iv[i])
    hwork        = lini(hgrid, hwmat[:,ib,jk,in], hp)

    xratio_base[i] = (1+τx)*xsim[i]   / ((1+τx)*xsim[i]*nsim[i]   + consim[i])
    #xratio_iv[i]   = (1+τx)*xsim_iv[i]/ ((1+τx)*xsim_iv[i]*nsim_iv[i] + consim_iv[i])
    denom_iv = (1 + τx) * xsim_iv[i] * nsim_iv[i] + consim_iv[i]
    if denom_iv != 0
        xratio_iv[i] = (1 + τx) * xsim_iv[i] / denom_iv
    else
        xratio_iv[i] = NaN
    end

    #println("xratio_base = ", xratio_base[i], "  xratio_iv = ", xratio_iv[i])
    incsim_iv[i]  = w * hp * hwork
end

# ── 5.4  IV coefficient ----------------------------------------------------
mask = (nsim .> 0) .& (incsim .< med_inc)
#println("mask = ", mask)
#Δratio = mean(xratio_iv[mask]) - mean(xratio_base[mask])

valid_mask = mask .& .!isnan.(xratio_iv) .& .!isnan.(xratio_base)
Δratio = mean(xratio_iv[valid_mask]) - mean(xratio_base[valid_mask])

println("Δratio = ", Δratio)
iv_coeff = 100*Δratio / xmean_p_chg         # identical to MATLAB formula
println("iv_coeff = ", iv_coeff)

println(elast0)

    ########################  Objective  ######################################
    objval = ((0.053-AvgN05[1])/0.053)^2 + ((0.196-prn[2])/0.196)^2 +
              ((0.631-prn[3])/0.631)^2 + ((0.144-prn[4])/0.144)^2 +
              ((0.039-iv_coeff)/0.039)^2 +
              ((avg_educinv_inc-0.092)/0.092)^2 + ((elast0-0.082)/0.082)^2 +
              ((0.301-avg_LS)/0.301)^2 + ((0.263-Giniinc)/0.263)^2 +
              ((0.698-elastx)/0.698)^2 + ((ige-0.32)/0.32)^2 + ((Y-Yss)/Yss)^2

println(elast0)

            
    ########################################################################
    # 6.  Display “Table 4” and “Table 5” (only when option == 1)           #
    ########################################################################
 
    if option == 1
        println(elast0)

        println()                                        # blank line
        println("----------------------------  Table 4  ----------------------------")
        println("****************** Parameters ******************")
        @printf("   phi_1=%6.2f   phi_2=%6.2f   phi_3=%6.2f\n",
                φgrid[2], φgrid[3], φgrid[4])
        @printf("   sigma_kappa=%6.3f   nu=%6.2f   sigma_b=%6.2f   chi=%6.2f\n",
                κσ, psil, psicsd, χ)
        @printf("   theta=%6.2f   alpha=%6.3f   rho_kappa=%6.3f   A_h=%6.2f\n",
                θ0, α1, κρ, Ah)

        println("****************** Target moments ******************")
        @printf("   Model=%6.3f   Data=%6.3f   Pr(n=1)\n", prn[2], 0.196)
        @printf("   Model=%6.3f   Data=%6.3f   Pr(n=2)\n", prn[3], 0.631)
        @printf("   Model=%6.3f   Data=%6.3f   Pr(n>=3)\n", prn[4], 0.144)
        @printf("   Model=%6.3f   Data=%6.3f   Gini income\n", Giniinc, 0.263)
        @printf("   Model=%6.3f   Data=%6.3f   Avg total hours worked\n", avg_LS, 0.301)
        @printf("   Model=%6.3f   Data=%6.3f   Income elasticity of fertility\n", elast0, 0.082)
        @printf("   Model=%6.3f   Data=%6.3f   Education inv. spillover estimate\n", iv_coeff, 0.039)
        @printf("   Model=%6.3f   Data=%6.3f   Childless in 1st income quintile\n", AvgN05[1], 0.053)
        @printf("   Model=%6.3f   Data=%6.3f   Avg investment-income ratio\n", avg_educinv_inc, 0.092)
        @printf("   Model=%6.3f   Data=%6.3f   Income elasticity of investment\n", elastx, 0.698)
        @printf("   Model=%6.3f   Data=%6.3f   Intergenerat. income elasticity\n", ige, 0.320)
        @printf("   Model=%6.3f   Data=%6.3f   Output per capita\n", Y, Yss)
        println("___________________________________________________________________________\n")

        println("----------------------------  Table 5  ----------------------------")
        println("************* Completed fertility *************")
        @printf("   Model=%6.2f   Data=%6.2f   All\n", avgfert, 1.91)
        for k in 1:5
            dataF = (1.80,1.91,1.87,1.93,2.03)[k]
            @printf("   Model=%6.2f   Data=%6.2f   %d%s income quintile\n",
                    AvgF5[k], dataF, k,
                    k==1 ? "st" : k==2 ? "nd" : k==3 ? "rd" : "th")
        end

        println("************* Childlessness rate (percent) *************")
        @printf("   Model=%6.1f   Data=%6.1f   All\n", 100*prn[1], 2.9)
        for k in 1:5
            dataC = (5.3,4.0,2.0,1.3,2.0)[k]
            @printf("   Model=%6.1f   Data=%6.1f   %d%s income quintile\n",
                    100*AvgN05[k], dataC, k,
                    k==1 ? "st" : k==2 ? "nd" : k==3 ? "rd" : "th")
        end

        println("************* Education spending per child rel. to income (percent) *************")
        @printf("   Model=%6.1f   Data=%6.1f   All\n", 100*avg_educinv_inc, 9.2)
        for k in 1:5
            dataX = (11.2,9.9,9.3,8.7,6.9)[k]
            @printf("   Model=%6.1f   Data=%6.1f   %d%s income quintile\n",
                    100*invinc5[k], dataX, k,
                    k==1 ? "st" : k==2 ? "nd" : k==3 ? "rd" : "th")
        end
        println("___________________________________________________________________________\n")

       #########################################################################TABLE 4
        # Parameters and their values
        param_names = [
            latexstring(" \\phi_1 = "), latexstring(" \\phi_2 = "), latexstring("\\phi_3 = "),
            latexstring("\\sigma_\\kappa = "), latexstring("\\nu = "), latexstring("\\sigma_b = "), latexstring("\\chi = "),
            latexstring("\\theta = "), latexstring("\\alpha = "), latexstring("\\rho_\\kappa = "), latexstring("A_h = ")
        ]

        param_vals = [
            string(round(φgrid[2], digits=2)), string(round(φgrid[3], digits=2)), string(round(φgrid[4], digits=2)),
            string(round(κσ, digits=3)), string(round(psil, digits=2)), string(round(psicsd, digits=3)), string(round(χ, digits=3)),
            string(round(θ0, digits=2)), string(round(α1, digits=3)), string(round(κρ, digits=3)), string(round(Ah, digits=2))
        ]

        # Moment descriptions and values
        moment_desc = [
            latexstring("\\mathrm{Pr}(\\#child=1)"), 
            latexstring("\\mathrm{Pr}(\\#child=2)"), 
            latexstring("\\mathrm{Pr}(\\#child\\geq 3)"), 
            "Gini income",
            "Avg total hours worked",
            "Income elasticity of fertility",
            "Education inv. spillover estimate",
            "Childless in first income quintile",
            "Avg. investment-income ratio",
            "Income elasticity of investment",
            "Intergenerat. income elasticity",
            "Output per capita (normalization)"
        ]

        model_vals = [
            string(round(prn[2], digits=3)), string(round(prn[3], digits=3)), string(round(prn[4], digits=3)),
            string(round(Giniinc, digits=3)), string(round(avg_LS, digits=3)), string(round(elast0, digits=3)),
            string(round(iv_coeff, digits=3)), string(round(AvgN05[1], digits=3)), string(round(avg_educinv_inc, digits=3)),
            string(round(elastx, digits=3)), string(round(ige, digits=3)), string(round(Y, digits=3))
        ]

        data_vals = [
            "0.196", "0.631", "0.144",
            "0.263", "0.301", "0.082",
            "0.039", "0.053", "0.092",
            "0.698", "0.320", string(round(Yss, digits=3))
        ]

        # Combine into rows
        nrows = max(length(param_names), length(moment_desc))
        data = Matrix{Any}(undef, nrows, 5)
        for i in 1:nrows
            data[i,1] = i <= length(param_names) ? param_names[i] : ""
            data[i,2] = i <= length(param_vals)  ? param_vals[i]  : ""
            data[i,3] = i <= length(moment_desc) ? moment_desc[i] : ""
            data[i,4] = i <= length(model_vals)  ? model_vals[i]  : ""
            data[i,5] = i <= length(data_vals)   ? data_vals[i]   : ""
        end

        # Column header
        header = ["Parameter and interpretation", "", "Moment", "Model", "Data"]

        latex_path = joinpath(m, "output/Table4.tex")

        # Write LaTeX table
        open(latex_path, "w") do io
            println(io, "\\begin{center}")
            println(io, "\\textbf{Table 4—Internally Calibrated Parameters and Moments: Model versus Data}")
            println(io, "\\end{center}")
            println(io, "")

            pretty_table(io, data;
                header=header,
                backend=Val(:latex),
                tf=tf_latex_booktabs,
                alignment=[:l, :r, :l, :r, :r],
                hlines=[1], # line under header
            )
        end


        

        ################TABLE 5

        header = (
            ["", "", "", "", "", "", ""],
            ["", "All", "1st", "2nd", "3rd", "4th", "5th"]
        )

        data = [
            "Completed fertility"  "" "" "" "" "" "";
            "Data"   string(round(1.91, digits=2))   string(round(1.80, digits=2))   string(round(1.91, digits=2))   string(round(1.87, digits=2))   string(round(1.93, digits=2))   string(round(2.03, digits=2));
            "Model"  string(round(avgfert, digits=2)) string(round(AvgF5[1], digits=2)) string(round(AvgF5[2], digits=2)) string(round(AvgF5[3], digits=2)) string(round(AvgF5[4], digits=2)) string(round(AvgF5[5], digits=2));
            "" "" "" "" "" "" "";
            "Childlessness rate (percent)" "" "" "" "" "" "";
            "Data"   string(round(2.9, digits=2))    string(round(5.3, digits=2))    string(round(4.0, digits=2))    string(round(2.0, digits=2))    string(round(1.3, digits=2))    string(round(2.0, digits=2));
            "Model"  string(round(100 * prn[1], digits=1))  string(round(100 * AvgN05[1], digits=1))  string(round(100 * AvgN05[2], digits=1))  string(round(100 * AvgN05[3], digits=1))  string(round(100 * AvgN05[4], digits=1))  string(round(100 * AvgN05[5], digits=1));
            "" "" "" "" "" "" "";
            "Education spending per child rel. to income (percent)" "" "" "" "" "" "";
            "Data"   string(round(9.2, digits=2))    string(round(11.2, digits=2))   string(round(9.9, digits=2))    string(round(9.3, digits=2))    string(round(8.7, digits=2))    string(round(6.9, digits=2));
            "Model"  string(round(100 * avg_educinv_inc, digits=1))  string(round(100 * invinc5[1], digits=1))  string(round(100 * invinc5[2], digits=1))  string(round(100 * invinc5[3], digits=1))  string(round(100 * invinc5[4], digits=1))  string(round(100 * invinc5[5], digits=1));
        ]

        hlines = [1, 4, 7, 10]

        latex_path_2 = joinpath(m, "output/Table5.tex")

        open(latex_path_2, "w") do io
            println(io, "\\begin{center}")
            println(io, "\\textbf{Table 5—Fertility and Education Spending across Income Quintiles}")
            println(io, "\\end{center}")
            println(io, "")

            pretty_table(io,
                data;
                header = header,
                backend = Val(:latex),
                tf = tf_latex_booktabs,
                alignment = :r,
                hlines = hlines
            )
        end
        

        ####################################################################
        # 7.  Save baseline results (MAT-style .mat file)                   #
        ####################################################################
        # Choose one of the two approaches below.  Comment out the other.  #
        ####################################################################

        # --- (a)  Save to a .mat file using MAT.jl ------------------------
        # using MAT                                            # (put at top of file)
        # matwrite("Baseline.mat"; objval, prn, Giniinc, AvgF5, AvgN05,
        #         invinc5, avg_LS, elast0, elastx, ige, avgfert, φgrid,
        #         parameters = (Ah, psicsd, psil, λ, φ1, φ2add, φ3add, γ,
        #                       κρ, κσ, α1, θ0, χ, τx, ω))

        # --- (b)  Save to a JLD2 file (native Julia) ----------------------
        # using JLD2                                           # (put at top of file)
        # @save "Baseline.jld2" objval prn Giniinc AvgF5 AvgN05 invinc5
        #                  avg_LS elast0 elastx ige avgfert φgrid
        #                  Ah psicsd psil λ φ1 φ2add φ3add γ κρ κσ α1 θ0 χ τx ω
    end

    return objval
end



################################################################################
# 5. Run the main function
################################################################################
gamma = 2.0;
lambda = 0.057*18/25; 

chi = 0.078211519;      
theta0 = 1.0787131;
Ah = 2.2254829;         
psil = 1.8678885;   
psicsd = 0.5865557;
phi1 = 1.3913291;
phi2add = 0.78402234;
phi3add = 0.39724406;
kappasd = 0.36132727;
alpha1 =  0.66838201;   
kapparho = 0.37556873;

taux= 0.0;
omega = 0.0;
eta = 0.0;

Prog_MainSS(1,Ah,psicsd,psil,lambda,phi1,phi2add,phi3add,gamma,kapparho,kappasd,alpha1,theta0,chi,taux,omega,0);


    

end

end