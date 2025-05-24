module Custom

#using Pkg
#Pkg.add("Distributions")

using Distributions, DataFrames

export HKprodfn, GSS, HES, qualutil10, Trn, lini, ols1, wls1, ginicoeff, tauchenmethod, gengrid, assign_cohort, mean_ci, age_to_cycle, _bin_index, rank_quintiles!, xtile_ado


    #----------------------------------------------------------------------------
    #For Table 4-5 Solve.m
    #----------------------------------------------------------------------------

    function HKprodfn(kappanow, hpnow, x, theta0, alpha1, alpha2)
        
        return hnext = kappanow .* (theta0 .+ (x .^ alpha1))
    end

    function GSS(fx, a::Float64, b::Float64, tol::Float64)
        φ = (sqrt(5) - 1) / 2      # golden ratio conjugate
        c = b - φ*(b - a)
        d = a + φ*(b - a)
        fc = fx(c)
        fd = fx(d)
        while (b - a) > tol
            if fc < fd
                b, fd = d, fc
                d = c
                c = b - φ*(b - a)
                fc = fx(c)
            else
                a, fc = c, fd
                c = d
                d = a + φ*(b - a)
                fd = fx(d)
            end
        end
        return (a + b) / 2
    end

    function HES(c, n)

        return sc_cons = c / ((1.5 + 0.3 * (n)) / 1.5)
    end

    function qualutil10(hnext, hbar, chi, sigma)

        return q = log(hnext .- chi .* hbar)
    end

    function Trn(omega, n) 

        return Tr = omega .* n
    end

    function lini(t, ft, z)

        # Linear interpolation evaluating a point z
        # t: grid (knots), assumed to be sorted ascending
        # ft: function values at those grid points
        # z: query point
        # Output: linearly interpolated value f(z)
        
        Nt = length(t)
        
        if z .<= t[1]
            return ft[1]
        elseif z .>= t[Nt]
            return ft[Nt]
        else
            ind = count(x -> x < z, t)
            w = (t[ind+1] - z) ./ (t[ind+1] .- t[ind])
            return w .* ft[ind] .+ (1 .- w) .* ft[ind+1]
        end
        
    end

    function ols1(yvec::Vector{<:Real}, xvec::Vector{<:Real})
        # Ordinary Least Squares (OLS) estimate of y on x
        # Returns the slope coefficient (betâ₁)
        
        NS = length(xvec)
        X = ones(NS, 2)
        X[:, 2] .= xvec
        
        # β = (X'X)⁻¹ X'y
        β = (X' * X) ./ (X' * yvec)
        
        return β[2]  # Return slope coefficient only
    end

    function wls1(yvec::Vector{<:Real}, xvec::Vector{<:Real}, Wmat::AbstractMatrix{<:Real})
        # Univariate weighted least squares regression: y = β₀ + β₁·x
        # yvec, xvec: column vectors of data
        # Wmat: NS×NS weight matrix (typically diagonal)
        # Returns: slope coefficient β₁
    
        NS = length(xvec)
        X = ones(NS, 2)
        X[:, 2] .= xvec
    
        β = (X' * Wmat * X) ./ (X' * Wmat * yvec)
    
        return β[2]  # Return slope coefficient
    end

    function ginicoeff(In::AbstractArray; dim::Int = 1, nosamplecorr::Bool = false)
        @assert dim == 1 || dim == 2 "dim must be 1 (columns) or 2 (rows)"
    
        # Copy and turn missing into NaN
        A = copy(In)
        if eltype(A) <: Union{Missing, Real}
            A = coalesce.(A, NaN)
        end
    
        # If vector, make it a column and fix dim
        if ndims(A) == 1
            A = reshape(A, :, 1)
            dim = 1
        end
    
        # NaN handling
        IDXnan = isnan.(A)
        notNaN = .!IDXnan
    
        # Invalid slices: any negative values or sum < 2
        invalid = mapslices(x -> any(x .< 0) || sum(x) < 2, A; dims=dim)
        IDX     = vec(invalid)
    
        # Zero-out invalid and NaN entries
        if dim == 1
            A[:, IDX] .= 0
        else
            A[IDX, :] .= 0
        end
        A[IDXnan] .= 0
    
        # Sort along each slice
        A = mapslices(x -> sort(x), A; dims=dim)
    
        # Build the “rank mask” by reversing, doing a cumsum, and reversing back
        rankmask = mapslices(x -> reverse(cumsum(reverse(x))), notNaN; dims=dim)
    
        # Totals per slice
        totNum   = dim == 1 ? rankmask[1, :] : rankmask[:, 1]
        tot      = sum(A; dims=dim)
        weighted = sum(A .* rankmask; dims=dim)
    
        # Compute Gini numerator/denominator
        coeff = totNum .+ 1 .- 2 .* (weighted ./ tot)
        if nosamplecorr
            coeff ./= totNum
        else
            coeff ./= (totNum .- 1)
        end
    
        return vec(coeff), IDX
    end

    function gengrid(xlow::Real, xhigh::Real, xnum::Int; option::Int = 0, degree::Real = 1)
        @assert xhigh ≥ xlow        "xhigh ≥ xlow"
        @assert xnum  ≥ 2           "Need two points"
        @assert option in (0, 1, 2) "option must be 0, 1 or 2"
    
        if option == 0                       # 1) linear
            return LinRange(xlow, xhigh, xnum)
    
        elseif option == 1                   # 2) log-spaced (shift & stretch)
            span = xhigh - xlow + 1          
            a, b = 0.0, log10(span)          
            return 10.0 .^ LinRange(a, b, xnum) .+ (xlow - 1)
    
        else                                 # 3) log-spaced with custom degree
            raw = 10.0 .^ LinRange(0, degree, xnum) .- 1
            raw ./= (10.0^degree - 1)        # normalize to [0, 1]
            return raw .* (xhigh - xlow) .+ xlow
        end
    end

    function tauchenmethod(mew, sigma, rho, znum, q; Parallel=0, Verbose=0)
        znum   = Int(znum)                         # force integer
        sigmaz = sigma / sqrt(1 - rho^2)
        # state grid
        z      = mew .+ range(-q*sigmaz, stop=q*sigmaz, length=znum)
        ω      = step(z)

        # transition matrix
        zi = repeat(z,  1, znum)
        zj = repeat(z', znum, 1)
        P1 = cdf.(Normal(), ((zj .+ ω/2 .- rho .* zi) .- mew) ./ sigma)
        P2 = cdf.(Normal(), ((zj .- ω/2 .- rho .* zi) .- mew) ./ sigma)
        P  = P1 .- P2

        P[:, 1]   .= P1[:, 1]
        P[:, end] .= 1 .- P2[:, end]

        return z, P
    end

    #----------------------------------------------------------------------------
    #For Table 2, 3_educost_stats.do
    #----------------------------------------------------------------------------

    """
    Function to do line 73 Figure 1 and line __ Table 2
    """
    function assign_cohort(
        birth_year::Real,
        cohort0min::Real, cohort0max::Real,
        cohort1min::Real, cohort1max::Real)

    # Propagate missing
    ismissing(birth_year) && return missing

    # Accept only four-digit years
    birth_year < 1000 && return missing         

        if cohort0min ≤ birth_year ≤ cohort0max
            return 0
        elseif cohort1min ≤ birth_year ≤ cohort1max
            return 1
        else
            return missing          
        end
    end


    # if birth_year is missing, just return missing
    assign_cohort(::Missing,
            ::Real, ::Real,
            ::Real, ::Real) = missing


    """
    Funtion to obtain confidence mean and confidence interval around it (as in line 134) 
    """
    function mean_ci(x; alpha = 0.05)
        data = collect(skipmissing(x))          # <- materialise the iterator
        n = length(data)
        if n == 0
            return (mean = missing, lb = missing, ub = missing)
        end
        μ  = mean(data)
        se = std(data) / sqrt(n)
        z  = quantile(Normal(), 1 - alpha/2)
        (mean = μ,
         lb   = μ - z*se,
        ub   = μ + z*se)
    end 

    #cycle_stats(df, stage; quintile_col = :inc_both1, cohort = 1)

    

    function age_to_cycle(x)
        ismissing(x) && return missing
        x ≤ 6  && return 1
        x ≤ 12 && return 2
        x ≤ 15 && return 3
        x ≤ 18 && return 4
        return missing      # 19 or older or missing
    end

    function _bin_index(x::Real, edges::AbstractVector{<:Real}) #tolot @inline
        idx = searchsortedlast(edges, x)         
        return min(idx, length(edges)-1)         
    end

    function rank_quintiles!(vec::AbstractVector, ng::Int = 5)
        codes = Vector{Union{Missing,Int}}(missing, length(vec))
    
        idx_nonmissing = findall(!ismissing, vec)      
        k = length(idx_nonmissing)
        k == 0 && return codes
    
        vals = float.(vec[idx_nonmissing])            
        local_ord = sortperm(vals)                    
        ord = idx_nonmissing[local_ord]               
    
        step   = fld(k, ng)
        extras = k - step * ng
        cut = cumsum(step .+ (1:ng .<= extras))       
    
        for (rank, idx) in enumerate(ord)
            codes[idx] = searchsortedfirst(cut, rank)   
        end
        codes
    end

    """
    xtile_ado(v::AbstractVector{<:Union{Missing,Real}}, ng::Int)

    Exactly replicates Stata’s `xtile … nq(ng)` (no weights).  Missing entries stay missing;
    non-missings are cut into ng groups of (almost) equal size, ties at boundaries are
    assigned to the *higher* bin, just as in Stata’s ado file.  Returns a Vector{Union{Missing,Int}}.
    """
    function xtile_ado(v::AbstractVector{<:Union{Missing,Real}}, ng::Int)
        # 1. Extract nonmissing
        idx = findall(!ismissing, v)
        N = length(idx)
        codes = Vector{Union{Missing,Int}}(missing, length(v))
        if N == 0
            return codes
        end

        # 2. Determine group sizes: first `extra` bins get one extra observation
        base  = fld(N, ng)
        extra = N - base*ng
        sizes = [ base + (i <= extra ? 1 : 0) for i in 1:ng ]

        # 3. Sort the nonmissing values
        vals = Float64.(v[idx])
        ord  = sortperm(vals)         # positions in idx, in ascending order
        sorted = vals[ord]

        # 4. Compute the actual boundary values
        boundaries = Float64[]
        pos = 0
        for sz in sizes[1:end-1]      # we need ng-1 cutpoints
            pos += sz
            push!(boundaries, sorted[pos])
        end

        # 5. Build the edges vector: [-Inf, b1, b2, …, b_{ng-1}, Inf]
        edges = [-Inf; boundaries...; Inf]

        # 6. For each original nonmissing, find its bin g:
        #    edges[g] ≤ x < edges[g+1], but if x == a boundary, assign to the *higher* bin.
        for ii in idx
            x = float(v[ii])
            g = searchsortedlast(edges, x)
            # tie on boundary? bump up
            if g > 1 && x == edges[g]
                g += 1
            end
            codes[ii] = min(g, ng)
        end

        return codes
    end

    

end #module




