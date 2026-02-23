# Estimating the spark motion in two dimensional elliptical polar cap
#
# Julia translation of elipsDrift.c
#
# Usage:
#   julia elipsDrift.jl  ntime  incl_angle(deg)  maj_axis(m)  ellip(b/a)  corot_angle(deg)
#
# Dependencies:
#   julia> import Pkg; Pkg.add("GLMakie")

using GLMakie


# ---------------------------------------------------------------------------
# 2D coordinate rotation (inverse / transpose rotation matrix)
#   [cos θ   sin θ] [x]
#   [-sin θ  cos θ] [y]
# This is equivalent to rotating the axes by +θ (or the point by -θ).
# ---------------------------------------------------------------------------
function coortrans(x_in::Float64, y_in::Float64, theta::Float64)
    x_out =  x_in * cos(theta) + y_in * sin(theta)
    y_out = -x_in * sin(theta) + y_in * cos(theta)
    return x_out, y_out
end


# ---------------------------------------------------------------------------
# Build spark positions and rasterize onto a fine grid.
#
# Returns (x_pts, y_pts): coordinates of every grid cell that falls inside
# at least one spark ellipse — used directly for plotting.
#
# Layout of th_sprk_u / th_sprk_d (flat 1-D arrays, 1-indexed):
#   ring ii (1-based) occupies indices  (ii-1)*trk_max+1 .. (ii-1)*trk_max+trk_max
# ---------------------------------------------------------------------------
function sparkconfig(th_sprk_u, th_sprk_d, N_up, N_dn, theta_sp,
                     h_sprk, h_drft, a_cap, b_cap, th_cap,
                     co_angl, x_cent, y_cent, N_trk, trk_max)

    x_sprk   = Float64[]
    y_sprk   = Float64[]
    a_sprk_v = Float64[]   # semi-axis of each spark (in rotated frame)

    a_out = a_cap
    a_in  = a_out - 2.0 * h_sprk
    b_out = b_cap
    b_in  = b_out - b_cap / a_cap * 2.0 * h_sprk

    for ii in 1:N_trk

        a_trk = 0.5 * (a_out + a_in)
        b_trk = 0.5 * (b_out + b_in)

        # Flat-array offset for ring ii (0-based shift, Julia array is 1-based)
        u_off = (ii - 1) * trk_max
        d_off = (ii - 1) * trk_max

        # ---- Upper half: clockwise track (θ decreasing from π toward 0) ----
        for jj in 1:N_up[ii]

            ang = th_sprk_u[u_off + jj]
            xi  = a_trk * cos(ang - th_cap)
            yi  = b_trk * sin(ang - th_cap)
            xs, ys = coortrans(xi, yi, -th_cap)
            xs += x_cent;  ys += y_cent
            a_s = h_sprk

            # Spark near the leading edge (θ ≈ π)
            if π - co_angl - ang <= theta_sp[ii] / 2
                theta_fh = 0.5 * (π - co_angl - ang) + theta_sp[ii] / 4
                a_s   = a_trk * sin(theta_fh)
                trk_a = a_trk + h_sprk - a_s
                xi = trk_a * cos(π - co_angl - theta_fh - th_cap)
                yi = trk_a * b_trk / a_trk * sin(π - co_angl - theta_fh - th_cap)
                xs, ys = coortrans(xi, yi, -th_cap)
                xs += x_cent;  ys += y_cent
            end

            # Spark near the trailing edge (θ ≈ 0)
            if ang <= theta_sp[ii] / 2 - co_angl
                theta_fh = 0.5 * (ang + theta_sp[ii] / 2 + co_angl)
                a_s   = a_trk * sin(theta_fh)
                trk_a = a_trk + h_sprk - a_s
                xi = trk_a * cos(theta_fh - th_cap - co_angl)
                yi = trk_a * b_trk / a_trk * sin(theta_fh - th_cap - co_angl)
                xs, ys = coortrans(xi, yi, -th_cap)
                xs += x_cent;  ys += y_cent
            end

            push!(x_sprk, xs);  push!(y_sprk, ys);  push!(a_sprk_v, a_s)
        end

        # ---- Lower half: anti-clockwise track (θ increasing from π toward 2π) ----
        for jj in 1:N_dn[ii]

            ang = th_sprk_d[d_off + jj]
            xi  = a_trk * cos(ang - th_cap)
            yi  = b_trk * sin(ang - th_cap)
            xs, ys = coortrans(xi, yi, -th_cap)
            xs += x_cent;  ys += y_cent
            a_s = h_sprk

            # Spark near the leading edge (θ ≈ π)
            if ang - π + co_angl <= theta_sp[ii] / 2
                theta_fl = 0.5 * (ang - π + co_angl) + theta_sp[ii] / 4
                a_s   = 0.5 * (2 * a_trk * sin(theta_fl))
                trk_a = a_trk + h_sprk - a_s
                xi = trk_a * cos(π - co_angl + theta_fl - th_cap)
                yi = trk_a * b_trk / a_trk * sin(π - co_angl + theta_fl - th_cap)
                xs, ys = coortrans(xi, yi, -th_cap)
                xs += x_cent;  ys += y_cent
            end

            # Spark near the trailing edge (θ ≈ 2π)
            if 2π - ang - co_angl <= theta_sp[ii] / 2
                theta_fl = 0.5 * (2π - ang - co_angl) + theta_sp[ii] / 4
                a_s   = 0.5 * (2 * a_trk * sin(theta_fl))
                trk_a = a_trk + h_sprk - a_s
                xi = trk_a * cos(2π - theta_fl - th_cap - co_angl)
                yi = trk_a * b_trk / a_trk * sin(2π - theta_fl - th_cap - co_angl)
                xs, ys = coortrans(xi, yi, -th_cap)
                xs += x_cent;  ys += y_cent
            end

            push!(x_sprk, xs);  push!(y_sprk, ys);  push!(a_sprk_v, a_s)
        end

        # ---- Plug gap at the beginning (θ = π junction) ----
        if th_sprk_d[d_off + 1] - th_sprk_u[u_off + 1] >= theta_sp[ii]
            gap   = 0.5 * (th_sprk_d[d_off + 1] - th_sprk_u[u_off + 1])
            a_s   = 0.5 * (2 * a_trk * sin(gap) - 2.0 * h_sprk)
            e_trk = a_trk + h_sprk - a_s
            xi = e_trk * cos(π - th_cap - co_angl)
            yi = e_trk * b_trk / a_trk * sin(π - th_cap - co_angl)
            xs, ys = coortrans(xi, yi, -th_cap)
            push!(x_sprk, xs + x_cent);  push!(y_sprk, ys + y_cent);  push!(a_sprk_v, a_s)
        end

        # ---- Plug gap at the end (θ = 0 / 2π junction) ----
        theta_fh = th_sprk_u[u_off + N_up[ii]]          # last angle of upper track
        theta_fl = 2π - th_sprk_d[d_off + N_dn[ii]]     # last angle of lower track (reflected)
        if theta_sp[ii] <= theta_fl + theta_fh <= 2 * theta_sp[ii]
            a_s   = 0.5 * (2 * a_trk * sin(0.5 * (theta_fl + theta_fh)) - 2.0 * h_sprk)
            e_trk = a_trk + h_sprk - a_s
            xi = e_trk * cos(2π - th_cap - co_angl)
            yi = e_trk * b_trk / a_trk * sin(2π - th_cap - co_angl)
            xs, ys = coortrans(xi, yi, -th_cap)
            push!(x_sprk, xs + x_cent);  push!(y_sprk, ys + y_cent);  push!(a_sprk_v, a_s)
        end

        a_out -= 2.0 * h_sprk;  a_in -= 2.0 * h_sprk
        b_out -= b_cap / a_cap * 2.0 * h_sprk
        b_in  -= b_cap / a_cap * 2.0 * h_sprk
    end

    # ---- Core spark (central fill) ----
    push!(x_sprk, x_cent);  push!(y_sprk, y_cent)
    push!(a_sprk_v, a_cap - N_trk * 2.0 * h_sprk - 1.5 * h_drft)

    ncnt = length(x_sprk)

    # ---- Rasterize: collect every grid cell that falls inside any spark ----
    # Grid step h_drft, covering [h_drft/2, 2·x_cent) × [h_drft/2, 2·y_cent)
    x_pts = Float64[]
    y_pts = Float64[]

    x_val = h_drft / 2
    while x_val < 2 * x_cent
        y_val = h_drft / 2
        while y_val < 2 * y_cent
            for k in 1:ncnt
                a_s = a_sprk_v[k]
                b_s = a_s * b_cap / a_cap
                if a_s > 0 && b_s > 0
                    dx, dy = x_val - x_sprk[k], y_val - y_sprk[k]
                    xt, yt = coortrans(dx, dy, th_cap)
                    if sqrt(xt^2 / a_s^2 + yt^2 / b_s^2) < 1.0
                        push!(x_pts, x_val)
                        push!(y_pts, y_val)
                        # Note: no break — matches C original (duplicates harmless for plot)
                    end
                end
            end
            y_val += h_drft
        end
        x_val += h_drft
    end

    return x_pts, y_pts
end


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
function main()
    #=
    if length(ARGS) < 5
        println(stderr,
            "usage: julia elipsDrift.jl  ntime  incl_angle(deg)  " *
            "maj_axis(m)  ellip(b/a)  corot_angle(deg)")
        exit(1)
    end

    ntime   = parse(Int,     ARGS[1])
    th_cap  = parse(Float64, ARGS[2]) * π / 180   # inclination angle [rad]
    a_cap   = parse(Float64, ARGS[3])              # polar cap major semi-axis [m]
    b_cap   = a_cap * parse(Float64, ARGS[4])      # minor semi-axis [m]
    co_angl = parse(Float64, ARGS[5]) * π / 180   # corotation angle [rad]
    =#
    #200 0 20 1.0 0
    ntime = 200
    th_cap = 0.0
    a_cap = 20.0
    b_cap = 20.0
    co_angl = 0.0




    # Physical parameters (hard-coded as in C original)
    h_sprk = 2.6   # spark diameter [m]
    h_drft = 0.1   # drift step size per iteration [m]

    a_sprk = h_sprk
    b_sprk = a_sprk * b_cap / a_cap

    # Centre of the polar cap in rotated display coordinates
    x_cent = sqrt((a_cap * cos(th_cap))^2 + (b_cap * sin(th_cap))^2)
    y_cent = sqrt((a_cap * sin(th_cap))^2 + (b_cap * cos(th_cap))^2)

    # Number of concentric spark rings
    N_trk = floor(Int, b_cap / (2 * b_sprk))   # = floor(a_cap / (2*h_sprk))

    println(stderr, "N_trk = $N_trk   a_cap = $a_cap   b_cap = $b_cap")

    # Angular spacing between sparks for each ring
    theta_sp = zeros(Float64, N_trk)
    let a_o = a_cap, a_i = a_o - 2*a_sprk,
        b_o = b_cap, b_i = b_o - 2*b_sprk
        for ii in 1:N_trk
            N_s = floor(Int, 0.75 * (a_o*b_o - a_i*b_i) / (a_sprk*b_sprk))
            theta_sp[ii] = 2π / N_s
            a_o -= 2*a_sprk;  a_i -= 2*a_sprk
            b_o -= 2*b_sprk;  b_i -= 2*b_sprk
        end
    end

    # Maximum sparks per half-ring (stride in flat angle arrays)
    trk_max = floor(Int,
        0.75 * (a_cap*b_cap - (a_cap - 2*a_sprk)*(b_cap - 2*b_sprk)) /
        (a_sprk * b_sprk) / 2) + 1

    # Flat arrays: ring ii (1-based) uses indices (ii-1)*trk_max+1 .. (ii-1)*trk_max+trk_max
    sz = 2 * trk_max * N_trk + 2
    th_sprk_u = zeros(Float64, sz)   # upper-half spark angles
    th_sprk_d = zeros(Float64, sz)   # lower-half spark angles
    N_up = zeros(Int, N_trk)         # number of sparks on upper track per ring
    N_dn = zeros(Int, N_trk)         # number of sparks on lower track per ring

    # Initialise spark positions
    for ii in 1:N_trk
        u_off = (ii - 1) * trk_max
        d_off = (ii - 1) * trk_max

        th_sprk_u[u_off + 1] = π - theta_sp[ii] / 2 - co_angl
        N_up[ii] = 1
        while th_sprk_u[u_off + N_up[ii]] >= theta_sp[ii] - co_angl
            th_sprk_u[u_off + N_up[ii] + 1] = th_sprk_u[u_off + N_up[ii]] - theta_sp[ii]
            N_up[ii] += 1
        end

        th_sprk_d[d_off + 1] = π + theta_sp[ii] / 2 - co_angl
        N_dn[ii] = 1
        while th_sprk_d[d_off + N_dn[ii]] <= 2π - theta_sp[ii] - co_angl
            th_sprk_d[d_off + N_dn[ii] + 1] = th_sprk_d[d_off + N_dn[ii]] + theta_sp[ii]
            N_dn[ii] += 1
        end
    end

    # Ellipse outline for the polar cap boundary
    ts      = range(0.0, 2π; length = 101)
    x_elips = [coortrans(a_cap*cos(t), b_cap*sin(t), -th_cap)[1] + x_cent for t in ts]
    y_elips = [coortrans(a_cap*cos(t), b_cap*sin(t), -th_cap)[2] + y_cent for t in ts]

    # ---- Set up GLMakie interactive window ----
    fig = Figure(size = (600, 600))
    title_str = Observable("Iteration # 1")
    ax = Axis(fig[1, 1];
              xlabel  = "X (m)",
              ylabel  = "Y (m)",
              title   = title_str,
              aspect  = DataAspect())
    xlims!(ax, 0.0, 2 * x_cent)
    ylims!(ax, 0.0, 2 * y_cent)

    lines!(ax, x_elips, y_elips; color = :black, linewidth = 2)

    pts = Observable(Point2f[])
    scatter!(ax, pts; color = :steelblue, markersize = 3, marker = :circle)

    display(fig)

    # Drift step (constant across all rings — faithful to C original where
    # a_out/a_in are not updated inside the jj loop)
    del_theta_drift = h_drft / (0.5 * (a_cap + (a_cap - 2 * a_sprk)))

    # ---- Time evolution loop ----
    for step in 1:ntime

        # Compute spark layout and rasterized grid points
        x_pts, y_pts = sparkconfig(th_sprk_u, th_sprk_d, N_up, N_dn, theta_sp,
                                    h_sprk, h_drft, a_cap, b_cap, th_cap,
                                    co_angl, x_cent, y_cent, N_trk, trk_max)

        # Update plot
        pts[]       = Point2f.(x_pts, y_pts)
        title_str[] = "Iteration # $step"
        sleep(0.05)   # ≈ 50 ms frame delay (matches C's delay(50))

        # Advance spark angles by one drift step
        for jj in 1:N_trk
            u_off = (jj - 1) * trk_max
            d_off = (jj - 1) * trk_max

            # Upper track drifts clockwise (angle decreases)
            th_sprk_u[u_off + 1] -= del_theta_drift
            if th_sprk_u[u_off + 1] < π - theta_sp[jj] - co_angl
                th_sprk_u[u_off + 1] += theta_sp[jj]
            end
            N_up[jj] = 1
            while th_sprk_u[u_off + N_up[jj]] >= theta_sp[jj] - co_angl
                th_sprk_u[u_off + N_up[jj] + 1] = th_sprk_u[u_off + N_up[jj]] - theta_sp[jj]
                N_up[jj] += 1
            end

            # Lower track drifts anti-clockwise (angle increases)
            th_sprk_d[d_off + 1] += del_theta_drift
            if th_sprk_d[d_off + 1] > π + theta_sp[jj] - co_angl
                th_sprk_d[d_off + 1] -= theta_sp[jj]
            end
            N_dn[jj] = 1
            while th_sprk_d[d_off + N_dn[jj]] <= 2π - theta_sp[jj] - co_angl
                th_sprk_d[d_off + N_dn[jj] + 1] = th_sprk_d[d_off + N_dn[jj]] + theta_sp[jj]
                N_dn[jj] += 1
            end
        end
    end

    println("Animation complete. Close the window to exit.")
    while events(fig).window_open[]
        sleep(0.1)
    end
end

main()
