# =============================================================================
# sparkfunc.jl
#
# Julia translation of sparkfunc.h (and the helper routines from usefunc.h
# that sparkfunc.h depends on).
#
# Provides:
#   Geometry helpers
#     rotate_coordinates_2d       — 2D coordinate rotation
#     interpolate_between         — linear interpolation of a mapped value
#     find_edge_position          — inverse linear interpolation (find x given y)
#
#   3D vector algebra (neutron-star surface geometry)
#     normalize_vector            — scale a vector to unit length
#     cross_product               — 3D cross product
#     dot_product                 — 3D dot product
#     rotate_vector               — Rodrigues' rotation around an arbitrary axis
#
#   Polar-cap spark physics
#     compute_spark_configuration — place sparks on concentric elliptical rings
#     project_los_to_polar_cap    — map a pulsar phase angle to a 2D cap position
#     compute_spark_amplitude     — evaluate the Gaussian signal at a cap point
# =============================================================================


# =============================================================================
# Geometry helpers
# (translated from coortrans, findbet, findedge in usefunc.h)
# =============================================================================

"""
    rotate_coordinates_2d(x, y, angle) -> (x_rotated, y_rotated)

Apply a 2D coordinate rotation by `angle` (radians) using the transpose
(inverse) rotation matrix:

    [ cos(angle)   sin(angle) ] [ x ]
    [-sin(angle)   cos(angle) ] [ y ]

This is equivalent to rotating the coordinate axes forward by +angle, which
transforms a point from a frame tilted by angle back to the standard frame.

Used throughout to switch between the tilted polar-cap frame and the flat
display frame.
"""
function rotate_coordinates_2d(x::Float64, y::Float64, angle::Float64)
    x_rotated =  x * cos(angle) + y * sin(angle)
    y_rotated = -x * sin(angle) + y * cos(angle)
    return x_rotated, y_rotated
end


"""
    interpolate_between(input_start, input_end, input_query,
                        output_start, output_end) -> Float64

Linear interpolation: given that `input_query` lies between `input_start`
and `input_end`, return the value at the same fractional position between
`output_start` and `output_end`.

    fraction = (input_query - input_start) / (input_end - input_start)
    result   = output_start + fraction * (output_end - output_start)

Used in `project_los_to_polar_cap` to interpolate the azimuthal cap angle
along the tabulated line-of-sight track.
"""
function interpolate_between(input_start::Float64, input_end::Float64,
                              input_query::Float64,
                              output_start::Float64, output_end::Float64)
    fraction = (input_query - input_start) / (input_end - input_start)
    return output_start + fraction * (output_end - output_start)
end


"""
    find_edge_position(x1, y1, x2, y2, y_query) -> Float64

Given two points (x1, y1) and (x2, y2) on a line, return the x-value
that corresponds to `y_query`.

Solves the line equation for x:
    slope     = (y2 - y1) / (x2 - x1)
    intercept = y2 - slope * x2
    result    = (y_query - intercept) / slope

Used in `project_los_to_polar_cap` to recover the polar colatitude (theta)
of the line-of-sight boundary at a given azimuthal angle.
"""
function find_edge_position(x1::Float64, y1::Float64,
                             x2::Float64, y2::Float64,
                             y_query::Float64)
    slope     = (y2 - y1) / (x2 - x1)
    intercept = y2 - slope * x2
    return (y_query - intercept) / slope
end


# =============================================================================
# 3D vector algebra
# (translated from unitvect, crossprod, dotprod, rotvect in sparkfunc.h)
# =============================================================================

"""
    normalize_vector(input_vector) -> unit_vector

Return a unit vector (length 1) pointing in the same direction as
`input_vector`.  `input_vector` must have exactly 3 elements.

Used by `rotate_vector` and `project_los_to_polar_cap` to construct
a properly normalized rotation axis.
"""
function normalize_vector(input_vector::Vector{Float64})
    magnitude = sqrt(sum(component^2 for component in input_vector))
    return input_vector ./ magnitude
end


"""
    cross_product(vector_a, vector_b) -> result_vector

Return the cross product of two 3D vectors, i.e. a vector perpendicular
to both inputs whose magnitude equals the area of the parallelogram they span.

    result[1] = a[2]*b[3] - a[3]*b[2]
    result[2] = a[3]*b[1] - a[1]*b[3]
    result[3] = a[1]*b[2] - a[2]*b[1]

Used inside `rotate_vector` as part of Rodrigues' rotation formula.
"""
function cross_product(vector_a::Vector{Float64}, vector_b::Vector{Float64})
    return [
        vector_a[2] * vector_b[3] - vector_a[3] * vector_b[2],
        vector_a[3] * vector_b[1] - vector_a[1] * vector_b[3],
        vector_a[1] * vector_b[2] - vector_a[2] * vector_b[1],
    ]
end


"""
    dot_product(vector_a, vector_b) -> scalar

Return the scalar (dot) product of two 3D vectors:
    result = a[1]*b[1] + a[2]*b[2] + a[3]*b[3]

Used inside `rotate_vector` to extract the component of the input vector
that is parallel to the rotation axis.
"""
function dot_product(vector_a::Vector{Float64}, vector_b::Vector{Float64})
    return sum(vector_a[i] * vector_b[i] for i in 1:3)
end


"""
    rotate_vector(input_vector, rotation_axis, rotation_angle) -> rotated_vector

Rotate `input_vector` by `rotation_angle` (radians) around `rotation_axis`
using Rodrigues' rotation formula:

    v_rot = v · cos(θ)
          + (axis × v) · sin(θ)
          + axis · (axis · v) · (1 − cos(θ))

`rotation_axis` must be a unit vector (use `normalize_vector` beforehand).
All vectors are 3-element Float64 arrays.

Used in `project_los_to_polar_cap` to tilt a 3D surface position from the
neutron-star spherical coordinate frame into the flat 2D polar-cap frame.
"""
function rotate_vector(input_vector::Vector{Float64},
                       rotation_axis::Vector{Float64},
                       rotation_angle::Float64)
    perpendicular_part = cross_product(rotation_axis, input_vector)
    parallel_length    = dot_product(input_vector, rotation_axis)

    rotated_vector = [
        input_vector[i] * cos(rotation_angle) +
        perpendicular_part[i] * sin(rotation_angle) +
        rotation_axis[i] * parallel_length * (1.0 - cos(rotation_angle))
        for i in 1:3
    ]
    return rotated_vector
end


# =============================================================================
# Spark configuration on the elliptical polar cap
# (translated from sparkconfig in sparkfunc.h)
# =============================================================================

"""
    compute_spark_configuration(
        upper_track_angles, lower_track_angles,
        sparks_on_upper_track, sparks_on_lower_track,
        angular_spacing_per_ring,
        spark_diameter, drift_step,
        cap_major_semiaxis, cap_minor_semiaxis,
        cap_tilt_angle, corotation_angle,
        cap_centre_x, cap_centre_y,
        number_of_rings, track_stride
    ) -> (spark_x, spark_y, spark_semiaxis, number_of_sparks)

Compute the 2D centre positions and effective semi-axes of all sparks
on the elliptical polar cap for one time step.

# Geometry overview
Sparks fill the cap in `number_of_rings` concentric elliptical rings.
Each ring is one spark-diameter wide.  Within a ring, sparks are arranged
in two half-tracks:

  - Upper half-track (angles decreasing from π toward 0): clockwise drift
  - Lower half-track (angles increasing from π toward 2π): anti-clockwise

The polar cap has major semi-axis `cap_major_semiaxis` and minor semi-axis
`cap_minor_semiaxis`.  The ellipse is tilted by `cap_tilt_angle` in the
display frame.  `corotation_angle` shifts the whole drift pattern.

# Flat angle array layout
`upper_track_angles` and `lower_track_angles` are flat 1D arrays.
Ring `r` (1-based) occupies indices:
    (r - 1) * track_stride + 1  …  (r - 1) * track_stride + track_stride
`sparks_on_upper_track[r]` gives the number of active sparks on ring r.

# Edge clipping
When a spark is within half an angular spacing of the leading (θ ≈ π) or
trailing (θ ≈ 0 / 2π) edge of its half-track, its semi-axis is reduced
and its centre is shifted outward so it fits snugly against the cap boundary.

# Returns
  - `spark_x`          : x-coordinates of spark centres in the display frame
  - `spark_y`          : y-coordinates of spark centres in the display frame
  - `spark_semiaxis`   : effective major semi-axis of each spark
  - `number_of_sparks` : total number of sparks (length of the above vectors)
"""
function compute_spark_configuration(
        upper_track_angles::Vector{Float64},
        lower_track_angles::Vector{Float64},
        sparks_on_upper_track::Vector{Int},
        sparks_on_lower_track::Vector{Int},
        angular_spacing_per_ring::Vector{Float64},
        spark_diameter::Float64,
        drift_step::Float64,
        cap_major_semiaxis::Float64,
        cap_minor_semiaxis::Float64,
        cap_tilt_angle::Float64,
        corotation_angle::Float64,
        cap_centre_x::Float64,
        cap_centre_y::Float64,
        number_of_rings::Int,
        track_stride::Int)

    spark_x        = Float64[]
    spark_y        = Float64[]
    spark_semiaxis = Float64[]

    # The default (nominal) semi-axis of a spark is half the spark diameter.
    nominal_spark_semiaxis = spark_diameter / 2.0

    # Ratio of minor to major axis — same for every spark as for the cap itself.
    ellipticity = cap_minor_semiaxis / cap_major_semiaxis

    # Start from the outermost ring and step inward ring by ring.
    # Each ring is one spark_diameter wide in the radial direction.
    outer_major = cap_major_semiaxis
    inner_major = outer_major - spark_diameter

    outer_minor = cap_minor_semiaxis
    inner_minor = outer_minor - ellipticity * spark_diameter

    for ring in 1:number_of_rings

        # Semi-axes of the mid-line of the current annular ring.
        track_major = 0.5 * (outer_major + inner_major)
        track_minor = 0.5 * (outer_minor + inner_minor)

        # Offset into the flat angle arrays for this ring (Julia is 1-based).
        upper_offset = (ring - 1) * track_stride
        lower_offset = (ring - 1) * track_stride

        # ------------------------------------------------------------------
        # Upper half-ring: spark angles decrease from π toward 0.
        # ------------------------------------------------------------------
        for spark_index in 1:sparks_on_upper_track[ring]

            angle              = upper_track_angles[upper_offset + spark_index]
            effective_semiaxis = nominal_spark_semiaxis

            # Default spark centre on the track mid-line, in the cap frame.
            xi = track_major * cos(angle - cap_tilt_angle)
            yi = track_minor * sin(angle - cap_tilt_angle)

            centre_x, centre_y = rotate_coordinates_2d(xi, yi, -cap_tilt_angle)
            centre_x += cap_centre_x
            centre_y += cap_centre_y

            # ---- Overlap with the leading edge (angle near π) ----
            # When the spark is within half a spacing of the leading boundary,
            # shrink its semi-axis and shift its centre outward.
            leading_gap = π - corotation_angle - angle
            if leading_gap <= angular_spacing_per_ring[ring] / 2
                half_leading_gap   = 0.5 * leading_gap + angular_spacing_per_ring[ring] / 4
                effective_semiaxis = track_major * sin(half_leading_gap)
                shifted_radius     = track_major + nominal_spark_semiaxis - effective_semiaxis

                xi = shifted_radius * cos(π - corotation_angle - half_leading_gap - cap_tilt_angle)
                yi = shifted_radius * track_minor / track_major *
                         sin(π - corotation_angle - half_leading_gap - cap_tilt_angle)

                centre_x, centre_y = rotate_coordinates_2d(xi, yi, -cap_tilt_angle)
                centre_x += cap_centre_x
                centre_y += cap_centre_y
            end

            # ---- Overlap with the trailing edge (angle near 0) ----
            if angle <= angular_spacing_per_ring[ring] / 2 - corotation_angle
                half_trailing_gap  = 0.5 * (angle + angular_spacing_per_ring[ring] / 2 +
                                            corotation_angle)
                effective_semiaxis = track_major * sin(half_trailing_gap)
                shifted_radius     = track_major + nominal_spark_semiaxis - effective_semiaxis

                xi = shifted_radius * cos(half_trailing_gap - cap_tilt_angle - corotation_angle)
                yi = shifted_radius * track_minor / track_major *
                         sin(half_trailing_gap - cap_tilt_angle - corotation_angle)

                centre_x, centre_y = rotate_coordinates_2d(xi, yi, -cap_tilt_angle)
                centre_x += cap_centre_x
                centre_y += cap_centre_y
            end

            push!(spark_x,        centre_x)
            push!(spark_y,        centre_y)
            push!(spark_semiaxis, effective_semiaxis)
        end

        # ------------------------------------------------------------------
        # Lower half-ring: spark angles increase from π toward 2π.
        # ------------------------------------------------------------------
        for spark_index in 1:sparks_on_lower_track[ring]

            angle              = lower_track_angles[lower_offset + spark_index]
            effective_semiaxis = nominal_spark_semiaxis

            # Default spark centre on the track mid-line.
            xi = track_major * cos(angle - cap_tilt_angle)
            yi = track_minor * sin(angle - cap_tilt_angle)

            centre_x, centre_y = rotate_coordinates_2d(xi, yi, -cap_tilt_angle)
            centre_x += cap_centre_x
            centre_y += cap_centre_y

            # ---- Overlap with the leading edge (angle near π) ----
            leading_gap = angle - π + corotation_angle
            if leading_gap <= angular_spacing_per_ring[ring] / 2
                half_leading_gap   = 0.5 * leading_gap + angular_spacing_per_ring[ring] / 4
                effective_semiaxis = track_major * sin(half_leading_gap)
                shifted_radius     = track_major + nominal_spark_semiaxis - effective_semiaxis

                xi = shifted_radius * cos(π - corotation_angle + half_leading_gap - cap_tilt_angle)
                yi = shifted_radius * track_minor / track_major *
                         sin(π - corotation_angle + half_leading_gap - cap_tilt_angle)

                centre_x, centre_y = rotate_coordinates_2d(xi, yi, -cap_tilt_angle)
                centre_x += cap_centre_x
                centre_y += cap_centre_y
            end

            # ---- Overlap with the trailing edge (angle near 2π) ----
            trailing_gap = 2π - angle - corotation_angle
            if trailing_gap <= angular_spacing_per_ring[ring] / 2
                half_trailing_gap  = 0.5 * trailing_gap + angular_spacing_per_ring[ring] / 4
                effective_semiaxis = track_major * sin(half_trailing_gap)
                shifted_radius     = track_major + nominal_spark_semiaxis - effective_semiaxis

                xi = shifted_radius * cos(2π - half_trailing_gap - cap_tilt_angle - corotation_angle)
                yi = shifted_radius * track_minor / track_major *
                         sin(2π - half_trailing_gap - cap_tilt_angle - corotation_angle)

                centre_x, centre_y = rotate_coordinates_2d(xi, yi, -cap_tilt_angle)
                centre_x += cap_centre_x
                centre_y += cap_centre_y
            end

            push!(spark_x,        centre_x)
            push!(spark_y,        centre_y)
            push!(spark_semiaxis, effective_semiaxis)
        end

        # ------------------------------------------------------------------
        # Gap filling at the θ = π junction (between upper and lower tracks).
        # If the angular gap between the first upper-track spark and the first
        # lower-track spark is too wide, a bridging spark is inserted there.
        #
        # NOTE: this block is commented out in the original sparkfunc.h.
        #       Included here for completeness and future use.
        # ------------------------------------------------------------------
        #=
        junction_angular_gap = lower_track_angles[lower_offset + 1] -
                               upper_track_angles[upper_offset + 1]
        if junction_angular_gap >= angular_spacing_per_ring[ring]
            half_junction_gap   = 0.5 * junction_angular_gap
            bridge_semiaxis     = 0.5 * (2 * track_major * sin(half_junction_gap) -
                                         spark_diameter)
            bridge_track_radius = track_major + nominal_spark_semiaxis - bridge_semiaxis

            xi = bridge_track_radius * cos(π - cap_tilt_angle - corotation_angle)
            yi = bridge_track_radius * track_minor / track_major *
                     sin(π - cap_tilt_angle - corotation_angle)
            bx, by = rotate_coordinates_2d(xi, yi, -cap_tilt_angle)
            push!(spark_x,        bx + cap_centre_x)
            push!(spark_y,        by + cap_centre_y)
            push!(spark_semiaxis, bridge_semiaxis)
        end

        # Gap filling at the θ = 0 / 2π junction (end of the ring).
        last_upper_angle = upper_track_angles[upper_offset + sparks_on_upper_track[ring]]
        last_lower_angle = 2π - lower_track_angles[lower_offset + sparks_on_lower_track[ring]]
        end_gap_total    = last_upper_angle + last_lower_angle

        if angular_spacing_per_ring[ring] <= end_gap_total <= 2 * angular_spacing_per_ring[ring]
            bridge_semiaxis     = 0.5 * (2 * track_major * sin(0.5 * end_gap_total) -
                                         spark_diameter)
            bridge_track_radius = track_major + nominal_spark_semiaxis - bridge_semiaxis

            xi = bridge_track_radius * cos(2π - cap_tilt_angle - corotation_angle)
            yi = bridge_track_radius * track_minor / track_major *
                     sin(2π - cap_tilt_angle - corotation_angle)
            bx, by = rotate_coordinates_2d(xi, yi, -cap_tilt_angle)
            push!(spark_x,        bx + cap_centre_x)
            push!(spark_y,        by + cap_centre_y)
            push!(spark_semiaxis, bridge_semiaxis)
        end
        =#

        # Step inward to the next ring.
        outer_major -= spark_diameter
        inner_major -= spark_diameter
        outer_minor -= ellipticity * spark_diameter
        inner_minor -= ellipticity * spark_diameter
    end

    # ------------------------------------------------------------------
    # Central (core) spark that fills the innermost region of the cap.
    # Its semi-axis is whatever space remains after all rings are laid down.
    # ------------------------------------------------------------------
    core_semiaxis = cap_major_semiaxis - number_of_rings * spark_diameter - 1.5 * drift_step
    push!(spark_x,        cap_centre_x)
    push!(spark_y,        cap_centre_y)
    push!(spark_semiaxis, core_semiaxis)

    number_of_sparks = length(spark_x)
    return spark_x, spark_y, spark_semiaxis, number_of_sparks
end


# =============================================================================
# Line-of-sight projection onto the polar cap
# (translated from pospolcap in sparkfunc.h)
# =============================================================================

"""
    project_los_to_polar_cap(
        los_track, phase_angle,
        cap_centre_theta, cap_centre_phi,
        max_boundary_phi, min_boundary_phi,
        number_of_los_points
    ) -> (x_on_cap, y_on_cap)

Given the current pulsar rotation phase angle, find the 2D coordinates of
the line-of-sight footprint on the flat polar-cap coordinate plane.

# How it works
1. `los_track` is a flat array storing a tabulated path of the line-of-sight
   across the open field-line boundary.  Each row has three values:
       phase angle [rad],  boundary colatitude theta [rad],  boundary azimuth phi [rad]
   Row `i` (1-based) occupies indices 3i-2, 3i-1, 3i.

2. For the given `phase_angle`, the two bracketing rows are found and used to
   interpolate the boundary angles (phi_cap, theta_cap) at that phase.

3. The 3D surface point at (theta_cap, phi_cap) is projected into the flat
   2D polar-cap frame by rotating it around the axis defined by the chord
   connecting the two extreme longitudes of the open field-line boundary.
   The rotation angle equals the cap centre colatitude, which tilts the
   curved surface into a local tangent plane.

# Arguments
- `los_track`           : flat array, row i → [phase, theta_boundary, phi_boundary]
- `phase_angle`         : current pulsar rotation phase [rad]
- `cap_centre_theta`    : colatitude of the polar cap centre [rad]
- `cap_centre_phi`      : longitude  of the polar cap centre [rad]
- `max_boundary_phi`    : maximum longitude of the open field-line boundary [rad]
- `min_boundary_phi`    : minimum longitude of the open field-line boundary [rad]
- `number_of_los_points`: number of rows in `los_track`

Returns `(x_on_cap, y_on_cap)` — the position in the 2D polar-cap frame.
"""
function project_los_to_polar_cap(
        los_track::Vector{Float64},
        phase_angle::Float64,
        cap_centre_theta::Float64,
        cap_centre_phi::Float64,
        max_boundary_phi::Float64,
        min_boundary_phi::Float64,
        number_of_los_points::Int)

    # Convenience accessors for the flat array layout (1-based indexing).
    phase_at(row) = los_track[3 * (row - 1) + 1]
    theta_at(row) = los_track[3 * (row - 1) + 2]
    phi_at(row)   = los_track[3 * (row - 1) + 3]

    # --- Step 1: interpolate the cap boundary position at this phase ---
    phi_at_boundary   = 0.0
    theta_at_boundary = 0.0

    for row in 1:(number_of_los_points - 1)
        if phase_angle >= phase_at(row) && phase_angle < phase_at(row + 1)
            # Interpolate the azimuthal cap angle across the bracket interval.
            phi_at_boundary = interpolate_between(
                phase_at(row), phase_at(row + 1), phase_angle,
                phi_at(row),   phi_at(row + 1)
            )

            # Recover the colatitude (theta) that corresponds to phi_at_boundary
            # by inverting the linear relationship between phi and theta.
            theta_at_boundary = find_edge_position(
                phi_at(row),   theta_at(row),
                phi_at(row + 1), theta_at(row + 1),
                phi_at_boundary
            )
        end
    end

    # --- Step 2: convert the 3D surface point to flat 2D cap coordinates ---
    # Notional neutron-star radius (same constant as in the original C code).
    neutron_star_radius = 10_000.0

    # Cartesian position of the polar cap centre on the stellar surface.
    cap_centre_cartesian = [
        neutron_star_radius * sin(cap_centre_theta) * cos(cap_centre_phi),
        neutron_star_radius * sin(cap_centre_theta) * sin(cap_centre_phi),
        neutron_star_radius * cos(cap_centre_theta),
    ]

    # Rotation axis: the chord connecting the two extreme azimuthal limits of
    # the open field-line boundary, at the cap centre colatitude.
    # This is tangent to the stellar surface and perpendicular to the
    # meridian plane of the cap centre.
    axis_raw = [
        sin(cap_centre_theta) * cos(max_boundary_phi) -
        sin(cap_centre_theta) * cos(min_boundary_phi),
        sin(cap_centre_theta) * sin(max_boundary_phi) -
        sin(cap_centre_theta) * sin(min_boundary_phi),
        0.0,
    ]
    rotation_axis = normalize_vector(axis_raw)

    # Vector from the cap centre to the interpolated boundary point.
    boundary_point_cartesian = [
        neutron_star_radius * sin(theta_at_boundary) * cos(phi_at_boundary) -
            cap_centre_cartesian[1],
        neutron_star_radius * sin(theta_at_boundary) * sin(phi_at_boundary) -
            cap_centre_cartesian[2],
        neutron_star_radius * cos(theta_at_boundary) -
            cap_centre_cartesian[3],
    ]

    # Rotate by -cap_centre_theta to flatten the curved cap surface into the
    # local tangent plane (the 2D polar-cap coordinate frame).
    flat_point = rotate_vector(
        boundary_point_cartesian, rotation_axis, -cap_centre_theta
    )

    x_on_cap = flat_point[1] + cap_centre_cartesian[1]
    y_on_cap = flat_point[2] + cap_centre_cartesian[2]

    return x_on_cap, y_on_cap
end


# =============================================================================
# Signal amplitude at a polar-cap point
# (translated from sparkamp in sparkfunc.h)
# =============================================================================

"""
    compute_spark_amplitude(
        spark_x, spark_y, spark_semiaxis,
        cap_major_semiaxis, cap_minor_semiaxis, cap_tilt_angle,
        query_x, query_y
    ) -> amplitude

Compute the total signal amplitude at the polar-cap point `(query_x, query_y)`
by summing Gaussian intensity profiles over all sparks.

# Spark emission profile
Each spark is an ellipse with the same axis ratio as the polar cap, tilted
by `cap_tilt_angle`.  Inside the ellipse the intensity follows a 2D Gaussian:

    intensity = exp( -x_rotated² / (2 σ_major²)
                     -y_rotated² / (2 σ_minor²) )

where the Gaussian widths are σ = semiaxis / 1.75 (the profile reaches 1/e
at about 57 % of the ellipse boundary).

The function tests sparks in order and exits as soon as the first spark
containing the query point is found (behaviour identical to the original C).

# Arguments
- `spark_x`, `spark_y`  : centre coordinates of each spark (display frame)
- `spark_semiaxis`      : major semi-axis of each spark
- `cap_major_semiaxis`  : major semi-axis of the polar cap (sets the ellipticity ratio)
- `cap_minor_semiaxis`  : minor semi-axis of the polar cap
- `cap_tilt_angle`      : tilt of the cap ellipse [rad]
- `query_x`, `query_y`  : point on the cap where the amplitude is evaluated

Returns the total Gaussian amplitude (0 if the point is outside all sparks).
"""
function compute_spark_amplitude(
        spark_x::Vector{Float64},
        spark_y::Vector{Float64},
        spark_semiaxis::Vector{Float64},
        cap_major_semiaxis::Float64,
        cap_minor_semiaxis::Float64,
        cap_tilt_angle::Float64,
        query_x::Float64,
        query_y::Float64)

    # Ratio of minor to major axis — same for every spark as for the cap.
    ellipticity = cap_minor_semiaxis / cap_major_semiaxis

    total_amplitude  = 0.0
    number_of_sparks = length(spark_x)

    for spark_index in 1:number_of_sparks

        # Translate query point into the coordinate frame of this spark.
        delta_x = query_x - spark_x[spark_index]
        delta_y = query_y - spark_y[spark_index]

        # Rotate into the cap-tilted (spark-aligned) frame so the ellipse
        # axes align with the coordinate axes.
        x_in_spark_frame, y_in_spark_frame =
            rotate_coordinates_2d(delta_x, delta_y, cap_tilt_angle)

        # Spark semi-axes (same ellipticity as the polar cap).
        spark_major_semiaxis = spark_semiaxis[spark_index]
        spark_minor_semiaxis = spark_major_semiaxis * ellipticity

        # Standard ellipse membership test: value < 1 means the point is inside.
        ellipse_value = sqrt(
            x_in_spark_frame^2 / spark_major_semiaxis^2 +
            y_in_spark_frame^2 / spark_minor_semiaxis^2
        )

        if ellipse_value < 1.0
            # Gaussian sigma chosen so the profile drops to 1/e at ~57 % of the edge.
            sigma_major = spark_major_semiaxis / 1.75
            sigma_minor = spark_minor_semiaxis / 1.75

            intensity = exp(
                -x_in_spark_frame^2 / (2.0 * sigma_major^2) -
                 y_in_spark_frame^2 / (2.0 * sigma_minor^2)
            )
            total_amplitude += intensity

            # Exit after the first containing spark (matches the C original).
            break
        end
    end

    return total_amplitude
end
