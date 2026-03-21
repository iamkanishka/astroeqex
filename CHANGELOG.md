# Changelog

All notable changes are documented here.
Format: [Keep a Changelog](https://keepachangelog.com/en/1.0.0/).
Versioning: [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

---

## [Unreleased]

---

## [0.2.0] — 2025-03-19

Major expansion of the library. The entire `AstroEquations.Physics.*` namespace is new. Astrophysics, mathematics, and statistics modules are substantially deepened. Total public function count grows from ~60 to **628**.

### Added — Physics (all new)

#### `AstroEquations.Physics.Electromagnetism` (78 functions)
Maxwell's equations (`gauss_law`, `gauss_law_differential`, `gauss_law_magnetism`, `faraday_law`, `faraday_emf`, `ampere_law`), Lorentz force (`lorentz_force_point`, `cyclotron_radius`, `cyclotron_frequency`), electric fields and potentials (`electric_field_point`, `electric_field_ring`, `dipole_field_parallel`, `dipole_field_perpendicular`, `dipole_moment`, `electric_potential_point`, `potential_difference_point`, `potential_energy`, `field_energy_density`, `field_energy`), charge and current densities (`surface_charge_density`, `linear_charge_density`, `current_density`, `drift_velocity`, `current_from_properties`, `continuity`), circuit theory (Ohm's law, power, resistance, series/parallel R/C/L, RC/RL/LC time constants, capacitor and inductor impedance), capacitors (`capacitance`, `capacitance_from_geometry`, `capacitor_energy`, `capacitor_field`), magnetic fields (`biot_savart`, `moving_charge_field`, `wire_magnetic_field`, `solenoid_field`, `toroid_field`, `parallel_wire_force`, `magnetic_flux`), inductors (`solenoid_inductance`, `inductor_emf`, `inductor_energy`, `magnetic_energy_density`, `mutual_inductance_emf`, `vector_potential`), dielectric and magnetic materials (`relative_permittivity`, `electric_susceptibility`, `absolute_permittivity`, `polarization`, `surface_bound_charge`, `volume_bound_charge`, `total_bound_charge`, `electric_displacement`, `magnetic_field_strength`, `magnetic_susceptibility`, `relative_permeability`, `magnetic_dipole_moment`, `bound_surface_current`, `gauss_law_polarization`, `gauss_law_displacement`), EM waves (`wave_speed`, `free_space_impedance`, `poynting_magnitude`, `radiation_pressure`, `skin_depth`, `plasma_frequency`, `brewsters_angle`).

#### `AstroEquations.Physics.Energy` (18 functions)
`work/2`, `work_angle/3`, `kinetic_energy/2`, `kinetic_energy_from_momentum/2`, `rotational_kinetic_energy/2`, `gravitational_pe/3`, `gravitational_pe_general/3`, `spring_pe/2`, `electric_pe/4`, `power_from_work/2`, `mechanical_power/3`, `power_vectors/2`, `work_energy_theorem/2`, `escape_velocity/2`, `rest_energy/2`, `binding_energy/2`, `shm_total_energy/2`, `shm_kinetic_energy/4`.

#### `AstroEquations.Physics.Forces` (22 functions)
`newtons_second_law/2`, `gravitational_force/3`, `weight/2`, `buoyancy/2`, `buoyancy_from_density/3`, `kinetic_friction/2`, `static_friction/2`, `rolling_friction/2`, `spring_force/2`, `centripetal_force/3`, `centripetal_acceleration/2`, `centripetal_force_angular/3`, `stokes_drag/2`, `stokes_drag_sphere/3`, `quadratic_drag/4`, `terminal_velocity/5`, `hydrostatic_pressure/3`, `pressure/2`, `normal_force_incline/3`, `incline_force_parallel/3`, `impulse/2`, `momentum/2`.

#### `AstroEquations.Physics.GeneralRelativity` (22 functions)
Metrics: `minkowski_metric/0`, `minkowski_interval/5`, `schwarzschild_metric/5`, `schwarzschild_interval/9`, `rindler_interval/5`, `flrw_interval/4`. Tensor algebra: `lower_index/4`, `raise_index/4`, `transform_tensor/3`, `four_vector_product/3`, `inverse_metric/1`. GR effects: `gravitational_time_dilation/2`, `gravitational_redshift/2`, `gravitational_wave_strain/2`, `orbital_precession/3`, `light_deflection/2`, `radial_freefall_velocity/2`. Cosmology: `hubble_from_density/1`, `critical_density/1`, `density_parameter/2`, `scale_factor_matter/3`, `lookback_time/2`.

#### `AstroEquations.Physics.Materials` (22 functions)
`density/2`, `continuous_density/3`, `number_density/2`, `normal_stress/2`, `normal_strain/2`, `youngs_modulus/2`, `extension/4`, `shear_stress/2`, `shear_modulus/2`, `bulk_modulus/3`, `elastic_energy_density/2`, `linear_expansion/3`, `volumetric_expansion/3`, `heat_conduction/4`, `thermal_resistance/3`, `viscous_shear_stress/2`, `kinematic_viscosity/2`, `reynolds_number/4`, `speed_of_sound_solid/2`, `speed_of_sound_gas/4`, `polytropic_pressure/3`.

#### `AstroEquations.Physics.Motion` (53 functions)
SUVAT kinematics (`velocity`, `acceleration`, `final_velocity`, `displacement_suvat`, `velocity_squared`, `displacement_average`), dynamics (`newtons_second_law`, `momentum`, `impulse`, `acceleration_from_force`, `kinetic_energy`, `rotational_kinetic_energy`, `total_kinetic_energy`, `centripetal_force`, `centripetal_acceleration`, `orbital_speed`), rotation (`angular_velocity`, `angular_acceleration`, `tangential_velocity`, `tangential_acceleration`, `angular_momentum`, `torque`, `rotational_newtons_law`, `angular_impulse`), moments of inertia (7 shapes + `parallel_axis`), projectile motion (`horizontal_displacement`, `vertical_displacement`, `vertical_velocity`, `time_of_flight`, `projectile_range`, `projectile_max_height`), SHM (`shm_displacement`, `shm_velocity`, `shm_acceleration`, `spring_angular_frequency`, `pendulum_angular_frequency`, `period_from_omega`), analytical mechanics (`lagrangian`, `generalized_momentum`, `hamiltonian`, `hamiltons_equations`, `poisson_bracket`), orbital (`keplers_third_law`, `vis_viva`).

#### `AstroEquations.Physics.NewtonGravity` (21 functions)
`force/3`, `field/2`, `potential/2`, `potential_energy/3`, `approximate_potential_energy/3`, `orbital_period/3`, `circular_orbit_speed/2`, `escape_velocity/2`, `vis_viva/3`, `specific_orbital_energy/2`, `specific_angular_momentum/2`, `semi_major_axis_from_period/2`, `surface_gravity/2`, `tidal_acceleration/3`, `roche_limit/3`, `hill_sphere/3`.

#### `AstroEquations.Physics.Oscillations` (24 functions)
`force/2`, `potential_energy/2`, `angular_frequency/2`, `pendulum_angular_frequency/2`, `physical_pendulum_frequency/4`, `period_from_omega/1`, `spring_period/2`, `pendulum_period/2`, `frequency_from_period/1`, `omega_from_frequency/1`, `shm_displacement/4`, `shm_velocity/4`, `shm_acceleration/4`, `shm_velocity_from_position/3`, `shm_total_energy/2`, `shm_kinetic_energy/4`, `damping_ratio/3`, `damped_frequency/2`, `underdamped_displacement/6`, `quality_factor/3`, `resonance_amplitude/5`, `resonant_frequency/3`, `lc_angular_frequency/2`, `beat_frequency/2`.

#### `AstroEquations.Physics.QuantumMechanics` (35 functions)
Fundamentals (`uncertainty_principle?`, `min_uncertainty_product`, `energy_time_uncertainty?`, `de_broglie_wavelength`, `de_broglie_wavelength_ke`, `photon_energy`, `photon_momentum`), energy levels (`hydrogen_energy_level`, `hydrogen_energy_joules`, `bohr_radius`, `rydberg_wavelength`, `infinite_well_energy`, `harmonic_oscillator_energy`), probability (`born_rule`, `expectation_braket`, `variance`, `standard_deviation`, `expectation_position`), operators (`pauli_x`, `pauli_z`, `identity`, `atomic_raise`, `atomic_lower`, `atomic_sigma_z`, `apply_raise`, `apply_lower`, `annihilate`, `create`), state (`density_matrix`, `purity`, `matrix_trace`, `hilbert_schmidt_norm`), scattering (`transmission_coefficient`, `reflection_coefficient`, `wave_vector`).

#### `AstroEquations.Physics.SpecialRelativity` (25 functions)
`gamma_factor/2`, `beta/2`, `rapidity/2`, `speed_from_rapidity/2`, `time_dilation/3`, `length_contraction/3`, `relative_velocity/3`, `relativistic_mass/3`, `relativistic_momentum/3`, `rest_energy/2`, `total_energy/3`, `kinetic_energy/3`, `energy_momentum_relation/3`, `speed_from_kinetic_energy/3`, `relativistic_doppler/3`, `transverse_doppler/3`, `four_vector/4`, `four_velocity/4`, `four_momentum/5`, `four_product/2`, `invariant_mass/3`, `galilean_transform/5`, `lorentz_boost/6`, `proper_time/5`, `spacetime_interval/5`.

#### `AstroEquations.Physics.Thermodynamics` (39 functions)
Ideal gas: `ideal_gas_law/1` (keyword-list polymorphic — solves for any nil variable), `molar_ideal_gas_law/5`, `mean_kinetic_energy/2`, `rms_speed/3`, `most_probable_speed/3`, `mean_speed/3`. Heat and work: `heat_energy/3`, `heat_capacity/2`, `specific_heat_capacity/2`, `first_law/2`, `isobaric_work/2`, `isothermal_work/5`, `adiabatic_work/5`, `adiabatic_pressure/4`, `adiabatic_temperature/4`. Entropy: `entropy/2`, `entropy_change/2`, `microstates/2`. Heat engines: `carnot_efficiency/2`, `cop_refrigerator/2`, `cop_heat_pump/2`. Heat capacity: `mayers_relation/2`, `heat_capacity_ratio/2`, `cv_monatomic/0`, `cv_diatomic/0`. Radiation: `photon_energy/1`, `wiens_displacement/2`, `stefan_boltzmann/2`, `stefan_boltzmann_total/3`, `planck_wavelength/5`, `planck_frequency/5`, `rayleigh_jeans/4`. Kinetics: `maxwell_boltzmann/4`, `newton_cooling/4`.

#### `AstroEquations.Physics.Waves` (34 functions)
`wave_number/1`, `wave_velocity/2`, `frequency/2`, `wavelength/2`, `period/1`, `angular_frequency/1`, `omega_from_frequency/1`, `wave_function/6`, `standing_wave/5`, `constructive_interference?/3`, `destructive_interference?/3`, `beat_frequency/2`, `single_slit_minimum/3`, `double_slit_bright/3`, `grating_angle/3`, `grating_resolving_power/2`, `doppler/4`, `relativistic_doppler/3`, `intensity_decibels/2`, `intensity_from_decibels/2`, `intensity_ratio_from_amplitude/2`, `string_harmonic/3`, `open_pipe_harmonic/3`, `closed_pipe_harmonic/3`, `phase_velocity/2`, `group_velocity/4`, `is_dispersionless?/2`, `wave_power/4`, `wave_intensity/2`, `spherical_wave_intensity/2`, `snells_law/3`, `critical_angle/2`, `refractive_index/2`, `malus_law/2`.

---

### Added — Astrophysics (expanded)

#### `AstroEquations.AstrophysicsAndAstronomy.Astrometry`
New: `relativistic_velocity/2`, `absolute_magnitude/2`, `absolute_magnitude_extinction/3`, `distance_from_modulus/1`, `flux_ratio_from_magnitudes/2`, `flux_from_magnitude/2`, `bolometric_magnitude/2`, `luminosity_from_bolometric/1`, `parallax_from_distance/1`, `proper_motion_total/2`, `transverse_velocity/2`, `space_velocity/2`, `angular_diameter/2`, `physical_size/2`, `hubble_distance/2`, `feh_to_z/1`.

#### `AstroEquations.AstrophysicsAndAstronomy.BlackHole`
New: `schwarzschild_radius_solar/1`, `evaporation_time/1`, `photon_sphere_radius/1`, `isco_radius/1`, `bekenstein_hawking_entropy/1`, `kerr_spin_parameter/2`, `ergosphere_radius/2`.

#### `AstroEquations.AstrophysicsAndAstronomy.Galaxies`
New: `de_vaucouleurs_profile/3`, `disk_density/5`, `nfw_profile/3`, `keplerian_rotation/2`, `flat_rotation_curve/1`, `tully_fisher/3`, `faber_jackson/2`, `specific_sfr/2`, `sfr_from_uv/1`, `sfr_from_halpha/1`, `schechter_function/4`, `virial_mass/2`.

#### `AstroEquations.AstrophysicsAndAstronomy.Instrumentation`
New: `image_scale/2`, `seeing_limit/2`, `total_resolution_limit/3`, `nyquist_sampling/4`, `plate_scale/1`, `fitting_error/2`, `adaptive_optics_error/9`, `strehl_ratio/1`, `signal_to_noise/6`, `limiting_magnitude/6`, `dynamic_range_db/2`, `resolving_power/2`, `grating_dispersion/3`, `blaze_wavelength/2`, `airmass/1`, `specific_impulse/2`.

#### `AstroEquations.AstrophysicsAndAstronomy.Stars`
New: `kelvin_helmholtz_timescale/3`, `nuclear_timescale/1`, `dynamical_timescale/2`, `gravitational_potential_energy/2`, `eddington_luminosity/1`, `eddington_mass/1`, `eddington_mass_loss_rate/1`, `radius_from_luminosity_temperature/2`, `planck_function/2`, `jeans_mass/3`, `jeans_radius/3`.

---

### Added — Mathematics (expanded)

#### `AstroEquations.Mathematics.Calculus`
New: full numerical calculus suite — `derivative/3`, `forward_difference/3`, `backward_difference/3`, `second_derivative/3`, `trapezoid/4`, `simpsons/4`, `rectangle/4`, `bisection/5`, `newton_raphson/4`, `euler_ode/5`, `rk4_ode/5`.

#### `AstroEquations.Mathematics.Geometry`
New: `angular_separation/4`, `angular_separation_deg/4`, `position_angle/4`, `solid_angle_cone/1`, `solid_angle_rectangle/3`, `sphere_solid_angle/0`, `semi_latus_rectum/2`, `orbital_radius/3`, `periapsis/2`, `apoapsis/2`, `eccentric_anomaly/3`, `true_anomaly/2`, `einstein_radius/4`, `microlensing_magnification/1`, `equatorial_to_cartesian/2`, `cartesian_to_equatorial/3`.

#### `AstroEquations.Mathematics.Notation`
New: 24 physical constants (`speed_of_light`, `gravitational_constant`, `planck_constant`, `hbar`, `boltzmann_constant`, `stefan_boltzmann`, `radiation_constant`, `proton_mass`, `electron_mass`, `elementary_charge`, `thomson_cross_section`, `wien_constant`, `avogadro`, `solar_mass`, `solar_radius`, `solar_luminosity`, `solar_temperature`, `astronomical_unit`, `parsec`, `light_year`, `earth_mass`, `earth_radius`, `jupiter_mass`, `hubble_constant`). Unit conversions for length, angle, energy, time, temperature. Sexagesimal (`hms_to_deg`, `dms_to_deg`, `deg_to_hms`, `deg_to_dms`). Julian Date (`calendar_to_jd`, `jd_to_mjd`, `mjd_to_jd`, `j2000`, `julian_centuries`).

#### `AstroEquations.Mathematics.Trigonometry`
New: `spherical_law_of_cosines/3`, `spherical_law_of_sines/3`, `altitude/3`, `azimuth/3`, `hour_angle/2`, `parallactic_angle/3`, `gmst/1`, `local_sidereal_time/2`, `sunrise_hour_angle/3`, `atmospheric_refraction/1`, `equatorial_to_ecliptic/3`, `equatorial_to_galactic/2`, `galactic_to_equatorial/2`.

---

### Added — Statistics (all new)

#### `AstroEquations.Statistics.Variance` (26 functions)
`population/1`, `sample/1`, `weighted/2`, `standard_error/1`, `welford_accumulate/1`, `welford_finalize/1`, `covariance/2`, `pearson_r/2`, `chi_squared/2`, `reduced_chi_squared/2`, `mad/1`, `iqr/1`.

#### `AstroEquations.Statistics.StandardDeviation` (16 functions)
`population/1`, `sample/1`, `weighted/2`, `standard_error/1`, `z_score/3`, `normalise/1`, `propagate_addition/2`, `propagate_product/4`, `propagate_power/3`, `propagate_log/2`, `photometric_snr/4`, `magnitude_uncertainty/1`, `limiting_flux/3`, `poisson_uncertainty/1`, `online/1`.

---

### Added — Package infrastructure

- `lib/astro_equations.ex` — package entry module with `version/0` and module index
- `test/test_helper.exs` — ExUnit bootstrap
- 15 test files covering all 24 modules, **675 tests** total
- `.credo.exs` — Credo strict configuration with `max_arity: 9` for scientific functions
- `mix.exs` — `dialyxir` and `excoveralls` dev/test dependencies

---

### Fixed

- `Electromagnetism.electric_susceptibility/1` — was returning `1 - εᵣ` (inverted sign); corrected to `εᵣ - 1`
- `Electromagnetism.bound_volume_current/2` — attempted struct field access on plain lists; rewritten with `Enum.at/2`
- `GeneralRelativity` — merged duplicate top-level modules into the single `AstroEquations.Physics.GeneralRelativity` namespace; `raise_index/4` now uses `inverse_metric/1` correctly
- `QuantumMechanics` — removed undeclared `Complex` library dependency; all matrix operations rewritten in pure Elixir; fixed `purity/1` arity; `matrix_multiply/2` now handles both matrix×vector and matrix×matrix; removed broken `diagonalize/1` stub
- `QuantumMechanics.matrix_multiply/2` — guard `not is_list(hd(b))` raised on empty list; fixed to `b != [] and not is_list(hd(b))`
- `Astrometry.apparent_magnitude_diff/4` — was ignoring its magnitude arguments; simplified to `apparent_magnitude_diff/2` using only the flux ratio
- `Astrometry.color_index/2` — was using `:math.log` (natural log) instead of `:math.log10`
- `Astrometry.metallicity/2` — domain error on negative input; corrected to accept absolute Fe/H ratios
- `Instrumentation.nyquist_sampling/4` — formula ignored `fsys` and `dt` entirely; replaced with the correct pixels-per-resolution-element formula
- `Instrumentation.atmospheric_extinction/3` — sign was inverted; corrected to `m + A·sec(z)`
- `Instrumentation` — module name typo `Instrumentaion` → `Instrumentation`
- `Trigonometry` — module name typo `Trignometry` → `Trigonometry`

### Fixed (documentation and type system)

- **`@spec` for 628 functions** — all public functions now carry complete type annotations
- **`@doc` for 628 functions** — all public functions now carry documentation with parameter descriptions
- **`## Examples` for 321 functions** — inline `iex>` doctestable examples (54% coverage)
- **8 `@spec` return type corrections** — functions that call `raise/1` on domain violations now declare `:: float | no_return()` rather than `:: float`
- **`welford_finalize/1` return type** — standardised to `{:ok, mean, pop_var, samp_var} | {:error, String.t()}` for consistency with `StandardDeviation.online/1`

### Changed

- `SpecialRelativity` — `four_vector/4` now returns a plain map `%{ct:, x:, y:, z:}` rather than a struct, consistent with the rest of the library
- `Thermodynamics.ideal_gas_law/1` — accepts a keyword list `[p:, v:, n:, k_b:, t:]` and solves for whichever key is `nil`, allowing the same function to solve all five variable forms

---

## [0.1.0] — 2024-01-01

Initial release.

### Added

- `AstroEquations.AstrophysicsAndAstronomy.Astrometry` — redshift, magnitude, flux, colour index, metallicity (~8 functions)
- `AstroEquations.AstrophysicsAndAstronomy.BlackHole` — Schwarzschild radius, Hawking temperature (~3 functions)
- `AstroEquations.AstrophysicsAndAstronomy.Galaxies` — Hubble classification, Sérsic profile, disk density (~5 functions)
- `AstroEquations.AstrophysicsAndAstronomy.Instrumentation` (as `Instrumentaion`) — lensmaker's equation, focal ratio, FOV, resolution limits, AO error, SNR, atmospheric extinction, rocket equation (~12 functions)
- `AstroEquations.AstrophysicsAndAstronomy.Stars` — hydrostatic equilibrium, mass conservation, energy equation, radiative transport, timescales, Eddington limits, mass-luminosity (~11 functions)
- `AstroEquations.Mathematics.Notation` — a small set of physical constants and unit conversions (~24 functions)
- `AstroEquations.Statistics.Variance` and `AstroEquations.Statistics.StandardDeviation` — stub implementations (~20 functions)

---

[Unreleased]: https://github.com/your-org/astroequations/compare/v0.2.0...HEAD
[0.2.0]: https://github.com/your-org/astroequations/compare/v0.1.0...v0.2.0
[0.1.0]: https://github.com/your-org/astroequations/releases/tag/v0.1.0
