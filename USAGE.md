# AstroEquations — Usage Guide

This guide walks through practical usage of every module with worked examples, pipe patterns, and notes on units and edge cases.

---

## Table of Contents

1. [Installation & Setup](#installation--setup)
2. [Astrophysics & Astronomy](#astrophysics--astronomy)
   - [Astrometry](#astrometry)
   - [Black Holes](#black-holes)
   - [Galaxies](#galaxies)
   - [Instrumentation](#instrumentation)
   - [Stars](#stars)
3. [Physics](#physics)
   - [Electromagnetism](#electromagnetism)
   - [Energy](#energy)
   - [Forces](#forces)
   - [General Relativity](#general-relativity)
   - [Materials](#materials)
   - [Motion](#motion)
   - [Newton Gravity](#newton-gravity)
   - [Oscillations](#oscillations)
   - [Quantum Mechanics](#quantum-mechanics)
   - [Special Relativity](#special-relativity)
   - [Thermodynamics](#thermodynamics)
   - [Waves](#waves)
4. [Mathematics](#mathematics)
   - [Calculus](#calculus)
   - [Geometry](#geometry)
   - [Notation](#notation)
   - [Trigonometry](#trigonometry)
5. [Statistics](#statistics)
   - [Variance](#variance)
   - [Standard Deviation](#standard-deviation)
6. [Patterns & Tips](#patterns--tips)

---

## Installation & Setup

```elixir
# mix.exs
def deps do
  [{:astroequations, "~> 0.2"}]
end
```

All modules are under the `AstroEquations` namespace. Use aliases for brevity:

```elixir
alias AstroEquations.AstrophysicsAndAstronomy.{Astrometry, BlackHole, Galaxies, Instrumentation, Stars}
alias AstroEquations.Physics.{Electromagnetism, Energy, Forces, GeneralRelativity,
                               Materials, Motion, NewtonGravity, Oscillations,
                               QuantumMechanics, SpecialRelativity, Thermodynamics, Waves}
alias AstroEquations.Mathematics.{Calculus, Geometry, Notation, Trigonometry}
alias AstroEquations.Statistics.{Variance, StandardDeviation}
```

---

## Astrophysics & Astronomy

### Astrometry

**Units:** wavelengths in metres, distances in parsecs (unless noted), velocities in km/s, flux in W/m², angles in arcseconds.

```elixir
alias AstroEquations.AstrophysicsAndAstronomy.Astrometry

# Redshift: shifted Hα from 656 nm to 700 nm
Astrometry.redshift(700, 656)
#=> 0.0671

# Recession velocity (non-relativistic)
Astrometry.recession_velocity(0.0671)
#=> 20_123.3  # km/s

# For high-z objects, use the relativistic form
Astrometry.relativistic_velocity(0.5)
#=> 120_067.5  # km/s

# Magnitude system: flux ratio 100 = 5 magnitudes
Astrometry.apparent_magnitude_diff(1.0, 100.0)
#=> -5.0

Astrometry.flux_ratio_from_magnitudes(10, 15)
#=> 100.0

# Distance modulus — for Andromeda (~770 kpc)
Astrometry.distance_modulus(770_000)
#=> 24.43

Astrometry.distance_from_modulus(24.43)
#=> 769_974.0  # parsecs

# Parallax — Proxima Centauri at 0.7687 arcsec
Astrometry.distance_from_parallax(0.7687)
#=> 1.3009  # pc

# Hubble distance from redshift
Astrometry.hubble_distance(0.1)
#=> 428.3  # Mpc

# Metallicity [Fe/H]
Astrometry.metallicity(0.017, 0.017)   #=> 0.0  (solar)
Astrometry.metallicity(0.0085, 0.017)  #=> -0.301  (half-solar)
Astrometry.feh_to_z(0.0)              #=> 0.017  (solar Z)
```

---

### Black Holes

**Units:** mass in kg, radii in metres, temperature in Kelvin, time in seconds.

```elixir
alias AstroEquations.AstrophysicsAndAstronomy.BlackHole

m_sun = 1.989e30  # kg

# Schwarzschild radius — the event horizon
BlackHole.schwarzschild_radius(m_sun)
#=> 2953.25  # m (~3 km for a solar-mass BH)

# Solar-mass convenience wrapper
BlackHole.schwarzschild_radius_solar(10)
#=> 29_532.5  # m for a 10 M☉ BH

# ISCO and photon sphere (multiples of r_s)
r_s = BlackHole.schwarzschild_radius(m_sun)
BlackHole.isco_radius(m_sun) / r_s     #=> 3.0  (innermost stable orbit = 3 r_s)
BlackHole.photon_sphere_radius(m_sun) / r_s  #=> 1.5

# Hawking temperature — stellar-mass BHs are extremely cold
BlackHole.hawking_temperature(m_sun)
#=> 6.17e-8  # K

# Kerr spin parameter (0 = non-rotating, 1 = maximally rotating)
BlackHole.kerr_spin_parameter(m_sun, 8.0e47)
#=> 0.427

# Gravitational time dilation at the solar surface radius
BlackHole.gravitational_time_dilation(m_sun, 6.957e8)
#=> 0.99999788  # clocks run very slightly slow
```

---

### Galaxies

```elixir
alias AstroEquations.AstrophysicsAndAstronomy.Galaxies

# Hubble morphological classification from axis ratio
Galaxies.hubble_classify(10, 7)   #=> "E3"
Galaxies.hubble_classify(10, 10)  #=> "E0"  (circular)

# Sérsic surface brightness profile
# At r = r_e: I(r_e) = I₀ for the standard Sérsic normalisation
Galaxies.sersic_profile(100, 1.0, 1.0, 4)
#=> 100.0

# NFW dark matter density profile
Galaxies.nfw_profile(1.0e7, 5.0, 3.0)
#=> density at r=5 kpc for r_s=3 kpc

# Rotation curve — Keplerian vs flat
r = 1.0e20   # ~3.2 kpc
Galaxies.keplerian_rotation(1.0e41, r)   #=> orbital speed (Keplerian)
Galaxies.flat_rotation_curve(220.0)      #=> 220.0  (asymptotic v_rot in km/s)

# Tully-Fisher: luminosity from rotation velocity
Galaxies.tully_fisher(220)   #=> L in solar luminosities

# Star formation rate from Hα luminosity
Galaxies.sfr_from_halpha(1.26e41)  #=> ~1.0  M☉/yr

# Virial mass from velocity dispersion and radius
Galaxies.virial_mass(200_000, 1.0e22)  #=> mass in kg
```

---

### Instrumentation

```elixir
alias AstroEquations.AstrophysicsAndAstronomy.Instrumentation

# Basic telescope optics
Instrumentation.focal_ratio(2.0, 0.4)         #=> 5.0   (f/5)
Instrumentation.field_of_view(0.036, 0.2, 5)  #=> FOV in radians

# Resolution limits
Instrumentation.diffraction_limit(500.0e-9, 0.4)   #=> 1.525e-6 rad (Rayleigh)
Instrumentation.seeing_limit(500.0e-9, 0.15)        #=> atmospheric limit

# Combined limit (quadrature)
Instrumentation.total_resolution_limit(500.0e-9, 8.2, 0.15)
#=> total seeing + diffraction in quadrature

# Adaptive optics performance
Instrumentation.strehl_ratio(0.1)    #=> 0.905  (10% Marechal wavefront RMS)
Instrumentation.fitting_error(0.1, 0.15)  #=> fitting error variance

# Signal-to-noise and detection
Instrumentation.signal_to_noise(1000, 50, 3, 900, 0.02, 2.0)
#=> CCD SNR

Instrumentation.limiting_magnitude(10, 3, 100, 0.02, 2.0, 5.0)
#=> sky-limited detection magnitude

# Spectroscopy
Instrumentation.resolving_power(550.0e-9, 0.05e-9)  #=> 11_000
Instrumentation.airmass(:math.pi() / 3)              #=> 2.0  (at 60° zenith angle)

# Atmospheric extinction
Instrumentation.atmospheric_extinction(5.0, 0.15, 1.5)
#=> 5.225  (magnitude after extinction at 1.5 airmasses)

# Rocketry
Instrumentation.tsiolkovsky_rocket_equation(4400, 100_000, 20_000)
#=> Δv in m/s
```

---

### Stars

```elixir
alias AstroEquations.AstrophysicsAndAstronomy.Stars

# Stellar structure differential equations
Stars.hydrostatic_equilibrium(1.989e30, 1408.0, 6.957e8)
#=> dP/dr  (always negative — inward)

Stars.mass_conservation(1408.0, 6.957e8)
#=> dM/dr  (always positive)

Stars.energy_equation(1408.0, 1.934e-7, 6.957e8)
#=> dL/dr

Stars.radiative_transport(0.2, 1408, 5778, 3.828e26, 6.957e8)
#=> dT/dr  (negative)

# Timescales (results in years or seconds depending on function)
Stars.kelvin_helmholtz_timescale(1.0, 1.0, 1.0)  #=> ~30 Myr
Stars.nuclear_timescale(1.0)                       #=> ~1 Gyr
Stars.dynamical_timescale(1.0, 1.0)               #=> seconds

# Mass-luminosity relation
Stars.mass_luminosity(1.0)    #=> 1.0   L☉ (by definition)
Stars.mass_luminosity(5.0)    #=> ~590  L☉
Stars.mass_luminosity(30.0)   #=> ~32_000 L☉ (Eddington limited)

# Eddington limit
Stars.eddington_luminosity(1.0)   #=> ~32_000 L☉

# Emission spectrum
Stars.wien_peak_wavelength(5778)   #=> 5.015e-7 m  (green light — Sun)
Stars.wien_peak_wavelength(30_000) #=> 9.66e-8 m   (UV — hot OB star)

# Chandrasekhar mass limit
Stars.chandrasekhar_mass()       #=> ~1.46 M☉ (default: C/O WD, μ_e=2)
Stars.chandrasekhar_mass(4.0)    #=> lower limit for He WD

# Jeans collapse
Stars.jeans_mass(10, 1.0e-17)    #=> Jeans mass in kg
Stars.jeans_radius(10, 1.0e-17)  #=> Jeans radius in m
```

---

## Physics

### Electromagnetism

```elixir
alias AstroEquations.Physics.Electromagnetism, as: EM

# Maxwell's equations
EM.gauss_law(1.0e-9)            #=> electric flux for 1 nC
EM.gauss_law_magnetism()        #=> 0  (no monopoles, always)
EM.faraday_emf(0.05)            #=> -0.05 V  (ε = -dΦ/dt)

# Point charge fields
EM.electric_field_point(1.0e-9, 0.1)         #=> E at 0.1 m from 1 nC
EM.electric_potential_point(1.0e-9, 0.1)     #=> V at 0.1 m

# Lorentz force: F = q(E + v×B)
EM.lorentz_force_point(1.602e-19, {1.0e4, 0, 0}, {0, 0, 0}, {0, 0, 1.0})
#=> {Fx, Fy, Fz}

# Cyclotron motion
EM.cyclotron_radius(9.109e-31, 1.0e6, 1.602e-19, 0.01)  #=> radius in m
EM.cyclotron_frequency(1.602e-19, 0.01, 9.109e-31)       #=> omega in rad/s

# Circuit theory
EM.ohms_law(2.0, 5.0)                          #=> 10.0 V
EM.series_resistance([100, 220, 330])           #=> 650 Ω
EM.parallel_resistance([100, 100])              #=> 50.0 Ω
EM.rc_time_constant(10_000, 100.0e-6)          #=> 1.0 s
EM.lc_resonant_frequency(1.0e-3, 10.0e-6)     #=> 10_000 rad/s

# Magnetic fields
EM.wire_magnetic_field(1.0, 0.1)              #=> 2.0e-6 T at 0.1 m
EM.solenoid_field(1000, 2.0)                  #=> B inside solenoid (n=1000/m, I=2A)

# EM waves
EM.wave_speed(1.257e-6, 8.854e-12)            #=> ~3.0e8 m/s (speed of light)
EM.free_space_impedance()                      #=> 376.73 Ω
EM.skin_depth(1.257e-6, 5.8e7, 6.283e6)      #=> ~66 μm (copper, 1 MHz)
EM.brewsters_angle(1.0, 1.5)                  #=> Brewster's angle for glass
```

---

### Energy

```elixir
alias AstroEquations.Physics.Energy

# Work
Energy.work([10, 0, 0], [5, 0, 0])    #=> 50.0 J  (dot product)
Energy.work_angle(10, 5, :math.pi()/3) #=> 25.0 J  (at 60°)

# Kinetic energy
Energy.kinetic_energy(2.0, 10.0)            #=> 100.0 J
Energy.kinetic_energy_from_momentum(20.0, 2.0) #=> 100.0 J
Energy.rotational_kinetic_energy(0.5, 20.0)    #=> 100.0 J

# Potential energy
Energy.gravitational_pe(10.0, 100.0)           #=> 9810.0 J (mgh)
Energy.gravitational_pe_general(5.972e24, 1000, 6.371e6)  #=> U = -Gm₁m₂/r
Energy.spring_pe(200.0, 0.1)                  #=> 1.0 J

# Power
Energy.power_from_work(1000.0, 5.0)   #=> 200.0 W
Energy.mechanical_power(50.0, 10.0)   #=> 500.0 W

# Special
Energy.escape_velocity(5.972e24, 6.371e6)  #=> 11_186 m/s
Energy.rest_energy(1.0)                     #=> 8.988e16 J  (E = mc²)
```

---

### Forces

```elixir
alias AstroEquations.Physics.Forces

Forces.newtons_second_law(70.0, 9.81)    #=> 686.7 N
Forces.gravitational_force(5.972e24, 7.348e22, 3.844e8)  #=> ~1.98e20 N
Forces.weight(70.0)                      #=> 686.7 N

# Friction
Forces.kinetic_friction(0.3, 100.0)    #=> 30.0 N
Forces.static_friction(0.5, 100.0)     #=> 50.0 N

# Drag
Forces.stokes_drag_sphere(1.81e-5, 0.001, 0.01)  #=> viscous drag on sphere
Forces.quadratic_drag(0.47, 1.225, 0.6, 30.0)    #=> aerodynamic drag
Forces.terminal_velocity(75.0, 0.47, 1.225, 0.6) #=> ~55 m/s for a person

# Centripetal
Forces.centripetal_force(1000.0, 30.0, 200.0)   #=> 4500 N
Forces.centripetal_acceleration(30.0, 200.0)     #=> 4.5 m/s²

# Inclines
Forces.normal_force_incline(10.0, :math.pi()/6)    #=> 85.0 N (30° slope)
Forces.incline_force_parallel(10.0, :math.pi()/6)  #=> 49.05 N along slope
```

---

### General Relativity

```elixir
alias AstroEquations.Physics.GeneralRelativity, as: GR

# Metric tensors
GR.minkowski_metric()
#=> [[-1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1]]

GR.schwarzschild_metric(1.989e30, 1.0e9, :math.pi()/2)
#=> 4×4 metric tensor at r=1000 km from solar-mass BH

# Spacetime intervals
GR.minkowski_interval(1, 0, 0, 0)        #=> -1  (timelike, using c=1)
GR.minkowski_interval(0, 1, 0, 0)        #=> +1  (spacelike)

# Tensor index gymnastics
m = GR.minkowski_metric()
GR.four_vector_product([1,0,0,0], [1,0,0,0], m)   #=> -1  (η_μν t^μ t^ν)

# Physical GR effects
GR.gravitational_time_dilation(1.989e30, 6.957e8)
#=> 0.999997879  (clocks run slow at solar surface)

GR.gravitational_redshift(1.989e30, 6.957e8)
#=> 2.12e-6  (z at solar surface)

# Solar limb light deflection: GR predicts 1.75 arcsec
alpha = GR.light_deflection(1.989e30, 6.957e8)
alpha * 180 / :math.pi() * 3600
#=> 1.75  arcsec

# Mercury perihelion precession: ~43 arcsec/century
delta_phi = GR.orbital_precession(1.989e30, 5.79e10, 0.206)
# Convert from rad/orbit to arcsec/century:
period_s = 88.0 * 24 * 3600
century_s = 100 * 365.25 * 24 * 3600
(delta_phi * 180 / :math.pi() * 3600) * (century_s / period_s)
#=> 42.98  arcsec/century

# Friedmann cosmology
h0_si = 67.36 * 1000 / 3.0856e22   # SI: s⁻¹
GR.critical_density(h0_si)          #=> 9.47e-27 kg/m³
GR.density_parameter(9.47e-27, 9.47e-27)  #=> 1.0  (flat universe)
GR.lookback_time(0.5, h0_si)        #=> lookback time in seconds
```

---

### Materials

```elixir
alias AstroEquations.Physics.Materials

# Density
Materials.density(7900.0, 1.0)        #=> 7900.0 kg/m³ (steel)

# Elasticity
Materials.normal_stress(1000.0, 0.01)    #=> 100_000 Pa
Materials.normal_strain(0.001, 1.0)      #=> 0.001
Materials.youngs_modulus(200.0e6, 0.001) #=> 200.0e9 Pa (steel)
Materials.extension(1000, 1.0, 200.0e9, 1.0e-4)  #=> 5.0e-5 m

# Thermal
Materials.linear_expansion(12.0e-6, 1.0, 100)   #=> 1.2e-3 m expansion
Materials.heat_conduction(401, 0.01, 100, 0.01)  #=> 40_100 W (copper rod)

# Fluid mechanics
Materials.reynolds_number(1.225, 30, 0.1, 1.81e-5)  #=> ~203_000 (turbulent!)
Materials.speed_of_sound_gas(1.4, 293.0, 0.029)     #=> ~343 m/s (air at 20°C)
Materials.speed_of_sound_solid(200.0e9, 7900.0)      #=> ~5032 m/s (steel)
```

---

### Motion

```elixir
alias AstroEquations.Physics.Motion

# SUVAT kinematics
Motion.final_velocity(0, 9.81, 2)        #=> 19.62 m/s (free fall, 2 s)
Motion.displacement_suvat(0, 3, 9.81)    #=> 44.145 m
Motion.velocity_squared(0, 9.81, 5)      #=> 9.905 m/s (from v² = 2as)

# Rotation
Motion.angular_momentum(5.0, 4.0)            #=> 20.0 kg·m²/s
Motion.torque(100.0, 0.5)                    #=> 50.0 N·m  (perpendicular)
Motion.thin_disk_moment_of_inertia(2.0, 0.5) #=> 0.25 kg·m²

# Projectile motion
Motion.projectile_range(20.0, :math.pi()/4)     #=> 40.775 m (optimal 45°)
Motion.time_of_flight(20.0, :math.pi()/4)        #=> 2.887 s
Motion.projectile_max_height(20.0, :math.pi()/4) #=> 10.194 m

# SHM
Motion.shm_displacement(0.1, 10.0, 0.0)   #=> 0.1 m (at t=0, x = A)
Motion.shm_velocity(0.1, 10.0, 0.0)       #=> 0.0 m/s (at max displacement)
Motion.spring_angular_frequency(100.0, 0.5) #=> ~14.14 rad/s

# Analytical mechanics
Motion.lagrangian(50.0, 20.0)   #=> 30.0 J (L = T - V)
Motion.keplers_third_law(1.496e11, 1.989e30)  #=> ~3.156e7 s ≈ 1 year

# Moments of inertia
Motion.point_mass_moment_of_inertia(1.0, 2.0)        #=> 4.0 kg·m²
Motion.solid_sphere_moment_of_inertia(1.0, 0.5)      #=> 0.1 kg·m²
Motion.parallel_axis(0.1, 1.0, 0.5)                  #=> 0.35 kg·m²
```

---

### Newton Gravity

```elixir
alias AstroEquations.Physics.NewtonGravity

m_earth = 5.972e24
r_earth = 6.371e6

# Basic gravitational functions
NewtonGravity.force(m_earth, 70.0, r_earth)          #=> 686.5 N
NewtonGravity.field(m_earth, r_earth)                 #=> 9.82 m/s²
NewtonGravity.surface_gravity(m_earth, r_earth)       #=> 9.82 m/s²
NewtonGravity.potential(m_earth, r_earth)             #=> -6.254e7 J/kg

# Orbital mechanics
NewtonGravity.circular_orbit_speed(m_earth, r_earth + 4.0e5)   #=> ~7660 m/s (ISS)
NewtonGravity.escape_velocity(m_earth, r_earth)                  #=> 11_186 m/s
NewtonGravity.orbital_period(3.844e8, m_earth, 7.348e22)        #=> ~2.36e6 s (Moon)

# Vis-viva: v² = GM(2/r - 1/a)
NewtonGravity.vis_viva(m_earth, r_earth + 4.0e5, r_earth + 4.0e5)
#=> circular orbit speed at ISS altitude (same as circular_orbit_speed)

# Exotic
NewtonGravity.roche_limit(1.737e6, m_earth, 7.348e22)   #=> ~9500 km
NewtonGravity.hill_sphere(1.496e11, m_earth, 1.989e30)   #=> ~1.5e9 m
NewtonGravity.tidal_acceleration(7.348e22, 3.844e8, r_earth)  #=> differential g
```

---

### Oscillations

```elixir
alias AstroEquations.Physics.Oscillations

# Spring-mass system
Oscillations.angular_frequency(100.0, 0.5)  #=> 14.14 rad/s
Oscillations.spring_period(0.5, 100.0)      #=> 0.444 s
Oscillations.pendulum_period(1.0)           #=> 2.007 s (1 m pendulum)

# SHM kinematics
Oscillations.shm_displacement(0.1, 14.14, 0.0)  #=> 0.1 m at t=0
Oscillations.shm_velocity_from_position(0.1, 14.14, 0.05)  #=> speed at x=0.05 m
Oscillations.shm_total_energy(100.0, 0.1)        #=> 0.5 J

# Damped oscillations
Oscillations.damping_ratio(2.0, 100.0, 0.5)   #=> ζ (< 1 = underdamped)
Oscillations.damped_frequency(14.14, 0.1)      #=> slightly less than omega₀
Oscillations.quality_factor(14.14, 0.5, 2.0)  #=> Q factor

# LC circuit analogy
Oscillations.lc_angular_frequency(1.0e-3, 10.0e-6)  #=> 10_000 rad/s

# Beats
Oscillations.beat_frequency(440, 444)  #=> 4.0 Hz
```

---

### Quantum Mechanics

```elixir
alias AstroEquations.Physics.QuantumMechanics, as: QM

# Fundamental relations
QM.uncertainty_principle?(1.0e-10, 1.0e-24)  #=> true
QM.de_broglie_wavelength(1.0e-24)            #=> 6.63e-10 m
QM.de_broglie_wavelength_ke(1.602e-19)       #=> ~1.23e-9 m (electron at 1 eV)

# Energy levels
QM.hydrogen_energy_level(1)  #=> -13.6 eV  (ground state)
QM.hydrogen_energy_level(2)  #=> -3.4 eV   (first excited state)
QM.rydberg_wavelength(2, 3)  #=> 6.563e-7 m  (Hα, Balmer series)
QM.infinite_well_energy(1, 9.109e-31, 1.0e-9)  #=> ~3.8e-19 J (electron in 1 nm box)

# Operators and states
QM.pauli_x()   #=> [[0,1],[1,0]]
QM.pauli_z()   #=> [[1,0],[0,-1]]

# Expectation values and uncertainty
state = [1.0, 0.0]  # |0⟩
QM.expectation_braket(QM.pauli_z(), state)   #=> 1.0   ⟨0|σ_z|0⟩
QM.variance(QM.pauli_z(), state)             #=> 0.0   (definite state)

# Two-level atom
QM.atomic_raise()   #=> [[0,1],[0,0]]   σ₊
QM.apply_raise([1.0, 0.0])   #=> [0, 1]  σ₊|g⟩ = |e⟩

# Photon ladder operators
QM.create([1.0, 0.0, 0.0], 0)    #=> [0.0, 1.0, 0.0]   â†|0⟩ = |1⟩
QM.annihilate([0.0, 1.0, 0.0], 1) #=> [1.0, 0.0, 0.0]   â|1⟩ = |0⟩

# Density matrix
rho = QM.density_matrix([[1.0, 0.0]], [1.0])
QM.purity(rho)   #=> 1.0  (pure state)

# Scattering
QM.transmission_coefficient(5.0e9, 3.0e9)   #=> T  (T + R = 1)
QM.reflection_coefficient(5.0e9, 3.0e9)     #=> R
```

---

### Special Relativity

```elixir
alias AstroEquations.Physics.SpecialRelativity, as: SR

c = 2.99792458e8

# Kinematics
SR.gamma_factor(0.9 * c)           #=> 2.2942
SR.beta(0.9 * c)                   #=> 0.9
SR.time_dilation(1.0, 0.6 * c)    #=> 1.25 s
SR.length_contraction(1.0, 0.6 * c)  #=> 0.8 m

# Velocity addition (always < c)
SR.relative_velocity(0.9 * c, -0.9 * c)  #=> 0.9945c

# Dynamics
SR.relativistic_momentum(1.0, 0.9 * c)   #=> γmv
SR.rest_energy(1.0)                        #=> 8.988e16 J
SR.kinetic_energy(1.0, 0.9 * c)          #=> (γ-1)mc²
SR.energy_momentum_relation(1.0e-22, 0)   #=> E = pc  (massless)

# Doppler
SR.relativistic_doppler(5.0e14, 0.5 * c)  #=> blueshift (~8.66e14 Hz approaching)

# Four-vectors
fv = SR.four_vector(1, 0, 0, 0)    #=> %{ct: 1, x: 0, y: 0, z: 0}
SR.spacetime_interval(1, 0, 0, 0)  #=> -1.0  (timelike)
SR.spacetime_interval(0, 1, 0, 0)  #=> +1.0  (spacelike)

# Lorentz boost (coordinate transformation)
{ct_p, x_p, y_p, z_p} = SR.lorentz_boost(0, 0, 0, 0, 0.6 * c)
```

---

### Thermodynamics

```elixir
alias AstroEquations.Physics.Thermodynamics

# Ideal gas — keyword-list interface (nil = unknown to solve for)
Thermodynamics.ideal_gas_law(p: nil, v: 0.0224, n: 6.022e23, k_b: 1.38e-23, t: 273.15)
#=> %{p: 101_325.0}

Thermodynamics.ideal_gas_law(p: 101_325, v: 0.0224, n: 6.022e23, k_b: 1.38e-23, t: nil)
#=> %{t: 273.15}

# Maxwell-Boltzmann distribution (N₂ at 300 K)
m_n2 = 4.65e-26
Thermodynamics.rms_speed(300, m_n2)           #=> 515 m/s
Thermodynamics.most_probable_speed(300, m_n2)  #=> 421 m/s
Thermodynamics.mean_speed(300, m_n2)           #=> 476 m/s

# Thermodynamic processes
Thermodynamics.first_law(500, 200)                         #=> 300 J  ΔU = Q - W
Thermodynamics.isothermal_work(6.022e23, 300, 0.001, 0.002) #=> work in isothermal expansion
Thermodynamics.adiabatic_pressure(101_325, 1.0, 0.5, 1.4)  #=> pressure after adiabatic compression

# Heat engines
Thermodynamics.carnot_efficiency(800, 300)   #=> 0.625  (62.5% theoretical max)
Thermodynamics.cop_refrigerator(300, 250)    #=> 5.0    (COP)
Thermodynamics.cop_heat_pump(300, 250)       #=> 6.0

# Radiation
Thermodynamics.wiens_displacement(5778)                    #=> 5.015e-7 m (Sun)
Thermodynamics.stefan_boltzmann(5778)                      #=> 6.316e7 W/m²
Thermodynamics.stefan_boltzmann_total(6.957e8, 5778)       #=> ~3.83e26 W  (≈ L☉)
Thermodynamics.planck_wavelength(500.0e-9, 5778)           #=> spectral radiance
```

---

### Waves

```elixir
alias AstroEquations.Physics.Waves

# Wave kinematics
Waves.wave_velocity(440, 0.780)   #=> 343.2 m/s (speed of sound)
Waves.wave_number(500.0e-9)       #=> 1.257e7 rad/m  (visible light)
Waves.angular_frequency(period: 0.002)    #=> 3141.6 rad/s
Waves.angular_frequency(frequency: 440)  #=> 2764.6 rad/s

# Interference
Waves.constructive_interference?(0.0, 500.0e-9, 0)    #=> true
Waves.destructive_interference?(250.0e-9, 500.0e-9, 0) #=> true
Waves.beat_frequency(440, 441)   #=> 1.0 Hz

# Diffraction
Waves.single_slit_minimum(1, 500.0e-9, 1.0e-4)   #=> first minimum angle
Waves.grating_angle(1, 500.0e-9, 1.0e-6)          #=> first-order maximum (~30°)
Waves.grating_resolving_power(1, 10_000)            #=> 10_000

# Sound
Waves.intensity_decibels(1.0e-6)               #=> 60.0 dB
Waves.intensity_from_decibels(60.0)            #=> 1.0e-6 W/m²
Waves.string_harmonic(1, 343, 0.686)           #=> 250 Hz (A♯4, 68.6 cm string)

# Optics
Waves.snells_law(1.0, :math.pi()/4, 1.5)      #=> refracted angle
Waves.critical_angle(1.5, 1.0)                 #=> ~41.8° (glass → air TIR)
Waves.malus_law(100.0, :math.pi()/4)           #=> 50.0 W/m² (at 45°)

# Doppler
Waves.doppler(440.0, 20.0, 0.0)               #=> ~466 Hz (approaching at 20 m/s)
Waves.relativistic_doppler(5.0e14, 0.5 * 2.998e8)  #=> relativistic blueshift
```

---

## Mathematics

### Calculus

```elixir
alias AstroEquations.Mathematics.Calculus

# Differentiation
Calculus.derivative(&:math.sin/1, 0.0)              #=> 1.0  (cos(0) = 1)
Calculus.derivative(&:math.exp/1, 1.0)              #=> 2.7183  (e)
Calculus.second_derivative(&:math.sin/1, 0.0)       #=> ~0.0  (−sin(0) = 0)
Calculus.forward_difference(&:math.sin/1, 0.0)      #=> ~1.0  (less accurate)

# Integration
Calculus.simpsons(&:math.sin/1, 0, :math.pi(), 1000)  #=> 2.0  (exact)
Calculus.trapezoid(fn x -> x * x end, 0, 1, 10_000)   #=> 0.3333  (1/3)

# Root finding
Calculus.bisection(fn x -> x * x - 2 end, 1.0, 2.0)  #=> 1.4142135623...
Calculus.newton_raphson(fn x -> x * x - 2 end, 1.5)   #=> 1.4142135623...

# ODE solvers: returns list of {x, y} pairs
Calculus.euler_ode(fn _x, y -> -y end, 0, 1.0, 1.0, 1000)
|> List.last()
#=> {1.0, 0.3679...}  (≈ e⁻¹)

Calculus.rk4_ode(fn _x, y -> -y end, 0, 1.0, 1.0, 100)
|> List.last()
#=> {1.0, 0.36787944...}  (much more accurate)
```

---

### Geometry

```elixir
alias AstroEquations.Mathematics.Geometry

# Unit conversions
Geometry.deg_to_rad(180)    #=> 3.14159...
Geometry.hours_to_deg(6)    #=> 90.0
Geometry.arcsec_to_rad(1)   #=> 4.848e-6

# Angular separation (Haversine — accurate for any angle)
Geometry.angular_separation(0.0, 0.0, :math.pi()/2, 0.0)
#=> 1.5708 rad  (90°)

Geometry.angular_separation_deg(10.0, 20.0, 10.0, 25.0)
#=> 5.0  degrees

# Orbital geometry
Geometry.periapsis(1.0, 0.5)   #=> 0.5  (a(1-e))
Geometry.apoapsis(1.0, 0.5)    #=> 1.5  (a(1+e))

# Kepler's equation: M = E - e sin(E)
e_anomaly = Geometry.eccentric_anomaly(1.5, 0.4)  #=> E given M=1.5, e=0.4
Geometry.true_anomaly(e_anomaly, 0.4)             #=> true anomaly ν

# Gravitational lensing
r_E = Geometry.einstein_radius(1.989e30, 1.0e22, 2.0e22, 1.0e22)
Geometry.microlensing_magnification(0.5)   #=> A ≈ 2.19 at u = r/r_E = 0.5

# Coordinate conversions (all in radians)
{x, y, z} = Geometry.equatorial_to_cartesian(1.5, 0.3)
Geometry.cartesian_to_equatorial(x, y, z)   #=> {1.5, 0.3}  (round-trip)
```

---

### Notation

```elixir
alias AstroEquations.Mathematics.Notation

# Physical constants
Notation.speed_of_light()        #=> 2.99792458e8   m/s
Notation.gravitational_constant()  #=> 6.6743e-11    m³ kg⁻¹ s⁻²
Notation.planck_constant()        #=> 6.62607015e-34 J·s
Notation.boltzmann_constant()     #=> 1.380649e-23   J/K
Notation.solar_mass()             #=> 1.98892e30     kg
Notation.parsec()                 #=> 3.086e16       m

# Unit conversions
Notation.ev_to_j(13.6)           #=> 2.179e-18 J  (ionisation energy)
Notation.pc_to_m(1)              #=> 3.086e16 m
Notation.au_to_m(1)              #=> 1.496e11 m
Notation.deg_to_rad(45)          #=> 0.7854 rad
Notation.k_to_celsius(5778)      #=> 5504.85 °C

# Sexagesimal conversions
Notation.hms_to_deg(5, 34, 32.0)   #=> 83.633°  (Orion belt, RA)
Notation.dms_to_deg(-5, 23, 28.0)  #=> -5.391°  (Orion belt, Dec)

{h, m, s} = Notation.deg_to_hms(83.633)
#=> {5, 34, 31.92}

# Julian Date
Notation.calendar_to_jd(2000, 1, 1, 12, 0, 0)   #=> 2_451_545.0  (J2000.0)
Notation.j2000()                                  #=> 2_451_545.0
Notation.julian_centuries(2_451_545.0)            #=> 0.0
```

---

### Trigonometry

```elixir
alias AstroEquations.Mathematics.Trigonometry

# Spherical triangle: cosine and sine rules
Trigonometry.spherical_law_of_cosines(:math.pi()/3, :math.pi()/3, :math.pi()/3)
#=> ~1.047 rad  (equilateral spherical triangle)

# Altitude and azimuth of an object
lat = 51.5 * :math.pi() / 180   # London
dec = 0.0                         # celestial equator
ha  = 0.0                         # on meridian
Trigonometry.altitude(dec, lat, ha)   #=> 0.669 rad  (38.5° altitude)
Trigonometry.azimuth(dec, lat, ha)    #=> π  (due south)

# Sidereal time
jd = Notation.calendar_to_jd(2024, 3, 20, 0, 0, 0)
Trigonometry.gmst(jd)                       #=> GMST in radians
Trigonometry.local_sidereal_time(jd, 0.0)  #=> LST at Greenwich

# Atmospheric refraction at 30° altitude
Trigonometry.atmospheric_refraction(30.0)  #=> ~1.67 arcmin

# Coordinate conversions (inputs and outputs in radians)
{lambda, beta} = Trigonometry.equatorial_to_ecliptic(1.5, 0.3)
{l, b}         = Trigonometry.equatorial_to_galactic(1.5, 0.3)
Trigonometry.galactic_to_equatorial(0.0, 0.0)   #=> RA/Dec of galactic centre
```

---

## Statistics

### Variance

```elixir
alias AstroEquations.Statistics.Variance

data = [2.0, 4.0, 4.0, 4.0, 5.0, 5.0, 7.0, 9.0]

# Basic statistics
Variance.population(data)   #=> 4.0   (σ² = Σ(x-μ)²/N)
Variance.sample(data)       #=> 4.571 (s² = Σ(x-x̄)²/(N-1), Bessel-corrected)

# Weighted variance
weights = [1, 1, 1, 2, 1, 1, 1, 1]  # double-weight the fourth observation
Variance.weighted(Enum.map(data, &(&1 * 1.0)), Enum.map(weights, &(&1 * 1.0)))

# Standard error of the mean
Variance.standard_error(data)   #=> σ/√N

# Welford single-pass streaming (memory efficient for large datasets)
acc = Variance.welford_accumulate(data)
{:ok, mean, pop_var, samp_var} = Variance.welford_finalize(acc)
# mean = 5.0, pop_var = 4.0, samp_var = 4.571

# Useful for streams or chunked data:
1..1_000_000
|> Variance.welford_accumulate()
|> Variance.welford_finalize()
#=> {:ok, 500_000.5, ..., ...}

# Correlation
Variance.pearson_r([1,2,3,4,5], [2,4,6,8,10])   #=> 1.0  (perfectly linear)
Variance.covariance([1,2,3], [1,2,3])            #=> 1.0

# Goodness of fit
obs = [10, 20, 30]
exp = [10, 20, 30]
Variance.chi_squared(obs, exp)          #=> 0.0   (perfect fit)
Variance.reduced_chi_squared(9.0, 3)   #=> 3.0

# Robust statistics (less sensitive to outliers)
Variance.mad([1, 2, 100])   #=> 1.0  (vs std dev of ~57)
Variance.iqr([1, 2, 3, 4, 5, 6, 7, 8])  #=> 4.0
```

---

### Standard Deviation

```elixir
alias AstroEquations.Statistics.StandardDeviation, as: SD

data = [2.0, 4.0, 4.0, 4.0, 5.0, 5.0, 7.0, 9.0]

SD.population(data)   #=> 2.0    (σ = √4)
SD.sample(data)       #=> 2.138  (s)
SD.standard_error(data)  #=> s/√N

# Z-scores and normalisation
SD.z_score(7.0, 5.0, 2.0)  #=> 1.0  (one sigma above mean)
normalised = SD.normalise(data)
# => list with mean=0, std=1

# Error propagation
SD.propagate_addition(0.05, 0.03)          #=> 0.0583  (quadrature sum)
SD.propagate_product(10.0, 0.5, 5.0, 0.25) #=> absolute uncertainty on 10×5
SD.propagate_power(3.0, 0.1, 2)            #=> σ_y = 2x σ_x = 0.6
SD.propagate_log(100.0, 5.0)              #=> σ_y = σ_x/x = 0.05

# Photometry
SD.photometric_snr(10_000, 500, 10, 5)  #=> CCD SNR
SD.magnitude_uncertainty(100.0)          #=> ~0.0109 mag  (1/SNR × 1.0857)
SD.poisson_uncertainty(250)              #=> 15.81  (√N)

# One-pass streaming mean and std dev
{:ok, mean, pop_std, samp_std} = SD.online([1.0, 2.0, 3.0, 4.0, 5.0])
# mean = 3.0, samp_std = 1.581
```

---

## Patterns & Tips

### Pipe-based survey pipelines

```elixir
alias AstroEquations.AstrophysicsAndAstronomy.Astrometry
alias AstroEquations.Statistics.{Variance, StandardDeviation}

# Compute weighted distance moduli for a set of standard candles
observations = [
  %{flux: 1.2e-15, flux_ref: 3.631e-23, parallax: nil, period: 10.2},
  %{flux: 0.8e-15, flux_ref: 3.631e-23, parallax: nil, period: 6.7},
]

distance_moduli =
  observations
  |> Enum.map(& Astrometry.apparent_magnitude_diff(&1.flux, &1.flux_ref))
  |> Enum.map(& Astrometry.distance_from_modulus(&1))

mean_distance =
  distance_moduli
  |> Variance.welford_accumulate()
  |> Variance.welford_finalize()
  |> then(fn {:ok, mean, _, _} -> mean end)
```

### Solving stellar structure ODEs with RK4

```elixir
alias AstroEquations.Mathematics.Calculus
alias AstroEquations.AstrophysicsAndAstronomy.Stars

# Simple polytropic star: dP/dr = hydrostatic equilibrium
g = 6.674e-11

# One step of the stellar structure ODE
dPdr = fn r, {p, m, l} ->
  rho = (p / 1.38e6) ** (3/5)  # polytropic P = K ρ^(5/3)
  dP = Stars.hydrostatic_equilibrium(m, rho, r)
  dM = Stars.mass_conservation(rho, r)
  dL = Stars.energy_equation(rho, 1.934e-7, r)
  {dP, dM, dL}
end
```

### Relativistic astronomy pipeline

```elixir
alias AstroEquations.Physics.SpecialRelativity, as: SR
alias AstroEquations.AstrophysicsAndAstronomy.Astrometry

# Convert observed redshift to recession velocity and Lorentz factor
z = 0.15
v_nr  = Astrometry.recession_velocity(z)          # non-relativistic
v_rel = Astrometry.relativistic_velocity(z)        # relativistic
gamma = SR.gamma_factor(v_rel * 1000)              # Lorentz factor

# Time dilation: how much younger is the source?
proper_time_ratio = 1.0 / gamma
```

### Error propagation for photometric pipeline

```elixir
alias AstroEquations.Statistics.StandardDeviation, as: SD

# Full photometric uncertainty budget
n_star   = 5_000.0   # source electrons
n_sky    = 200.0     # sky per pixel
n_pix    = 25        # pixels in aperture
n_dark   = 0.01      # dark current per pixel per second
read_noise = 4.0     # electrons RMS

snr = SD.photometric_snr(n_star, n_sky, n_pix, read_noise * read_noise)
sigma_mag = SD.magnitude_uncertainty(snr)

# Propagate through a colour index calculation
sigma_b  = 0.02  # mag uncertainty in B band
sigma_v  = 0.015 # mag uncertainty in V band
sigma_bv = SD.propagate_addition(sigma_b, sigma_v)  #=> 0.025 mag
```
