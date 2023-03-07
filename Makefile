#####################################################
###                                               ###
###             Makefile for GRAPUS               ###
###                                               ###
###         Duncan H. Forgan (24/01/2018)         ###
###       				          ###
###                                               ###
#####################################################

# Compiler variable:
FC     = gfortran
VPATH = src/main/ src/disc/ src/embryo/ src/eos/ src/io/ src/nbody/ src/planet/ src/wind/ src/collapse/


# For serial runs use these flags
FFLAGS = -O3 -frecord-marker=4 -fdefault-real-8 -fdefault-double-8 -fimplicit-none -fbounds-check -ffpe-trap=invalid,zero,overflow -Wunused -fbacktrace

# For OpenMP runs
##FFLAGS = -O0 -g -frecord-marker=4 -fdefault-real-8 -fdefault-double-8 -fimplicit-none -ffpe-trap=invalid,zero,overflow -fbounds-check -fcheck=all -fopenmp -Wunused -fbacktrace
##FFLAGS = -O0 -g -fdefault-real-8 -fdefault-double-8 -fimplicit-none -fopenmp -fbacktrace
#FFLAGS = -O3 -frecord-marker=4 -fdefault-real-8 -fdefault-double-8 -fimplicit-none -fbounds-check -fopenmp -Wunused

# Create object files:

%.o: %.f90
	$(FC) $(FFLAGS) -c $<

%.o : %.F90
	$(FC) -c $(FFLAGS) -c $< -o $@

SOURCESAF90 = unit_module.f90 star_module.f90 embryo_module.f90 eosmodule.f90 \
	planet_module.f90 wind_module.f90 cloud_module.f90 main.f90 \
	accrete_gas.f90 calc_core_formation.f90 \
	calc_core_radiative_feedback.f90 \
	calc_grain_growth.f90 calc_grain_sedimentation.f90 \
	calc_tidal_disruption.f90 calc_typeI_migration.f90 check_embryo_thermal_state.f90 \
	check_for_dead_embryos.f90 compute_planet_torques.f90 compute_wind.f90 \
	disc_properties.f90 Eacc_calc.f90 eosread.f90 eos_cs.f90 eos_T.f90 \
	evolve.f90  evolve_disc_interpmodel.f90 evolve_radius.f90 \
	evolve_temperature.f90 evolve_embryos.f90 find_planets_in_disc.f90 \
	generate_disc.f90 generate_embryos.f90 \
	generate_star.f90 initial.f90 \
	interpolate_1D.f90 interpolate_2D.f90 mdotcalc.f90 \
	move_embryos.f90 migration_timescales.f90 \
	midplane_properties.f90 midplane_properties_grav_fixedalpha_1.f90 \
        midplane_properties_grav_fixedalpha_2.f90 \
	migrate_planets.f90 nbody_integrate.f90 nbody_timestep.f90 nbody_acceleration.f90 \
	nbody_grav_acceleration.f90 nbody_drag_terms.f90 \
	nbody_deallocate_arrays.f90 nbody_orbits.f90 \
	nbody_output.f90 nbody_system_properties.f90  nbody_rk4.f90 \
	sample_gaussian.f90 setup_cloud.f90 setup_planets.f90 \
	setup_wind.f90 sigma_mdot.f90 timestep.f90 wind_profiles.f90 \
	write_population_snapshot.f90 write_dump.F90 xraywind.f90

OBJECTSA    = $(SOURCESAF90:.f90=.o)
OBJECTS     = $(OBJECTSA:.F90=.o)

# Create executable files:
build: grapusvisag

grapusvisag: $(OBJECTS) 
	$(FC) $(FFLAGS) -o $@ $(OBJECTS)
 
# Clean statements:
clean: 
	\rm *.o *.mod grapusvisag

# End Makefile
