# Compiler and Flags Configuration
FC       := gfortran
FCFLAGS  := -ffpe-summary=none

.PHONY: all
all: distance_distribution concentrations clean

distance_distribution: src/math_lib/inv_Langevin_function.f src/halfspace_WLC_simulations/WLC_simulation_posZ_endpoint.f src/halfspace_WLC_simulations/WLC_distance_distribution.f
	$(FC) $(FCFLAGS) -c src/math_lib/inv_Langevin_function.f
	$(FC) $(FCFLAGS) -c src/halfspace_WLC_simulations/WLC_simulation_posZ_endpoint.f
	$(FC) $(FCFLAGS) -c src/halfspace_WLC_simulations/WLC_distance_distribution.f
	$(FC) $(FCFLAGS) WLC_distance_distribution.o WLC_simulation_posZ_endpoint.o inv_Langevin_function.o -o WLC_distance_distribution

concentrations: src/math_lib/inv_Langevin_function.f  src/halfspace_WLC_simulations/WLC_simulation_posZ_stepwise.f src/halfspace_WLC_simulations/WLC_concentration_ABCD.f
	$(FC) $(FCFLAGS) -c src/math_lib/inv_Langevin_function.f
	$(FC) $(FCFLAGS) -c src/halfspace_WLC_simulations/WLC_simulation_posZ_stepwise.f
	$(FC) $(FCFLAGS) -c src/halfspace_WLC_simulations/WLC_concentration_ABCD.f
	$(FC) $(FCFLAGS) WLC_concentration_ABCD.o WLC_simulation_posZ_stepwise.o inv_Langevin_function.o -o WLC_concentration_ABCD

.PHONY: clean
clean:
	rm inv_Langevin_function.o WLC_simulation_posZ_endpoint.o WLC_simulation_posZ_stepwise.o WLC_distance_distribution.o WLC_concentration_ABCD.o
