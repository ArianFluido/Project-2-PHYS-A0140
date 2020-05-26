#! /usr/bin/env python
import friction_tools as ft
import numpy as np
import os

# global constants
lattice_constant_Au = 4.0782
room_temperature = 270
moving_force = 0.04

# global variables
filename = 'simulation'
global_interval = 5

def get_equilibrium_distance(xy_cells):
    return np.divide(np.sqrt(50), xy_cells).item()

def get_xy_cells(lattice_constant):
    return int(np.round(10.0/lattice_constant))

xy_cells_Au = get_xy_cells(lattice_constant_Au)
equilibrium_distance_Au = get_equilibrium_distance(get_xy_cells(lattice_constant_Au))

forces = [0, 0.0100, 0.0200, 0.0300, 0.0400, 0.0500]

for force in forces:
    simu = ft.FrictionSimulation()

    # create two slabs of gold
    simu.create_slab(element='Au',xy_cells=xy_cells_Au,z_cells=3,bottom_z=12.0)
    simu.create_slab(element='Au',xy_cells=xy_cells_Au,z_cells=3,top_z=0.0)

    # check the system
    simu.list_atoms()

    # create interaction between the two slabs
    simu.create_interaction(['Au','Au'], strength=1.0, equilibrium_distance=equilibrium_distance_Au)

    top_indices = simu.get_indices_z_more_than(12.0+6)
    bot_indices= simu.get_indices_z_less_than(-3.5)

    # add velocity to the moving slab
    simu.fix_velocities(indices=top_indices, velocity=[0, 0.005, 0], xyz=[True,True,False])

    # set moving force of the moving slab
    simu.add_constant_force(top_indices, [0,0,-force])
    
    # make the stationary slab stationary
    simu.fix_positions(bot_indices)

    simu.set_temperature(temperature=room_temperature)
    simu.create_dynamics(dt=global_interval, temperature=room_temperature, coupled_indices=bot_indices)

    # capture data
    simu.gather_energy_and_temperature_during_simulation(interval=global_interval, filename='{}/data/energy_{}.txt'.format(os.getcwd(), force))
    simu.save_trajectory_during_simulation(interval=global_interval, filename='{}/data/{}_{}.traj'.format(os.getcwd(), filename, force)) # 5 fs
    simu.gather_average_force_during_simulation(interval=global_interval,indices=top_indices,filename='{}/data/forces_{}.txt'.format(os.getcwd(),force))
    
    # - monitor the simulation by printing info to stdout
    simu.print_stats_during_simulation(interval=50)

    # time the simulation
    import time
    t0 = time.time()

    # run the simulation
    print "starting simulation"
    simu.run_simulation(time=3000.0*5)
    print "finished simulation"
    
    t1 = time.time()
    print "time taken {ti} s".format(ti=str(int(t1-t0)))
    
    # Save the results
    file=np.loadtxt('{}/data/forces_{}.txt'.format(os.getcwd(),force))
    ft.trajectory_to_xyz(filename_in='{}/data/{}_{}.traj'.format(os.getcwd(), filename, force), filename_out='{}/data/{}_{}.xyz'.format(os.getcwd(), filename, force))
