from Utils import *

def TwoCollide():
    snap_shot_t = 1
    final_t = 1
    dt = 1
    init_bodies = [
        Body((0,0,0), (0,0,0), 1.0),
        Body((2e-3,0,0), (0,0,0), 1.0)
    ]
    sim_setup = SimArgs(snap_shot_t, final_t, dt, init_bodies)

    final_bodies = [(1e-3, 0, 0)]
    expected_res = SimSolution(len(final_bodies), final_bodies, float(final_t)/float(dt), v_max=0, dx_min=2e-3)

    return (sim_setup, expected_res)
    
StationaryColide = TestCase("Stationary 2 Collision", "Test to see if 2 stationary points in collision distance are combined", TwoCollide)