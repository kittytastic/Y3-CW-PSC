import math
from Utils import *

def TwoCollide():
    snap_shot_t = 1
    final_t = 3
    dt = 0.001
    init_bodies = [
        Body((0,0,0), (0,0,0), 1.0),
        Body((2,0,0), (0,0,0), 1.0)
    ]
    sim_setup = SimArgs(snap_shot_t, final_t, dt, init_bodies)

    final_bodies = [(1, 0, 0)]
    expected_res = SimSolution(len(final_bodies), final_bodies, float(final_t)/float(dt), v_max=0, dx_min=math.inf)

    return {"setup": sim_setup, "expected":expected_res}
    
StationaryColide = TestCase("Stationary 2 Collision", "Test to see if 2 stationary points in collision distance are combined", TwoCollide)


def ThreeCollide():
    snap_shot_t = 1
    final_t = 2
    dt = 1
    init_bodies = [
        Body((0,0,0), (0,0,0), 1.0),
        Body((2e-3,0,0), (0,0,0), 1.0),
        Body((-2e-3,0,0), (0,0,0), 1.0)
    ]
    sim_setup = SimArgs(snap_shot_t, final_t, dt, init_bodies)

    final_bodies = [(0, 0, 0)]
    expected_res = SimSolution(len(final_bodies), final_bodies, float(final_t)/float(dt), v_max=0, dx_min=math.inf)

    return {"setup": sim_setup, "expected":expected_res}
    
StationaryThreeColide = TestCase("Stationary 3 Collision", "Test to see if 3 stationary points in collision distance are combined", ThreeCollide)

def CubeCollide():
    snap_shot_t = 1
    final_t = 2
    dt = 1
    init_bodies = [
        # Origin
        Body((0,0,0), (0,0,0), 1.0),
        # X
        Body((2e-3,0,0), (0,0,0), 1.0),
        Body((-2e-3,0,0), (0,0,0), 1.0),
        # Y
        Body((0,-2e-3,0), (0,0,0), 1.0),
        Body((0, 2e-3,0), (0,0,0), 1.0),
        #Z
        Body((0,0,-2e-3), (0,0,0), 1.0),
        Body((0,0, 2e-3), (0,0,0), 1.0),
    ]
    sim_setup = SimArgs(snap_shot_t, final_t, dt, init_bodies)

    final_bodies = [(0, 0, 0)]
    expected_res = SimSolution(len(final_bodies), final_bodies, float(final_t)/float(dt), v_max=0, dx_min=math.inf)

    return {"setup": sim_setup, "expected":expected_res}
    
StationaryCubeColide = TestCase("Stationary Cube Collision", "Test to see if 6 points around a origin combine", CubeCollide)