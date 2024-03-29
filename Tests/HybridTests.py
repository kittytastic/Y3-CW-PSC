import math
from Utils import *

def StandingCubeCollide():
    snap_shot_t = 0.1
    final_t = 1.4
    dt = 0.000001

    staring_d = 1.5
    init_bodies = [
        # Origin
        Body((0,0,0), (0,0,0), 1.0),
        # X
        Body((staring_d,0,0), (0,0,0), 1.0),
        Body((-staring_d,0,0), (0,0,0), 1.0),
        # Y
        Body((0,-staring_d,0), (0,0,0), 1.0),
        Body((0, staring_d,0), (0,0,0), 1.0),
        #Z
        Body((0,0,-staring_d), (0,0,0), 1.0),
        Body((0,0, staring_d), (0,0,0), 1.0),
    ]
    sim_setup = SimArgs(snap_shot_t, final_t, dt, init_bodies)

    final_bodies = [(0.0, 0.0, 0.0)]
    expected_res = SimSolution(len(final_bodies), final_bodies, float(final_t)/float(dt), v_max=0, dx_min=math.inf)

    return {"setup":sim_setup, "expected":expected_res, "error":1e-6}
    
StandingCubeColide = TestCase("Standing Cube Collision", "Test to see if 6 points around an origin point accelerate and combine", StandingCubeCollide)

def MovingThreeCollide():
    snap_shot_t = 1
    final_t = 10
    dt = 0.01
    init_bodies = [
        Body((0,0,0), (0,0,0), 1.0),
        Body((2,0,0), (0,0,0), 1.0),
        Body((-2,0,0), (0,0,0), 1.0)
    ]
    sim_setup = SimArgs(snap_shot_t, final_t, dt, init_bodies)

    final_bodies = [(0, 0, 0)]
    expected_res = SimSolution(len(final_bodies), final_bodies, float(final_t)/float(dt), v_max=0, dx_min=math.inf)

    return {"setup": sim_setup, "expected":expected_res}
    
MovingThreeColideTest = TestCase("Dynamic 3 Collision", "Test to see if 3 stationary points in collision distance are combined", MovingThreeCollide)

def OrbitThreeCollide():
    snap_shot_t = 1
    final_t = 10
    dt = 0.01
    init_bodies = [
        Body((0,0,0), (0,0,0), 1.0),
        Body((2,0,0), (0,0.01,0), 1.0),
        Body((-2,0,0), (0,-0.01,0), 1.0)
    ]
    sim_setup = SimArgs(snap_shot_t, final_t, dt, init_bodies)

    final_bodies = [(0, 0, 0)]
    expected_res = SimSolution(len(final_bodies), final_bodies, float(final_t)/float(dt), v_max=0, dx_min=math.inf)

    return {"setup": sim_setup, "expected":expected_res}
    
OrbitThreeColideTest = TestCase("Orbit 3 Collision", "Test to see if 3 stationary points in collision distance are combined", OrbitThreeCollide)


def MovingCubeCollide():
    snap_shot_t = 1
    final_t = 10
    dt = 0.001

    s = 2 
    init_bodies = [
        # Origin
        Body((0,0,0), (0,0,0), 1.0),
        # X
        Body((s,0,0), (0,0,0), 1.0),
        Body((-s,0,0), (0,0,0), 1.0),
        # Y
        Body((0,-s,0), (0,0,0), 1.0),
        Body((0, s,0), (0,0,0), 1.0),
        #Z
        Body((0,0,-s), (0,0,0), 1.0),
        Body((0,0, s), (0,0,0), 1.0),
    ]
    sim_setup = SimArgs(snap_shot_t, final_t, dt, init_bodies)

    final_bodies = [(0, 0, 0)]
    expected_res = SimSolution(len(final_bodies), final_bodies, float(final_t)/float(dt), v_max=0, dx_min=math.inf)

    return {"setup": sim_setup, "expected":expected_res}
    
MovingCubeColideTest = TestCase("Dynamic Cube Collision", "Test to see if 3 stationary points in collision distance are combined", MovingCubeCollide)