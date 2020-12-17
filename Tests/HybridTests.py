import math
from Utils import *

def StandingCubeCollide():
    snap_shot_t = 0.1
    final_t = 1.35
    dt = 0.0000001

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

    return (sim_setup, expected_res)
    
StandingCubeColide = TestCase("Standing Cube Collision", "Test to see if 6 points around an origin point accelerate and combine", StandingCubeCollide)