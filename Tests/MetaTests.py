import math
from Utils import *

def Meta():
    snap_shot_t = 0.5
    final_t = 3
    dt = 0.0000001

    staring_d = 1.5
    init_bodies = [
        Body((0,0, 0), (0,-1,0), 1.0),
        Body((0,0, 5), (0,1,0), 1.0),
    ]
    sim_setup = SimArgs(snap_shot_t, final_t, dt, init_bodies)

    final_bodies = [(0.0, -2.95259, 0.142064), (0.0, 2.95259, 4.85794)]
    expected_res = SimSolution(len(final_bodies), final_bodies, float(final_t)/float(dt), v_max=0.965569703992783, dx_min=7.55715938141388)

    return {"setup":sim_setup, "expected":expected_res, "error": 1e-8}
    
MetaTest = TestCase("Meta Test", "Check v_max and dx_min", Meta)