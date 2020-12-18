from Utils import *
import math

def TwoBodyOrbit():
    snap_shot_t = 0.01
    final_t = 10
    dt = 0.001
    init_bodies = [Body((1,0,0), (0,0.5,0), 2.0), Body((-1,0,0), (0,-0.5,0), 2.0)]
    sim_setup = SimArgs(snap_shot_t, final_t, dt, init_bodies)

    final_bodies = [(1.06395, -0.0442575, 0), (-1.06647, 0.04427, 0)]
    expected_res = SimSolution(len(final_bodies), final_bodies, float(final_t)/float(dt), v_max=0.478111503159813, dx_min=2.13225348945739)

    return {"setup": sim_setup, "expected":expected_res}
    
TwoBodyOrbitTest = TestCase("Two Body Orbit", "Test to see if 2 body orbits match experimental data", TwoBodyOrbit)