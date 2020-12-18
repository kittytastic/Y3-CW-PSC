from Utils import *
import math

import random

def TwoBodyOrbit():
    snap_shot_t = 0.01
    final_t = 10
    dt = 0.001
    init_bodies = [Body((1,0,0), (0,0.5,0), 2.0), Body((-1,0,0), (0,-0.5,0), 2.0)]
    sim_setup = SimArgs(snap_shot_t, final_t, dt, init_bodies)

    final_bodies = [(1.06395, -0.0442575, 0), (-1.06647, 0.04427, 0)]
    expected_res = SimSolution(len(final_bodies), final_bodies, float(final_t)/float(dt), v_max=0.478111503159813, dx_min=2.13225348945739)

    return {"setup": sim_setup, "expected":expected_res, "error": 10e-3}
    
TwoBodyOrbitTest = TestCase("Two Body Orbit", "Test to see if 2 body orbits match experimental data", TwoBodyOrbit)


def FiftyBody():
    snap_shot_t = 0.01
    final_t = 10
    dt = 0.00001

    world_bounds = 50
    v_bounds = 4
    m_bounds = 5

    random.seed(42)
    scaleX = lambda x: (x-0.5)*world_bounds
    scaleV = lambda x: (x-0.5)*v_bounds
    scaleM = lambda x: x*m_bounds

    init_bodies = [Body(
        (scaleX(random.random()), scaleX(random.random()), scaleX(random.random())),
        (scaleV(random.random()), scaleV(random.random()), scaleV(random.random())),
        scaleM(random.random())
        ) for i in range(50)]
    sim_setup = SimArgs(snap_shot_t, final_t, dt, init_bodies)

    final_bodies = [(-3.52304, -8.66556, -5.996), (-28.6653, -3.98223, -38.5197), (8.84364, 12.2392, -28.6022), (25.1154, -15.1033, -27.5478), (22.2413, 7.28363, 30.8867), (12.1514, 13.0284, 8.13657), (-21.6512, -23.3154, -25.5679), (-10.9414, 11.6706, -12.0444), (-18.3298, 27.7895, -8.48581), (0.102802, -2.20614, 4.06237), (-20.7538, 24.3291, 10.8322), (2.99047, -16.7246, -3.06008), (-2.20349, -27.3313, 4.73003), (-16.6738, 5.0203, 13.9605), (28.6928, 23.0105, -20.0179), (4.60876, 2.60021, 9.14095), (16.0983, -21.519, 23.5757), (-2.31183, 11.6816, -16.4337), (-5.91751, 13.7704, 0.568219), (0.269422, -6.40897, 7.62847), (12.5987, -2.99819, -19.6148), (20.7174, -10.9924, -7.2039), (-11.9787, -20.2078, -18.4765), (16.7696, -27.4787, 7.74498), (-8.52826, -16.9069, 29.987), (-15.0437, -18.0874, 5.06758), (6.67359, -9.47221, -8.17144), (-11.1883, 11.6313, -6.89509), (14.8003, -1.70881, 36.1247), (-16.2949, -18.1914, -16.8649), (-14.0668, 27.8652, 17.6112), (5.08903, -16.4742, 20.8782), (-3.65021, -9.30188, 17.1325), (-22.6914, -10.3198, 22.2009), (6.35693, 5.92147, 13.7081), (24.0684, -15.7025, -14.9453), (2.49565, -19.3187, 11.033), (13.1389, 2.4667, 28.773), (-14.6202, 31.4493, -0.669322), (3.60516, -13.5777, -23.194), (0.523728, -13.5393, 21.3898), (-33.7216, 6.95172, -6.41617), (-7.89119, -14.2877, 10.8081), (-6.91134, 0.662791, -2.8471), (1.91237, 3.63421, 1.85648), (1.50896, 12.2915, 17.2987), (23.2332, 9.31823, -3.82763), (25.4499, -2.08479, 3.52438), (-15.877, 21.5877, -9.46673)]
    expected_res = SimSolution(len(final_bodies), final_bodies, float(final_t)/float(dt), v_max=3.40921332558821, dx_min=5.02924144197184)

    return {"setup": sim_setup, "expected":expected_res, "error": 10e-5}
    
FiftyBodyTest = TestCase("50 Body", "Test to see if 50 bodies match experimental data", FiftyBody)