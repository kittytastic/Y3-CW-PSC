from Utils import *
import random

def GalaxySpeedTest():
    snap_shot_t = 0.1
    final_t = 10
    dt = 0.1

    world_bounds = 10000
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
        ) for i in range(10000)]
    sim_setup = SimArgs(snap_shot_t, final_t, dt, init_bodies)

    return {"setup":sim_setup}

GalaxySpeedTest = SpeedTestCase("Galaxy Speed", "Loads of bodies", GalaxySpeedTest)