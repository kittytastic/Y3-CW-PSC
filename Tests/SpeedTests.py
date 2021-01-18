from Utils import *
import random
def randomBodies(num_bodies, world_bounds, v_bounds, m_bounds):

    random.seed(42)
    scaleX = lambda x: (x-0.5)*world_bounds
    scaleV = lambda x: (x-0.5)*v_bounds
    scaleM = lambda x: x*m_bounds

    init_bodies = [Body(
        (scaleX(random.random()), scaleX(random.random()), scaleX(random.random())),
        (scaleV(random.random()), scaleV(random.random()), scaleV(random.random())),
        scaleM(random.random())
        ) for i in range(num_bodies)]

    return init_bodies


def HugeGalaxy():
    snap_shot_t = 0.1
    final_t = 10
    dt = 0.1

    bodies = randomBodies(10000, 10000, 4, 5)
    sim_setup = SimArgs(snap_shot_t, final_t, dt, bodies)

    return {"setup":sim_setup}

HugeGalaxySpeedTest = SpeedTestCase("Galaxy Speed", "Loads of bodies", HugeGalaxy)


def QuickGalaxy():
    snap_shot_t = 0.1
    final_t = 10
    dt = 0.1

    bodies = randomBodies(1000, 10000, 4, 5)
    sim_setup = SimArgs(snap_shot_t, final_t, dt, bodies)

    return {"setup":sim_setup}

QuickGalaxySpeedTest = SpeedTestCase("Quick Galaxy Speed", "Many bodies, few timesteps", QuickGalaxy)

def MediumGalaxy():
    snap_shot_t = 0.1
    final_t = 10
    dt = 0.1

    bodies = randomBodies(3000, 10000, 4, 5)
    sim_setup = SimArgs(snap_shot_t, final_t, dt, bodies)

    return {"setup":sim_setup}

MediumGalaxySpeedTest = SpeedTestCase("Medium", "Many bodies, few timesteps", MediumGalaxy)

def FineGalaxy():
    snap_shot_t = 0.1
    final_t = 10
    dt = 0.0025

    bodies = randomBodies(500, 10000, 4, 5)
    sim_setup = SimArgs(snap_shot_t, final_t, dt, bodies)

    return {"setup":sim_setup}

FineGalaxySpeedTest = SpeedTestCase("Fine grain speed", "Medium bodies at small ts", FineGalaxy)

def UltraFineGalaxy():
    snap_shot_t = 0.1
    final_t = 10
    dt = 0.00001

    bodies = randomBodies(20, 10000, 4, 5)
    sim_setup = SimArgs(snap_shot_t, final_t, dt, bodies)

    return {"setup":sim_setup}

UltraFineGalaxySpeedTest = SpeedTestCase("Ultra Fine grain speed", "Low bodies at small ts", UltraFineGalaxy)

def ClusterF():
    snap_shot_t = 0.1
    final_t = 10
    dt = 0.002

    bodies = randomBodies(500, 60, 4, 5)
    sim_setup = SimArgs(snap_shot_t, final_t, dt, bodies)

    return {"setup":sim_setup, "cluster":True}

ClusterFSpeedTest = SpeedTestCase("Cluster F speed", "Low bodies at small ts", ClusterF)