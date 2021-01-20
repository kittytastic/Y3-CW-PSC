import math
from Utils import *
import random

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


def PartialCollapse():
    snap_shot_t = 1
    final_t = 3
    dt = 0.000001

    wb = 12
    vb = 10
    mb = 10
    bodies = 15

    random.seed(0)
    scaleX = lambda x: (x-0.5)*wb
    scaleV = lambda x: (x-0.5)*vb
    scaleM = lambda x: x*mb

    init_bodies = [Body(
        (scaleX(random.random()), scaleX(random.random()), scaleX(random.random())),
        (scaleV(random.random()), scaleV(random.random()), scaleV(random.random())),
        scaleM(random.random())
        ) for i in range(bodies)]

    sim_setup = SimArgs(snap_shot_t, final_t, dt, init_bodies)

    final_bodies = [('3.8639', '0.694393', '1.47927'), ('3.13977', '-1.21242', '-9.26418'), ('15.6299', '8.22849', '10.2423'), ('0.666698', '8.25475', '2.05358'), ('14.5213', '15.0794', '3.34475'), ('5.70571', '6.2821', '-13.5011'), ('13.6096', '-9.67883', '1.0752')]
    expected_res = SimSolution(len(final_bodies), final_bodies, float(final_t)/float(dt), v_max=5.79448857476936, dx_min=8.22866500731216)

    return {"setup":sim_setup, "expected":expected_res, "error":1e-5}

PartialCollapseTest = TestCase("Partial Collapse", "Only some bodies merge", PartialCollapse)