from Utils import *


def lateralMovement():
    snap_shot_t = 0.01
    final_t = 10
    dt = 0.001
    init_bodies = [Body((0,0,0), (1.0,0,0), 1.0)]
    sim_setup = SimArgs(snap_shot_t, final_t, dt, init_bodies)

    final_bodies = [(10, 0, 0)]
    expected_res = SimSolution(len(final_bodies), final_bodies, float(final_t)/float(dt))

    return (sim_setup, expected_res)
    

lateralMovementTest = TestCase("Lateral Movement", "Test to see if a body moves left to right", lateralMovement)