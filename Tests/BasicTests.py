from Utils import *


def MoveRight():
    snap_shot_t = 0.01
    final_t = 10
    dt = 0.001
    init_bodies = [Body((0,0,0), (1.0,0,0), 1.0)]
    sim_setup = SimArgs(snap_shot_t, final_t, dt, init_bodies)

    final_bodies = [(10, 0, 0)]
    expected_res = SimSolution(len(final_bodies), final_bodies, float(final_t)/float(dt))

    return (sim_setup, expected_res)
    
MoveRightTest = TestCase("Move Right", "Test to see if a body moves right", MoveRight)

def MoveDiag():
    snap_shot_t = 0.01
    final_t = 10
    dt = 0.001
    init_bodies = [Body((0,0,0), (1.0, -1.0, 2.0), 1.0)]
    sim_setup = SimArgs(snap_shot_t, final_t, dt, init_bodies)

    final_bodies = [(10, -10, 20)]
    expected_res = SimSolution(len(final_bodies), final_bodies, float(final_t)/float(dt))

    return (sim_setup, expected_res)
    
MoveDiagTest = TestCase("Move Diagonally", "Test to see if a body moves diagonally", MoveDiag)