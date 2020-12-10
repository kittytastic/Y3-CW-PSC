from Utils import *


def lateralMovement(bin_name):
    #res = run_bin_with_args(bin_name, "0.01  100.0  0.001  0.0 0.0 0.0  1.0 0.0 0.0  1.0")
    res = run_bin_with_args(bin_name, "0.01 100.0 0.0001 1.0 0.0 0.0 0.0 1.0 0.0 0.4 0.0 0.0 0.0 -0.1 0.0 0.0 2.0 2.0 0.0 0.0 0.1 0.0 0.0 2.0")
    if not res[0]:
        return (False, res[1])
    else:
        return (True, None)

lateralMovementTest = TestCase("Lateral Movement", "Test to see if a body moves left to right", lateralMovement)