import subprocess
import shlex
import xml.etree.ElementTree as ET
import math
import time
import os

def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

def softEprint(*args, **kwargs):
    print(*args, file=sys.stderr, end="", flush=True, **kwargs)

def softPrint(*args, **kwargs):
    print(*args, end="", flush=True, **kwargs)

def overPrint(*args, **kwargs):
    print("\r", *args, **kwargs)

def run_bin_with_args(binary, args):
    try:
        processOut = subprocess.run(["./"+binary, *shlex.split(args)], capture_output=True)
        if processOut.returncode!=0:
            return (False, processOut.stderr.decode("utf-8"))
        else:
            return (True, processOut.stdout.decode("utf-8"))
    except Exception as e:
        return (False, str(e))


class TestEqError(Exception):
    pass


class GeneralTestCase():
    def __init__(self, name, desc, test_func):
        self.name = name
        self.description = desc
        self.test_func = test_func

    def run_test(self, target_bin):
        self.print_prefix()
        res = self.test_func(target_bin)
        if res[0]:
            if res[1]:
                self.print_success_custom(res[1])
            else:
                self.print_success()
            return True
        else:
            self.print_error(res[1])
            return False


    def print_prefix(self):
        softPrint("Running %s test... "% self.name)

    def print_success(self):
        print("   ✔️")
    
    def print_success_custom(self, msg):
        print(msg)

    def print_error(self, error_msg):
        print("   ❌")
        print("Error: Test Failed")
        print("Test Description: %s"%self.description)
        print()
        print(error_msg)


class TestCase(GeneralTestCase):
    def __init__(self, name, desc, test_description):
        super(TestCase, self).__init__(name, desc, self.test)
        desc = test_description()
        self.setup = desc["setup"]
        self.expected  = desc["expected"]
        self.error = desc["error"] if "error" in desc else 1e-9

    def test(self, bin_name):
        
        args = self.setup.genArg()
        
        if os.environ.get("VERBOSE")=="T":
            print()
            print(args)    

        t0 = time.time()
        runtime_res = run_bin_with_args(bin_name, args)
        t1 = time.time()

        # Would be nice to pipe output instread of output at end... 
        verbose_adjustment = ""
        if os.environ.get("VERBOSE")=="T":
            verbose_adjustment = runtime_res[1]
            verbose_adjustment += "\n\n[TEST MESSAGE]: "

        elapsed = t1-t0
        if not runtime_res[0]:
            return (False, "%s Test ran for %.3fs\n%s"%(verbose_adjustment, elapsed, runtime_res[1]))
        
        # Get final v_max and dx_max from stdout 
        runtime_out = runtime_res[1].split("\n")
        final_out = runtime_out[-4]
        final_out = final_out.split(",")
        v_max = float(final_out[4].split("=")[1])
        dx_max = float(final_out[5].split("=")[1])
        
        results = parseResultFiles()
        results.setOutPart(v_max, dx_max) 

        try: 
            AssertSimSolutionEq(self.expected, results, self.error)
        except TestEqError as e:
            return (False, "%s(%.2fs) %s"%(verbose_adjustment, elapsed, str(e)))
        
        return (True, "%s   ✔️   (%.2fs)"%(verbose_adjustment, elapsed))

class SpeedTestCase(GeneralTestCase):
    def __init__(self, name, desc, test_description):
        super(SpeedTestCase, self).__init__(name, desc, self.test)
        desc = test_description()
        self.setup = desc["setup"]
        self.cluster = "cluster" in desc
        

    def test(self, bin_name):
        
        args = self.setup.genArg()
        
        if os.environ.get("VERBOSE")=="T":
            print()
            print(args)
        

        t0 = time.time()
        runtime_res = run_bin_with_args(bin_name, args)
        t1 = time.time()

        elapsed = t1-t0

        if not runtime_res[0]:
            return (False, "Test ran for %.3fs\n%s"%(elapsed, runtime_res[1]))        

        results = parseResultFiles()

        requirement = False
        if self.cluster:
            requirement = len(self.setup.bodies) < results.num_bodies * 0.5
        else:
            requirement = results.num_bodies < len(self.setup.bodies) * 0.8

        if requirement:
            return (False, "(%.2fs) Speed test collapsed to %d bodies"%(elapsed, results.num_bodies))
    
        return (True, "   ✔️   (%.2fs)"%(elapsed))
        


#####################################
#              Bodies               #
#####################################

def smart_is_close(a, b, tol=1e-9):
    if math.isclose(a,0, abs_tol=1e-9) or math.isclose(b,0, abs_tol=1e-9):
        # If we are comparing to 0
        return math.isclose(a,b, abs_tol=tol)
    else:
        return math.isclose(a,b,rel_tol=tol)

def turp_is_approx(turp_1, turp_2, allowed_error):
    if len(turp_1) != len(turp_2):
        return False
    
    for i in range(len(turp_1)):
        if not smart_is_close(float(turp_1[i]), float(turp_2[i]), tol=allowed_error):
            return False

    return True

def AssertPointArraysEqual(expected, observed, allowed_error):
    if len(expected) != len(observed):
        raise TestEqError("Point arrays have dirrent length; Expected: %d  Observed: %d"%(len(expected), len(observed)))

    matched = [False] * len(observed)

    for i in range(len(expected)):
        has_partner = False
        j = 0
        while j < len(observed) and not has_partner:
            if not matched[j] and turp_is_approx(expected[i], observed[j], allowed_error):
                matched[j] = True
                has_partner = True
            
            j+=1
        
        if not has_partner:
            raise TestEqError("Expected point is unmatched: %s"%str(expected[i]))
    
    return True


def AssertSimSolutionEq(expected, observed, allowed_error):
        if expected.num_bodies != observed.num_bodies:
            raise TestEqError("num_bodies != num_bodies; Expected: %d  Observed:  %d"%(expected.num_bodies, observed.num_bodies))

        AssertPointArraysEqual(expected.points, observed.points, allowed_error)

        if not smart_is_close(expected.v_max, observed.v_max, tol=allowed_error):
            raise TestEqError("v_max != v_max; Expected: %s  Observed:  %s"%(expected.v_max, observed.v_max))

        if not smart_is_close(expected.dx_min, observed.dx_min, tol=allowed_error):
            raise TestEqError("dx_min != dx_min; Expected: %s  Observed:  %s"%(expected.dx_min, observed.dx_min))
        
        return True 


class SimSolution():
    def __init__(self, num_bodies, points, steps, v_max = None, dx_min = None):
        self.num_bodies = num_bodies
        self.points = points
        self.steps = steps
        self.v_max = v_max
        self.dx_min = dx_min

    def setOutPart(self, v_max, dx_max):
        self.v_max = v_max
        self.dx_min = dx_max

    def __str__(self):
        outS = "------ SimSolution ------\n"
        outS += "Time Steps: %s\n"%self.steps
        outS += "Num Bodies: %s\n"%self.num_bodies
        outS += "Bodies: %s"%self.points
        return outS


class Body():
    def __init__(self, pos_turp, v_turp, mass):
        if len(pos_turp)!=3:
            raise Error("A body's position should be a 3-tuple")

        if len(pos_turp)!=3:
            raise Error("A body's velocity should be a 3-tuple")

        self.pos = pos_turp
        self.v = v_turp
        self.mass = mass

    def asArg(self):
        return "{self.pos[0]:f} {self.pos[1]:f} {self.pos[2]:f} {self.v[0]:f} {self.v[1]:f} {self.v[2]:f} {self.mass:f}".format(self=self)


class SimArgs():
    def __init__(self, snap_frequency, final_time, dt, bodies):
        self.snap_frequency = snap_frequency
        self.final_time = final_time
        self.dt = dt
        self.bodies = bodies

    def genArg(self):
        bodies_arg = " ".join([body.asArg() for body in self.bodies])
        return "{self.snap_frequency} {self.final_time} {self.dt} {bodies_arg}".format(self=self, bodies_arg=bodies_arg)



#####################################
#           Result Parsing          #
#####################################
def parseResultFiles():
    # Parse PVD
    root = ET.parse("result.pvd").getroot()
    *_, last_step = root.iter('DataSet')
   
    num_steps = int(last_step.get("timestep"))
    last_file = last_step.get("file")
    
    # Parse last VTP
    root = ET.parse(last_file).getroot()

    piece = next(root.iter('Piece'))
    num_bodies = int(piece.get("NumberOfPoints"))
   
    dataArray = next(root.iter('DataArray'))
    all_points = dataArray.text.split()
    all_points = [float(x) for x in all_points]
    points = [(all_points[i], all_points[i+1], all_points[i+2]) for i in range(0, len(all_points), 3)]
   

    return SimSolution(num_bodies, points, num_steps)
