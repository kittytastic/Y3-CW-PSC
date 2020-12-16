import subprocess
import shlex
import xml.etree.ElementTree as ET
import math
import time

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
        self.setup, self.expected  = test_description()

    def test(self, bin_name):
        args = self.setup.genArg()

        t0 = time.time()
        runtime_res = run_bin_with_args(bin_name, args)
        t1 = time.time()

        elapsed = t1-t0
        if not runtime_res[0]:
            return (False, ("Test ran for %.3fs\n"%elapsed)+res[1])
        

        file_out_res = parseResultFiles() 
        if not self.expected == file_out_res:
            return (False, "Paraview files output didn't match expected result")
        
        return (True, "   ✔️   (%.2fs)"%elapsed)
        


#####################################
#              Bodies               #
#####################################
def turp_is_approx(turp_1, turp_2):
    if len(turp_1) != len(turp_2):
        return False
    
    for i in range(len(turp_1)):
        if not math.isclose(float(turp_1[i]), float(turp_2[i])):
            return False

    return True

def point_arrays_equal(points_1, points_2):
    if len(points_1) != len(points_2):
        return False

    matched = [False] * len(points_2)

    for i in range(len(points_1)):
        has_partner = False
        j = 0
        while j < len(points_2) and not has_partner:
            if not matched[j] and turp_is_approx(points_1[i], points_2[j]):
                matched[j] = True
                has_partner = True
            
            j+=1
        
        if not has_partner:
            return False
    
    return True

class SimSolution():
    def __init__(self, num_bodies, points, steps):
        self.num_bodies = num_bodies
        self.points = points
        self.steps = steps

    def __str__(self):
        outS = "------ SimSolution ------\n"
        outS += "Time Steps: %s\n"%self.steps
        outS += "Num Bodies: %s\n"%self.num_bodies
        outS += "Bodies: %s"%self.points
        return outS

    def __eq__(self, other):
        if self.num_bodies != other.num_bodies:
            return False
        
        return point_arrays_equal(self.points, other.points)

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
