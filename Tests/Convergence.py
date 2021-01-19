import sys
import time
import os
import subprocess
import shlex
import argparse
import matplotlib.pyplot as plt

from Utils import *
from SpeedTests import randomBodies

bin_file_name = "test.out"

def compile(file_name, destination="test.out", compiler_call="g++ -O3 -fopenmp"):
    softPrint("Compiling...")
    success = True    
    try:
        invocation = compiler_call + " -o " + destination + " " + file_name
        processOut = subprocess.run(shlex.split(invocation), capture_output=True)
        if processOut.returncode!=0:
            print("  ❌")
            print("ERROR: Compilation failed with the following errors:")
            print()
            print(processOut.stderr.decode("utf-8"))
            return False
        else:
            print("  ✔️")
            return True
    except Exception as e:
        print("  ❌")
        print("ERROR: Compilation failed with the following python exception 🐍: ")
        print()
        print(str(e))
        return False

def runConvergenceTest(sim_args):
    args = sim_args.genArg()
    dt = sim_args.dt

    t0 = time.time()
    runtime_res = run_bin_with_args(bin_file_name, args)
    t1 = time.time()

    elapsed = t1-t0
    if not runtime_res[0]:
        print("ERROR: Test ran for %.3fs\n%s"%(dt,runtime_res[1]))
        raise Exception("There was a runtime error")
        
    # Get final v_max and dx_max from stdout 
    runtime_out = runtime_res[1].split("\n")
    #print(runtime_out)
    last_obj = [ float(i) for i in runtime_out[-2].split(":")[1].split(',')]
    final_out = runtime_out[-4]
    final_out = final_out.split(",")
    v_max = float(final_out[4].split("=")[1])
    dx_max = float(final_out[5].split("=")[1])
        
    results = parseResultFiles()
    results.setOutPart(v_max, dx_max) 

    return (last_obj, results, elapsed)

def CreateBodies(seed):
    snap_shot_t = 1
    final_t = 30
    dt = None

    init_bodies = randomBodies(50, 100, 10, 10, r_seed=seed)

    sim_setup = SimArgs(snap_shot_t, final_t, dt, init_bodies)
    

    

    return sim_setup

out_file_name = "Convergence.csv"

if __name__ =="__main__":

    parser = argparse.ArgumentParser(description='Particles - Test Suite')
    parser.add_argument("--intel", dest="intel", help="Use Intel compiler", action='store_true')

    args = parser.parse_args()
    

    print("-----------------------------")
    print("|    Running Convergence    |")
    print("-----------------------------")
    

    sf_name = ["../Solution/step-1.cpp", "../Solution/step-3.cpp"]
    sf_col = ["r", "g"]

    for sf in range(2):
        print()
        print("Target File: %s"%sf_name[sf])

        compiler_call = "g++ -O3 -fopenmp"
        if args.intel:
            compiler_call = "icpx -O3 -fopenmp --std=c++0x"

        if not compile(sf_name[sf], bin_file_name, compiler_call=compiler_call): 
            exit()

        for i in range(4):
            print("------ Seed %d ------"%(i))
            sim_args = CreateBodies(i)
            dt = 1
            k=0
            dt_min = 1e-4

            dt_plot = []
            diff_plot = []

            last_x, last_y, last_z = 0,0,0
            
            while dt >= dt_min:
                softPrint("dt: %.2e  "%dt)
                sim_args.dt = dt
                exact_obj, results, elapsed = runConvergenceTest(sim_args)
                assert(results.num_bodies>0)
                x,y,z = exact_obj
                
                dx = x - last_x
                dy = y - last_y
                dz = z - last_z

                mag = math.sqrt(dx*dx + dy*dy + dz*dz)
                
                dt_plot.append(k)
                diff_plot.append(mag)

                last_x, last_y, last_z = x,y,z

                print("  ✔️  (%.2f)  bodies: %d  mag: %.2e"%(elapsed, results.num_bodies, mag))
                dt = dt/2
                k+=1

            dt_plot.pop(0)
            diff_plot.pop(0)
            plt.plot(dt_plot, diff_plot, color=sf_col[sf])
    plt.xlabel('k where $h=1/2^k$')
    plt.ylabel('$f_h(T)$')
    plt.yscale('log')
    plt.title("Convergence") 

    plt.savefig('convergence.png')


    

    


    
