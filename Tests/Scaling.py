import sys
import time
import os
import subprocess
import shlex
import argparse
import matplotlib.pyplot as plt

from Utils import *
from SpeedTests import randomBodies


def run_bin(binary, split_args, env_args):
    #print(["./"+binary, *split_args])
    my_env = os.environ.copy()
    processOut = subprocess.run(["./"+binary, *split_args], env={**my_env, **env_args}, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if processOut.returncode!=0:
        return (False, processOut.stderr.decode("utf-8"))
    else:
        return (True, processOut.stdout.decode("utf-8"))

   

def runConvergenceTest(sim_args, bin_name, core_count):
    args = sim_args.genArg()

    split_args = shlex.split(args)
    #print(split_args)
    t0 = time.time()
    runtime_res = run_bin(bin_name, split_args, {"OMP_NUM_THREADS":str(core_count)})
    t1 = time.time()

    elapsed = t1-t0
    if not runtime_res[0]:
        print("ERROR: Test ran for %.3fs\n%s"%(elapsed, runtime_res[1]))
        raise Exception("There was a runtime error")

    num_bodies = int(runtime_res[1].split('\n')[-3].split(":")[1])
    return (elapsed, num_bodies)

def CreateBodies(seed, scale):
    snap_shot_t = 0
    final_t = 30
    dt = 0.001

    init_bodies = randomBodies(int(500*scale), 100, 10, 10, r_seed=seed)

    sim_setup = SimArgs(snap_shot_t, final_t, dt, init_bodies)
    

    return sim_setup

def CreateBodies2(seed, scale):
    snap_shot_t = 0
    final_t = 3
    dt = 0.001

    init_bodies = randomBodies(int(1500*scale), 100, 10, 10, r_seed=seed)

    sim_setup = SimArgs(snap_shot_t, final_t, dt, init_bodies)
    

    return sim_setup

out_file_name = "Convergence.csv"

if __name__ =="__main__":

    parser = argparse.ArgumentParser(description='Test Suite - Scaling')
    parser.add_argument("file_name")
    parser.add_argument("--scaling", dest="scaling", help="Scaling value, default 1", default=1 )
    parser.add_argument("--max_core", dest="max_core", help="Max Cores", default=1, required=True )
    parser.add_argument("--core_step", dest="core_step", help="Core Step", default=1, required=True )
    

    args = parser.parse_args()
    file_name = str(args.file_name)
    max_cores = int(args.max_core)
    step = int(args.core_step)
    scaling = float(args.scaling)
    

    print("-----------------------------")
    print("|      Running Scaling      |")
    print("-----------------------------")
    
    
    print("Using Binary: %s"%file_name)
    seed_col = ['r', 'g'] 
    scenario = [CreateBodies2, CreateBodies]
    speed_1 = None

    core_options = [2**i for i in range(int(math.ceil(math.log(step, 2))))] + list(range(step, max_cores+1, step))
    print("Conducting %d tests per seed"%(len(core_options)))
    print("With cores counts of: %s"%core_options)
    
    start_time = time.time()
    for i in range(2):
        print("------ Seed %d ------"%(i))
        sim_args = scenario[i](i, scaling)

        time_plot = []
        cores_plot = []

        
        
        for cc in core_options:
            softPrint("Using %d cores... "%cc)
            elapsed, num_bodies = runConvergenceTest(sim_args, file_name, cc)

            if speed_1 == None:
                speed_1 = elapsed
            
            speed_up = speed_1 / elapsed
            time_plot.append(speed_up)
            cores_plot.append(cc)
            
            print("  DONE  (%.2f)  remaining bodies: %d"%(elapsed, num_bodies))
            

        plt.plot(cores_plot, time_plot, color=seed_col[i])
    end_time = time.time()
    total_time = end_time-start_time
    print("Total runtime %.2fs"%(total_time))
    
    plt.xlabel('cores')
    plt.ylabel('scaling')
    plt.title("Scaling") 

    plt.savefig('scaling.png')


    

    


    
