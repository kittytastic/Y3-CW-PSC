import matplotlib.pyplot as plt
import numpy as np

with open("./scaling.data", "r") as f:
    raw = f.readlines()
    raw = ''.join(raw)
    raw_runs = raw.split("----------\n")
    partial_runs = [r.split("\n") for r in raw_runs]
    partial_runs = partial_runs[1:]
    print(partial_runs)
    run_bc = [int(pr[0]) for pr in partial_runs]
    print(run_bc) 

    runs = []
    for pr in partial_runs:
        this_run = [] 
        for i in range(1, len(pr)-1):
            splt = pr[i].split()
            c = int(splt[0])
            t = float(splt[1])
            this_run.append((c,t))
        runs.append(this_run)
    
    print(runs)

    all_f = []
    for i in range(len(runs)):
        print("---- %d ----"%(run_bc[i]))
        fs = []
        ccs = []
        t1 = runs[i][0][1]

        for p, tp in runs[i][1:]:
            nume = tp - t1/p
            denom = t1 - t1/p
            f = nume/denom
            fs.append(f)
            ccs.append(p)
        
        plt.plot(ccs, fs, ".", color="C%d"%i)
        plt.plot(ccs, np.zeros(len(fs))+np.mean(fs), "--", color="C%d"%i)

        all_f.append(fs)

    plt.xlabel("cores")
    plt.ylabel("f (from strong scaling model)")
    plt.title("Analysis of f values assuming strong scaling")
    plt.savefig("f.png")
    
    print(all_f)