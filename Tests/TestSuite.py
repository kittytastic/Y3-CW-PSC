import sys
import time
import os
import subprocess
import shlex
import argparse

from Utils import *
from BasicTests import *
from CollisionTests import *
from ExperimentalTests import *
from HybridTests import *


testSuites = {
    "full": [MoveRightTest, ShortMoveDiagTest, StationaryColide, TwoBodyOrbitTest, StandingCubeColide],
    "live": [ShortMoveDiagTest, StationaryColide],
    "collide": [StandingCubeColide],
    "gen": [StandingCubeColide]
}

def compile(file_name, destination="test.out", compiler_call="g++ -O3"):
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


if __name__ =="__main__":

    parser = argparse.ArgumentParser(description='Particles - Test Suite')
    parser.add_argument("file_name")
    parser.add_argument("--type", dest="testType", help="Select a test type from: %s"%(",".join(testSuites.keys())),default="live" )
    parser.add_argument("-v", dest="verbose", help="Verbose", action='store_true')

    args = parser.parse_args()
    file_name = str(args.file_name)
    testType = str(args.testType)
    os.environ['VERBOSE'] = "T" if args.verbose else "F"
    bin_file_name = "test.out"

    print("----------------------------")
    print("|    Running Test Suite    |")
    print("----------------------------")
    print("Target File: %s"%str(args.file_name))

    if not compile(str(args.file_name), bin_file_name): 
        exit()

    if str(args.testType) not in testSuites:
        print("Error: No suite called \"%s\""%str(args.testType))
        exit()

    for test in testSuites[str(args.testType)]:
        test.run_test(bin_file_name)

    
