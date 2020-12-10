import sys
import time
import os
import subprocess
import shlex
import argparse


def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

def softEprint(*args, **kwargs):
    print(*args, file=sys.stderr, end="", flush=True, **kwargs)

def softPrint(*args, **kwargs):
    print(*args, end="", flush=True, **kwargs)

def overPrint(*args, **kwargs):
    print("\r", *args, **kwargs)


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

def runTest(args):
    print("Running Test")
    time.sleep(3)
    eprint("Test 1.... ❌")



if __name__ =="__main__":
    parser = argparse.ArgumentParser(description='Particles - Test Suite')
    parser.add_argument("file_name")
    args = parser.parse_args()

    print(args)

    print("----------------------------")
    print("|    Running Test Suite    |")
    print("----------------------------")
    print("Target File: %s"%str(args.file_name))

    if not compile(str(args.file_name)): 
        exit()
