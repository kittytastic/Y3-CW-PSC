import subprocess
import shlex

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


class TestCase():
    def __init__(self, name, desc, test_func):
        self.name = name
        self.description = desc
        self.test_func = test_func

    def run_test(self, target_bin):
        self.print_prefix()
        res = self.test_func(target_bin)
        if res[0]:
            if res[1]:
                self.print_success_verbose(res[1])
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
    
    def print_success_verbose(self):
        print("   ✔️")

    def print_error(self, error_msg):
        print("   ❌")
        print("Error: Test Failed")
        print("Test Description: %s"%self.description)
        print()
        print(error_msg)
    