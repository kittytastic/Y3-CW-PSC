import sys
import subprocess
import os
import re
import time



class Test(object):
  """
   Global database of all of the answers I did/want give students
  """
  
  
  """
  This is the scaling I use whenever I apply a timeout, i.e. if the student's
  solution needs more than TimeoutScalingFactor the time of my simple code, 
  then I consider the submission to be wrong.
  """
  TimeoutScalingFactor = 4
    
  def __init__(self, submitted_file, benchmark_file = ""):
    """
     submitted_file: String
       File submitted by student
    
     benchmark_file: String
       File to benchmark/validate against. Can be empty if no such thing is 
       available or shall be provided. 
    """
    self.submitted_file                         = submitted_file
    self.benchmark_file                         = benchmark_file
    
    self.compilation_success                    = True 
    self.yields_correct_result                  = True
    
    self.last_output                            = ""
    
    self.timeout                                = 4000
    
    self.runtime                                = 0
    pass


  def ignore_wrong_output(self):
    if self.compilation_success:
      self.yields_correct_result                = True


  def __get_submitted_executable_name(self):
    return "submitted.out"


  def __get_benchmark_executable_name(self):
    return "benchmark.out"


  def compile(self,compiler_call="g++ -O3"):
    """
    
      Compiles the code (and also compiles the solution). The routine needs the 
      report file to write errors to, and it wants to know what the solution file
      is called. That's the file submitted by the student.
      
      Returns 1 if successful. Otherwise 0
      
    """
    success = True    
    try:
      invocation = compiler_call + " -o " + self.__get_submitted_executable_name() + " " + self.submitted_file
      returnCode = subprocess.run(invocation, shell=True)
      if returnCode.returncode!=0:
        self.last_output = "ERROR: tried to translate with " + invocation + " and got " + str(returnCode.stderr)
        return 0
      else:
        return 1
    except Exception as e:
      self.last_output = "ERROR: wanted to translate but got an exception " + str(e)
      return 0
      

  def run(self,arguments="",environment_variables={}):
    """
    
      Compiles the code (and also compiles the solution). The routine needs the 
      report file to write errors to, and it wants to know what the solution file
      is called. That's the file submitted by the student.
      
      Returns 1 if successful. Otherwise 0
      
    """
    success = True    
    encoding = 'utf-8'
    try:
      invocation = "./" + self.__get_submitted_executable_name() + " " + arguments
      my_environment = os.environ
      for i in environment_variables:
        my_environment[ i[0] ] = i[1]

      start_solution = time.time()
      returnCode = subprocess.run(invocation, timeout=self.timeout, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True, env=my_environment )
      if returnCode.returncode!=0:
        self.last_output = "ERROR: tried to run code with " + arguments + " and got " + str(returnCode.stderr)
        return 0
      else:
        self.last_output = returnCode.stdout.decode(encoding).splitlines()[-1]
        self.runtime     = time.time() - start_solution
        return 1
    except Exception as e:
      self.last_output = "ERROR: wanted to run code but got an exception " + str(e)
      return 0

  
  def search_for_pattern_in_output(self,search_pattern):
    result = ""
    m = re.findall( search_pattern, self.last_output )
    if m:
      result = str(m[-1])
    return result
        

    
      

  #    if self.benchmark_file!="":
  #      invocation = compiler_call + " -o " + self.__get_benchmark_executable_name() + " " + self.benchmark_file
  #      returnCode = subprocess.run(invocation, shell=True)
      
      
    