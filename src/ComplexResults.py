import os
tests = ["compcld316k", "hook1498k", "pflow742k", "fault742k", "ice245k"]
bash_command = "lrun -n 1 ./test/ij -fromfile test/coofile.coo -solver 0 -max_iter 10000 > out.txt"

for test in tests:

   bc2 = bash_command.replace("coofile", test).replace("-solver 0", "-solver 31")

   print(bc2)
   Iters = ""
   WCTSetup = ""
   WCTSolve = ""
   Density = ""
   there = False
   there2 = False
   os.system(bc2)
   f = open("out.txt", "r")
   lines = f.readlines()
   f.close()
   for line in lines:
      if(line.find("PCG Setup") != -1):
         there = True
      elif(line.find("PCG Solve") != -1):
         there2 = True
      elif(line.find("wall clock time") != -1 and there == True and there2 == False):
         temp = line.split()
         WCTSetup = temp[4]
      elif(line.find("wall clock time") != -1 and there == True and there2 == True):
         temp = line.split()
         WCTSolve = temp[4]
      elif(line.find("Iterations") != -1):
         temp = line.split()
         Iters = temp[2]
      elif(line.find("Prec. density:") != -1):
         temp = line.split()
         Density = temp[3]
   print(Iters + ", " + WCTSetup + ", " + WCTSolve + ", " + Density)
 
   bc2 = bc2.replace("31", "8")     
   print(bc2)
   Iters = ""
   WCTSetup = ""
   WCTSolve = ""
   Density = ""
   there = False
   there2 = False
   os.system(bc2)
   f = open("out.txt", "r")
   lines = f.readlines()
   f.close()
   for line in lines:
      if(line.find("PCG Setup") != -1):
         there = True
      elif(line.find("PCG Solve") != -1):
         there2 = True
      elif(line.find("wall clock time") != -1 and there == True and there2 == False):
         temp = line.split()
         WCTSetup = temp[4]
      elif(line.find("wall clock time") != -1 and there == True and there2 == True):
         temp = line.split()
         WCTSolve = temp[4]
      elif(line.find("Iterations") != -1):
         temp = line.split()
         Iters = temp[2]
      elif(line.find("Prec. density:") != -1):
         temp = line.split()
         Density = temp[3]
   print(Iters + ", " + WCTSetup + ", " + WCTSolve + ", " + Density)
  
