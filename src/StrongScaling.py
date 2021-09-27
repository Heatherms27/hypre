import os
normal = "-fs_max_steps 10 -fs_max_step_size 1 -fs_kap_tol 1e-3"
tests = ["-fs_max_steps 10 -fs_max_step_size 1 -fs_kap_tol 1e-1", "-fs_max_steps 3 -fs_max_step_size 3 -fs_kap_tol 1e-3"]
bash_command = "lrun -n 1 -N 1 ./test/ij -solver 31 -fs_max_steps 10 -fs_max_step_size 1 -fs_kap_tol 1e-3 -9pt -n 1 1 1 -P 1 1 1 -max_iter 5000 > out.txt"

for i in range(2):

   bash_command_t = bash_command.replace(normal, tests[i])

   print(bash_command_t)
   P1 = [1, 2, 2, 4, 4, 8, 8, 16, 16]
   P2 = [1, 1, 2, 2, 4, 4, 8, 8, 16]
   P3 = [1, 1, 1, 1, 1, 1, 1, 1, 1]
   N = [1, 1, 1, 1, 1, 1, 2, 4, 8]
   Iters = ""
   WCTSetup = ""
   WCTSolve = ""
   Density = ""
   print("-9pt")
   for j in range(9):
      j1 = 2**j
      bc2 = bash_command_t.replace("-n 1 -N 1", "-n " + str(j1) + " -N " + str(N[j])).replace("-n 1 1 1", "-n 1264 1264 1").replace("-P 1 1 1", "-P " + str(P1[j]) + " " + str(P2[j]) + " " + str(P3[j]))
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
 
   P1 = [1, 2, 2, 2, 4, 4, 4, 8, 8]
   P2 = [1, 1, 2, 2, 2, 4, 4, 4, 8]
   P3 = [1, 1, 1, 2, 2, 2, 4, 4, 4]
     
   print("-27pt")
      
   for j in range(9):
      j1 = 2**j
      bc2 = bash_command_t.replace("-n 1 -N 1", "-n " + str(j1) + " -N " + str(N[j])).replace("-n 1 1 1", "-n 256 256 256").replace("-P 1 1 1", "-P " + str(P1[j]) + " " + str(P2[j]) + " " + str(P3[j])).replace("9pt", "27pt")
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
    
   
