f1 = open("BLBLNew.mol", "r")
f2 = open("vesta.xyz", "w")


for line in f1:
   a = line.split()
   if (len(a) == 10):
      
#      f2.write("C {} {} {} \n".format(a[4], a[5], a[6]))
      if int(a[2]) == 1: # Hydrogen  
          f2.write("H {} {} {} \n".format(a[4], a[5], a[6]))
#      if int(a[2]) == 1: # Carbon 
#          f2.write("C {} {} {} \n".format(a[4], a[5], a[6]))
      if int(a[2]) == 2: # Boron 
          f2.write("C {} {} {} \n".format(a[4], a[5], a[6]))
#      elif int(a[2]) == 3: # Nitrogen 
#          f2.write("N {} {} {} \n".format(a[4], a[5], a[6]))
      #elif int(a[2]) == 4: # Hydrogen 
      #    f2.write("H {} {} {} \n".format(a[4], a[5], a[6]))

   elif (len(a)==2) and (a[1]=="atoms"): 
      f2.write(f"{a[0]}\n\n")

f2.close()

