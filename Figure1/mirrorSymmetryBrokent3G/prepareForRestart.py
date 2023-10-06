f1 = open("BLBLNew.mol", "r")
f2 = open("generate.xyz", "r")
g = open("BL.mol.restart", "w")

for i in range(4):
   f2.readline()

for line in f1:
   a = line.split()
   if (len(a) == 10):
      b = f2.readline().split()
      g.write("{} {} {} {} {} {} {} {} {} {}\n".format(a[0], a[1], a[2], a[3], b[1], b[2], b[3], a[7], a[8], a[9]))
   else:
      g.write(line)
