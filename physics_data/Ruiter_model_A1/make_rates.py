import numpy as np


sntypes = ['DDS','Sch','SDS','HeRS']

data_tables = {}
for sn in sntypes:

    data_tables[sn] = np.genfromtxt(sn + '.dat')


# now interpolate each and write out points

points = np.array(list(np.arange(0.0,2.07,0.1)) + list(np.arange(2.0,14.1,0.5)))
rtotal = np.zeros(np.size(points))
for i,p in enumerate(points):
    for sn in sntypes:
        rtotal[i] += 10.0**(np.interp(p, data_tables[sn][:,0], data_tables[sn][:,1]))

    rtotal[i] = np.log10(rtotal[i])

print(np.size(points))

print("{",end='')
for p in points:
    print(" %3.1f,"%(p),end='')
print("};")

print("{",end='')
for p in points:
    r = np.interp(p, points, rtotal)
    print(" %3.3f,"%(r),end='')
print("},")

for sn in sntypes:
    print(sn)

    print("{",end='')
    for p in points:

        r = np.interp(p, data_tables[sn][:,0], data_tables[sn][:,1])

        print(" %3.3f,"%(r),end='')

    print("},")
