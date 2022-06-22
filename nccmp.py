#!/bin/env python
import sys

def nccmp(filename1,filename2):
    import os.path
    import netCDF4 as nc
    import hashlib as hsh

    try:
        ncid1 = nc.Dataset(filename1,"r")
    except RuntimeError:
        sys.stderr.write("%s does not exist.\n" % filename1)
        return 1

    try:
        ncid2 = nc.Dataset(filename2,"r")
    except RuntimeError:
        sys.stderr.write("%s does not exist.\n" % filename2)
        return 1
     
    vdiffs = []
    for v in ncid1.variables.iterkeys():
        md51 = hsh.md5(ncid1.variables[v][:]).hexdigest()
        md52 = hsh.md5(ncid2.variables[v][:]).hexdigest()
        if md51 != md52:
            vdiffs.append(v)
    if len(vdiffs) > 0:
        sys.stdout.write("Files %s and %s have differences in fields:\n" % (filename1,filename2))
        for v in vdiffs:
            sys.stdout.write("%s\n" % v)
    else:
        sys.stdout.write("No differences found in files %s and %s \n" % (filename1,filename2))

    return 0

if __name__ == "__main__":
    filename1 = sys.argv[1]
    filename2 = sys.argv[2]
    nccmp(filename1,filename2)

    


##
## End of script.
##
