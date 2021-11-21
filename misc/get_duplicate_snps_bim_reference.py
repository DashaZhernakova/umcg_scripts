import sys
from collections import defaultdict


bim_f = open(sys.argv[1])
ref_bim_f  = open(sys.argv[2])

ref_dict = {}
for l in ref_bim_f:
    spl = l.strip().split()
    ref_dict[spl[1]] = set([spl[4], spl[5]])
sys.stderr.write("Read reference alleles\n")

duplicates = defaultdict(list)
visited = set()
cur_spl = bim_f.readline().strip().split()

#for l in bim_f:
while True:
    l = bim_f.readline()
    if not l:
        print(cur_spl[1] + "_" + cur_spl[4] + "_" + cur_spl[5])
        break
    spl = l.strip().split()
    snp = spl[1]
    if snp == cur_spl[1]:
        sys.stderr.write("Duplicate! " + snp + "\n")
        if snp in ref_dict:
            ref_alls = ref_dict[snp]
            
            if  set([cur_spl[4], cur_spl[5]]) == ref_alls:
                print(cur_spl[1] + "_" + cur_spl[4] + "_" + cur_spl[5])
                l = bim_f.readline()
                cur_spl = l.strip().split()
            elif set([spl[4], spl[5]]) == ref_alls:
                print(spl[1] + "_" + spl[4] + "_" + spl[5])
                l = bim_f.readline()
                cur_spl = l.strip().split()
            else:
                print(cur_spl[1] + "_" + cur_spl[4] + "_" + cur_spl[5])
                l = bim_f.readline()
                cur_spl = l.strip().split()


        else:
            print(cur_spl[1] + "_" + cur_spl[4] + "_" + cur_spl[5])
            l = bim_f.readline()
            cur_spl = l.strip().split()

    else:
        print(cur_spl[1] + "_" + cur_spl[4] + "_" + cur_spl[5])
        cur_spl = spl
