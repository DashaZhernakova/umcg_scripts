import sys

def calculate_FRS(age, sex, smk, cho, hdc, sbp, bp_med, diab):
    
    cho = 38.67 * cho
    hdc = 38.67 * hdc
    
    if sex == 1:
        if not bp_med:
            sbp_beta = 2.76157
        else:
            sbp_beta = 2.82263
        risk_factors = (math.log(age) * 2.32888) + (math.log(cho) * 1.20904) + (math.log(hdc) * -0.70833) + (math.log(sbp) * sbp_beta) + smk * 0.52873 + 0.69154 * diab - 26.1931
        risk = 100* (1 - pow(0.95012, math.exp(risk_factors)))
    else: # men
        if not bp_med:
            sbp_beta = 1.93303
        else:
            sbp_beta = 1.99881
        risk_factors = (math.log(age) * 3.06117) + (math.log(cho) * 1.12370) + (math.log(hdc) * -0.93263) + (math.log(sbp) * sbp_beta) + smk * 0.65451 + 0.57367 * diab - 23.9802
        risk = 100* (1 - pow(0.88936, math.exp(risk_factors)))
    #print (risk_factors)
    return(risk)


agesex_fname = sys.argv[1]
pheno_fname = sys.argv[2]
antihyper_fname = sys.argv[3]
diab_fname = sys.argv[4]

cnt = 1
pheno_dict = {}
with open(pheno_fname) as f:
    header =  {col: col_num for col_num, col in enumerate(f.readline().rstrip().split("\t"))}
    cho_col = header['CHO']
    hdl_col = header['HDC']
    sbp_col = header['SBP']
    for l in f:
        spl = l.rstrip().split("\t")
        if not (spl[cho_col] == "NA" or spl[hdl_col] == "NA" or spl[sbp_col] == "NA"):
            pheno_dict[spl[0]] = [float(spl[cho_col]), float(spl[hdl_col]), float(spl[sbp_col])]


antihyper_users = set()
diab_samples = set()
with open(antihyper_fname) as f: 
    for l in f:
        antihyper_users.add(l.strip())

with open(diab_fname) as f:
    f.readline()   
    for l in f:
        diab_samples.add(l.strip())

print ("id\tFRS")
with open(agesex_fname) as f:
    header =  {col: col_num for col_num, col in enumerate(f.readline().rstrip().split("\t"))}
    age_col = header['age']
    sex_col = header['gender_F1M2']
    smk_col = header['SMK3']
    for l in f:
        spl = l.rstrip().split("\t")
        
        antihyper = True if spl[0] in antihyper_users else antihyper = False
        diab = True if spl[0] in diab_samples else diab = False
                
        pheno = pheno_dict.get(spl[0])
        if pheno:
            frs = calculate_FRS(int(spl[age_col]), spl[sex_col], spl[smk_col], pheno[0], pheno[1], pheno[2], antihyper, diab)
        else:
            frs = "NA"
        
        print (spl[0] + "\t" + str(frs))