import sys

def calculate_FRS(age, sex, smoker, total_cholesterol, hdl_cholesterol, systolic_blood_pressure, blood_pressure_med_treatment, diabetes = False):
    points = 0
    percent_risk = 0
    #Process males -----------------------------------------------------------
    if sex == "2":
        # Age - male        
        if  30 <= age <= 34:
            points-=0
        if  35 <= age <= 39:
            points+=2
        if  40 <= age <= 44:
            points += 5
        if  45 <= age <= 49:
            points += 6
        if  50 <= age <= 54:
            points += 8
        if  55 <= age <= 59:
            points+= 10
        if  60 <= age <= 64:
            points+=11
        if  65 <= age <= 69:
            points+=12
        if  70 <= age <= 74:
            points+=14
        if  75 <= age:
            points+=15

        #Total cholesterol, mg/dL - Male ------------------------
        
        if total_cholesterol < 160:
            points +=0
        if 160 <= total_cholesterol <= 199:
            points+=1
        if 200 <= total_cholesterol <= 239:
            points+=2
        if 240 <= total_cholesterol <= 279:
            points+=3
        if total_cholesterol > 280:
            points+=4
        
        #smoking - male
        if smoker:
            points+=4 
        else: # nonsmoker
            points += 0
            
        #hdl cholesterol
        if hdl_cholesterol > 60:
            points-=2
        if 50 <= hdl_cholesterol <= 59:
            points-=1
        if 45 <= hdl_cholesterol <= 49:
            points+=0
        if 35 <= hdl_cholesterol <= 44:
            points+=1
        if hdl_cholesterol < 35:
            points+=2
            
        #systolic blood pressure
        if not blood_pressure_med_treatment:
            if systolic_blood_pressure < 120:
                points -= 2
            if 120 <= systolic_blood_pressure <= 129:
                points+=0
            if 130 <= systolic_blood_pressure <= 139:
                points+=1           
            if 140 <= systolic_blood_pressure <= 159:
                points+=2
            if systolic_blood_pressure >= 160:
                points +=3
        else: #if the patient is on blood pressure meds
            if systolic_blood_pressure < 120:
                points+=0
            if 120 <= systolic_blood_pressure <= 129:
                points+=1
            if 130 <= systolic_blood_pressure <= 139:
                points+=1           
            if 140 <= systolic_blood_pressure <= 159:
                points+=2
            if systolic_blood_pressure >= 160:
                points +=3
        if diabetes:
            points += 3
        
        #calulate % risk for males
        if points <= -3:
            percent_risk = 1
        elif points == -2:
            percent_risk = 1.1
        elif points == -1:
            percent_risk = 1.4
        elif points == 0:
            percent_risk = 1.6
        elif points == 1:
            percent_risk = 1.9
        elif points == 2:
            percent_risk = 2.3
        elif points == 3:
            percent_risk = 2.8
        elif points == 4:
            percent_risk = 3.3
        elif points == 5:
            percent_risk = 3.9
        elif points == 6:
            percent_risk = 4.7
        elif points == 7:
            percent_risk = 5.6
        elif points == 8:
            percent_risk = 6.7
        elif points == 9:
            percent_risk = 7.9
        elif points == 10:
            percent_risk = 9.4
        elif points == 11:
            percent_risk = 11.2
        elif points == 12:
            percent_risk = 13.2
        elif points == 13:
            percent_risk = 15.6
        elif points == 14:
            percent_risk = 18.4
        elif points == 15:
            percent_risk = 21.6
        elif points == 16:
            percent_risk = 25.3
        elif points == 17:
            percent_risk = 29.4
        elif points >= 18:
            percent_risk = 30
            
    #process females ----------------------------------------------------------
    elif sex == "1":
        # Age - female        
        if  30 <= age <= 34:
            points-=0
        if  35 <= age <= 39:
            points-=2
        if  40 <= age <= 44:
            points+=4
        if  45 <= age <= 49:
            points+=5
        if  50 <= age <= 54:
            points+=7
        if  55 <= age <= 59:
            points+=8
        if  60 <= age <= 64:
            points+=9
        if  65 <= age <= 69:
            points+=10
        if  70 <= age <= 74:
            points+=11
        if  75 <= age:
            points+=12

        #Total cholesterol, mg/dL - Female ------------------------
        
        if total_cholesterol < 160:
            points +=0
        if 160 <= total_cholesterol <= 199:
            points+=1
        if 200 <= total_cholesterol <= 239:
            points+=3
        if 240 <= total_cholesterol <= 279:
            points+=4
        if total_cholesterol > 280:
            points+=5
        
        #smoking - female
        if smoker == '1':
            points+=3  
        else: #nonsmoker
            points += 0
            
        #hdl cholesterol - female
        if hdl_cholesterol > 60:
            points-=2
        if 50 <= hdl_cholesterol <= 59:
            points-=1
        if 45 <= hdl_cholesterol <= 49:
            points+=0
        if 35 <= hdl_cholesterol <= 44:
            points+=1
        if hdl_cholesterol < 35:
            points+=2
            
        #systolic blood pressure
        if not blood_pressure_med_treatment: #untreated
            if systolic_blood_pressure < 120:
                points-=3
            if 120 <= systolic_blood_pressure <= 129:
                points+=0
            if 130 <= systolic_blood_pressure <= 139:
                points+=1           
            if 140 <= systolic_blood_pressure <= 149:
                points+=2
            if 150 <= systolic_blood_pressure <= 159:
                points+=4
            if systolic_blood_pressure >= 160:
                points +=5
        else: #if the patient is on blood pressure meds
            if systolic_blood_pressure < 120:
                points-=1
            if 120 <= systolic_blood_pressure <= 129:
                points+=2
            if 130 <= systolic_blood_pressure <= 139:
                points+=3           
            if 140 <= systolic_blood_pressure <= 149:
                points+=5
            if 150 <= systolic_blood_pressure <= 159:
                points+=6
            if systolic_blood_pressure >= 160:
                points +=7
        
        #calulate % risk for females
        if points <= -2:
            percent_risk = 1
        elif points == -1:
            percent_risk = 1.0
        elif points == 0:
            percent_risk = 1.2
        elif points == 1:
            percent_risk = 1.5
        elif points == 2:
            percent_risk = 1.7
        elif points == 3:
            percent_risk = 2.0
        elif points == 4:
            percent_risk = 2.4
        elif points == 5:
            percent_risk = 2.8
        elif points == 6:
            percent_risk = 3.3
        elif points == 7:
            percent_risk = 3.9
        elif points == 8:
            percent_risk = 4.5
        elif points == 9:
            percent_risk = 5.3
        elif points == 10:
            percent_risk = 6.3
        elif points == 11:
            percent_risk = 7.3
        elif points == 12:
            percent_risk = 8.6
        elif points == 13:
            percent_risk = 10.0
        elif points == 14:
            percent_risk = 11.7
        elif points == 15:
            percent_risk = 13.7
        elif points == 16:
            percent_risk = 15.9
        elif points == 17:
            percent_risk = 18.5
        elif points == 18:
            percent_risk = 21.5
        elif points == 19:
            percent_risk = 24.8
        elif points == 20:
            percent_risk = 28.5
        elif points >= 21:
            percent_risk = 30
    return percent_risk


agesex_fname = sys.argv[1]
pheno_fname = sys.argv[2]
antihyper_fname = sys.argv[3]

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

with open(antihyper_fname) as f:
    f.readline()
    
    for l in f:
        antihyper_users.add(l.strip())
print ("id\tFRS")
with open(agesex_fname) as f:
    header =  {col: col_num for col_num, col in enumerate(f.readline().rstrip().split("\t"))}
    age_col = header['age']
    sex_col = header['gender_F1M2']
    smk_col = header['SMK3']
    for l in f:
        spl = l.rstrip().split("\t")
        if spl[0] in antihyper_users:
            antihyper = True
        else:
            antihyper = False
        pheno = pheno_dict.get(spl[0])
        if pheno:
            frs = calculate_FRS(int(spl[age_col]), spl[sex_col], spl[smk_col], pheno[0], pheno[1], pheno[2], antihyper)
        else:
            frs = "NA"
        
        print (spl[0] + "\t" + str(frs))