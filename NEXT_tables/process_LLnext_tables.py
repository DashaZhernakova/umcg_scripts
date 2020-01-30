
# coding: utf-8

# In[4]:


from collections import defaultdict
import glob
import re


# In[164]:


f = open("/Users/dashazhernakova/Documents/UMCG/data/NEXT_tables/Registration_LL_NEXT_samples_NEXTnumberAdded_PartMerge_v4.0_comments.txt")


# In[165]:


out = open("/Users/dashazhernakova/Documents/UMCG/data/NEXT_tables/Registration_LL_NEXT_samples_NEXTnumberAdded_PartMerge_v4.0_comments_mvd.txt", "w")


# In[166]:


# first read the 4 header lines
spl_h1 = f.readline().rstrip('\r\n').split("\t")
spl_h2 = f.readline().rstrip('\r\n').split("\t")
spl_h3 = f.readline().rstrip('\r\n').split("\t")
spl_h4 = f.readline().rstrip('\r\n').split("\t")


# In[167]:


#
# Process the header, link LL ids to column numbers
#
colnames = ["Barcode tube", "Barcode swab tube", "Date withdrawal", "Time withdrawal", "Date processing", "Time processing", "comment"]
header_dict = defaultdict(list)
id_colns = []

comm_list = ['']*3

cur_coln = 3 #current LL id column number

coln = cur_coln + 1
while coln < len(spl_h1):
    col = spl_h1[coln]
    
    #print coln, col, spl_h2[coln], spl_h3[coln], spl_h4[coln]
    if col == "Lifelines + NEXT number":
        ll_id = spl_h2[cur_coln]
        
        #fill in the lld_id -> col names dict
        header_dict[ll_id].append(cur_coln)
        #print "\t", ll_id, cur_coln
        # keep the col number where new id starts
        id_colns.append(cur_coln)
        
        #write comments
        spl_h1[coln - 1] = comm_list[0].replace(";", "", 1)
        spl_h2[coln - 1] = comm_list[1].replace(";", "", 1)
        spl_h3[coln - 1] = comm_list[2].replace(";", "", 1) 
        comm_list = ['']*3

        
        # check that the number of fields per id are the same
        assert coln - cur_coln == 7
        
        cur_coln = coln
        
    #add all text to comments string
    elif (not spl_h4[coln] == "comment"):
        if len(spl_h1[coln]) > 1:
            comm_list[0] += ";" + spl_h1[coln]
            spl_h1[coln] = ""
        if len(spl_h2[coln]) > 1:
            comm_list[1] += ";" + spl_h2[coln]
            spl_h2[coln] = ""
        if len(spl_h3[coln]) > 1:
            comm_list[2] += ";" + spl_h3[coln]
            spl_h3[coln] = ""

    
    coln += 1

# Finish with the last id
ll_id = spl_h2[cur_coln]
header_dict[ll_id].append(cur_coln)
id_colns.append(cur_coln)

#write comments
spl_h1[coln - 1] = comm_list[0].replace(";", "", 1)
spl_h2[coln - 1] = comm_list[1].replace(";", "", 1)
spl_h3[coln - 1] = comm_list[2].replace(";", "", 1) 

#print "\t", ll_id, cur_coln
out.write("\t".join(spl_h1) + "\n")
out.write("\t".join(spl_h2) + "\n")
out.write("\t".join(spl_h3) + "\n")
out.write("\t".join(spl_h4) + "\n")
    


# In[168]:


# define how sample type corresponds to sample group
stype2group = {}
with open("/Users/dashazhernakova/Documents/UMCG/data/NEXT_tables/stype2sgroup.txt") as stype_f:
    stype2group = {l.rstrip().split("\t")[0] : l.rstrip().split("\t")[1] for l in stype_f}
#stype2group


# In[169]:


# process the main table
res_dict = defaultdict(dict)
cnt = 0
n_fields = 7
#cur_when = ''
#cur_who = ''
for l in f:
    cnt += 1
    spl = l.rstrip('\r\n').split("\t")
    if len(spl[0]) > 0:
        when = spl[0].replace('Birth','B')
    if len(spl[1]) > 0:
        who = spl[1][0]
    if len(spl[2]) > 0:
        what = spl[2]
    #print what
    what = re.sub(r" %s$" % when, "", what)
    sgroup = stype2group.get(what)
       
    # loop over all ids for this sample type
    for ll_id, colnums in header_dict.items():
        for coln in colnums:
            # TODO don't forget Bristol stool samples later
            
            #move comments to the last field
            spl[coln:coln + n_fields] = move_comments(spl[coln:coln + n_fields])

            #fill the dict, skip empty samples and irrelevant sample types
            res_list = spl[coln:coln + n_fields]
            
            if sgroup and any(res_list) and not res_list[-1] == 'o':
                #what_stripped = re.sub(r" %s$" % when, "", what)
                #res_list2 = res_list
                #print ll_id, who, when, what, res_list
                res_list2 = search_fill_missing_dates(ll_id, who, when, what, res_list, res_dict)
                
                spl[coln:coln + n_fields] = res_list2
                res_dict[sgroup][ll_id + ":" + who + ":" + when + ":" + what] = res_list2
                #if ll_id == '007274':
                #    print ll_id, who, when, what, res_list, spl[coln:coln + n_fields]
            #check colnames
            assert spl_h5[coln:coln + n_fields] == colnames, spl_h5[coln:coln + n_fields]     
    out.write("\t".join(spl) + "\n")
        
print "Processed lines:", cnt
f.close()
out.close()


# In[170]:


#
# process and fill storage files
#
st_header_pattern = ["Box number", "Position", "Lifelines number", "Lifelines Next number", "Time point", "Aliquot number", "Barcode", "Sample issued (yes/No)", "Date Sample issued ", "Remarks"]
storage_fdir = "/Users/dashazhernakova/Documents/UMCG/data/NEXT_tables/Fixed_tables/"
out_storage_fdir = "/Users/dashazhernakova/Documents/UMCG/data/NEXT_tables/Fixed_tables_DZ/"
for filepath in glob.iglob(storage_fdir + '* (edited Justin).txt'):
    fname = filepath.split("/")[-1]
    m = re.match(r"(.*)_([BMF]) \(edited Justin\)\.txt", fname)
    file_sgroup = m.group(1)
    file_who_abbr = m.group(1)
    print filepath, file_sgroup, file_who_abbr
    out_fpath = out_storage_fdir + fname.replace(" (edited Justin).txt", ".edited_DZ.txt", 1)

    st_f = open(filepath)
    out_st_f = open(out_fpath, 'w')
    header_spl = st_f.readline().rstrip().split('\t')

    #TODO: Comments also side by side
    out_st_f.write("\t".join(header_spl[:7]) + "\tReg_barcode\tReg_barcode_swab\tMerged barcode\tComment barcode\t" + "\t".join(colnames[2:-1]) + "\t" + "\t".join(header_spl[7:10]) + "\tComment reg file" + "\n")
    #out_st_f.write("\t".join(header_spl[:7]) + "\t" + "\t".join(colnames) + "\t" + header_spl[9] + "\n")
    assert header_spl == st_header_pattern, header_spl
    #cnt = 0

    for l in st_f:
        spl = l.rstrip('\r\n').split('\t')
        lst_to_write = [""] * 19
        lst_to_write[:7] = spl[:7]

        if spl[3] == 'null' or len(spl[3]) == 0:
            out_st_f.write("\t".join(lst_to_write) + "\n")
            continue

        ll_id = "%06d" % (int(spl[3]),)  
        res = res_dict[file_sgroup].get(ll_id + ":" + file_who + ":" + spl[4] + ":" + spl[5])

        if not res:
            lst_to_write[10] = "WARNING: no info in the registration file"
            if ll_id not in header_dict:
                lst_to_write[10] = "ERROR: no such id in the registration file"
            out_st_f.write("\t".join(lst_to_write) + "\n")
            continue

        #compare barcodes, write a merged column
        #error_barcode_mismatch
        #error_no_barcode_in_reg_file
        lst_to_write[7:9] = res[0:2]
        lst_to_write[9] = "NA"
        bc = spl[6]
        if len(res[1]) > 0: # if swap barcode present
            if len(bc) > 0 and not bc == res[1] and not bc.lower() == 'no barcode':
                    lst_to_write[10] = "ERROR: barcodes don't match"
                    #lst_to_write[9] = "NA"
            else:
                lst_to_write[9] = res[1]
        elif len(res[0]) > 0 and not res[0] == 'NA':
            if len(bc) > 0 and not bc == res[0] and not bc.lower() == 'no barcode':
                    lst_to_write[10] = "ERROR: barcodes don't match"
                    #lst_to_write[9] = "NA"
            else:
                lst_to_write[9] = res[0]
        elif len(bc) > 0 and not bc.lower() == 'no barcode':
            lst_to_write[9] = bc
            lst_to_write[10] = "WARNING: no barcode in registration file"

        lst_to_write[11:15] = res[2:-1] # time and date from reg file
        lst_to_write[15:18] = spl[7:10] # the rest of the storage file
        lst_to_write[18] = res[-1] # comments from reg file
        out_st_f.write("\t".join(lst_to_write) + "\n")
        #cnt += 1
        #if cnt > 10:
        #    break
    st_f.close()
    out_st_f.close()


# In[163]:


def move_comments(id_spl):
    exceptions = ['na', 'n/a', 'morning', 'evening', 'afternoon', 'no barcode', 'no  barcode']
    if id_spl[0].lower() == 'no barcode' or id_spl[0].lower() == 'no  barcode':
        id_spl[0] = 'NA'
        id_spl[-1] += ';' + 'no barcode'
    
    # unknown time
    if id_spl[3].lower() in ['unknown', 'no time']:
        id_spl[3] = 'NA'
        id_spl[-1] += ';' + 'Unknown time withdrawal'
    if id_spl[5].lower() in ['unknown', 'no time']:
        id_spl[5] = 'NA'
        id_spl[-1] += ';' + 'Unknown time processing'
    #unknown date
    if id_spl[2].lower() in ['unknown', 'no date']:
        id_spl[2] = 'NA'
        id_spl[-1] += ';' + 'Unknown date withdrawal'
    if id_spl[4].lower() in ['unknown', 'no date']:
        id_spl[4] = 'NA'
        id_spl[-1] += ';' + 'Unknown date processing'
    
    for i,val in enumerate(id_spl[:-1]):
        if val[:1].isalpha() and not val[:-1].isdigit() and not val.lower() in exceptions and not val.startswith('ID:') and not val.startswith('SF'):
            id_spl[-1] += ";" + val
            id_spl[i] = ""
            
    id_spl[-1] = id_spl[-1].replace(";", "", 1)
    return id_spl


# In[149]:


def search_fill_missing_dates(ll_id, who, when, what, res_list, res_dict):
    stype_list = []
    res_to_use = []
    #if not any(res_list[:6]):
    #    print "Missing values", ll_id + ":" + who + ":" + when, res_list
    #    return res_list
    if not (res_list[1] or res_list[0]) or res_list[0] == 'NA':
            return res_list
    if not all(res_list[2:6]):
        if what.startswith("Plasma-") and who == "Baby" and when.startswith('M'):
            res_to_use = res_dict['temp_group'].get(ll_id + ":" + who + ":" + when + ":" + 'EDTA 0.5mL-1')
            return fill_missing_dates(res_list, res_to_use)
        if what.startswith("Serum CA-") and who == "Baby" and when.startswith('M'):
            res_to_use = res_dict['temp_group'].get(ll_id + ":" + who + ":" + when + ":" + 'Serum CA 0.5mL-1')
            return fill_missing_dates(res_list, res_to_use)
        if what.startswith("Plasma-"):
            res_to_use = res_dict['temp_group'].get(ll_id + ":" + who + ":" + when + ":" + '1x K2EDTA 10 mL')
            return fill_missing_dates(res_list, res_to_use)
        if what.startswith("Serum SSTII-"):
            res_to_use = res_dict['temp_group'].get(ll_id + ":" + who + ":" + when + ":" + '1x Serum SSTII 8.5 mL')
            return fill_missing_dates(res_list, res_to_use)
        if what.startswith("Serum CA-"):
            res_to_use = res_dict['temp_group'].get(ll_id + ":" + who + ":" + when + ":" + '1x Serum CA 6 mL')
            return fill_missing_dates(res_list, res_to_use)
        if what.startswith("Citraat-"):
            res_to_use = res_dict['temp_group'].get(ll_id + ":" + who + ":" + when + ":" + '1x Citraat 6 mL')
            return fill_missing_dates(res_list, res_to_use)
        if what.startswith("Heparine-"):
            res_to_use = res_dict['temp_group'].get(ll_id + ":" + who + ":" + when + ":" + 'Heparine gel 0.6mL-1')
            return fill_missing_dates(res_list, res_to_use)
        if what.startswith("Urine-"):
            res_to_use = res_dict['temp_group'].get(ll_id + ":" + who + ":" + when + ":" + '1x Urine')
            return fill_missing_dates(res_list, res_to_use)
        if what.startswith("Placenta mother-"):
            res_to_use = res_dict['temp_group'].get(ll_id + ":" + who + ":" + when + ":" + '1x Placenta mother')
            return fill_missing_dates(res_list, res_to_use)
    return res_list
def fill_missing_dates(res_list, res_to_use):
    if not res_to_use:
        return res_list
    for i in range(2,6):
        if len(res_list[i]) == 0 and len(res_to_use[i]) > 0:
            res_list[i] = res_to_use[i]
    return res_list

