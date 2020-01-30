
# coding: utf-8

# In[743]:


from collections import defaultdict
import glob
import re


# In[850]:


f = open("/Users/dashazhernakova/Documents/UMCG/data/NEXT_tables/Registration_LL_NEXT_samples_edSJ_NEXTnumbersAdded_comments.txt")


# In[851]:


out = open("/Users/dashazhernakova/Documents/UMCG/data/NEXT_tables/Registration_LL_NEXT_samples_edSJ_NEXTnumbersAdded_comments_mvd.txt", "w")


# In[852]:


# first read the 4 header lines
spl_h1 = f.readline().rstrip('\r\n').split("\t")
spl_h2 = f.readline().rstrip('\r\n').split("\t")
spl_h3 = f.readline().rstrip('\r\n').split("\t")
spl_h4 = f.readline().rstrip('\r\n').split("\t")
spl_h5 = f.readline().rstrip('\r\n').split("\t")


# In[853]:


#
# Process the header, link LL ids to column numbers
#
colnames = ["Barcode tube", "Barcode swab tube", "Date withdrawal", "Time withdrawal", "Date processing", "Time processing", "comment"]
header_dict = defaultdict(list)
id_colns = []

comm_list = ['']*4

cur_coln = 3 #current LL id column number

coln = cur_coln + 1
while coln < len(spl_h1):
    col = spl_h1[coln]
    
    #print coln, col, spl_h2[coln], spl_h3[coln], spl_h4[coln],spl_h5[coln]
    if col == "Lifelines + NEXT number":
        ll_id = spl_h4[cur_coln]
        
        #fill in the lld_id -> col names dict
        header_dict[ll_id].append(cur_coln)
        #print "\t", ll_id, cur_coln
        # keep the col number where new id starts
        id_colns.append(cur_coln)
        
        #write comments
        spl_h1[coln - 1] = comm_list[0].replace(";", "", 1)
        spl_h2[coln - 1] = comm_list[1].replace(";", "", 1)
        spl_h3[coln - 1] = comm_list[2].replace(";", "", 1) 
        spl_h4[coln - 1] = comm_list[3].replace(";", "", 1) 
        comm_list = ['']*4

        
        # check that the number of fields per id are the same
        assert coln - cur_coln == 7
        
        cur_coln = coln
        
    #add all text to comments string
    elif (not spl_h5[coln] == "comment"):
        if len(spl_h1[coln]) > 1:
            comm_list[0] += ";" + spl_h1[coln]
            spl_h1[coln] = ""
        if len(spl_h2[coln]) > 1:
            comm_list[1] += ";" + spl_h2[coln]
            spl_h2[coln] = ""
        if len(spl_h3[coln]) > 1:
            comm_list[2] += ";" + spl_h3[coln]
            spl_h3[coln] = ""
        if len(spl_h4[coln]) > 1:
            comm_list[3] += ";" + spl_h4[coln]
            spl_h4[coln] = ""
    
    coln += 1

# Finish with the last id
ll_id = spl_h4[cur_coln]
header_dict[ll_id].append(cur_coln)
id_colns.append(cur_coln)

#write comments
spl_h1[coln - 1] = comm_list[0].replace(";", "", 1)
spl_h2[coln - 1] = comm_list[1].replace(";", "", 1)
spl_h3[coln - 1] = comm_list[2].replace(";", "", 1) 
spl_h4[coln - 1] = comm_list[3].replace(";", "", 1) 

#print "\t", ll_id, cur_coln
out.write("\t".join(spl_h1) + "\n")
out.write("\t".join(spl_h2) + "\n")
out.write("\t".join(spl_h3) + "\n")
out.write("\t".join(spl_h4) + "\n")
out.write("\t".join(spl_h5) + "\n")
    


# In[854]:


# define how sample type corresponds to sample group
stype2group = {}
with open("/Users/dashazhernakova/Documents/UMCG/data/NEXT_tables/stype2sgroup.txt") as stype_f:
    stype2group = {l.rstrip().split("\t")[0] : l.rstrip().split("\t")[1] for l in stype_f}
#stype2group


# In[855]:


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
        who = spl[1]
    if len(spl[2]) > 0:
        what = spl[2]
    #print what
    what = re.sub(r" %s$" % when, "", what)
    #print what
    sgroup = stype2group.get(what)
       
    # loop over all ids for this sample type
    for ll_id, colnums in header_dict.items():
        for coln in colnums:
            # TODO don't forget Bristol stool samples later
            
            #move comments to the last field
            spl[coln:coln + n_fields] = move_comments(spl[coln:coln + n_fields])

            #fill the dict, skip empty samples and irrelevant sample types
            res_list = spl[coln:coln + n_fields]
            #print ll_id, who, when, what, res_list 
            if sgroup and any(res_list) and not res_list[-1] == 'o':
                #what_stripped = re.sub(r" %s$" % when, "", what)
                #res_list2 = res_list
                #print ll_id, who, when, what, res_list
                res_list2 = search_fill_missing_dates(ll_id, who, when, what, res_list, res_dict)
                
                spl[coln:coln + n_fields] = res_list2
                res_dict[sgroup][ll_id + ":" + who + ":" + when + ":" + what] = res_list2

            #check colnames
            assert spl_h5[coln:coln + n_fields] == colnames, spl_h5[coln:coln + n_fields]     
    out.write("\t".join(spl) + "\n")
        
print "Processed lines:", cnt
f.close()
out.close()


# In[789]:


#
# process and fill storage files
#
st_header = ["Box number", "Position", "Lifelines number", "Lifelines Next number", "Time point", "Aliquot number", "Barcode", "Sample issued (yes/No)", "Date Sample issued ", "Remarks"]
storage_fdir = "/Users/dashazhernakova/Documents/UMCG/data/NEXT_tables/Fixed_tables/"
#for filepath in glob.iglob(storage_fdir + '*.txt'):
filepath = storage_fdir + "Breast_milk_M (edited justin).txt"
out_fpath = filepath[:-4] + ".DZprocessed.txt"
file_sgroup = "Breast_milk"
file_who = "Mother"
st_f = open(filepath)
out_st_f = open(out_fpath, 'w')
first_line = st_f.readline().rstrip().split('\t')

#TODO: Comments also side by side
out_st_f.write("\t".join(first_line[:6]) + "\t" + "\t".join(colnames) + "\t" + first_line[9] + "\n")
assert first_line == st_header, first_line
#cnt = 0
lst_to_write = spl[:7]
for l in st_f:
    spl = l.rstrip('\r\n').split('\t')
    ll_id = "%06d" % (int(spl[3]),)  
    res = res_dict[file_sgroup][ll_id + ":" + file_who + ":" + spl[4] + ":" + spl[5]]
    
    lst_to_write.append(res[0])
    lst_to_write.append(res[1])
    
    #compare barcodes, write a merged column
    lst_to_write.append("")
    bc = spl[6]
    if len(res[1]) > 1:
        if len(bc) > 0:
            if not bc == res[0]:
                lst_to_write[-1] = "ERROR: barcodes don't match")
            else:
                lst_to_write[-1] = res[1]
    elif len(bc) > 0:
        if len(res[0]) > 0:
            if not bc == res[0]:
                lst_to_write[-1] = "ERROR: barcodes don't match")
            else:
                lst_to_write[-1] = bc
        else:
            lst_to_write[-1] = bc
    elif len(res[0]) > 0:
        lst_to_write[-1] = res[0]
        
    
    out_st_f.write("\t".join(spl[:6]) + "\t" + "\t".join(res) + "\t" + spl[9] + "\n")
    #cnt += 1
    #if cnt > 10:
    #    break
st_f.close()
out_st_f.close()


# In[611]:


filepath = storage_fdir + "Breast_milk_M (edited justin).txt"
out_fpath = filepath[:-4] + ".processed.txt"
filepath[:-4]


# In[837]:


res_dict.get('Fces')


# In[595]:


s = res_dict['PAXgene']['011774:Mother:P28:1x PAXgene 2.5 mL']
id_spl = s[:]
#move_comments(id_spl)
#print id_spl
print s


# In[523]:


def move_comments(id_spl):
    for i,val in enumerate(id_spl[:-1]):
        if val[:1].isalpha() and not val[:-1].isdigit() and not val == 'NA' and not val == 'n/a' and not val.startswith('ID:') and not val.startswith('SF'):
            id_spl[-1] += ";" + val
            id_spl[i] = ""
    id_spl[-1] = id_spl[-1].replace(";", "", 1)
    return id_spl


# In[515]:


val='SF04317376'
print val[:1].isalpha()
print not val[:-1].isdigit()


# In[804]:


not any(res_list[:6])


# In[849]:


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


# In[624]:


when = 'B'
r" %s$" % when
#re.sub(r" B$", "", 'Plasma-2B B')
re.search(r"Placenta mother-", "Placenta mother-2", 'Plasma-2 B')


# In[705]:


search_fill_missing_dates('012359','Mother', 'P12', 'Plasma-7', ['80564482', '', '', '', '', '', ''], res_dict)
#res_dict['temp_group']['012359:Mother:P12:1x K2EDTA 10 mL']


# In[840]:


res_dict['Plasma_samples'].get('002576:Mother:P12:Plasma-1')


# In[643]:


range(2,6)

