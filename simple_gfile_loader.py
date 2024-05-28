# -*- coding: utf-8 -*-
"""
Created on Mon May 27 10:47:48 2024

@author: user
"""

import load_coord_method as lcm
from colorama import Fore, Style
import re




simu_loc = 'C:/Users/ychuang/Documents/SOLPS_data/simulation_data/mast/027205'
simu_case = 'gfile_lab'
g_loc = '{}/{}/prac_g.00275_efitpp'.format(simu_loc, simu_case)


shift_value = 0.55

with open(g_loc) as f:
     lines = f.readlines()
     

aline = lines[1]
    
x = aline.split()
# new_x = x[1:]


print(x)
rcentr = x[2]
rleft = x[3]

rcentr_n = float(rcentr) + shift_value
rleft_n = float(rleft) + shift_value
R_rcentr = '{:.9E}'.format(rcentr_n)
R_rleft = '{:.9E}'.format(rleft_n)
x[2] = R_rcentr
x[3] = R_rleft


writelist = ' '.join(str(a) for a in x)
lines[1] = ' ' + writelist + "\n"


bline = lines[2]


sp_b = re.findall('[-+]?\d+\.?\d+[eE]?[-+]\d+', lines[2])
mylist =[]
for ii, dat in enumerate(sp_b):
    
    if ii == 0:
        a = float(dat) + shift_value
        b = "{:.9E}".format(a)
    elif ii!= 0 and float(dat) > 0:
        b = ' '+ dat
    else:
        b = dat
    mylist.append(b)



print(mylist)


writelist = ''.join(str(a) for a in mylist)
lines[2] = ' ' + writelist + "\n"



# p = Fore.RED + lines[916]
# print(p)

def line_index_finder(file_lines):
    
    index_dic = {}
    
    for j, string in enumerate(file_lines):
    
        if '361   37' in string:
            # print(j)
            index_dic['bd'] = j
        
        elif 'RCENTR' in string:
            print('rcentr is on line {}'.format(str(j)))
            index_dic['rcentr'] = j
        
        elif 'RBDRY' in string:
            print('rbdry is on line {}'.format(str(j)))
            index_dic['rbdry'] = j
        
        elif 'XLIM' in string:
            print('xlim is on line {}'.format(str(j)))
            index_dic['xlim'] = j
            
        
        
        
        
    
    return index_dic

index_dic = line_index_finder(file_lines = lines)
bd_st = index_dic['bd'] + 1
print(bd_st) 

def hop_corrector(line_total, line_start, lines):
    for ln in range(line_total):
        
        aline = lines[line_start + ln]
        line_end = line_start + ln
        
        spl = re.findall('[-+]?\d+\.?\d+[eE]?[-+]\d+', aline)
        mylist =[]
        for ii, dat in enumerate(spl):
            
            if ln % 2 == 0:
             
                if ii % 2 == 0:
                    a = float(dat) + shift_value
                    if ii!= 0 and float(dat) >= 0:
                        
                        b = ' ' + "{:.9E}".format(a)
                    else:
                        
                        b = "{:.9E}".format(a)
                        
                else:
                    if ii!= 0 and float(dat) > 0:
                        
                        b = ' '+ dat
                    else:
                        b = dat
            else:
                if ii % 2 != 0:
                    a = float(dat) + shift_value
                    if ii!= 0 and float(dat) >= 0:
                        
                        b = ' ' + "{:.9E}".format(a)
                    else:
                        
                        b = "{:.9E}".format(a)
                        
                else:
                    if ii!= 0 and float(dat) > 0:
                        
                        b = ' '+ dat
                    else:
                        b = dat
            
            mylist.append(b)
        
        # print(mylist)
        
        writelist = ''.join(str(a) for a in mylist)
        
        if float(spl[0]) < 0:
            lines[line_end] = writelist + "\n"
        else:
            lines[line_end] = ' ' + writelist + "\n"
    
    return line_end, lines

lk, lines = hop_corrector(line_total = 145, line_start = bd_st, lines = lines)       
    
print('check the end line: {}'.format(lines[lk]))
lim_st = lk + 1


lk, lines = hop_corrector(line_total = 15, line_start = lim_st, lines = lines)
print('check the end line: {}'.format(lines[lk]))




in_rcentr = index_dic['rcentr']
number = re.findall('[-+]?\d+\.\d+', lines[in_rcentr])


search_text = number[0]
print(number[0])

shift_axis = float(number[0]) + shift_value
print(shift_axis)


exp_sp = lines[in_rcentr].split()
print(exp_sp)

exp_sp[2] = '{:.8f}'.format(shift_axis)

writelist = ' '.join(str(a) for a in exp_sp)
lines[in_rcentr] = writelist + "\n"


in_rbd = index_dic['rbdry']


aline = lines[in_rbd]
tk = aline.split()
x = tk[:2]

sp_b = re.findall('[-+]?\d+\.\d+', aline)
mylist =[]
for ii, dat in enumerate(sp_b):
    
    if float(dat) <= 0:
        print(dat)
        b = "{:.6f}".format(float(dat))
    
    else:
        a = float(dat) + shift_value
        b = "{:.6f}".format(a)
        
        
    mylist.append(b)


print(mylist)

write_text = ''.join(str(a) + ' ' for a in x)
writelist = ''.join("     " + str(a)  for a in mylist)
lines[in_rbd] = write_text + writelist + "\n"



rbd_st = index_dic['rbdry'] + 1


for ln in range(15):
    
    aline = lines[rbd_st + ln]
    line_end = rbd_st + ln
    
    
    sp_b = re.findall('[-+]?\d+\.\d+', aline)
    mylist =[]
    for ii, dat in enumerate(sp_b):
        
        if float(dat) <= 0:
            print(dat)
            b = "{:.6f}".format(float(dat))
        
        else:
            a = float(dat) + shift_value
            b = "{:.6f}".format(a)
            
            
        mylist.append(b)


    print(mylist)


    writelist = ''.join(str(a) + "     " for a in mylist)
    lines[line_end] = "      " + writelist + "\n"



'======== xlimit ==========='


in_xlim = index_dic['xlim']


aline = lines[in_xlim]
tk = aline.split()
x = tk[:2]

sp_b = re.findall('[-+]?\d+\.\d+', aline)
mylist =[]
for ii, dat in enumerate(sp_b):
    
    if float(dat) <= 0:
        print(dat)
        b = "{:.6f}".format(float(dat))
    
    else:
        a = float(dat) + shift_value
        b = "{:.6f}".format(a)
        
        
    mylist.append(b)


print(mylist)

write_text = ''.join(str(a) + ' ' for a in x)
writelist = ''.join("     " + str(a)  for a in mylist)
lines[in_xlim] = write_text + writelist + "\n"



xlim_st = index_dic['xlim'] + 1


for ln in range(6):
    
    aline = lines[xlim_st + ln]
    line_end = xlim_st + ln
    
    
    sp_b = re.findall('[-+]?\d+\.\d+', aline)
    mylist =[]
    for ii, dat in enumerate(sp_b):
        
        if float(dat) <= 0:
            print(dat)
            b = "{:.6f}".format(float(dat))
        
        else:
            a = float(dat) + shift_value
            b = "{:.6f}".format(a)
            
            
        mylist.append(b)


    print(mylist)


    writelist = ''.join(str(a) + "     " for a in mylist)
    lines[line_end] = "      " + writelist + "\n"





m_gfile = '{}/{}/m_prac_g.00275_efitpp'.format(simu_loc, simu_case)


with open(m_gfile,'w') as g:
    for i, line in enumerate(lines):         ## STARTS THE NUMBERING FROM 1 (by default it begins with 0)    
        g.writelines(line)



"""

j= index_dic['r']

while j< index_dic['z']:
    string2 = datalist[j]
    gdata = re.findall('[-+]?\d+\.?\d+[eE]?[-+]\d+', string2)
    mylist =[]
    for ii in gdata:
        a = np.float64(ii) +np.float64(shift_value)
        b = "{:.8E}".format(a)
        mylist.append(b)
    writelist = ' '.join(str(x)+"   " for x in mylist)
    datalist[j] = '    ' + writelist + "\n"
    j= j+ 1

    
    def replacetext(search_text, replace_text):
        with open(g_loc, "r+") as f:
            file_content = f.read()
            file_content = re.sub(search_text, replace_text, file_content)
            f.seek(0)
            f.write(file_content)
            f.truncate()
            return "Text replaced"

    search_text = rcentr

    shift_axis = float(rcentr) + shift_value
    print(shift_axis)

    replace_text = '{:.9E}'.format(shift_axis)
    print(replacetext(search_text = search_text, replace_text = replace_text))
"""



mod_gfile = lcm.loadg(m_gfile)
g_file = lcm.loadg(g_loc)