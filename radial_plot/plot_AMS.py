# -*- coding: utf-8 -*-
"""
Created on Sun Apr 13 20:16:19 2025

@author: ychuang
"""


class midnd_plot:
    
    def __init__(self, DF, data):
        
        self.DF = DF
        self.data = data
    
    
    
    def AM_NTplot_method(self, iterlist, cl_dic, A_dic, AM_flag, scandetail, scan_style):
        
        
        fig, axs = plt.subplots(2, 1)
        
        for aa in iterlist:
            
            psi_coord, mid_neu_pro, mid_atom_temp_pro = self.AM_NT_midprof(itername = aa, 
                                                                AM_flag= AM_flag)
            
            if AM_flag == 'atom':
                
                anchored_text = AnchoredText('(a){}'.format('atomic density [$m^{-3}$]'), loc='upper center')
                anchored_text2 = AnchoredText('(b){}'.format('atomic temperature [eV]'), loc= 'upper center')
            
            elif AM_flag == 'mol':
                
                anchored_text = AnchoredText('(a){}'.format('molecular density [$m^{-3}$]'), loc= 'upper center')
                anchored_text2 = AnchoredText('(b){}'.format('molecular temperature [eV]'), loc= 'upper center')
                

            
            if self.series_flag == 'twin_scan':
                
                if scan_style == 'tempscan':
                    
                    ad = aa[1]
                
                elif scan_style == 'denscan':
                    
                    ad = aa[0]
                
                else:
                    print('AM_NTplot_method, please check scan_style')
                
                
            
            else:
                ad = aa

            if scan_style == 'denscan':
                
                axs[0].set_title('Density scan with Te = {} eV'.format(scandetail))
                axs[0].plot(psi_coord, mid_neu_pro, color = cl_dic[ad], 
                            label= '{}'.format(A_dic[ad]))
                axs[1].plot(psi_coord, mid_atom_temp_pro, color = cl_dic[ad])
                axs[0].add_artist(anchored_text)
                axs[1].add_artist(anchored_text2)
                axs[1].set_xlabel('$\psi_N$')
                axs[0].legend()

                            
            elif scan_style == 'tempscan':
                
                
                axs[0].set_title('Temperature scan with Ne = {}'.format(scandetail))
                axs[0].plot(psi_coord, mid_neu_pro, color = cl_dic[ad])
                axs[1].plot(psi_coord, mid_atom_temp_pro, color = cl_dic[ad], 
                            label= '{}'.format(A_dic[ad]))
                axs[0].add_artist(anchored_text)
                axs[1].add_artist(anchored_text2)
                axs[1].set_xlabel('$\psi_N$')
                axs[1].legend()
            
            else:
                print('AM_NTplot_method, please check the scan parameter')

            plt.subplots_adjust(hspace=.0)


    def twinscan_plot_method(self, dircomp, scan_style, AM_flag):
        
        
        if scan_style == 'tempscan':
            
            key_a = 'denscan_list'
            key_b = 'tempscan_list'
        
        elif scan_style == 'denscan':
            
            key_a = 'tempscan_list'
            key_b = 'denscan_list'
        
        else:
            print('twinscan_plot_method, please check the scan_style!')
        
        keylist_a = [str(x) for x in dircomp[key_a]]
        
        for ta in keylist_a:
            
            keylist_b = [str(x) for x in dircomp[key_b]]
            color_list = ['red', 'orange', 'green', 'blue', 'purple']
            
            color_dic = self.pair_dic(keys = keylist_b, values = color_list)
            
            scan_list = []
            iter_key = []
            
            for tb in keylist_b:
                
                
                if scan_style == 'tempscan':
                    
                    it_in = (ta, tb)
                
                elif scan_style == 'denscan':
                    
                    it_in = (tb, ta)
                
                else:
                    print('twinscan_plot_method, please check the scan_style!')
                    
                    
                psi_coord, mid_ne_pro, mid_te_pro = self.nete_midprof(itername = it_in)
                
                if scan_style == 'tempscan':
                    
                    scan_add = '{:.1f} eV'.format(mid_te_pro[0])
                
                elif scan_style == 'denscan':
                    
                    scan_add = '{:.2E}'.format(mid_ne_pro[0])
                
                else:
                    print('twinscan_plot_method, please check the scan_style!')
                    
                
                
                scan_list.append(scan_add)
                iter_key.append(it_in)
            
            
            if scan_style == 'tempscan':
                psi_coord, mid_ne_pro, mid_te_pro = self.nete_midprof(itername = (ta, '3.73'))
                scan_title = '{:.2E}'.format(mid_ne_pro[0])
            
            elif scan_style == 'denscan':
                psi_coord, mid_ne_pro, mid_te_pro = self.nete_midprof(itername = ('5.02', ta))
                scan_title = '{:.1f}'.format(mid_te_pro[0])
            
            else:
                print('twinscan_plot_method, please check the scan_style!')
            
            label_dic = self.pair_dic(keys = keylist_b, values = scan_list)
            
            self.AM_NTplot_method(iterlist = iter_key, AM_flag = AM_flag,
                    cl_dic = color_dic, A_dic = label_dic,
                    scandetail = scan_title, scan_style = scan_style)
                
            print(scan_list)




    def AtomNT_plot(self, AM_flag, scan_style):
        
        if self.withshift == True and self.withseries == False:
            
            color_dic = {'org': 'red', 'dot3': 'orange', 'dot5': 'green',
                         'dot7': 'blue', 'one': 'purple'}
            
            label_dic = {'org': '1.4', 'dot3': '2.0', 'dot5': '2.4',
                      'dot7': '2.8', 'one': '3.4'}
            
            asp_ch = self.data['dircomp']['multi_shift']
            
            self.AtomNTplot_method(iterlist = asp_ch, 
                                      cl_dic = color_dic, A_dic = label_dic, scan = 'not')
        
        elif self.withshift == False and self.withseries == True:
            
            
            if self.series_flag == 'twin_scan':
                
                dircomp = self.data['dircomp']
                
                self.twinscan_plot_method(dircomp = dircomp, scan_style = scan_style, 
                                    AM_flag = AM_flag)
            
             
            else:
                print('AtomNT_plot, please check the series flag')



