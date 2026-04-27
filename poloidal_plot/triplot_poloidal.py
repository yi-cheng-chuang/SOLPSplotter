# -*- coding: utf-8 -*-
"""
Created on Mon Nov  3 05:01:08 2025

@author: ychuang
"""

from SOLPS_input.header import *
import matplotlib.pylab as pylab


class series_triplots:

    def __init__(self, DF, data):

        self.DF = DF
        self.data = data

    def twscan_opacity_polplot_method(self, iterlist, cl_dic, ang_list, format_option, space_option):

        if format_option == '3x1':

            fig, axs = plt.subplots(3, 1)

        elif format_option == '2x1':

            fig, axs = plt.subplots(2, 1)

        elif format_option == '1x1':

            fig, axs = plt.subplots()

        nx = self.data['b2fgeo']['nx']
        ny = self.data['b2fgeo']['ny']

        plt.subplots_adjust(hspace=.0)
        anchored_text_1 = AnchoredText('{}'.format('Density pedestal width [mm]'),
                                       loc='upper left')
        anchored_text_2 = AnchoredText('{}'.format('Neutral penetration length [mm]'),
                                       loc='upper center')
        anchored_text_3 = AnchoredText('{}'.format('Neutral opaqueness'),
                                       loc='upper left')
        anchored_text_4 = AnchoredText('{}'.format('Neutral density'),
                                       loc='upper center')

        marker = ['o', 's', 'x', '^', 'v', '*', '.', 'o', 's']
        lins_map = dict(zip(iterlist, marker))

        for aa in iterlist:

            """
            label= 'core density {} $10^{19}$'.format(aa)

            """

            efold = self.data['opacity_poloidal'][aa]['efold_length'] * \
                pow(10, 3)
            width = self.data['opacity_poloidal'][aa]['pedestal_width'] * \
                2 * pow(10, 3)
            opq = self.data['opacity_poloidal'][aa]['dimensionless_opaqueness']
            nd_sep = self.data['opacity_poloidal'][aa]['neutral_density']
            efold_psiN = self.data['opacity_poloidal'][aa]['efold_length_psiN']
            width_psiN = self.data['opacity_poloidal'][aa]['pedestal_width_psiN']

            nd_sep_std = self.data['opacity_poloidal']['1.0']['neutral_density']

            nd_change_percent = np.divide(
                nd_sep, nd_sep_std, out=np.zeros_like(nd_sep, dtype=float), where=nd_sep_std != 0) * 100

            if format_option == '3x1':

                # axs[0].set_title('opacity triplot')

                if space_option == 'psiN':
                    axs[0].plot(ang_list, width_psiN, lins_map[aa] + '-', color=cl_dic[aa], lw=1.5,
                                label='{}'.format(aa))
                    axs[1].plot(ang_list, efold_psiN,
                                lins_map[aa] + '-', color=cl_dic[aa], lw=1.5)

                elif space_option == 'real':

                    axs[0].plot(ang_list, width, lins_map[aa] + '-', color=cl_dic[aa], lw=1.5,
                                label='{}'.format(aa))
                    axs[1].plot(ang_list, efold, lins_map[aa] +
                                '-', color=cl_dic[aa], lw=1.5)
                    axs[0].set_ylim(0, 20)
                    axs[1].set_ylim(0, 20)
                    axs[2].set_ylim(0, 4)

                axs[2].plot(ang_list, opq, lins_map[aa] +
                            '-', color=cl_dic[aa], lw=1.5)
                axs[2].set_xlabel('poloidal angle')
                axs[0].legend(loc='best')

            elif format_option == '2x1':

                # axs[0].set_title('opacity twoplot')
                axs[0].plot(ang_list, width, lins_map[aa] +
                            '-', lw=1.5, color=cl_dic[aa], label='{}'.format(aa))
                axs[1].plot(ang_list, opq, lins_map[aa] +
                            '-', lw=1.5, color=cl_dic[aa])
                axs[0].set_ylim(0, 20)
                axs[1].set_ylim(0.5, 2.0)
                axs[1].set_xlabel('poloidal angle')
                axs[0].legend(loc='best')

            elif format_option == '1x1':

                # axs.set_title('atomic neutral density change percentage')
                # axs.plot(ang_list, nd_change_percent, lins_map[aa] + '-', lw=1.5,
                #          color=cl_dic[aa], label='{}'.format(aa))
                # axs.plot(ang_list, nd_sep, lins_map[aa] + '-', lw=1.5,
                #          color=cl_dic[aa], label='{}'.format(aa))
                axs.plot(ang_list, efold, lins_map[aa] + '-', lw=1.5,
                         color=cl_dic[aa], label='{}'.format(aa))
                axs.set_xlabel('poloidal angle')
                # axs.set_yscale('log')
                axs.set_ylim(4, 20)
                axs.legend(loc='best')

        if format_option == '3x1':

            axs[0].add_artist(anchored_text_1)
            axs[1].add_artist(anchored_text_2)
            axs[2].add_artist(anchored_text_3)
            axs[0].axvline(x=0, color='black', lw=3, ls='-')
            axs[1].axvline(x=0, color='black', lw=3, ls='-', label='LFS')
            axs[2].axvline(x=0, color='black', lw=3, ls='-')
            axs[0].axvline(x=180, color='brown', lw=3, ls='-')
            axs[1].axvline(x=180, color='brown', lw=3, ls='-', label='HFS')
            axs[2].axvline(x=180, color='brown', lw=3, ls='-')
            axs[1].legend(loc='best')

        elif format_option == '2x1':

            axs[0].add_artist(anchored_text_1)
            axs[1].add_artist(anchored_text_3)
            axs[0].axvline(x=0, color='black', lw=3, ls='-')
            axs[1].axvline(x=0, color='black', lw=3, ls='-', label='LFS')
            axs[0].axvline(x=180, color='brown', lw=3, ls='-')
            axs[1].axvline(x=180, color='brown', lw=3, ls='-', label='HFS')
            axs[1].legend(loc='best')
            axs[0].grid(True)
            axs[1].grid(True)

        elif format_option == '1x1':

            axs.add_artist(anchored_text_2)
            # axs.add_artist(anchored_text_2)
            axs.axvline(x=0, color='black', lw=3, ls='-', label='LFS')
            axs.axvline(x=180, color='brown', lw=3, ls='-', label='HFS')
            axs.legend(loc='best')
            axs.grid(True)

    def tri_opacity_polplot(self, format_option, space_option):

        # print(adj_list)
        # for i in itemname:
        #     plt.figure()
        # if log_flag:
        #     plt.yscale('log')
        # else:
        #     pass

        if self.DF.withshift == False and self.DF.withseries == True:

            density_dic = {}
            for k in self.data['dircomp']['Attempt'].keys():
                density_dic[k] = k

            labels = list(self.data['dircomp']['Attempt'].keys())
            # li = len(itemname)
            # colors = plt.cm.tab10(np.linspace(0, 1, len(labels)))
            colors = pylab.cm.gnuplot2(np.linspace(
                0, 1, len(labels) + 1))
            color_map = dict(zip(labels, colors))
            # marker = ['o', 's', 'x', '^', 'v', '*', '.', 'o', 's']
            # lins_map = dict(zip(labels, marker))
            # print(color_map)

            # Example usage:
            # plt.scatter(x, y, color=color_map['C'])

            # result = self.data['opacity_poloidal'][aa]
            pol_loc = self.data['angle']['angle_list']

            self.twscan_opacity_polplot_method(
                iterlist=labels, cl_dic=color_map, ang_list=pol_loc, format_option=format_option, space_option=space_option)

        elif self.DF.withshift == True and self.DF.withseries == True:
            print('poloidal_plot function is not there yet!')

        else:
            print('poloidal_plot function has a bug.')

        """
        
        
            if self.DF.withshift == False and self.DF.withseries == False:

                result = self.data['opacity_poloidal']

                unit = self.opacity_study_unit()
                pol_loc = self.data['angle']['angle_list']
                xpoint = self.data['angle']['xpoint_angle']
                a_shift = self.data['dircomp']['a_shift']
                A_val = A_dic[a_shift]
                color = color_dic[a_shift]

                self.opacity_poloidal_plot_method(item=i, pol_angle=pol_loc,
                                                  result_dic=result, color_code=color, A_value=A_val, unit_dic=unit)

                self.poloidal_label(angle_fix=pol_loc,
                                    item=i, xpoint_fix=xpoint)
        
        """


"""   
    
    def twscan_opacity_polplot(self, scan_style, plot_option, format_option, space_option):
        
        
        withshift = self.DF.withshift
        withseries = self.DF.withseries
        
        
        if withshift == False and withseries == True:
            
            # series_flag = self.DefaultSettings['series_flag']
            
            
            if self.DF.series_flag == 'twin_scan':
                
                dircomp = self.data['dircomp']
                
                if scan_style == 'tempscan':
                    
                    key_a = 'denscan_list'
                    key_b = 'tempscan_list'
                
                elif scan_style == 'denscan':
                    
                    key_a = 'tempscan_list'
                    key_b = 'denscan_list'
                
                else:
                    print('twinscan_plot_method, please check the scan_style!')
                
                keylist_a = []
                
                # pol_list_a = []
                # for i in range(48):
                #     pol_list_a.append('{}'.format(25 + i))
                
                
                # self.calc_pol_angle(pol_list = pol_list_a, plot_angle= False)
                ang_list = self.data['angle']['angle_list']
                
                for x in dircomp[key_a]:
                    keylist_a.append('{:.3f}'.format(x))
                
                for ta in keylist_a:
                    
                    keylist_b = []
                    
                    for x in dircomp[key_b]:
                        keylist_b.append('{:.3f}'.format(x))
                    
                    
                    iter_key, color_dic= self.twa.twscan_plot_prep(ta = ta, 
                    keylist_b = keylist_b, scan_style = scan_style)
                    
                    
                    print('check:')
                    print(iter_key)
                    print(color_dic)
                    # print(label_dic)
                    
                    self.twscan_opacity_polplot_method(iterlist = iter_key, cl_dic = color_dic, 
                                scan_style = scan_style, ang_list = ang_list, 
                    plot_option = plot_option, format_option = format_option, space_option = space_option)
                    
             
            else:
                print('neteTS_plot, please check the series flag')


"""
