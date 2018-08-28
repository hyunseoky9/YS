import os
import numpy as np
import pandas
import matplotlib.pyplot as plt
from scipy.interpolate import spline
def display(direc,win,kmplot,p2print,sortpar,repcap,histrep):
    strings = os.listdir(direc)
    sortinfo = []
    sortamount = len(sortpar)
    dcm = {'timestep':0,'krecord':1,'rep':2,'s':3,'N0':4,'K0':5,'u':6,'gen_num':7,'c':8,'r':9,'kmax':10,'host_num':11,
    'mig':12,'tr':13,'kf':14, 'kfvar':15, 'evolrate':16,'evolq':17}
    dcc = {'type':0,'q':1,'a':2,'b':3,'back':4,'rep':5,'L':6,'s':7,'N0':8,'K':9,'mu':10,
    'gen_num':11,'cost':12,'r':13,'N1r':14,'timestep':15,'krecord':16,'gen_num':17}
    for i in range(len(strings)):
        s = strings[i].split('_')
        model = (s[0][0],s[0][1:-1]) # ver = (model type,version)
        params = s[1].split('(')[0]
        params = params.split(',')
        if model[0] == 'm':
            dc = dcm
            timestep = int(params[0])
            params[0] = int(params[0])
            params[1] = int(params[1])
            params[2] = int(params[2])
            params[3] = float(params[3])
            params[4] = int(params[4])
            params[5] = int(params[5])
            params[6] = float(params[6])
            params[7] = int(params[7])
            params[8] = float(params[8])
            params[9] = float(params[9])
            params[10] = int(params[10])
            params[11] = int(params[11])
            if model[1] == '1.1.3':
                params[12] = float(params[12])
                params[13] = float(params[13])
        elif model[0] == 'c':
            dc = dcc
            timestep = int(params[0])
            params[0] = int(params[0])
            params[1] = float(params[1])
            params[2] = float(params[2])
            params[3] = float(params[3])
            params[4] = int(params[4])
            params[5] = int(params[5])
            params[6] = int(params[6])
            params[7] = float(params[7])
            params[8] = int(params[8])
            params[9] = int(params[9])
            params[10] = float(params[10])
            params[11] = int(params[11])
            params[12] = float(params[12])
            params[13] = float(params[13])
            params[14] = float(params[14])
            params[15] = int(params[15])
            params[16] = int(params[16])
            params[17] = int(params[17])

        print(params)
        # for differently labeled head row in meta1.1.3 data
        sortup = []
        for j in range(len(sortpar)):
            sortup.append(params[dc[sortpar[j]]])
        sortup.append(i)
        sortinfo.append(sortup)

        #if model[0] == 'm':
        #    timestep = int(params[0])
        #    krecord = int(params[1])
        #   rep = int(params[2])
        #    s = float(params[3])
        #    N0 = int(params[4])
        #    K = int(params[5])
        #    u = float(params[6])
        #    gen_num = int(params[7])
        #    c = float(params[8])
        #    r = float(params[9])
        #    kmax = int(params[10])
        #    host_num = int(params[11])
        #    if model[1] == '1.1.3':
        #        mig = float(params[12])
        #        tr = float(params[13])
        # for differently labeled head row in meta1.1.3 data
    sortinfo.sort()
    for i in range(len(sortinfo)):
        print(strings[sortinfo[i][sortamount]])
        s = strings[sortinfo[i][sortamount]].split('_')
        model = (s[0][0],s[0][1:-1]) # ver = (model type,version)
        params = s[1].split('(')[0]
        params = params.split(',')
        if model[0] == 'm':
            dc = dcm
            timestep = int(params[0])
            params[0] = int(params[0])
            params[1] = int(params[1])
            params[2] = int(params[2])
            params[3] = float(params[3])
            params[4] = int(params[4])
            params[5] = int(params[5])
            params[6] = float(params[6])
            params[7] = int(params[7])
            params[8] = float(params[8])
            params[9] = float(params[9])
            params[10] = int(params[10])
            params[11] = int(params[11])
            if model[1] == '1.1.3':
                params[12] = float(params[12])
                params[13] = float(params[13])
        elif model[0] == 'c':
            dc = dcc
            timestep = int(params[0])
            params[0] = int(params[0])
            params[1] = float(params[1])
            params[2] = float(params[2])
            params[3] = float(params[3])
            params[4] = int(params[4])
            params[5] = int(params[5])
            params[6] = int(params[6])
            params[7] = float(params[7])
            params[8] = int(params[8])
            params[9] = int(params[9])
            params[10] = float(params[10])
            params[11] = int(params[11])
            params[12] = float(params[12])
            params[13] = float(params[13])
            params[14] = float(params[14])
            params[15] = int(params[15])
            params[16] = int(params[16])
            params[17] = int(params[17])
        if model[1] == '1.1.3':
            k1 = 'k1.0'
            k2 = 'k2.0'
            pop1 = 'pop1.0'
            pop2 = 'pop2.0'
        else:
            k1 = 'k1'
            k2 = 'k2'
            pop1 = 'pop1'
            pop2 = 'pop2'
        if win:
            pop1dead = 0
            pop2dead = 0
        data_a = pandas.read_csv(direc+'/'+strings[sortinfo[i][sortamount]])
        data_a2 = data_a.loc[data_a['rep']==1]
        if len(repcap) == 0:
            repcap = [0,int(params[dc['rep']])]
        else:
            if repcap[1] > int(params[dc['rep']]):
                repcap[1] = int(params[dc['rep']])
        endmeank1 = 0
        endmeank2 = 0
        ### graphing and calculating win ratio
        if params[dc['timestep']]:
            for i in range(len(p2print)):
                print('%s=%s'%(p2print[i],params[dc[p2print[i]]]))
            if params[dc['krecord']] == 0 or params[dc['krecord']] == 1: # for mean k and min k plot
                if repcap[0] <= 1 and repcap[1] >= 1: # graph only if the rep is in the repcap range  //and -1 in list(data_a2[k1])
                    plt.figure(figsize=(18, 16), dpi= 80, facecolor='w', edgecolor='k')
                    plt.subplot(131)
                    if kmplot == 3:
                        plt.plot(list(data_a2[k2]),'r',label='k2')
                    elif kmplot == 2:
                        plt.plot(list(data_a2[k1]),'g',label='k1')
                    else:
                        plt.plot(list(data_a2[k1]),'g',label='k1')
                        plt.plot(list(data_a2[k2]),'r',label='k2')
                    plt.legend()
                    plt.ylabel('mean k')
                    plt.xlabel('generations')
                    plt.title('k mean')
                    plt.ylim(0,25)
                    plt.xlim(0,params[dc['gen_num']])

                    plt.subplot(132)
                    if kmplot ==3:
                        plt.plot(list(data_a2[pop2]),'r',label='pop2')
                    elif kmplot == 2:
                        plt.plot(list(data_a2[pop1]),'g',label='pop1')
                    else:
                        plt.plot(list(data_a2[pop1]),'g',label='pop1')
                        plt.plot(list(data_a2[pop2]),'r',label='pop2')
                    plt.legend()
                    plt.ylabel('pop#')
                    plt.xlabel('generations')
                    plt.title('pop')

                    if kmplot == 1:
                        plt.subplot(133)
                        ratio = np.array(data_a2[pop2])/(np.array(data_a2[pop2]) + np.array(data_a2[pop1]))
                        plt.plot(ratio,'r')
                        plt.ylabel('seg2ratio')
                        plt.xlabel('generations')
                        plt.title('seg2 ratio')
                if win:
                    if -1 in list(data_a2[k1]) and -1 not in list(data_a2[k2]):
                        pop1dead += 1
                    if -1 in list(data_a2[k2]) and -1 not in list(data_a2[k1]):
                        pop2dead += 1
                    endmeank2 += np.mean(list(data_a2[k2])[100::])
                    endmeank1 += np.mean(list(data_a2[k1])[100::])
            else: # histogram
                for k in histrep:
                    data_a2 = data_a.loc[data_a['rep']==k]
                    sidec = np.linspace(204,0,params[dc['gen_num']]) #colors for the 2 other in rgb value for pop1(green) and pop2(red) color
                    fig=plt.figure(figsize=(18, 16), dpi= 80, facecolor='w', edgecolor='k')
                    #print('sidec=',sidec)
                    for index, row in data_a2.iterrows():
                        gen = int(list(row)[1])
                        if gen == 1900: # print certain generation or generations
                            print(list(row))
                            row = list(row)[4::]
                            #print(row)
                            #print(index )
                            pop1 = []
                            pop2 = []
                            print(row)
                            for i in range(len(row)):
                                if i % 2 == 0:
                                    pop1.append(row[i])
                                else:
                                    pop2.append(row[i])
                            # get mean and variance.
                            p1var = 0
                            p2var = 0
                            p2mean = 0
                            p1mean = 0 # [mean,var]
                            p2stat = []
                            for i in range(len(pop1)):
                                p2mean += pop2[i]*i
                                p1mean += pop1[i]*i
                            p2mean = p2mean/np.sum(pop2)
                            p1mean = p1mean/np.sum(pop1)
                            for i in range(len(pop1)):
                                p2var += pop2[i]*(i - p2mean)**2
                                p1var += pop1[i]*(i - p1mean)**2
                            p2var = p2var / np.sum(pop2)
                            p1var = p1var / np.sum(pop1)
                            print('seg1 var: %f, mean: %f'%(p1var,p1mean))
                            print('seg2 var: %f, mean: %f'%(p2var,p2mean))
                            x = np.linspace(0,2*params[dc['kmax']],2*params[dc['kmax']] + 1)
                            x_smooth = np.linspace(0,2*params[dc['kmax']],200)
                            p1_smooth = spline(x,pop1,x_smooth)
                            p2_smooth = spline(x,pop2,x_smooth)

                            c1 = (sidec[gen]/255,255/255,sidec[gen]/255)
                            c2 = (255/255,sidec[gen]/255,sidec[gen]/255)
                            if kmplot == 1:
                                plt.plot(x_smooth,p1_smooth,color=c1,label='pop1')
                                plt.plot(x_smooth,p2_smooth,color=c2,label='pop2')
                            elif kmplot == 2:
                                plt.plot(x_smooth,p1_smooth,color=c1,label='pop1')
                            else:
                                plt.plot(x_smooth,p2_smooth,color=c2,label='pop2')
                    plt.legend()
                    plt.ylabel('count')
                    plt.xlim(0,10)
                    plt.xlabel('generations')
                    plt.title('rep=%d'%(1))
                    plt.show()
            print('reps=',int(params[dc['rep']]))
            num = 0
            for i in range(int(params[dc['rep']])):
                data_a2 = data_a.loc[data_a['rep']==(i+1)]
                if params[dc['krecord']] == 0 or params[dc['krecord']] == 1:
                    if win:
                        if -1 in list(data_a2[k1]) and -1 not in list(data_a2[k2]):
                            pop1dead += 1
                        if -1 in list(data_a2[k2]) and -1 not in list(data_a2[k1]):
                            pop2dead += 1
                        endmeank1 += np.mean(list(data_a2[k2])[100::])
                        endmeank2 += np.mean(list(data_a2[k1])[100::])

                    ##########DELETE -1 IN LIST(...)!!!!!######
                    if repcap[0] <= i and repcap[1] >= i: #and -1 in list(data_a2[k2]): # graph only if the rep is in the repcap range ///
                        num+=1
                        plt.subplot(131)
                        if kmplot == 3:
                            plt.plot(list(data_a2[k2]),'r')
                        elif kmplot == 2:
                            plt.plot(list(data_a2[k1]),'g')
                        else:
                            plt.plot(list(data_a2[k1]),'g')
                            plt.plot(list(data_a2[k2]),'r')

                        plt.subplot(132)
                        if kmplot == 3:
                            plt.plot(list(data_a2[pop2]),'r')
                        elif kmplot == 2:
                            plt.plot(list(data_a2[pop1]),'g')
                        else:
                            plt.plot(list(data_a2[pop1]),'g')
                            plt.plot(list(data_a2[pop2]),'r')
                        if kmplot == 1:
                            plt.subplot(133)
                            ratio = np.array(data_a2[pop2])/(np.array(data_a2[pop2]) + np.array(data_a2[pop1]))
                            plt.plot(ratio,'r')
                '''
                if params[dc['krecord']] == 2:
                    if repcap[0] <= i and repcap[1] >= i: # graph only if the rep is in the repcap range  //and -1 in list(data_a2[k1])
                        sidec = np.linspace(204,0,params[dc['gen_num']]) #colors for the 2 other in rgb value for pop1(green) and pop2(red) color
                        fig=plt.figure(figsize=(18, 16), dpi= 80, facecolor='w', edgecolor='k')
                        #print('sidec=',sidec)
                        for index, row in data_a2.iterrows():
                            if index % 100 == 0:
                                gen = int(list(row)[1])
                                row = list(row)[4::]
                                #print(row)
                                #print(index )
                                pop1 = []
                                pop2 = []
                                for j in range(len(row)):
                                    if i % 2 == 1:
                                        pop1.append(row[j])
                                    else:
                                        pop2.append(row[j])
                                x = np.linspace(0,2*params[dc['kmax']],2*params[dc['kmax']] + 1)
                                x_smooth = np.linspace(0,2*params[dc['kmax']],200)
                                p1_smooth = spline(x,pop1,x_smooth)
                                p2_smooth = spline(x,pop2,x_smooth)
                                c1 = (sidec[gen]/255,255/255,sidec[gen]/255)
                                c2 = (255/255,sidec[gen]/255,sidec[gen]/255)
                                if kmplot == 1:
                                    plt.plot(x_smooth,p1_smooth,color=c1,label='pop1')
                                    plt.plot(x_smooth,p2_smooth,color=c2,label='pop2')
                                elif kmplot == 2:
                                    plt.plot(x_smooth,p1_smooth,color=c1,label='pop1')
                                else:
                                    plt.plot(x_smooth,p2_smooth,color=c2,label='pop2')
                        plt.ylabel('count')
                        plt.xlabel('generations')
                        plt.title('rep=%d'%(params[dc['rep']]))
                        plt.show()
            '''
            if params[dc['krecord']] != 2:
                print("num=%d"%(num))
                if win:
                    if pop1dead == 0 and pop2dead == 0:
                        print('2seg winning prob: nan')
                    else:
                        print('2 seg winning prob: %.3f'%(pop1dead/(pop1dead + pop2dead)))
                        print('endmeank1=',endmeank1/int(params[dc['rep']]))
                        print('endmeank2=',endmeank2/int(params[dc['rep']]))
                plt.show()
        else:   
            print("data needs to have timestep to display plot!")
            break
