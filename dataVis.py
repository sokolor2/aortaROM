import pandas as pd
import os
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import seaborn as sns
import numpy as np


def pulseHist(resistanceID: int, resistances: list):
    #Plots the pulse histograms for all flows
    fig, ax = plt.subplots(nrows=4, ncols=5, figsize=(50,5))

    for i in range(0,4):
        for j in range(0,5):
            file_name = f'/home/sokolor2/Data/Flow{j}.csv'

            df = pd.read_csv(file_name)
            filteredFlow = df[df['Resistance'].astype(int) == resistances[i*2]].copy()
            timeCols = [f'TimePoint:{k}' for k in range(0,19)]
            filteredFlow.loc[:, timeCols] = filteredFlow[timeCols].apply(pd.to_numeric, errors='coerce')
    
            result = (
                filteredFlow[['Calico ID'] + timeCols]
                .assign(Pulse_Pressure=lambda x: x[timeCols].max(axis=1) - x[timeCols].min(axis=1))
                .reset_index()
            )
    
            ax[i,j].hist(result['Pulse_Pressure'], bins=30, range=(10,50))
            #ax[i,j].set_title(f'Flow {i}: Resistance {resistanceID}')
            ax[i,j].set_ylim(0, 250)

            if j != 0:
                ax[i, j].set_yticklabels([])
            

    result.to_csv("oneCond.csv")
    fig.suptitle("Pulse Pressure Histogram For Thoracic Aorta Geometries")
    fig.text(0.04, 0.5, 'Density', va='center', rotation='vertical')
    fig.text(0.5, 0.04, 'Pulse Pressure (mmHg)', va='center', rotation='horizontal')
    fig.text(0.96, 0.5, 'Decreasing Resistance', va='center', rotation='vertical')
    fig.text(0.5, 0.90, 'Increasing Peak Flow', va='center', rotation='horizontal')
    plt.show()

    
def getChestFlows(flow: int):
    #Seperates out the flows by chest and abdominal CT scans
    df1 = pd.read_csv('/home/sokolor2/FlowData/Flow' + str(flow) + '.csv')
    df1['Calico ID'] = df1['Calico ID'].str.replace("_TS_part1", "")
    df1['Calico ID'] = df1['Calico ID'].str.strip()

    df2 = pd.read_csv('/home/sokolor2/Data/aortaDataChest.csv')
    df2['Calico ID'] = df2['Calico ID'].str.strip()

    #mergedDf = pd.merge(df1,df2, on='Calico ID', how = 'left', indicator=True)


    #filteredDf = mergedDf[mergedDf['_merge'] == 'left_only'].drop(columns=['_merge'] + list(df2.columns.difference(['Calico ID'])))
    selected_rows = df1[df1['Calico ID'].isin(df2['Calico ID'])]

    print(selected_rows)
    selected_rows.to_csv('/home/sokolor2/Data/flow' + str(flow) + 'Chest.csv', index=False)

def pulsePheComp(resistanceID: int, resistances, pheCodeTrue, pheCodeFalse, pheName):
    #Plots the pulse pressure histograms for different PHEcodes (density plot)

    fig, ax = plt.subplots(nrows=4, ncols=5, figsize=(50,5))
    pd.set_option('display.max_columns', None)

    for i in range(0,4):
        for j in range(0,5):
            dfFlow = pd.read_csv(f'/home/sokolor2/Data/Flow{j}.csv')
            dfFlow['Calico ID'] = dfFlow['Calico ID'].str.replace("_TS", "")

            dfTrue = pd.read_csv('/home/sokolor2/Data/' + pheCodeTrue)


            mapping_dict = dfTrue.set_index('Calico ID')['Length'].to_dict()
            dfTrue = dfFlow[dfFlow['Calico ID'].isin(dfTrue['Calico ID'])]
            dfTrue['Length'] = dfFlow['Calico ID'].map(mapping_dict)
            dfTrue.insert(0, pheName, True)

            dfFalse = pd.read_csv('/home/sokolor2/Data/' + pheCodeFalse)

            mapping_dict = dfFalse.set_index('Calico ID')['Length'].to_dict()
            dfFalse = dfFlow[dfFlow['Calico ID'].isin(dfFalse['Calico ID'])]
            dfFalse['Length'] = dfFlow['Calico ID'].map(mapping_dict)
            dfFalse.insert(0, pheName, False)

            flows = pd.concat([dfTrue, dfFalse], axis = 0)
            filteredFlow = flows[flows['Resistance'].astype(int) == resistances[i*2]].copy()



            timeCols = [f'TimePoint:{k}' for k in range(0,19)]
        
            filteredFlow.loc[:, timeCols] = filteredFlow[timeCols].apply(pd.to_numeric, errors='coerce')

            result = (
                filteredFlow[['Calico ID'] + [pheName] + timeCols + ['Length']]
                .assign(Pulse_Pressure=lambda x: x[timeCols].max(axis=1) - x[timeCols].min(axis=1))
                .reset_index()
            )
            result['Ratio'] = result['Pulse_Pressure']/result['Length']


            resultTrue = result[result[pheName] == True]['Pulse_Pressure']
            resultFalse = result[result[pheName] == False]['Pulse_Pressure']
            #sns.histplot(data=pd.DataFrame(result[result[pheName] == True]['Pulse_Pressure']), x = "Pulse_Pressure", stat="density", kde=True,ax=ax[i,j])
            #sns.histplot(data=pd.DataFrame(result[result[pheName] == False]['Pulse_Pressure']), x = "Pulse_Pressure", stat="density", kde=True,ax=ax[i,j])

            ax[i,j].hist(resultTrue, bins=30, density=True, color='blue', alpha=0.7, label=pheName + ' True', range=(10,50))
            ax[i,j].hist(resultFalse, bins=30, density=True, color='orange', alpha=0.7, label=pheName + ' False', range=(10,50))
            ax[i,j].set_ylim(0, 0.15)

            if j != 0:
                ax[i, j].set_yticklabels([])

    #result.to_csv("/home/sokolor2/Data/aValveDisRatio.csv", index=False)
    fig.legend(labels=[pheName + ' True', pheName + ' False'], loc='lower right')
    fig.suptitle("Pulse Pressure Histogram For Thoracic Aorta Geometries With and Without " + pheName + " PHECode")
    fig.text(0.04, 0.5, 'Density', va='center', rotation='vertical')
    fig.text(0.5, 0.04, 'Pulse Pressure (mmHg)', va='center', rotation='horizontal')
    fig.text(0.96, 0.5, 'Decreasing Resistance', va='center', rotation='vertical')
    fig.text(0.5, 0.90, 'Increasing Peak Flow', va='center', rotation='horizontal')
    
    plt.show()


def lengthDiamCor(pathToCSV: str):
    #Correlation plot between aortic length and max diameter
    df = pd.read_csv(pathToCSV)
    length = df['Average Diameter (mm)']
    maxDiam = df['Max Diameter (mm)']

    correlationCoef = np.corrcoef(length, maxDiam)[0,1]
    plt.scatter(length, maxDiam, color='blue', label=f'Data Points (Correlation: {correlationCoef: .2f})')

    slope, intercept = np.polyfit(length, maxDiam, 1)
    plt.plot(length, slope*length+intercept, color='red', label='Best Fit Line')
    
    plt.xlabel("Average Aortic Diameter")
    plt.ylabel("Max Aortic Diameter")
    plt.title("Correlation Between Average Diameter and Max Diameter")

    plt.legend()
    
    plt.show()

def oneHistEach(flowID: int, resistanceId: int, resistances: list):
    pd.set_option('display.max_columns', None)
    timeCols = [f'TimePoint:{k}' for k in range(0,19)]
    
    fig, ax = plt.subplots(nrows=1, ncols=4, figsize=(20,5.5))
    dfFlow = pd.read_csv(f'/home/sokolor2/Data/flow{flowID}Chest.csv')

    filteredFlow = dfFlow[dfFlow['Resistance'].astype(int) == resistances[resistanceID]].copy()
    
    result = (
                filteredFlow[['Calico ID'] + timeCols]
                .assign(Pulse_Pressure=lambda x: x[timeCols].max(axis=1) - x[timeCols].min(axis=1))
                .reset_index()
            )

##    plt.hist(result['Pulse_Pressure'], bins=30, range=(10,50))
##    plt.xlabel("Pulse Pressure (mmHg)")
##    plt.ylabel("Quantity")
##    plt.title("Pulse Pressure Distribution", fontweight='bold')
##    plt.savefig("pulseDist.png", dpi=600)
##    plt.show()

    
    aneurysm = pd.read_csv('/home/sokolor2/Data/aneurysmRatio.csv')
    otherAns = pd.read_csv('/home/sokolor2/Data/otherAnsRatio.csv')
    valveDis = pd.read_csv('/home/sokolor2/Data/valveDisRatio.csv')
    aortaValveDis = pd.read_csv('/home/sokolor2/Data/aValveDisRatio.csv')
    
    ax[0].hist(aneurysm.loc[aneurysm['Aortic Aneurysm'], 'Pulse_Pressure'],density=True, bins=30, range=(10,50), alpha=0.7)
    ax[0].hist(aneurysm.loc[aneurysm['Aortic Aneurysm']==False, 'Pulse_Pressure'],density=True, bins=30, range=(10,50), alpha=0.7)
    ax[0].set_xlabel("Pulse Pressure (mmHg)", fontweight='bold', fontsize=14)

    ax[1].hist(otherAns.loc[otherAns['Other Aneurysm'], 'Pulse_Pressure'],density=True, bins=30, range=(10,50), alpha=0.7)
    ax[1].hist(otherAns.loc[otherAns['Other Aneurysm']==False, 'Pulse_Pressure'],density=True, bins=30, range=(10,50), alpha=0.7)
    ax[1].set_xlabel("Pulse Pressure (mmHg)", fontweight='bold', fontsize=14)

    ax[2].hist(valveDis.loc[valveDis['Heart Valve Disorders'], 'Pulse_Pressure'],density=True, bins=30, range=(10,50), alpha=0.7)
    ax[2].hist(valveDis.loc[valveDis['Heart Valve Disorders']==False, 'Pulse_Pressure'],density=True, bins=30, range=(10,50), alpha=0.7)
    ax[2].set_xlabel("Pulse Pressure (mmHg)", fontweight='bold', fontsize=14)

    ax[3].hist(aortaValveDis.loc[aortaValveDis['Aortic Valve Disorders'], 'Pulse_Pressure'],density=True, bins=30, range=(10,50), alpha=0.7)
    ax[3].hist(aortaValveDis.loc[aortaValveDis['Aortic Valve Disorders']==False, 'Pulse_Pressure'],density=True, bins=30, range=(10,50), alpha=0.7)
    ax[3].set_xlabel("Pulse Pressure (mmHg)", fontweight='bold', fontsize=14)

    for i in range(0,4):
        ax[i].set_ylim(0,0.12)
    ax[0].set_ylabel("Density", fontweight='bold', fontsize=14)
    fig.suptitle("Pulse Pressure Histograms for Aortic Aneurysm, Other Aneurysm, Heart Valve Disorder, and Aortic Valve Disorder", fontweight='bold', fontsize=20)
    fig.legend(labels=['PheCode True',  'PheCode False'], loc='lower right')
    fig.text(0.17, 0.02, 'Aortic Aneurysm', va='center', rotation='horizontal', fontweight='bold', fontsize=12)
    fig.text(0.37, 0.02, 'Other Aneurysm', va='center', rotation='horizontal', fontweight='bold', fontsize=12)
    fig.text(0.56, 0.02, 'Heart Valve Disorders', va='center', rotation='horizontal', fontweight='bold',fontsize=12)
    fig.text(0.76, 0.02, 'Aortic Valve Disorders', va='center', rotation='horizontal', fontweight='bold', fontsize=12)
    plt.savefig("oneHistEach.png", dpi=600, bbox_inches='tight')
    
    plt.show()

def pressureTimeCurve(resistances:list, flowID: int):
    timeCols = [f'TimePoint:{k}' for k in range(0,19)]
    df = pd.read_csv(f"/home/sokolor2/Data/Flow{flowID}.csv")
    plt.rcParams['axes.prop_cycle'] = plt.cycler(color=['teal', 'black', 'maroon', 'darkgreen'])
    plt.figure(figsize=(12,7))
    pulses = []
    maxes = []
    mins = []
    for j in range(0,4):
        filteredFlow = df[df['Resistance'].astype(int) == resistances[j*2]].copy()
        filteredFlow = filteredFlow[timeCols]
        avgFlow = np.mean(filteredFlow, axis=0)
        stdFlow = np.std(filteredFlow, axis=0)
        pulses.append(np.max(filteredFlow,axis=0)-np.min(filteredFlow,axis=0))
        maxes.append(np.max(filteredFlow,axis = 0))
        mins.append(np.min(filteredFlow,axis=0))
        plt.errorbar(np.linspace(0,0.75,19), avgFlow, stdFlow, label=round(resistances[j*2]//1000000)*10, capsize=5, elinewidth=3)


    legend = plt.legend(title='Resistance (dynes-s/cm$^{5}$)')
    legend.get_title().set_fontweight('bold')
    legend.get_frame().set_edgecolor('black')
    legend.get_title().set_fontsize(14)
    
    plt.xlabel("Time (s)", fontweight='bold', fontsize=14)
    plt.ylabel("Pressure (mmHg)", fontweight='bold', fontsize=14)
    plt.title("Average Pressure vs Time for Varying Resistances", fontweight='bold', fontsize=16)
    fig_manager = plt.get_current_fig_manager()

        
    print(np.mean(pulses))
    print(np.std(pulses))
    print(np.max(maxes))
    print(np.min(mins))
    plt.savefig("PvT.png", dpi=600, bbox_inches='tight')
    plt.show()

def surfPlot(resistances):
    timeCols = [f'TimePoint:{k}' for k in range(0,19)]
    avgPulses = np.zeros((8,5))
    plt.figure(figsize=(12,7))
    ax = plt.gca()
    for i in range(5):
        df = pd.read_csv(f"/home/sokolor2/Data/Flow{i}.csv")
        for j in range(8):
            pulses = []
            filteredFlow = df[df['Resistance'].astype(int) == resistances[j]].copy()
            
            result = (
                filteredFlow[['Calico ID'] + timeCols]
                .assign(Pulse_Pressure=lambda x: x[timeCols].max(axis=1) - x[timeCols].min(axis=1))
                .reset_index()
            )
            avgPulses[j,i] = result['Pulse_Pressure'].mean()
    
    img = plt.imshow(avgPulses, cmap='viridis', interpolation='bilinear')
    plt.gca().invert_yaxis()
    cbar = plt.colorbar(img)
    cbar.set_label('Pulse Pressure (mmHg)', rotation=270, labelpad=15, fontsize=16)
    for i in range(len(resistances)):
        resistances[i] = round((resistances[i]//1000000)*10)

    plt.yticks(range(len(resistances)),resistances)
    xRange = [180,190,200,210,220]
    plt.xticks(range(len(xRange)), xRange)

    plt.ylabel("Resistance (dynes-s/cm$^{5}$)", fontweight='bold', fontsize=18)
    plt.xlabel("Peak Flow Rate (cm$^{3}$/s)", fontweight='bold', fontsize=18)
    plt.title("Average Pulse Pressure vs Resistance and Flow Rate", fontweight='bold', fontsize=20)

    fig_manager = plt.get_current_fig_manager()    
    plt.savefig("surfPlot.png", dpi=600)
    plt.imshow(avgPulses, cmap='viridis', interpolation='bilinear', aspect='auto')

    plt.show()

def correlationMatrix():
    df = pd.read_csv("aortaIDs.csv")
    selectedCols = ['Length (mm)', 'Max Diameter (mm)', 'Average Diameter (mm)', 'Average Curvature (mm^-1)', 'Max Curvature (mm^-1)', 'Average Torsion (mm^-1)', 'Max Torsion Magnitude (mm^-1)', 'Sex', 'Age']
    correlation_matrix = df[selectedCols].corr()
    #plt.figure(figsize(10,8))
    sns.heatmap(correlation_matrix, annot=True, cmap='coolwarm', fmt='.2f', linewidths=.5)
    plt.show()


        

resistances = [92281250, 106468800, 149031200, 163218800, 205781200, 219968800, 262531200, 276718800]
#getChestFlows(0)
resistanceID = 7
#pulseHist(resistanceID, resistances)
#pulsePheComp(0, resistances, 'aorticValveDisTrue.csv', 'aorticValveDisFalse.csv', 'Aortic Valve Disorders')
#lengthDiamCor("/home/sokolor2/aortaGeoms(V.2).csv")
#oneHistEach(2, 4, resistances)
pressureTimeCurve(resistances, 2)
#surfPlot(resistances)
#correlationMatrix()
            

