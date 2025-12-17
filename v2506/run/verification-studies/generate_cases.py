import numpy as np
import pandas as pd
from pathlib import Path
import shutil

templateCase='template'
overviewCases='overview'
prefixFolders='case'

pathTemplate = Path.cwd() / templateCase
pathOverview = Path.cwd() / overviewCases

df = pd.read_csv(pathOverview, sep="\t")

nCases = df.shape[0]
print(f"Generating {nCases} cases from template case '{templateCase}' using overview file '{overviewCases}'")

for i in range(nCases):
    if (Path.cwd() / f"{prefixFolders}_{str(i).zfill(4)}").exists():
        shutil.rmtree(f"{prefixFolders}_{str(i).zfill(4)}")
    caseNumber = str(i).zfill(4)
    caseName = f"{prefixFolders}_{caseNumber}"
    if (Path.cwd() / caseName).exists():
        print(f"Case {caseName} already exists - doing nothing.")
    else:
        print(f"Creating case {caseName}")
        pathCase = Path.cwd() / caseName

        shutil.copytree(pathTemplate, pathCase)
        pathSettingsFile = pathCase / 'settings'
        
        parameter = {}
        parameter['R'] = df['R'][i]
        parameter['Bi'] = df['Bi'][i]
        parameter['Da'] = df['Da'][i]
        parameter['theta'] = df['theta'][i]
        parameter['zeta'] = df['zeta'][i]
        parameter['epsilon'] = df['epsilon'][i]

        parameter['n_R'] = 250
        parameter['n_z'] = 250
        parameter['n_procs'] = 4

        parameter['D_CO2_l'] = 2e-10
        parameter['H_cc'] = 3.25
        parameter['C_S'] = 1.0
        
        parameter['nu'] = 1e-6
        parameter['Uavg'] = 3.84e-5
        parameter['K_ext'] = 0.0013
        parameter['D_S'] = 4.8e-11
        parameter['C_CO2_g'] = 1.6e-05
        parameter['C_CO2_l'] = 5.2e-05
        parameter['k_rxn'] = 96.0

        parameter['z'] = 3.0*(parameter['R']/(250e-6))**2
        #parameter['Uavg'] = parameter['z']*parameter['D_CO2_l']/(parameter['R']*parameter['R']*parameter['zeta'])
        #parameter['K_ext'] = parameter['Bi']*parameter['H_cc']*parameter['D_CO2_l']/parameter['R']
        #parameter['D_S'] = parameter['theta']*parameter['D_CO2_l']
        #parameter['C_CO2_g'] = parameter['epsilon']*parameter['C_S']/parameter['H_cc']
        #parameter['C_CO2_l'] = parameter['H_cc']*parameter['C_CO2_g']
        #parameter['k_rxn'] = parameter['Da']*parameter['D_CO2_l']/(parameter['C_S']*parameter['R']*parameter['R'])


        with open(pathSettingsFile, 'a') as file:
            for key, value in parameter.items():
                file.write(f"{key} {value};\n")

                #find and replace all occurrences of {key} in all files in pathCase
                for path in pathCase.rglob('*'):
                    if path.is_file():
                        with open(path, 'r') as f:
                            filedata = f.read()
                        filedata = filedata.replace(f"{{{key}}}", str(value))
                        with open(path, 'w') as f:
                            f.write(filedata)

        