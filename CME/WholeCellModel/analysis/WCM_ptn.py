"""

Description: functions to analyze ptns for CMEODE WCM
"""
import pandas as pd
import numpy as np


def get_categories_ptn():
    """

    Return: the list of locusNums that with experimental protemomic counts larger than 10 (non-ribosomal ptns)
        ptn_list contains 338 ptns, ribosomal_ptn_list contains 52, ten_list 65
    """


    ptn_list = []
    
    ribosomal_ptn_list = []

    ten_ptn_list = []

    proteomics = pd.read_excel('/home/enguang/CMEODE/CMEODE_Hook/input_data/initial_concentrations.xlsx', sheet_name='Comparative Proteomics')
    
    for row, protein in proteomics.iterrows():  
        if row != 0:
            locusTag = protein['Locus Tag']
            locusNum = locusTag.split('_')[1]
            PtnName = protein['Gene Product']
    #         print(PtnName)
            if 'S ribosomal' not in PtnName: #exclude the ribosomal proteins
                initial_count = protein['Sim. Initial Ptn Cnt']
                if initial_count > 10:
                    ptn_list.append(locusNum)
                else:
                    ten_ptn_list.append(locusNum)
            else:
                ribosomal_ptn_list.append(locusNum)


    return ptn_list, ribosomal_ptn_list, ten_ptn_list

