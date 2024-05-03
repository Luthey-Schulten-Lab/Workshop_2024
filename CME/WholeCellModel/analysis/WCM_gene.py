"""

Description: functions to analyze gene for CMEODE WCM
"""

import matplotlib.pyplot as plt

import numpy as np


def mapDNA(genome):
    """
    
    Return: Multi-layer dictionary containing all the genome information
    """
    DNAmap = {}
    chromosome_length = int((genome.features[0].location.end+1)/10)
    ori_ter_rotation_factor = int(chromosome_length/2)
    print(ori_ter_rotation_factor)   
    for feature in genome.features:
        strand = feature.strand     
        if strand == 1:    
            start = int(feature.location.start/10)        
            if start<=ori_ter_rotation_factor:
                start = start + ori_ter_rotation_factor           
            elif start>ori_ter_rotation_factor:               
                start = start - ori_ter_rotation_factor
            end = int(feature.location.end/10)
            if end<=ori_ter_rotation_factor:             
                end = end + ori_ter_rotation_factor
            elif end>ori_ter_rotation_factor:
                end = end - ori_ter_rotation_factor
        elif strand == -1:      
            start = int(feature.location.end/10)        
            if start<=ori_ter_rotation_factor:
                start = start + ori_ter_rotation_factor
            elif start>ori_ter_rotation_factor:
                start = start - ori_ter_rotation_factor
            end = int(feature.location.start/10)
            if end<=ori_ter_rotation_factor:            
                end = end + ori_ter_rotation_factor        
            elif end>ori_ter_rotation_factor:
                end = end - ori_ter_rotation_factor
        if feature.type == 'CDS':
            
            if('protein_id' in feature.qualifiers.keys()):
                # Excluding 3 pseudo genes
                # Totoal coding genes are 455
                locusTag = feature.qualifiers['locus_tag'][0]
                DNAmap[locusTag] = {}
                DNAmap[locusTag]['Type'] = 'protein'
                DNAmap[locusTag]['startIndex'] = [int(start)]
                DNAmap[locusTag]['originalStart'] = int(start)
                DNAmap[locusTag]['endIndex'] = [int(end)]           
                DNAmap[locusTag]['originalEnd'] = int(end)
                DNAmap[locusTag]['RNAsequence'] = str(feature.location.extract(genome.seq).transcribe())
                DNAmap[locusTag]['AAsequence'] = str(feature.location.extract(genome.seq).transcribe().translate(table=4))            
                DNAmap[locusTag]['GeneName'] = str(feature.qualifiers['product'][0])     
        elif feature.type == 'tRNA':
            locusTag = feature.qualifiers['locus_tag'][0]       
            DNAmap[locusTag] = {} 
            DNAmap[locusTag]['Type'] = 'tRNA'
            DNAmap[locusTag]['startIndex'] = [int(start)]          
            DNAmap[locusTag]['originalStart'] = int(start)
            DNAmap[locusTag]['endIndex'] = [int(end)]
            DNAmap[locusTag]['originalEnd'] = int(end) 
            DNAmap[locusTag]['RNAsequence'] = str(feature.location.extract(genome.seq).transcribe())
            DNAmap[locusTag]['GeneName'] = str(feature.qualifiers['product'][0])      
        elif feature.type == 'rRNA':           
            locusTag = feature.qualifiers['locus_tag'][0]           
            DNAmap[locusTag] = {}          
            DNAmap[locusTag]['Type'] = 'rRNA'
            DNAmap[locusTag]['startIndex'] = [int(start)]
            DNAmap[locusTag]['originalStart'] = int(start)
            DNAmap[locusTag]['endIndex'] = [int(end)]
            DNAmap[locusTag]['originalEnd'] = int(end)          
            DNAmap[locusTag]['RNAsequence'] = str(feature.location.extract(genome.seq).transcribe())
            DNAmap[locusTag]['GeneName'] = str(feature.qualifiers['product'][0])
            
        else:
            # 2 misc binding features, 2 ncRNA genes, 1 tm RNA gene, 
            try:
                locusTag = feature.qualifiers['locus_tag'][0]           
                DNAmap[locusTag] = {}          
                DNAmap[locusTag]['Type'] = feature.type
                DNAmap[locusTag]['startIndex'] = [int(start)]
                DNAmap[locusTag]['originalStart'] = int(start)
                DNAmap[locusTag]['endIndex'] = [int(end)]
                DNAmap[locusTag]['originalEnd'] = int(end)          
                DNAmap[locusTag]['RNAsequence'] = str(feature.location.extract(genome.seq).transcribe())
                DNAmap[locusTag]['GeneName'] = str(feature.qualifiers['product'][0])
            except:
                continue
            
    return DNAmap

def categorizeGenes(genome):

    geneTypes = []
    geneNumbers = []
    LocusNumstoTypes = {}

    for locusTag, locusDict in genome.items():
        locusNum = locusTag.split('_')[1]

        if locusDict['Type'] not in geneTypes:
            geneTypes.append(locusDict['Type'])
            LocusNumstoTypes[locusDict['Type']] = []

        LocusNumstoTypes[locusDict['Type']].append(locusNum)
    
    for type, locusNums in LocusNumstoTypes.items():
        geneNumbers.append(len(locusNums))

    # map(str, list) will convert the list's elements into strings.
    print('Six types of genes in Syn3A are {0} with respective numbers {1}.'.format(', '.join(geneTypes), ', '.join(map(str, geneNumbers))))

    return LocusNumstoTypes, geneTypes

def get_filamentlength_new(w):
    """
    Input: 
    
    Description: Get the trajectories of filament length of dnaA on mother and two daughter chromosomes for all replicates
    """

    time = w.get_t()
    
    reps = np.arange(1, w.N_reps+1)
    
    filas_mother = np.zeros((len(time), len(reps)))
    filas_d1 = np.zeros((len(time), len(reps)))
    filas_d2 = np.zeros((len(time), len(reps)))


    for i_rep in range(len(reps)):
        fila_mother = np.zeros(len(time))
        fila_d1 = np.zeros(len(time))
        fila_d2 = np.zeros(len(time))
        
        rep = reps[i_rep]

        for i in range(1,31):
            unbound_specie = 'ssdnaAFila_' + str(i)
            unbound2_specie = 'ssdnaAFila2_' + str(i)
            unbound3_specie = 'ssdnaAFila3_' + str(i)
            fila_mother += w.get_specie_trace_single_rep(unbound_specie, rep)*i
            fila_d1 += w.get_specie_trace_single_rep(unbound2_specie, rep)*i
            fila_d2 += w.get_specie_trace_single_rep(unbound3_specie, rep)*i
        
        filas_mother[:,i_rep] = fila_mother
        filas_d1[:,i_rep] = fila_d1
        filas_d2[:,i_rep] = fila_d2

    return filas_mother, filas_d1, filas_d2


def get_filamentlength(w):
    """
    Input: 

    Description: Get the trajectories of filament length of dnaA on mother and two daughter chromosomes for all replicates
    """

    time = w.get_t()
    
    reps = np.arange(1, w.N_reps+1)
    
    filas_mother = np.zeros((len(time), len(reps)))
    filas_d1 = np.zeros((len(time), len(reps)))
    filas_d2 = np.zeros((len(time), len(reps)))


    for i_rep in range(len(reps)):
        fila_mother = np.zeros(len(time))
        fila_d1 = np.zeros(len(time))
        fila_d2 = np.zeros(len(time))
        
        rep = reps[i_rep]

        for i in range(1,31):
            unbound_specie = 'ssDNAunboundSite_' + str(i)
            unbound2_specie = 'ssDNAunboundSite2_' + str(i)
            unbound3_specie = 'ssDNAunboundSite3_' + str(i)
            fila_mother += w.get_specie_trace_single_rep(unbound_specie, rep)*i
            fila_d1 += w.get_specie_trace_single_rep(unbound2_specie, rep)*i
            fila_d2 += w.get_specie_trace_single_rep(unbound3_specie, rep)*i
        
        filas_mother[:,i_rep] = fila_mother
        filas_d1[:,i_rep] = fila_d1
        filas_d2[:,i_rep] = fila_d2

    return filas_mother, filas_d1, filas_d2

def get_filamentlength_mother(w):
    """
    Input: 

    Description: Get the trajectories of filament length of dnaA on mother and two daughter chromosomes for all replicates
    """

    time = w.get_t()
    
    reps = np.arange(1, w.N_reps+1)
    
    filas_mother = np.zeros((len(time), len(reps)))


    for i_rep in range(len(reps)):
        fila_mother = np.zeros(len(time))
        
        rep = reps[i_rep]

        for i in range(1,31):
            unbound_specie = 'ssDNAunboundSite_' + str(i)
            unbound2_specie = 'ssDNAunboundSite2_' + str(i)
            unbound3_specie = 'ssDNAunboundSite3_' + str(i)
            fila_mother += w.get_specie_trace_single_rep(unbound_specie, rep)*i
        
        filas_mother[:,i_rep] = fila_mother

    return filas_mother


def analyze_initiation(w, reps):
    """
    Input: w
            reps: serial numbers of replicates, e.g. [1,4,9,..,12]
    Return:
            ini_rounds: [1,0,2,..,3]
            ini_times: [[800], [], [700, 5000], ..., [500, 4000, 8000]]
    Description:
    """
    
    reps = [rep - 1 for rep in reps]

    rep_num = len(reps)

    rep_str = ', '.join(str(rep) for rep in reps)

    print('Following is the initiation analysis of {0} replicates {1}.'.format(rep_num, rep_str))

    ini_check = w.get_specie_trace('RepInitCheck')
    
    ini_rounds = []
    
    ini_times = []

    for rep in reps:

        ini_check_rep = ini_check[:,rep]

        diff = ini_check_rep[1:] - ini_check_rep[:-1]
        
        locs = np.where(diff == 1)[0]
        
        ini_rounds.append(len(locs))

        ini_times.append(locs)

        time_str = ''

        for loc in locs:
            time_str = time_str + str(loc/60) + ','

        
        print('Replicate {0} finished {1} rounds of replication initiation at {2} minutes'.format(rep, len(locs), time_str))
    
    print('Average rounds of initiation is {0:.2f} over the whole cell cycle'.format(sum(ini_rounds)/len(ini_rounds)))
    
    print('Minimum rounds of replication initiation is {0}, maximum rounds is {1}'.format(min(ini_rounds), max(ini_rounds)))

    print('Distribution of initiation: 0 round {0:.2f}, 1 round {1:.2f}, 2 rounds {2:.2f}, 3 rounds {3:.2f}, 4 rounds {4:.2f}, 5 rounds {5:.2f}'
          .format(ini_rounds.count(0)/len(ini_rounds), ini_rounds.count(1)/len(ini_rounds), 
                  ini_rounds.count(2)/len(ini_rounds), ini_rounds.count(3)/len(ini_rounds),
                  ini_rounds.count(4)/len(ini_rounds),ini_rounds.count(5)/len(ini_rounds)))
    
    ini_mother_times, ini_daughter1_times, ini_daughter2_times = convert_ini_times(ini_rounds, ini_times)

    # print('Average time to finish first, second, and third round initiation are {0}, {1} and {2} seconds'.format(sum(ini_mother_times)/len(ini_mother_times), 
    #                                             sum(ini_daughter1_times)/len(ini_daughter1_times), sum(ini_daughter2_times)/len(ini_daughter2_times)))

    return ini_rounds, ini_times, ini_mother_times, ini_daughter1_times, ini_daughter2_times


def convert_ini_times(ini_rounds, ini_times):
    """
    
    Description: Convert replication initiation information into a dictionary
    """
    
    max_ini = max(ini_rounds)

    ini_mother_times = []
    ini_daughter1_times = []
    ini_daughter2_times = []

    for ini_time in ini_times:
        try:
            ini_mother_times.append(ini_time[0])
            ini_daughter1_times.append(ini_time[1])
            ini_daughter2_times.append(ini_time[2])
        except:
            None
            # print('Not enough initiation')
    
    return ini_mother_times, ini_daughter1_times, ini_daughter2_times



def analyze_filament(w, reps):
    
    filas_mother, fila_d1, filas_d2 = get_filamentlength(w)



    return None


