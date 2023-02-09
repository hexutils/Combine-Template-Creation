import re
import tqdm
import uproot
import argparse
import numpy as np
import mplhep as hep
import Template_creator
import matplotlib.pyplot as plt

def place_that_list(filename):
    global insertionList
    if '_pure' in filename:
        if 'BW1' in filename:
            insertionList[0] = filename
        elif 'BW2' in filename:
            insertionList[1] = filename
        else: #BW3
            insertionList[2] = filename
    elif '_phi_pi_over_2' in filename:
        if 'BW2' in filename and 'BW1' in filename:
            insertionList[4] = filename
        elif 'BW2' in filename and 'BW3' in filename:
            insertionList[7] = filename
        else: #BW1 and BW3
            insertionList[6] = filename
    elif "_phi_0_" in filename and "_phi_0" in filename:
        if 'BW2' in filename and 'BW1' in filename:
            insertionList[3] = filename
        elif 'BW2' in filename and 'BW3' in filename:
            insertionList[8] = filename
        else: #BW1 and BW3
            insertionList[5] = filename
    else:
        print("whoops!")
        
if __name__ == "__main__":
    plt.style.use(hep.style.ROOT)

    parser = argparse.ArgumentParser()
    # parser.add_argument('filename')
    parser.add_argument('-n', '--nbins', default=40)
    parser.add_argument('-o', '--outFolder', default='./')
    parser.add_argument('-c', '--crossSection', required=True)
    parser.add_argument('-b', '--backgrounds', nargs='+', required=True)
    parser.add_argument('-a', '--areas', nargs=3, type=float)
    args = parser.parse_args()

    coupling_hunter = re.compile(r'\w+_ghzpzp(\d)_?\S+')

    data_samples = {}
    cross_section_samples = {}
    bkg_samples = {}
    with open(args.crossSection) as f:
        f.readline()
        for line in tqdm.tqdm(f):
            line = line.strip().split(',')
            with uproot.open(line[0]) as dataFile:
                dataFile = dataFile[dataFile.keys()[0]]
                sample = dataFile["M4L"].array(library="np")
                
                data_samples[line[0].split('/')[-1]] = sample
                cross_section_samples[line[0].split('/')[-1]] = float(line[1])
    
    insertionList = [None]*len(data_samples)
    
    for sampleName in data_samples.keys():
        place_that_list(sampleName)
    
    for bkg in args.backgrounds:
        with uproot.open(bkg) as dataFile:
            dataFile = dataFile[dataFile.keys()[0]]
            sample = dataFile["M4L"].array(library="np")
            bkg_samples[bkg.split('/')[-1].split('.')[0].split('_')[0]] = sample
    
    
    
    Three_BW_Creation = Template_creator.Interf_Reso_template_creator_1D(args.outFolder, "Mass_Template",
                                                     bkg_samples.values(), bkg_samples.keys(), [823, 5515], 6, 9,
                                                     *list(map(data_samples.get,insertionList)),
                                                     *list(map(cross_section_samples.get, insertionList)),
                                                     args.nbins, *args.areas)
    Three_BW_Creation.create_datacards()
    Three_BW_Creation.stackPlot()
    Three_BW_Creation.plot_overall_interference()