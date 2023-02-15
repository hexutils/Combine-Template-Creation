import re
import tqdm
import uproot
import argparse
import numpy as np
import mplhep as hep
import Template_creator
import matplotlib.pyplot as plt

def place_that_list(filename):
    """Places different filenames according to their name when input. 
    Filenames should be either of form:
    BW<Reso Num>_pure if a pure sample
    BW<Reso Num>_phi_<phase>_BW<Reso Num>_phi_<phase>
    where phase is either 0 or pi_over_2

    Parameters
    ----------
    filename : string
        the filename that is being given
    """
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
    elif "_phi_0_" in filename and "_phi_0" in filename[filename.find('_phi_0')+6:]:
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
    parser.add_argument('-n', '--nbins', default=40,
                        help="The number of bins you want")
    parser.add_argument('-o', '--outFolder', default='./',
                        help="The directory you'd like to output to")
    parser.add_argument('-c', '--crossSection', required=True,
                        help="The cross section file with all your sample names inside")
    parser.add_argument('-b', '--backgrounds', nargs='+', required=True,
                        help="Your ROOT files containing the background samples")
    parser.add_argument('-ba', '--bkgAreas', nargs='+', required=True, type=float,
                        help="The areas for each background")
    parser.add_argument('-a', '--areas', nargs=3, type=float,
                        help="The areas for the three signals")
    args = parser.parse_args()
    
    """
    Cross section files should be arranged as follows:
    <absolute file path>, <cross section>, <uncertainty in cross section> (last one is optional)
    The easiest to way to generate these is with the get_cross_section_from_LHE_file function in the lhe2root repo:
    https://github.com/hexutils/lhe2root/blob/main/lhe2root_methods.py
    """
    
    if len(args.bkgAreas) != len(args.backgrounds):
        raise argparse.ArgumentError("Background argument and bkgArea arguments should be of the same length!")

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
                                                     bkg_samples.values(), bkg_samples.keys(), args.bkgAreas, 6, 9,
                                                     *list(map(data_samples.get,insertionList)),
                                                     *list(map(cross_section_samples.get, insertionList)),
                                                     args.nbins, *args.areas)
    Three_BW_Creation.create_datacards()
    Three_BW_Creation.stackPlot()
    Three_BW_Creation.plot_overall_interference()