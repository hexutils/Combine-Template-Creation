import Mass_interference_helper_methods as mihm
import matplotlib.pyplot as plt
import Template_helper_methods
import mplhep as hep
import pandas as pd
import numpy as np
import warnings
import shutil
import uproot
import tqdm
import ROOT
import time
import copy
import os

plt.style.use(hep.style.ROOT)

class Template_creator(object):
    def __init__(self, output_directory, fname, bkgs, bkgNames, bkg_areas, lowerlim, upperlim):
        """This serves as a parent class for all other templates made

        Arguments:
            output_directory -- The directory you would like to output all your data in
            fname -- The filenames contained in all your outputs
            bkgs -- A list of background iterables. These can be lists of mass distriutions, lists of angle distributions, etc.
                NOTE: bkgs is essentially a list of lists
            bkgNames -- A list of names for each background. This should be the same dimension as bkgs.
            bkg_areas -- A list of areas for each background that you would like to scale to
            lowerlim -- The lower limit of your attribute's range (i.e. if it was phi, lowerlim would be -pi)
            upperlim -- The upper limit of your attribute's range (i.e. if it was phi, upperlim would be pi)
        """
        self.running_location = os.getcwd()
        if 'HexUtils' in self.running_location:
            self.HexUtils_path = self.running_location[:self.running_location.find('HexUtils') + 8] #8 is the number of characters in "HexUtils"
        else:
            self.HexUtils_path = ""
            warnings.warn("Not in an instance of HexUtils! This will lead to undefined behavior!")
        
        self.dimension = 0 #This attribute is set when creating a certain dimension of template
        
        
        self.output_directory = os.path.abspath(output_directory)
        if self.output_directory[-1] != '/':
            self.output_directory += '/' #make sure there is a slash at the end to ensure the shell knows it is a directory
        self.fname = fname.split('.')[0] #make sure there are no extensions in the filename
        
        self.lowerlim = lowerlim
        self.upperlim = upperlim
        self.discr_range = (0,1) #This is (0,1) for most templates, but (-1,1) for hypothesis interference templates
        
        self.signals = {} #this is a dictionary with values that look like (<raw_data>, <area desired>)
        self.scaled_signals = {} #this is a dictionary that should contain scaled histograms of your raw data.
        #so the values should look like (<counts>, <bins>)
        self.signal_weights = {} #this is a dictionary that should hold the probabilities from MELA should you choose to use them
        self.discr_signals = {} #this is a dictionary that should hold the discriminants that your signals are using
        #both signal weights and discr_weights contain ITERABLES! Each value should be an iterable of values!!
        
        ################### The dictionaries below are the versions of the above dictionaries but for background ###########
        self.bkgs = {}
        self.scaled_bkgs = {}
        self.bkg_weights = {}
        self.discr_bkgs = {}
        
        for name, bkg_sample, bkg_area in zip(bkgNames, bkgs, bkg_areas):
            self.bkgs[name] = (bkg_sample, bkg_area) #preprocessing the background samples
            
    def create_datacards(self, verbose=False, clean=True):
        """This function uses DatacardMaker_OnShell and MakeInputRoot_OnShell to generate the datacards and input for Higgs Combine

        Keyword Arguments:
            verbose -- Whether you would like verbosity (default: {False})
            clean -- Whether you would like to wipe your output folders before putting anything else in them (default: {False})
        """
        filename = self.output_directory + self.fname + ".root"
        in_folder = self.output_directory + self.fname
        out_folder = in_folder+"_out"
        
        in_folder += '/'
        out_folder += '/'
        
        if not os.path.isdir(in_folder):
            os.mkdir(in_folder)
        if clean:
            os.system("rm " + in_folder + '*')
        shutil.move(filename, in_folder)
        
        if not os.path.isdir(out_folder):
            os.mkdir(out_folder)
        if clean:
            os.system("rm " + out_folder + '*')  
        
        runstr = "python3 MakeInputRoot_OnShell.py "
        runstr += in_folder + " " + out_folder
        if not verbose:
            runstr += " > /dev/null"
        os.system(runstr)
        # print(runstr)
        
        runstr = "python3 DatacardMaker_OnShell.py "
        runstr += out_folder
        if not verbose:
            runstr += "> /dev/null"
        os.system(runstr)
        # print(runstr)
    
    def scale_and_add_bkgs(self):
        """This is a placeholder for the same function in the 1D and 2D template cases
        """
        return 
    
    def stackPlot(self, nbins=40):
        """Generates a stack plot of your samples

        Keyword Arguments:
            nbins -- The number of plots you want (default: {40})
        """
        bins = []
        weights = []
        labels = []
        for bkg_name, (bkg, area) in self.bkgs.items():
            if any(bins):
                bkg, _ = np.histogram(bkg, bins=bins, range=(self.lowerlim, self.upperlim))
            else:
                bkg, bins = np.histogram(bkg, bins=nbins, range=(self.lowerlim, self.upperlim))
            bkg = Template_helper_methods.scale(bkg, area)
            labels.append(bkg_name)
            weights.append(bkg)
        
        for sig_name, (sig, area) in self.signals.items():
            sig, _ = np.histogram(sig, bins=bins, range=(self.lowerlim, self.upperlim))
            sig = Template_helper_methods.scale(bkg, area)
            plt.cla()
            hep.histplot(weights + [sig], bins=bins, label=labels + [sig_name], stack=True, lw=3)
            plt.legend()
            plt.savefig(self.output_directory + sig_name + '_stack.png')
        
class Template_Creator_1D(Template_creator):
    """This is a class that is a parent for all 1-dimensional templates for Higgs Combine

    Arguments:
        Template_creator -- The parent class for all templates
    """
    def __init__(self, output_directory, fname, bkgs, bkgNames, bkg_areas, lowerlim, upperlim):
        """This initialization takes in all the same inputs as the parent class.

        Arguments:
            output_directory -- The directory you would like to output all your data in
            fname -- The filenames contained in all your outputs
            bkgs -- A list of background iterables. These can be lists of mass distriutions, lists of angle distributions, etc.
                NOTE: bkgs is essentially a list of lists
            bkgNames -- A list of names for each background. This should be the same dimension as bkgs.
            bkg_areas -- A list of areas for each background that you would like to scale to
            lowerlim -- The lower limit of your attribute's range (i.e. if it was phi, lowerlim would be -pi)
            upperlim -- The upper limit of your attribute's range (i.e. if it was phi, upperlim would be pi)
        """
        super().__init__(output_directory, fname, bkgs, bkgNames, bkg_areas, lowerlim, upperlim)
        self.dimension = 1 #resets the dimension to 1
        
    def scale_and_add_bkgs(self, bins=40, scaleTo=True):
        """This is the 1-dimensional version of the function. 
        It serves to bin and scale the backgrounds given to their respective areas, then add them into one histogram.

        Keyword Arguments:
            bins -- Either the number of bins you want, or a list of the bins you want (default: {40})
            scaleTo -- If true, this function will scale the backgrounds before adding them (default: {True})

        Returns:
            an overall histogram pair of (counts, bins) a la a numpy histogram
        """
        names_samples_and_areas = list(self.bkgs.items())
            
        name, (sample, area) = names_samples_and_areas[0]
        bkg_sample, bins = np.histogram(sample, range=(self.lowerlim, self.upperlim), bins=bins)
        if scaleTo:
            bkg_sample = Template_helper_methods.scale(bkg_sample, area)
            self.scaled_bkgs[name] = (bkg_sample, bins)
        overall = bkg_sample
        
        for name, (sample, area) in names_samples_and_areas[1:]:
            bkg_sample, _ = np.histogram(sample, range=(self.lowerlim, self.upperlim), bins=bins)
            if scaleTo:
                bkg_sample = Template_helper_methods.scale(bkg_sample, area)
                self.scaled_bkgs[name] = (bkg_sample, bins)
            overall += bkg_sample
        
        return overall, bins

class Template_Creator_2D(Template_creator):
    def __init__(self, output_directory, fname, bkgs, bkgNames, bkg_areas, lowerlim, upperlim):
        """This initialization takes in all the same inputs as the parent class.

        Arguments:
            output_directory -- The directory you would like to output all your data in
            fname -- The filenames contained in all your outputs
            bkgs -- A list of background iterables. These can be lists of mass distriutions, lists of angle distributions, etc.
                NOTE: bkgs is essentially a list of lists
            bkgNames -- A list of names for each background. This should be the same dimension as bkgs.
            bkg_areas -- A list of areas for each background that you would like to scale to
            lowerlim -- The lower limit of your attribute's range (i.e. if it was phi, lowerlim would be -pi)
            upperlim -- The upper limit of your attribute's range (i.e. if it was phi, upperlim would be pi)
        """
        super().__init__(output_directory, fname, bkgs, bkgNames, bkg_areas, lowerlim, upperlim)
        self.dimension = 2   
        
    def scale_and_add_bkgs(self, bins=40, scaleTo=True):
        """This is the 2-dimensional version of the function. 
        It serves to bin and scale the backgrounds given to their respective areas, then add them into one histogram.

        Keyword Arguments:
            bins -- Either the number of bins you want, or a list of the bins you want (default: {40})
            scaleTo -- If true, this function will scale the backgrounds before adding them (default: {True})

        Returns:
            an overall histogram pair of (counts, binsx, binsy) a la a numpy histogram
        """
        overall = np.ndarray((2,2))
        keys_to_follow = list(self.bkgs.keys())
        
        current_key = keys_to_follow[0]
        sample, area = self.bkgs[current_key]
        discr = self.discr_bkgs[current_key]
        
        bkg_sample, binsx, binsy = np.histogram2d(sample, discr, bins, range=[(self.lowerlim, self.upperlim), self.discr_range])
        if scaleTo:
            bkg_sample = Template_creator.scale(bkg_sample, area)
            self.scaled_bkgs[current_key] = (bkg_sample, binsx, binsy)
        overall = bkg_sample
        
        for current_key in keys_to_follow[1:]:
            sample, area = self.bkgs[current_key]
            discr = self.discr_bkgs[current_key]
            
            bkg_sample, _, _ = np.histogram2d(sample, discr, (binsx, binsy), range=[(self.lowerlim, self.upperlim), self.discr_range])
            if scaleTo:
                bkg_sample = Template_creator.scale(bkg_sample, area)
                self.scaled_bkgs[current_key] = (bkg_sample, binsx, binsy)
            overall += bkg_sample
            
        return overall, binsx, binsy 

class Interf_Coupling_template_creator(Template_Creator_2D): #WIP
    def __init__(self, output_directory, fname, bkgs, bkgNames, bkg_areas, lowerlim, upperlim,
                 pure1_weights, pure2_weights, interf_weights, weight_of_generation_hypothesis, mass_iterable,
                 interf_name, pure1_name, pure2_name,
                 bkg_pure1_weights, bkg_pure2_weights, bkg_interf_weights):
        """CURRENTLY A WORK IN PROGRESS. Designed to be a 2d template between different hypotheses

        Arguments:
            output_directory -- The directory you would like to output all your data in
            fname -- The filenames contained in all your outputs
            bkgs -- A list of background iterables. These can be lists of mass distriutions, lists of angle distributions, etc.
                NOTE: bkgs is essentially a list of lists
            bkgNames -- A list of names for each background. This should be the same dimension as bkgs.
            bkg_areas -- A list of areas for each background that you would like to scale to
            lowerlim -- The lower limit of your attribute's range (i.e. if it was phi, lowerlim would be -pi)
            upperlim -- The upper limit of your attribute's range (i.e. if it was phi, upperlim would be pi)
            pure1_weights -- The weights for the first pure sample. Should be an iterable!
            pure2_weights -- The weights for the second pure sample. Should be an iterable!
            interf_weights -- The weights for the sample of the interference between the 2. Should be an iterable!
            weight_of_generation_hypothesis -- The weight for what the sample was generated at. This is used to normalize weights.
            mass_iterable -- This is the mass iterable that you are using for one axis of the template
            interf_name -- The name of your interference sample
            pure1_name -- The name of your first pure sample
            pure2_name -- The name of your second pure sample
            bkg_pure1_weights -- The weights of your background to your first pure sample's hypothesis
            bkg_pure2_weights -- The weights of your background to your second pure sample's hypothesis
            bkg_interf_weights -- The weights of your background to your interference sample's hypothesis
        """
        super().__init__(output_directory, fname, bkgs, bkgNames, bkg_areas, lowerlim, upperlim)
        
        self.discr_range = (-1,1)
        
        weight_of_generation_hypothesis = np.array(weight_of_generation_hypothesis)
        self.signal_weights[pure1_name] = np.array(pure1_weights)/weight_of_generation_hypothesis
        self.signal_weights[pure2_name] = np.array(pure2_weights)/weight_of_generation_hypothesis
        
        self.signal_weights[interf_name] = np.array(interf_weights)/weight_of_generation_hypothesis
        self.signal_weights[interf_name] -= self.weights[pure1_name] + self.weights[pure2_name]
        
        mass_iterable = np.array(mass_iterable)
        self.signals.update(dict.fromkeys([interf_name, pure1_name, pure2_name], mass_iterable ))
        
        d_interference_signal = self.signal_weights[interf_name]
        d_interference_signal /= 2*np.sqrt(self.signal_weights[pure1_name]*self.signal_weights[pure2_name])
        
        d_interference_bkg = bkg_interf_weights
        d_interference_bkg /= 2*np.sqrt(bkg_pure1_weights*bkg_pure2_weights)
        
class Interf_Reso_template_creator_1D(Template_Creator_1D):
    def __init__(self, output_directory, fname, bkgs, bkgNames, bkg_areas, lowerlim, upperlim,
                 BW1_0_0, BW2_0_0, BW3_0_0, BW12_0_0, BW12_05_0, BW13_0_0, BW13_0_05, BW23_0_0, BW23_0_05,
                 CS_BW1, CS_BW2, CS_BW3, CS_BW12_0_0, CS_BW12_05_0, CS_BW13_0_0, CS_BW13_0_05, CS_BW23_0_0, CS_BW23_0_05,
                 nbins, area1, area2, area3):
        """_summary_

        Arguments:
            output_directory -- The directory you would like to output all your data in
            fname -- The filenames contained in all your outputs
            bkgs -- A list of background iterables. These can be lists of mass distriutions, lists of angle distributions, etc.
                NOTE: bkgs is essentially a list of lists
            bkgNames -- A list of names for each background. This should be the same dimension as bkgs.
            bkg_areas -- A list of areas for each background that you would like to scale to
            lowerlim -- The lower limit of your attribute's range (i.e. if it was phi, lowerlim would be -pi)
            upperlim -- The upper limit of your attribute's range (i.e. if it was phi, upperlim would be pi)
            BW1_0_0 -- The first pure sample
            BW2_0_0 -- The second pure sample
            BW3_0_0 -- The third pure sample
            BW12_0_0 -- The interference sample with 0 phase between 1 & 2
            BW12_05_0 -- The interference sample with pi/2 phase between 1 & 2 (hence the "05" indicating 0.5pi)
            BW13_0_0 -- The interference sample with 0 phase between 1 & 3
            BW13_0_05 -- The interference sample with pi/2 phase between 1 & 3
            BW23_0_0 -- The interference sample with 0 phase between 2 & 3
            BW23_0_05 -- The interference sample with pi/2 phase between 2 & 3
            CS_BW1 -- The cross section of the first pure sample
            CS_BW2 -- The cross section of the second pure sample
            CS_BW3 -- The cross section of the third pure sample
            CS_BW12_0_0 -- The cross section of the variable's name
            CS_BW12_05_0 -- The cross section of the variable's name
            CS_BW13_0_0 -- The cross section of the variable's name
            CS_BW13_0_05 -- The cross section of the variable's name
            CS_BW23_0_0 -- The cross section of the variable's name
            CS_BW23_0_05 -- The cross section of the variable's name
            nbins -- The number of bins you would like to use
            area1 -- The area desired of the first signal sample
            area2 -- The area desired of the second signal sample
            area3 -- The area desired of the third signal sample
        """
        super().__init__(output_directory, fname, bkgs, bkgNames, bkg_areas, lowerlim, upperlim)
        string_forms = ["BW1_0_0", "BW2_0_0", "BW3_0_0", 
                        "BW12_0_0", "BW12_05_0", "BW13_0_0", "BW13_0_05", "BW23_0_0", "BW23_0_05"]
        
        with uproot.recreate(self.output_directory + self.fname + ".root") as f:
            
            self.signals["BW1_0_0"] = (np.array(BW1_0_0), area1)
            
            BW1_0_0, bins = np.histogram(BW1_0_0, bins=nbins, range=(lowerlim, upperlim))            
            BW1_0_0 = Template_helper_methods.scale(BW1_0_0, CS_BW1)
            if np.any(BW1_0_0): #checks if the array is nonzero at any point
                temp = Template_helper_methods.scale(BW1_0_0, area1)
                f["ggH_0PM_BW1_0_0"] = (temp, bins)
                self.scaled_signals["BW1_0_0"] = (temp, bins)
            
            self.signals["BW2_0_0"] = (np.array(BW2_0_0), area2)
            BW2_0_0, _ = np.histogram(BW2_0_0, bins=bins, range=(lowerlim, upperlim))
            BW2_0_0 = Template_helper_methods.scale(BW2_0_0, CS_BW2)
            if np.any(BW2_0_0):
                temp = Template_helper_methods.scale(BW2_0_0, area2)
                f["ggH_0PM_BW2_0_0"] = (temp, bins)
                self.scaled_signals["BW2_0_0"] = (temp, bins)
            
            self.signals["BW3_0_0"] = (np.array(BW3_0_0), area3)
            BW3_0_0, _ = np.histogram(BW3_0_0, bins=bins, range=(lowerlim, upperlim))
            BW3_0_0 = Template_helper_methods.scale(BW3_0_0, CS_BW3)
            if np.any(BW3_0_0):
                temp = Template_helper_methods.scale(BW3_0_0, area3)
                f["ggH_0PM_BW3_0_0"] = (temp, bins)
                self.scaled_signals["BW3_0_0"] = (temp, bins)
            
            interfList = [BW12_0_0, BW12_05_0, BW13_0_0, BW13_0_05, BW23_0_0, BW23_0_05] #list of all the interference terms
            interfCSList = [CS_BW12_0_0, CS_BW12_05_0, CS_BW13_0_0, CS_BW13_0_05, CS_BW23_0_0, CS_BW23_0_05]
            for n, interference_term in enumerate(interfList):
                
                interference_term, _ = np.histogram(interference_term, bins=bins, range=(lowerlim, upperlim))
                interference_term = Template_helper_methods.scale(interference_term, interfCSList[n])
                
                interf_area = 0
                
                if "12" in string_forms[n + 3]:
                    interf_area = np.sqrt(area1*area2)
                    self.signals[string_forms[n+3]] = (np.array(interference_term), interf_area)
                    
                    interference_term -= BW1_0_0 + BW2_0_0
                elif "13" in string_forms[n + 3]:
                    interf_area = np.sqrt(area1*area3)
                    self.signals[string_forms[n+3]] = (np.array(interference_term), interf_area)
                    
                    interference_term -= BW1_0_0 + BW3_0_0
                else:
                    interf_area = np.sqrt(area2*area3)
                    self.signals[string_forms[n+3]] = (np.array(interference_term), interf_area)
                    
                    interference_term -= BW2_0_0 + BW3_0_0
                
                interference_term = Template_helper_methods.scale(interference_term, interf_area)
                self.scaled_signals[string_forms[n+3]] = (interference_term, bins)
                
                
                pos = np.maximum(interference_term.copy(),0) #splits the template up into positive and negative as you're supposed to
                neg = -1*np.minimum(interference_term.copy(),0)
                
                if np.any(pos):
                    f["ggH_0PM_" + string_forms[n+3] + "_pos"] = (pos, bins)
                if np.any(neg):
                    f["ggH_0PM_" + string_forms[n+3] + "_neg"] = (neg, bins)

            f["bkg_ggzz"] = self.scale_and_add_bkgs(self.bkgs, bins, scaleTo=True)
            
    def plot_overall_interference(self):
        # print(self.scaled_signals)
        pures = ([key for key in self.scaled_signals.keys() if ("13" not in key and "12" not in key and "23" not in key)], 
                 [value for key, value in self.scaled_signals.items() if ("13" not in key and "12" not in key and "23" not in key)])
        for interf12 in tqdm.tqdm(["BW12_0_0", "BW12_05_0"], desc="Top Level of interference loop"):
            for interf13 in tqdm.tqdm(["BW13_0_0", "BW13_0_05"], leave=False, desc="Second Level of interference loop"):
                for interf23 in tqdm.tqdm(["BW23_0_0", "BW23_0_05"], leave=False, desc="Bottom Level of interference loop"):
                    names, terms = copy.deepcopy(pures)
                    names += [interf12, interf13, interf23]
                    # print(names)
                    terms += [self.scaled_signals[interf12], self.scaled_signals[interf13], self.scaled_signals[interf23]]
        
                    mihm.plot_overall_interference(terms, names, 
                                                   self.output_directory, interf12+"_"+interf13+"_"+interf23+"_"+self.fname)

class Significance_Hypothesis_template_creator_1D(Template_Creator_1D):
    def __init__(self, output_directory, fname, bkgs, bkgNames, bkg_areas, lowerlim, upperlim,
                 signal1, signal1_name, signal2, signal2_name, signal_area, nbins):
        super().__init__(output_directory, fname, bkgs, bkgNames, bkg_areas, lowerlim, upperlim)

        with uproot.recreate(self.output_directory + self.fname + ".root") as f:
            self.signals[signal1_name] = signal1
            signal1, bins = np.histogram(signal1, bins=nbins, range=(lowerlim, upperlim))
            signal1 = Template_helper_methods.scale(signal1, signal_area)
            f["ggH_0PM"] = (signal1, bins)
            
            self.signals[signal2_name] = signal2
            signal2, _ = np.histogram(signal2, bins=bins, range=(lowerlim, upperlim))
            signal2 = Template_helper_methods.scale(signal2, signal_area)
            f["ggH_0M"] = (signal2, bins)
            
            f["bkg_ggzz"] = self.scale_and_add_bkgs(self.bkgs, bins, scaleTo=True)