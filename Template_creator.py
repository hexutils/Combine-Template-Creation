import matplotlib.pyplot as plt
import mplhep as hep
import pandas as pd
import numpy as np
import shutil
import uproot
import ROOT
import time
import os

plt.style.use(hep.style.ROOT)

def normalize_to_area(sample1, area, nbins, lowerlim, upperlim, dimension=1, sample2=None, lowerlim2=None, upperlim2=None):
    """normalizes a sample to be histogrammed to a certain area (This is a different version of the scale function)

    Arguments:
        sample1 -- the first sample
        area -- the area you want the histogram to be
        nbins -- the number of bins you want. This can either be a number or a list of bins
        lowerlim -- the lower limit of your histogram range
        upperlim -- the upper limit of your histogram range

    Keyword Arguments:
        dimension -- If this is 2 then do a 2d histogram plot (default: {1})
        sample2 -- If the dimension is 2 then you need 2 samples! (default: {None})
        
    Raises:
        ValueError: If you have a 2d histogram you need 2 samples!

    Returns:
        A numpy histogram tuple that has been properly normalized to a specific area. 
    """
    if dimension == 1:
        normed_counts, normed_bins = np.histogram(sample1, bins=nbins, range=(lowerlim,upperlim))
        normed_counts = normed_counts.astype(float)
        signs = np.sign(normed_counts) #keep the signs for future reference
        normed_counts = np.abs(normed_counts) #take the absolute value in case your histogram goes into the negatives
        
        normed_counts = signs*normed_counts*area/np.sum(normed_counts)
        return normed_counts, normed_bins
    else:
        if sample2 == None or lowerlim2 == None or upperlim2 == None:
            raise ValueError("Second sample variables 'sample2', 'lowerlim2', and 'upperlim2' required for 2d histogram!")
        
        normed_counts, *normed_bins = np.histogram2d(sample1, sample2, bins=nbins, range=( (lowerlim,upperlim),(lowerlim2,upperlim2) ))
        normed_counts = normed_counts.astype(float)
        signs = np.sign(normed_counts) #keep the signs for future reference
        normed_counts = np.abs(normed_counts) #take the absolute value in case your histogram goes into the negatives
        normed_counts = signs*normed_counts*area/np.sum(normed_counts)
        return normed_counts, normed_bins

class Template_creator(object):
    def __init__(self, g1sample, g4sample, signal_area, fname, d0_g1, d0_g4, d0_bkg, bkgs, bkgNames, discr_name,
                 lowerlim=6, upperlim=9, reweightFactors = None, signal_sample_region=None):
        """This class is a class to generate templates of both 1 dimension and 2 dimensions to use in Higgs Combine

        Arguments:
            g1sample -- This is the "first" signal sample. It is known as the "g1" sample simply because that is what it was known during development
                NOTE: The format of this input should be (sample, name)
                NOTE: The sample should be an iterable of masses
                
            g4sample -- This is the "second" signal sample. It is known as the "g4" sample simply because that is what it was known during development
                NOTE: The format of this input should be (sample, name)
                NOTE: The sample should be an iterable of masses
                
            signal_area -- The desired area of the signal m4L histogram.
            fname -- The base filename you want to use when dumping data into the templates
            d0_g1 -- The discriminant for your "first" signal sample. 
                NOTE: Should be an iterable of values
            d0_g4 -- The discriminant for your "second" signal sample. 
                NOTE: Should be an iterable of values
            d0_bkg -- The discriminant list for your background sample. Should be the same length as your backgrounds. 
                NOTE: Should be a list of value iterables
            bkgs -- The list of background samples for your template.
                NOTE: This should be a list of mass iterables
            bkgNames -- The names of each background for easy naming
            discr_name -- The name of the discriminant you are using in a way that looks pretty for matplotlib (i.e. D_{0-})

        Keyword Arguments:
            lowerlim -- The lower limit you want to use for your histograms (default: {6})
            upperlim -- The upper limit you want to use for your histograms (default: {9})
            reweightFactors -- currently a placeholder for reweighting in the future (default: {None})
            signal_sample_region -- The region that you want to "oversample," or more precisely sample compared to the overall (default: {None})
                NOTE: if this is kept at the default, it will distribute the bins evenly
                NOTE: The format of this argument should be in the form (<lower end of region>, <upper end of the region>)

        Raises:
            ValueError: If the background arguments are not the same length then it cannot proceed!
        """
        
        self.g1Sample, self.g1Name = g1sample
        self.g4Sample, self.g4Name = g4sample
        self.signal_area = signal_area
        
        if '.' in fname:
            fname = fname.split('.')[0] #DO NOT ACCEPT FILE EXTENSIONS
        self.fname = fname
        
        self.d0_g1 = d0_g1
        self.d0_g4 = d0_g4
        
        self.lowerlim = lowerlim
        self.upperlim = upperlim
        self.no_g1 = True
        self.no_g4 = True
        
        self.discr_name = discr_name
        
        for i in range(len(bkgs)):
            bkgs[i] = (bkgs[i][0].resonance_mags.to_numpy(), bkgs[i][1])
        
        self.bkgs_and_their_areas = bkgs #backgrounds should be in form [sample, area]
        
        self.bkgNames = bkgNames #Name your goddamn backgrounds!!!
        
        self.d0_bkg = d0_bkg #d0 for backgrounds should be a list in the same order
        
        if len(bkgs) != len(bkgNames) or len(bkgs) != len(d0_bkg):
            raise ValueError("Background elements should be the same length!")
        
        
        self.g4_is_reweighted = False
        
        self.bkg_1d_hist = None
        self.g1_1d_hist = None
        self.g4_1d_hist = None
        
        self.bkg_2d_hist = None
        self.g1_2d_hist = None
        self.g4_2d_hist = None
        
        # if g4sample == None:
        #     if reweightFactors == None:
        #         raise ValueError("Needs a set of reweight factors in order to reweight!")
        #     self.g4_is_reweighted = True
        
        # if g1sample == None:
        #     if reweightFactors == None:
        #         raise ValueError("Needs a set of reweight factors in order to reweight")
        #     if g4sample == None:
        #         raise ValueError("One of these samples must be a JPsiRoot class in order to reweight properly!")
        #     self.g1_is_reweighted = True
        
        # if signal_sample_region != None:
        #     if signal_sample_region[0] < self.lowerlim or signal_sample_region[1] > self.upperlim:
        #         raise ValueError("Signal region must be a subset of the overall region!")
        #     if not isinstance(signal_sample_region, Iterable) or len(signal_sample_region) != 2:
        #         raise ValueError("Signal region must be an iterable of form (lowerlim, upperlim)!")
        
        self.signal_region = signal_sample_region
    
    def plot_histograms(self, fname='', show=True):
        """Plots all the histograms that are created when you do rooFit. 
        If you haven't run the 1d input yet it wil do that for you using default parameters

        Keyword Arguments:
            fname -- the filename. If it is not '' then it will save the plots at that path (default: {''})
        """
        if '.' in fname:
            fname = fname.split('.')[0] #DO NOT ACCEPT FILE EXTENSIONS
            
        if self.bkg_1d_hist == None:
            print("Running rooFit_Input_1d with default parameters...")
            self.rooFit_Input_1d(save=False)
            
        if self.bkg_2d_hist == None:
            print("Running rooFit_Input_2d with default parameters...")
            self.rooFit_Input_2d(save=False)
            
        plotList = [self.g1_2d_hist, self.g4_2d_hist]
        plotNames = ['g1', 'g4']
        for n, i in enumerate(self.bkg_2d_hist):
            plotList.append(i)
            plotNames.append(self.bkgNames[n])
        
        plotList.append(self.bkg_overall_2d)
        plotNames.append('overall_bkg')
        
        for sample, name in zip(plotList, plotNames):
            plt.figure()
            
            hep.hist2dplot(sample)
            plt.xlabel(r'$m_{4\mu} (GeV)$', horizontalalignment='center', fontsize=40)
            plt.ylabel(self.discr_name, horizontalalignment='center', fontsize=40)
            plt.tight_layout()
            if fname != '':
                plt.savefig('plots/'+'2d_'+fname+'_'+name+'.png')
                plt.savefig('plots/'+'2d_'+fname+'_'+name+'.pdf')
            if show:
                plt.show()
            
        plotList = [self.g1_1d_hist, self.g4_1d_hist]
        plotNames = ['g1', 'g4']
        for n, i in enumerate(self.bkg_1d_hist):
            plotList.append(i)
            plotNames.append(self.bkgNames[n])
                
        plotList.append(self.bkg_overall_1d)
        plotNames.append('overall_bkg')
        
        for sample, name in zip(plotList, plotNames):
            plt.figure()

            hep.histplot(sample)
            plt.xlabel(r'$m_{4\mu} (GeV)$', horizontalalignment='center', fontsize=40)
            plt.xlim(6,9)
            plt.tight_layout()
            if fname != '':
                plt.savefig('plots/'+fname+'_'+name+'.png')
                plt.savefig('plots/'+fname+'_'+name+'.pdf')
            if show:
                plt.show()
    
    def stackPlot(self, nbins=100, fname='', show=True):
        """_summary_

        Keyword Arguments:
            nbins -- The number of bins in the stackplot. This is greater than the actual number of
                bins to make it prettier (default: {100})
            fname -- The filename. If it is not '' the plot will be saved with that name (default: {''})
        """
        
        if '.' in fname:
            fname = fname.split('.')[0] #DO NOT ACCEPT FILE EXTENSIONS
            
        bkgs = []
        for bkg, area in self.bkgs_and_their_areas:
            count, _ = normalize_to_area(bkg, area, nbins, self.lowerlim, self.upperlim)
            bkgs.append(count)
            
        plt.figure()
        sig_counts, bins = normalize_to_area(self.g1Sample.resonance_mags.to_numpy().reshape(-1),
                                             self.signal_area, nbins, self.lowerlim, self.upperlim)
        
        
        plt.stackplot(list(bins[:-2]) + [bins[-1]], *bkgs, sig_counts, labels=list(self.bkgNames) + [self.g1Name + ' Signal'])
        plt.xlim([self.lowerlim, self.upperlim])
        plt.xlabel(r'$m_{4\mu}$' + ' (GeV)', fontsize=40, horizontalalignment='center')
        plt.legend(loc='upper right')
        
        if fname:
            plt.savefig('plots/'+fname+'_g1.png')
            plt.savefig('plots/'+fname+'_g1.pdf')
            print('plots dumped to', 'plots/'+fname+'_g1(.pdf/.png)')
        if show:
            plt.show()
        
        plt.figure()
        sig_counts, bins = normalize_to_area(self.g4Sample.resonance_mags.to_numpy().reshape(-1),
                                             self.signal_area, nbins, self.lowerlim, self.upperlim)
        plt.stackplot(list(bins[:-2]) + [bins[-1]], *bkgs, sig_counts, labels=list(self.bkgNames) + [self.g4Name + ' Signal'])
        plt.xlim([self.lowerlim, self.upperlim])
        plt.xlabel(r'$m_{4\mu}$' + ' (GeV)', fontsize=40, horizontalalignment='center')
        plt.legend(loc='upper right')
        
        if fname:
            plt.savefig('plots/'+fname+'_g4.png')
            plt.savefig('plots/'+fname+'_g4.pdf')
            print('plots dumped to', 'plots/'+fname+'_g4(.pdf/.png)')
        if show:
            plt.show()
        
        # plt.figure()
        # plotList = [np.histogram(self.d0_g1, bins=10, range=(0,1)),
        #             np.histogram(self.d0_g4, bins=10, range=(0,1))]
        # for i in self.d0_bkg:
        #     plotList = [np.histogram(i, bins=10, range=(0,1))] + plotList
            
        # names = self.bkgNames + [r'$D_{0-}(0^+)$', r'$D_{0-}(0^-)$']
        
        # hep.histplot(plotList, stack=True, label=names)
        # plt.legend()
        
        plt.figure()
        bkgs = []
        sig_counts_g1, bins = normalize_to_area(self.d0_g1, 
                                                self.signal_area, nbins, 0, 1)
        for n, bkg in enumerate(self.d0_bkg):
            count, _ = normalize_to_area(bkg, 
                                         self.bkgs_and_their_areas[n][1],
                                         bins, 0, 1)
            bkgs.append(count)
        
        sig_counts_g4, _ = normalize_to_area(self.d0_g4, 
                                            self.signal_area, bins, 0, 1)
        
        # plt.stackplot(list(bins[:-2]) + [bins[-1]], *bkgs, sig_counts_g1, colors=['blue', 'orange', 'green'])
        # plt.stackplot(list(bins[:-2]) + [bins[-1]], *bkgs, sig_counts_g4, colors=['blue', 'orange', 'red'])
        # bkg1 = mpl.patches.Patch(color='blue', label=self.bkgNames[0])
        # bkg2 = mpl.patches.Patch(color='orange', label=self.bkgNames[1])
        # sig_g1 = mpl.patches.Patch(color='green', label=r'$D_{0-} 0^+$')
        # sig_g4 = mpl.patches.Patch(color='red', label=r'$D_{0-} 0^-$')
        # plt.legend(handles=[bkg1, bkg2, sig_g4, sig_g1])
        
        plt.stackplot(list(bins[:-2]) + [bins[-1]], *bkgs, sig_counts_g1, sig_counts_g4, 
                      labels=list(self.bkgNames) + [self.g1Name, self.g4Name])
        plt.xlim(0,1)
        plt.xlabel(r'$D_{0-}$')
        plt.legend()
        if fname:
            plt.savefig('plots/'+fname+'_d0_stack.png')
            plt.savefig('plots/'+fname+'_d0_stack.pdf')
            print('plots dumped to', 'plots/'+fname+'_g4(.pdf/.png)')
        if show:
            plt.show()
    
    def sample_signal_region(self, data, area, nbins_in_signal_region, nbins_overall):
        """Takes the region that is indicated to be sampled more heavily and creates a histogram with unequal binning

        Arguments:
            data -- The data to bin
            area -- The overall area of the signal desired
            nbins_in_signal_region -- How many bins you would like in the signal region
            nbins_overall -- How many bins you would like overall

        Returns:
            A numpy histogram of counts and bins matching the given area requested sampling the signal region more heavily
        """
        nbins_in_signal_region = int(nbins_in_signal_region)
        n_non_signal_bins = (nbins_overall - nbins_in_signal_region)
        
        signal_low, signal_high = self.signal_region
        overall_low, overall_high = self.lowerlim, self.upperlim
        
        prefactor = n_non_signal_bins//((signal_low - overall_low) + (overall_high - signal_high))
        
        n_non_signal_low_bins = int((signal_low - overall_low)*prefactor)
        n_non_signal_high_bins = int((overall_high - signal_high)*prefactor)
        
        if nbins_in_signal_region + n_non_signal_low_bins + n_non_signal_high_bins != nbins_overall:
            nbins_in_signal_region += (nbins_overall - nbins_in_signal_region - n_non_signal_low_bins - n_non_signal_high_bins)
            nbins_in_signal_region = int(nbins_in_signal_region)
                        
        non_signal_region_low_hist = ([], [])
        signal_region_hist = ([], [])
        non_signal_region_high_hist = ([],[])
        
        if n_non_signal_low_bins != 0:
            non_signal_region_low_hist = np.histogram(data, bins=n_non_signal_bins, range=(overall_low, signal_low))
            
        signal_region_hist = np.histogram(data, bins=nbins_in_signal_region, range=self.signal_region)
        
        if n_non_signal_high_bins != 0:
            non_signal_region_high_hist = np.histogram(data, bins=n_non_signal_bins, range=(signal_high, overall_high))
        
        bins = np.array(
            list(non_signal_region_low_hist[1][:-1]) + list(signal_region_hist[1][:-1]) + list(non_signal_region_high_hist[1])
        )
        counts = np.array(
            list(non_signal_region_low_hist[0]) + list(signal_region_hist[0]) + list(non_signal_region_high_hist[0]),
            dtype=float
        )
        counts *= area/np.sum(counts)
        
        return counts, bins
    
    def sample_signal_region_2d(self, datax, datay, area, nbins_in_signal_region, nbins_overall):
        """A 2 dimensional version of sample_signal_region

        Arguments:
            datax -- The x data. This is the data that will be sampled at the range!
            datay -- The y data. This data's binning will remain unchanged
            area -- The overall area of the signal desired
            nbins_in_signal_region -- How many bins you would like in the signal region
            nbins_overall -- How many bins you would like overall

        Returns:
            A numpy histogram of counts and bins matching the given area requested sampling the signal region more heavily
        """
        if self.d0_bkg == [None] or self.d0_g1 == None or self.d0_g4 == None:
            raise ValueError("1 Dimensional input to 2 dimensional template!")
        
        nbins_in_signal_region = int(nbins_in_signal_region)
        n_non_signal_bins = (nbins_overall - nbins_in_signal_region)
        
        
        signal_low, signal_high = self.signal_region
        overall_low, overall_high = self.lowerlim, self.upperlim
        
        prefactor = n_non_signal_bins//((signal_low - overall_low) + (overall_high - signal_high))
        
        n_non_signal_low_bins = int((signal_low - overall_low)*prefactor)
        n_non_signal_high_bins = int((overall_high - signal_high)*prefactor)
        
        if nbins_in_signal_region + n_non_signal_low_bins + n_non_signal_high_bins != nbins_overall:
            nbins_in_signal_region += (nbins_overall - nbins_in_signal_region - n_non_signal_low_bins - n_non_signal_high_bins)
            nbins_in_signal_region = int(nbins_in_signal_region)
        
        non_signal_region_low_hist = ([], [])
        signal_region_hist = ([], [])
        non_signal_region_high_hist = ([],[])
        
        if n_non_signal_low_bins != 0:
            non_signal_region_low_hist = np.histogram2d(datax, datay, bins=(n_non_signal_bins, nbins_overall), range=((overall_low, signal_low), (0,1)))
            
        signal_region_hist = np.histogram2d(datax, datay, bins=(nbins_in_signal_region, nbins_overall), range=(self.signal_region, (0,1)))
        
        if n_non_signal_high_bins != 0:
            non_signal_region_high_hist = np.histogram2d(datax, datay, bins=(n_non_signal_bins, nbins_overall), range=((signal_high, overall_high), (0,1)))
        
        
        binsx = np.array(
            list(non_signal_region_low_hist[1][:-1]) + list(signal_region_hist[1][:-1]) + list(non_signal_region_high_hist[1])
        )
        
        counts = np.array(
            list(non_signal_region_low_hist[0]) + list(signal_region_hist[0]) + list(non_signal_region_high_hist[0])
        )
        
        counts *= area/np.sum(counts)
        
        return counts, binsx, signal_region_hist[2] #binsy is unchanged when sampling for the mass

    def rooFit_Input_1d(self, nbins=30, g1SampleName='ggH_0PM', g4SampleName='ggH_0M', save=True):
        """This is a function for getting a simple resonance mass template for two samples: 0+ and 0- as signals,
        and an overall background sample

        Keyword Arguments:
            nbins -- The number of bins you want each histogram to have (default: {30})
            g1SampleName -- The name that you want to name your g1 sample histogram (default: {'ggH_0PM'})
            g4SampleName -- The name that you want to name your g4 sample histogram (default: {'ggH_0M'})
            save -- Whether you would like to save these templates (default: {True})

        Returns:
            None
        """        
        
        bkg_counts, bkg_bins = normalize_to_area(*self.bkgs_and_their_areas[0],
                                                 nbins, self.lowerlim, self.upperlim)
        # bkg_counts *= self.bkgs_and_their_areas[0][1]*np.diff(bkg_bins)
        
        tot_bkg_list = [(bkg_counts, bkg_bins.copy())]
        if self.signal_region != None:
            
            bkg_counts, bkg_bins = self.sample_signal_region(*self.bkgs_and_their_areas[0],nbins*3//4,nbins )
            tot_bkg_list = [(bkg_counts, bkg_bins.copy())]
        
        for bkg, area in self.bkgs_and_their_areas[1:]:
            if self.signal_region != None:
                tot_bkg_list.append(self.sample_signal_region(bkg, area, nbins*3//4, nbins))
                continue
            
            temp_bkg_counts, _ = normalize_to_area(bkg, area, 
                                                   bkg_bins, self.lowerlim, self.upperlim)
            # temp_bkg_counts *= area*np.diff(bkg_bins)
            tot_bkg_list.append( (temp_bkg_counts, bkg_bins.copy()) )
            
            
        tot_bkg = sum(i for i,j in tot_bkg_list)
        
        self.bkg_1d_hist = tot_bkg_list
        bkg = (tot_bkg, bkg_bins)
        self.bkg_overall_1d = bkg
        
        for (colname, data_g1), (_, data_g4) in zip(self.g1Sample,
                                 self.g4Sample): #for best accuracy, loop over resonance mags and not tot_vec!
            
            with uproot.recreate('fit_files/' + self.fname + '.root') as f:
                if self.signal_region != None:
                    sig_proper = data_g1
                    sig_counts, _ = self.sample_signal_region(sig_proper, self.signal_area, nbins*3//4, nbins)
                    self.g1_1d_hist = (sig_counts, bkg_bins)
                    if save:
                        # print(self.g1_1d_hist, len(self.g1_1d_hist[0]), len(self.g1_1d_hist[1]))
                        f[g1SampleName] = self.g1_1d_hist
                    
                    sig_proper = data_g4
                    sig_counts, _ = self.sample_signal_region(sig_proper, self.signal_area, nbins*3//4, nbins)
                    self.g4_1d_hist = (sig_counts, bkg_bins)
                    if save:
                        f[g4SampleName] = self.g4_1d_hist
                        f['bkg_ggzz'] = bkg
                    continue
                    
                sig_proper = data_g1
                sig_counts, _ = normalize_to_area(sig_proper, self.signal_area,
                                                  bkg_bins, self.lowerlim, self.upperlim)
                
                # sig_counts = sig_counts*self.signal_area*np.diff(bkg_bins) 
                
                if save:
                    f[g1SampleName] = (sig_counts, bkg_bins)
                self.g1_1d_hist = (sig_counts, bkg_bins)
                
                sig_proper = data_g4
                sig_counts, _ = normalize_to_area(sig_proper, self.signal_area, 
                                                  bkg_bins, self.lowerlim, self.upperlim)
                # sig_counts = sig_counts*self.signal_area*np.diff(bkg_bins)
                
                if save:
                    f[g4SampleName] = (sig_counts, bkg_bins)
                self.g4_1d_hist = (sig_counts, bkg_bins)
                if save:
                    f['bkg_ggzz'] = bkg
                
    def Unroll_2D_OnShell(self):
        """Code written by Jeffrey Davis of happy hour cocktail fame to unroll a 2 dimensional histogram
        """
        
        histfile = ROOT.TFile.Open('fit_files/'+self.fname+'_2d.root', "READ")
        # print('fit_files/'+self.fname+'_2d_out.root')
        fout = ROOT.TFile('fit_files/'+self.fname+'_2d_out.root',"RECREATE")
        fout.cd()
        
        for keyname in ['ggH_0PM', 'ggH_0M', 'bkg_ggzz']:
            hist = histfile.Get(keyname)
            
            xbins = hist.GetNbinsX()
            ybins = hist.GetNbinsY()

            temp_pos = ROOT.TH1F("temp_pos","",xbins*ybins,0,xbins*ybins)
            temp_neg = ROOT.TH1F("temp_neg","dif",xbins*ybins,0,xbins*ybins)

            #Unroll Hists

            indk = 0
            has_negative = False 
            for y in range (1,ybins+1):
                for x in range (1,xbins+1):
                    binx_c = hist.GetXaxis().GetBinCenter(x)
                    biny_c = hist.GetYaxis().GetBinCenter(y)
                    ibin =  hist.FindBin(binx_c,biny_c)
                    cont  = hist.GetBinContent(ibin)
                    #put small values in empty background bins
                    if cont == 0 : 
                        if "bkg" in hist.GetName():
                            intt = hist.Integral()
                            nb = ybins*xbins
                            contt = 0.1*intt*1.0/nb
                            # print ("found empty bin",contt)
                            hist.SetBinContent(ibin,contt)
                            # print (cont)
                    if cont  < 0 :
                        has_negative = True
                        
            for y in range (1,ybins+1):
                for x in range (1,xbins+1):
                    binx_c = hist.GetXaxis().GetBinCenter(x)
                    biny_c = hist.GetYaxis().GetBinCenter(y)
                    ibin =  hist.FindBin(binx_c,biny_c)
                    cont  = hist.GetBinContent(ibin)
                    if cont  < 0 :
                        temp_neg.Fill(indk,-1*cont)
                    else :
                        temp_pos.Fill(indk,cont)
                        temp_pos.SetBinError(indk, np.sqrt(cont))
                    indk = indk +1

            temp_name = hist.GetName()
            
            tpname = temp_name
            tnname = temp_name

            if (has_negative and ( "bkg" in tnname or "Data" in tnname  or "0PH" in tnname or "0PM" in tnname or "L1" in tnname or "0M" in tnname)):
                for y in range (1,ybins+1):
                    for x in range (1,xbins+1):
                        binx_c = hist.GetXaxis().GetBinCenter(x)
                        biny_c = hist.GetYaxis().GetBinCenter(y)
                        ibin =  hist.FindBin(binx_c,biny_c)
                        cont  = hist.GetBinContent(ibin)

                    #put small values in negtative background bins
                    #Also put 0 in negative signal bins
                    if cont  < 0 :
                        hist.SetBinContent(ibin,0)
                        print ("found negative bin",cont)
                        cont = 0
                    if cont == 0 :
                        if "bkg" in hist.GetName():
                            intt = hist.Integral()
                            nb = ybins*xbins
                            contt = 0.1*intt*1.0/nb
                            print ("found empty bin",contt)
                            hist.SetBinContent(ibin,contt)
                            print (cont)
                
                temp_neg.SetName(tnname)
                temp_pos.SetName(tpname)

            elif (has_negative or not ( "bkg" in tnname or "Data" in tnname  or "0PH" in tnname or "0PM" in tnname or "L1" in tnname or "0M" in tnname) ):

                if "up" in tpname or "dn" in tpname :
                    tpnm = tpname.split("_")
                    tpnm.insert(2,"positive")
                    tpname= tpnm[0]
                    for ist in range(1,len(tpnm)):
                        tpname = tpname+"_"+tpnm[ist]  

                else :     
                    tpname = tpname+"_positive"


                if "up" in tnname or "dn" in tnname :
                    tnnm = tnname.split("_")
                    tnnm.insert(2,"negative")
                    tnname= tnnm[0]
                    for ist in range(1,len(tnnm)):
                        tnname = tnname+"_"+tnnm[ist]  
                else :     
                    tnname = tnname+"_negative"

                    
                temp_neg.SetName(tnname)
                temp_pos.SetName(tpname)

            else:
            
                tnname = tnname.replace("0Xff_","0Mff_")
                tpname = tpname.replace("0Xff_","0Mff_")   
                
                temp_neg.SetName(tnname)
                temp_pos.SetName(tpname)

            if "data" in  tnname or "Data" in tnname : 
                
                temp_neg.SetName("data_obs")
                temp_pos.SetName("data_obs")
            
            
            
            if temp_pos.Integral() > 0:
                temp_pos.Write()
            if temp_neg.Integral() > 0:
                temp_neg.Write()
                
            # print('Dumped Histogram into fit_files/'+self.fname+'_2d_out.root')

    def rooFit_Input_2d(self, nbins=10, g1SampleName='ggH_0PM', g4SampleName='ggH_0M', save=True):
        """This is a function for getting a 2d resonance mass/discriminant template for two samples: 0+ and 0- as signals,
        and an overall background sample

        Keyword Arguments:
            nbins -- The number of bins you want each histogram to have (default: {30})
            g1SampleName -- The name that you want to name your g1 sample histogram (default: {'ggH_0PM'})
            g4SampleName -- The name that you want to name your g4 sample histogram (default: {'ggH_0M'})
            save -- Whether you would like to save these templates (default: {True})

        Returns:
            None
        """    
        
        bkg_proper = self.bkgs_and_their_areas[0]
        
        bkg_d0 = self.d0_bkg[0]
        bkg_counts, *bkg_bins = np.histogram2d(bkg_proper[0].reshape(-1), 
                                               bkg_d0.reshape(-1), 
                                               density=True, bins=nbins, range=( (self.lowerlim,self.upperlim),(0,1) ))
        
        bkg_counts *= bkg_proper[1]/np.sum(bkg_counts)
        tot_bkg_list = [(bkg_counts, bkg_bins.copy())]
        
        if self.signal_region != None:
            
            bkg_counts, *bkg_bins = self.sample_signal_region_2d(bkg_proper[0].reshape(-1),
                                                                 bkg_d0.reshape(-1),
                                                                 bkg_proper[1],
                                                                 nbins*3//4,nbins )
            tot_bkg_list = [(bkg_counts, bkg_bins.copy())]


        for n, (bkg, area) in enumerate(self.bkgs_and_their_areas[1:]):
            if self.signal_region != None:
                tot_bkg_list.append(self.sample_signal_region_2d(bkg.reshape(-1), 
                                                                 self.d0_bkg[n+1].reshape(-1),
                                                                 area, 
                                                                 nbins*3//4, nbins))
                continue
            temp_bkg_counts, *_ = np.histogram2d(bkg.reshape(-1), 
                                                    self.d0_bkg[n+1].reshape(-1), 
                                                    density=True, bins=nbins, range=( (self.lowerlim,self.upperlim),(0,1) ))
            
            temp_bkg_counts *= area/np.sum(temp_bkg_counts)
            
            
            tot_bkg_list.append((temp_bkg_counts, bkg_bins.copy()))
        
        tot_bkg = sum(i for i,*j in tot_bkg_list)
        self.bkg_2d_hist = tot_bkg_list        
        bkg = (tot_bkg, bkg_bins)
        self.bkg_overall_2d = bkg
        
        for (colname, data_g1), (_, data_g4) in zip(self.g1Sample,
                                 self.g4Sample):
            
            # if os.path.isfile('fit_files/' + self.fname + '_2d.root'):
            #     print("File already exists!")
            #     return True
            
            with uproot.recreate('fit_files/' + self.fname + '_2d.root') as f:
                if self.signal_region != None:
                    sig_proper = data_g1#.apply(four_vector.extractMass)
                    sig_d0 = self.d0_g1
                    sig_counts, *_ = self.sample_signal_region_2d(sig_proper, 
                                                                   sig_d0,
                                                                   self.signal_area, 
                                                                   nbins*3//4, nbins)
                    self.g1_2d_hist = (sig_counts, bkg_bins)
                    if save:
                        f[g1SampleName] = self.g1_2d_hist
                    
                    sig_proper = data_g4#.apply(four_vector.extractMass)
                    sig_d0 = self.d0_g4
                    sig_counts, *_ = self.sample_signal_region_2d(sig_proper, 
                                                                   sig_d0,
                                                                   self.signal_area, 
                                                                   nbins*3//4, nbins)
                    self.g4_2d_hist = (sig_counts, bkg_bins)
                    if save:
                        f[g4SampleName] = self.g4_2d_hist
                        f['bkg_ggzz'] = bkg
                    continue
                
                sig_proper = data_g1#.apply(four_vector.extractMass)
                sig_d0 = self.d0_g1
                sig_counts, *sig_bins = np.histogram2d(sig_proper, sig_d0, density=True, bins=nbins, range=( (self.lowerlim,self.upperlim),(0,1) ) )
                sig_counts = self.signal_area/np.sum(sig_counts)
                
                if save:
                    f[g1SampleName] = (sig_counts, bkg_bins)
                
                self.g1_2d_hist = (sig_counts, bkg_bins)
                
                sig_proper = data_g4#.apply(four_vector.extractMass)
                
                sig_d0 = self.d0_g4
                sig_counts, *sig_bins = np.histogram2d(sig_proper, sig_d0, density=True, bins=nbins, range=( (self.lowerlim,self.upperlim),(0,1) ) )
                sig_counts = sig_counts*self.signal_area*np.diff(bkg_bins[0])*np.diff(bkg_bins[1])
                
                if save:
                    f[g4SampleName] = (sig_counts, bkg_bins)
                self.g4_2d_hist = (sig_counts, bkg_bins)
                
                if save:
                    f['bkg_ggzz'] = bkg
                    
        time.sleep(2) #this line is necessary for the operating system to register that the files have been created    
        self.Unroll_2D_OnShell() #unroll the 2d histogram to make sure the template works fine
        
        
def instantiate_template_helper(reso1, reso1_name, reso2, reso2_name, reso_area, bkgs, bkg_areas, bkg_names, fname, 
                                d_reso1=None, d_reso2=None, d_bkgs=[None], d_name=None,
                                lowerlim=6, upperlim=9, signal_sample_region=None,
                                nbins_1d=40, nbins_2d=20, verbose=False):
    """This function takes in a variety of different parameters in order to make template creation easy

    Arguments:
        reso1 -- The first signal mass resonance. This should be an iterable of masses
        reso1_name -- The name of the first signal sample
        reso2 -- The second signal mass resonance. This should be an iterable of masses
        reso2_name -- The name of the second signal sample
        reso_area -- The area for the signal histogram that you desire
        bkgs -- The background samples for the template. This should a list of mass iterables
        bkg_areas -- The areas for each of the background samples. This should be a list of areas with the same ordering as the bkgs list
        bkg_names -- The names for each of the background samples. This should be a list of names with the same ordering as the bkgs list
        fname -- The filenames you want to use

    Keyword Arguments:
        d_reso1 -- A discriminant for the first signal resonance (default: {None})
        d_reso2 -- A discriminant for the second signal resonance (default: {None})
        d_bkgs -- A list of discriminants for the background samples. This should be a list of discriminants with the same ordering as the bkgs list (default: {[None]})
        d_name -- The name of the discriminant you are using (default: {None})
        lowerlim -- The lower limit for M4L in GeV (default: {6})
        upperlim -- The upper limit for M4L in GeV (default: {9})
        signal_sample_region -- Whether you would like to oversample a region within the range. Provide a range if you do (default: {None})
        nbins_1d -- The number of bins for the 1-dimensional template (default: {40})
        nbins_2d -- The number of bins for the 2-dimensional template (default: {20})
        verbose -- Whether you would like for the creation to be verbose (default: {False})

    Raises:
        ValueError: This is raised if the background lists do not have the same size

    Returns:
        The template class that was created using these parameters after creating the histograms
    """
    
    check_len_list = [bkgs, bkg_areas, bkg_names]
    
    if d_bkgs:
        check_len_list += [d_bkgs]
    
    for i in check_len_list:
        for j in check_len_list:
            if len(i) != len(j):
                raise ValueError("backgrounds lists should have the same length!")
    
    is_2d = False
    if d_reso1 and d_reso2 and d_bkgs and d_name:
        is_2d = True
    
    bkgs_and_their_area = []
    for bkg_sample, bkg_sample_area in zip(bkgs, bkg_areas):
        bkgs_and_their_area.append((bkg_sample, bkg_sample_area))
    
    template_generator = Template_creator(
        (reso1, reso1_name),
        (reso2, reso2_name),
        reso_area,
        fname,
        d_reso1,
        d_reso2,
        d_bkgs,
        bkgs_and_their_area,
        bkg_names,
        d_name,
        lowerlim=lowerlim,
        upperlim=upperlim,
        signal_sample_region=signal_sample_region
    )
    
    template_generator.rooFit_Input_1d(nbins_1d)
    
    if is_2d:
        template_generator.rooFit_Input_2d(nbins_2d)
        
    template_generator.plot_histograms(fname=fname, show=False)
    
    one_d_filename = fname+'.root'
    two_d_filename = fname+'_2d_out.root'
    
    one_d_in_folder = '1d_'+fname
    two_d_in_folder = '2d_'+fname
    
    if not os.path.isdir(one_d_in_folder): #if the folder doesn't exist - make it!
        os.mkdir(one_d_in_folder)
    
    if not os.path.isdir(two_d_in_folder):
        os.mkdir(two_d_in_folder)
    
    shutil.move(one_d_filename, one_d_in_folder+'/'+one_d_filename)
    shutil.move(two_d_filename, two_d_in_folder+'/'+two_d_filename)
    
    one_d_out_folder = one_d_in_folder+'_out'
    two_d_out_folder = two_d_in_folder+'_out'
    
    if not os.path.isdir(one_d_out_folder):
        os.mkdir(one_d_out_folder)
    
    if not os.path.isdir(two_d_out_folder):
        os.mkdir(two_d_out_folder)
        
    for in_folder, out_folder in zip([one_d_in_folder, two_d_in_folder],[one_d_out_folder, two_d_out_folder]):
        runstr = 'python MakeInputRoot_OnShell.py '
        runstr += in_folder + ' ' + out_folder
        if not verbose:
            runstr += ' > /dev/null'
        os.system(runstr)
        
        runstr = 'python DatacardMaker_OnShell.py '
        runstr += out_folder
        if not verbose:
            runstr += ' > /dev/null'
        
        os.system(runstr)
        
    return template_generator