# CMS_Analysis
Analysis codes for the low pile up JEC/JER studies.

#jec_dijet_data_simu_diff.C 

		-> For different alpha conditions, this code plots asymmetry distribution for different pT ranges and within different eta region in data and simulation. Afterward extract response as function of eta for different pT ranges for each alpha conditions. The pT ranges are specified according to the HTL trigger cuts. The ratio of response between MC/Data vs eta is plotted for different pT range and save them as histograms in the root file for each alpha conditions.

#data_simu_ratio_vs_alpha_newMethod.C

   ->Use root files created using " jec_dijet_data_simu_diff.C " are now used here to extract the L2Residual correction factors. Here I used two ways of extraction of L2Residuals corrections. In the first method, which is referred to as oldMethod, plot the Mc/Data ratio Vs alpha condition as in TGraph for each ptRange and eta bins. For each eta bins, MC/Data vs alpha is plotted. The  graph is extrapolated to get the value when alpha=0. Once you have those values from extrapolation, The MC/Data vs eta is plotted for each pT ranges. Also potted is the MC/Data vs |eta| as well.  TO get the correction factor, a for each eta bins, MC/Data vs pT is plotted and is fitted by some dedicated fit function and those fitted parameters are extracted in the textual format. In first method, L2Residual_lowPU_pp13TeV.txt is generated. In the second method, we use different approach, where instead of extrapolating the MC/Data vs alpha plot to 0, for each eta range, the MC/Data(pT,eta,alpja) is normalized by the MC/Data(pT,eta,alpha<0.3). These noramlized values  of MC/Data ratio are plotted as function of alpha simultaneously for different pT ranges within each eta bins and a linar fitting is done, which gives the slope and this slope is the kFSR factor. Then the  extrapolation to alpha=0 is given by MC/Data(pT,eta,alpha<0.3) - ( 0.3 * KFSR* MC/Data(pT,eta,alpha<0.3)). The L2residuals correction is extracted by fiting the MC/Data(pt,eta,alpha=0) vs pT for each etabins and saved in textfile.
