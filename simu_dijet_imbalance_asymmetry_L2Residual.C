#include <vector>
#include <sstream>
#include "TF1.h"
#include "TF2.h"
#include "TF3.h"
#include <TFile.h>
#include <TDirectoryFile.h>
#include <TTree.h>
#include <TKey.h>
#include <iostream>
#include "JetCorrector.h"

void simu_dijet_imbalance_asymmetry_L2Residual() {

  using namespace std;
  gROOT->Reset();

  // TFile *HiFor_file = new TFile("/afs/cern.ch/user/u/uacharya/eos/phys_heavyions/2024ppref_hiforest/JRA/QCD_pThat-15_Dijet_TuneCP5_pp_13p6TeV_pythia8/uacharya-pp_ref_lowpileup_simu_Raw_files_Reco-83737e01b3c9e1d2088c1f0f2fc648a0/USER/QCD_pThat-15_Dijet_TuneCP5_pp_13p6TeV_pythia8/pp_ref_lowpileup_Simu_pp_Reconstructed/241018_155410/0000/JRA_miniTree/merged_JRA_all.root");
  TFile *HiFor_file = new TFile("/eos/cms/store/group/phys_heavyions/uacharya/simu/merged_JRA_all.root");


  //TFile *HiFor_file = new TFile("eos/phys_heavyions/2024ppref_hiforest/JRA/QCD_pThat-15_Dijet_TuneCP5_pp_13p6TeV_pythia8/uacharya-pp_ref_lowpileup_simu_Raw_files_Reco-83737e01b3c9e1d2088c1f0f2fc648a0/USER/QCD_pThat-15_Dijet_TuneCP5_pp_13p6TeV_pythia8/pp_ref_lowpileup_Simu_pp_Reconstructed/241018_155410/0000/JRA_miniTree/merged_JRA_all.root");

  // File and Histogram Setup
  //JetMet0_V2//
  //TFile *HiFor_file = new TFile(Form("/eos/cms/store/group/phys_heavyions/uacharya/JetMet0_v2/HiForestMiniAOD_%d.root", infilInd));


  //JetMet1_V2//
  // TFile *HiFor_file = new TFile(Form("/eos/cms/store/group/phys_heavyions/uacharya/JetMet1_v2/HiForestMiniAOD_%d.root",infilInd));

  TDirectory *dir = (TDirectory*)HiFor_file->Get(Form("ak4pf"));
  if (!dir) {
    cout << "Directory ak4PFJetAnalyzer not found!" << endl;
    return;
  }

  TTree *tree;
  dir->GetObject("t", tree);
  if (!tree) {
    cout << "Tree t not found!" << endl;
    return;
  }
 
  //return;
  
  
  vector<string> JEC_input_file;
  JEC_input_file.push_back("/afs/cern.ch/user/u/uacharya/eos/cms_analaysis/JEC_analysis/JRA_step3harvest_L2Relative_AK4PF.txt");
  JetCorrector JECCor(JEC_input_file);

  const int mx = 100000;

				
  //Declaration of leaves types
  vector<int>     npus;
  vector<float>   tnpus;
  vector<float>   zpositions;
  vector<int>     bxns;
  vector<float>   sumpt_lowpt;
  vector<float>   sumpt_highpt;
  vector<int>     ntrks_lowpt;
  vector<int>     ntrks_highpt;
  vector<float>   rhos;
  Float_t         rho;
  Float_t         pthat;
  Float_t         beta;
  Float_t         betaStar;
  Float_t         weight;
  Float_t         refpvz;
  Float_t         pudensity;
  Float_t         gpudensity;
  Long64_t        npv;
  Long64_t        run;
  Long64_t        lumi;
  Long64_t        evt;
  UChar_t         nref;
  vector<unsigned char> refrank;
  vector<int>     refpdgid;
  vector<float>   refe;
  vector<float>   refpt;
  vector<float>   refeta;
  vector<float>   refphi;
  vector<float>   refy;
  vector<float>   refdrjt;
  vector<float>   refarea;
  vector<float>   jte;
  vector<float>*   jtpt=nullptr;
  vector<float>*  jteta=nullptr;
  vector<float>*  jtphi=nullptr;
  vector<float>*  jty=nullptr;
  vector<float>   jtjec;
  vector<float>   jtarea;
  vector<float>   jtchf;
  vector<float>   jtnhf;
  vector<float>   jtnef;
  vector<float>   jtcef;
  vector<float>   jtmuf;
  vector<float>   jthfhf;
  vector<float>   jthfef;
  vector<int>     refnMult;
  vector<int>     refchMult;
  vector<int>     jtnMult;
  vector<int>     jtchMult;
  vector<float>   refdzvtx;
  // // Set branch addresses.
  // // tree->SetBranchAddress("npus",npus);
  // // tree->SetBranchAddress("tnpus",tnpus);
  // // tree->SetBranchAddress("zpositions",zpositions);
  // // tree->SetBranchAddress("bxns",bxns);
  // // tree->SetBranchAddress("sumpt_lowpt",sumpt_lowpt);
  // // tree->SetBranchAddress("sumpt_highpt",sumpt_highpt);
  // // tree->SetBranchAddress("ntrks_lowpt",ntrks_lowpt);
  // // tree->SetBranchAddress("ntrks_highpt",ntrks_highpt);
  // // tree->SetBranchAddress("rhos",rhos);
  // // tree->SetBranchAddress("rho",rho);
  // // tree->SetBranchAddress("pthat",pthat);
  // // tree->SetBranchAddress("beta",beta);
  // // tree->SetBranchAddress("betaStar",betaStar);
  // // tree->SetBranchAddress("weight",weight);
  // // tree->SetBranchAddress("refpvz",refpvz);
  // // tree->SetBranchAddress("pudensity",pudensity);
  // // tree->SetBranchAddress("gpudensity",gpudensity);
  // // tree->SetBranchAddress("npv",npv);
  // // tree->SetBranchAddress("run",run);
  // // tree->SetBranchAddress("lumi",lumi);
  // // tree->SetBranchAddress("evt",evt);
  tree->SetBranchAddress("nref",&nref);
  // // tree->SetBranchAddress("refrank",refrank);
  // // tree->SetBranchAddress("refpdgid",refpdgid);
  // // tree->SetBranchAddress("refe",refe);
  // // tree->SetBranchAddress("refpt",refpt);
  // // tree->SetBranchAddress("refeta",refeta);
  // // tree->SetBranchAddress("refphi",refphi);
  // // tree->SetBranchAddress("refy",refy);
  // // tree->SetBranchAddress("refdrjt",refdrjt);
  // // tree->SetBranchAddress("refarea",refarea);
  // // tree->SetBranchAddress("jte",jte);
  tree->SetBranchAddress("jtpt",&jtpt);
  tree->SetBranchAddress("jteta",&jteta);
  tree->SetBranchAddress("jtphi",&jtphi);
  tree->SetBranchAddress("jty",&jty);
  // // tree->SetBranchAddress("jtjec",jtjec);
  // // tree->SetBranchAddress("jtarea",jtarea);
  // // tree->SetBranchAddress("jtchf",jtchf);
  // // tree->SetBranchAddress("jtnhf",jtnhf);
  // // tree->SetBranchAddress("jtnef",jtnef);
  // // tree->SetBranchAddress("jtcef",jtcef);
  // // tree->SetBranchAddress("jtmuf",jtmuf);
  // // tree->SetBranchAddress("jthfhf",jthfhf);
  // // tree->SetBranchAddress("jthfef",jthfef);
  // // tree->SetBranchAddress("refnMult",refnMult);
  // // tree->SetBranchAddress("refchMult",refchMult);
  // // tree->SetBranchAddress("jtnMult",jtnMult);
  // // tree->SetBranchAddress("jtchMult",jtchMult);
  // // tree->SetBranchAddress("refdzvtx",refdzvtx);

 
  int nEta = 22;
  float etaBin[] = {-5.0, -3.5, -3.0, -2.7, -2.4, -2.0, -1.7, -1.3, -1.0, -0.7, -0.4, 0.0, 0.4, 0.7, 1.0, 1.3, 1.7, 2.0, 2.4, 2.7, 3.0, 3.5, 5.0};
  const int nptBin = 5;
  double avgPtBin[]={50.0,90.0,120.0,160.0,200.0,500.0};
  
  TH1F *hPt_eta = new TH1F("hPt_eta", "Leading Jet Pt by Eta Bin", nEta, etaBin);
  TH1F *hPt_eta13 = new TH1F("hPt_eta13", "Second Jet Pt by Eta Bin", nEta, etaBin);
  TH1F *hPt_eta_third = new TH1F("hPt_eta_third", "Third Jet Pt by Eta Bin", nEta, etaBin);
  TH1F *hPt_eta_ratio;
  Long64_t nentries = tree->GetEntries();
  Float_t jtPf_CHF,jtPf_MUF,jtPf_NHF,jtPf_NEF;
  Int_t jtPf_CHM;

  TH1F* hAssymBins[nptBin][nEta];

  for (int ipt = 0; ipt <nptBin; ipt++) {
    for (int i = 0; i < nEta; i++) {
      hAssymBins[ipt][i] = new TH1F(Form("hAssym_etaBin_%d_%d", ipt,i), Form("Asymmetry in Eta Bin %.1f to %.1f", etaBin[i], etaBin[i+1]), 200, -1, 1);
      hAssymBins[ipt][i] ->Sumw2();
    }
  }
  // Event Loop
  bool alternate=true;
  for (Long64_t i = 0; i < nentries; i++) {
  //for (Long64_t i = 0; i < 1000; i++) {
  //for (Long64_t i = 0; i < 500000; i++) {
  //for (Long64_t i = 500000; i < nentries; i++) {
    tree->GetEntry(i);

    float maxPt = -9999, maxEta = -9999, maxPhi = -9999;
    float secondMaxPt = -9999, secondMaxEta = -9999, secondMaxPhi = -9999;
    float thirdMaxPt = -9999, thirdMaxEta = -9999, thirdMaxPhi = -9999;
    float assym;
    float avg_pT;
    float alpha_Val;  
    for (int j = 0; j < nref; j++) {
      
      float pt = (*jtpt)[j];
      float eta = (*jteta)[j];
      float phi = (*jtphi)[j];

      //cout<<i<<'\t'<<j<<'\t'<<"raw: " << pt <<'\t'<< eta << '\t'<< phi <<endl;

      if (pt > maxPt) {
	// Update second and third jets before changing maxPt
	thirdMaxPt = secondMaxPt;
	thirdMaxEta = secondMaxEta;
	thirdMaxPhi = secondMaxPhi;
	//cout<<i<<'\t'<<j<<'\t'<<"third1: " << thirdMaxPt <<'\t'<< thirdMaxEta << '\t'<< thirdMaxPhi <<endl;

	secondMaxPt = maxPt;
	secondMaxEta = maxEta;
	secondMaxPhi = maxPhi;
	//cout<<i<<'\t'<<j<<'\t'<<"second1: " << secondMaxPt <<'\t'<< secondMaxEta << '\t'<< secondMaxPhi <<endl;

	maxPt = pt;
	maxEta = eta;
	maxPhi = phi;
	//cout<<i<<'\t'<<j<<'\t'<<"first: " << pt <<'\t'<< eta << '\t'<< phi <<endl;

      } else if (pt > secondMaxPt) {
	// Update third jet before changing secondMaxPt
	thirdMaxPt = secondMaxPt;
	thirdMaxEta = secondMaxEta;
	thirdMaxPhi = secondMaxPhi;
	//cout<<i<<'\t'<<j<<'\t'<<"third2: " << thirdMaxPt <<'\t'<< thirdMaxEta << '\t'<< thirdMaxPhi <<endl;

	secondMaxPt = pt;
	secondMaxEta = eta;
	secondMaxPhi = phi;
	//cout<<i<<'\t'<<j<<'\t'<<"second2: " << secondMaxPt <<'\t'<< secondMaxEta << '\t'<< secondMaxPhi <<endl;

      } else if (pt > thirdMaxPt) {
	// Only update third jet
	thirdMaxPt = pt;
	thirdMaxEta = eta;
	thirdMaxPhi = phi;
	//cout<<i<<'\t'<<j<<'\t'<<"third3: " << thirdMaxPt <<'\t'<< thirdMaxEta << '\t'<< thirdMaxPhi <<endl;

      }
      
    }
  

#if 1
    if(1){
      //if (phiCondition)// && etaCondition)// && alphaCondition)
      float alphaValue=0.075;
      for(int ipt=0;ipt<nptBin; ipt++){
	for (int ib = 0; ib < nEta; ib++) {
	  if( fabs(maxEta)<1.3 &&fabs(secondMaxEta)>1.3){
	    float delPhi = TMath::Abs(maxPhi - secondMaxPhi);
	    bool phiCondition = (delPhi >2.7);
	    if(!phiCondition) continue;
	    avg_pT=0.5*(maxPt+secondMaxPt);
	    alpha_Val=thirdMaxPt/avg_pT;
	    bool alphaCondition1101=(alpha_Val<alphaValue);
	    if(!alphaCondition1101) continue;
	    if((avg_pT >= avgPtBin[ipt] && avg_pT < avgPtBin[ipt+1])){
	      if(secondMaxEta >= etaBin[ib] && secondMaxEta < etaBin[ib+1]) {
		assym=(secondMaxPt-maxPt)/(secondMaxPt+maxPt);
		if(fabs(assym)>0.55) continue;	      
		hAssymBins[ipt][ib]->Fill(assym);
		cout<<"case40: "<<i<<'\t'<<etaBin[ib]<<'\t'<<etaBin[ib+1]<<'\t'<<secondMaxEta<<'\t'<<secondMaxPt<<'\t'<<maxEta<<'\t'<<maxPt<<'\t'<<avg_pT<<'\t'<<assym<<endl;
	      }	  
	    }
	  }
	  else if(fabs(secondMaxEta)<1.3 &&fabs(maxEta)>1.3){
	    float delPhi = TMath::Abs(maxPhi - secondMaxPhi);
	    bool phiCondition = (delPhi >2.7);
	    if(!phiCondition) continue;
	    avg_pT=0.5*(maxPt+secondMaxPt);
	    alpha_Val=thirdMaxPt/avg_pT;
	    bool alphaCondition1102=(alpha_Val<alphaValue);
	    if(!alphaCondition1102) continue;
	    //if((avg_pT>=avgPtBin[0] &&avg_pT<avgPtBin[1])){
	    if(avg_pT>=avgPtBin[ipt] && avg_pT<avgPtBin[ipt+1]){
	      if( maxEta>= etaBin[ib] &&  maxEta< etaBin[ib+1]){
		assym=(maxPt-secondMaxPt)/(maxPt+secondMaxPt);
		if(fabs(assym)>0.55) continue;	      
		hAssymBins[ipt][ib]->Fill(assym);
		cout<<"case80: "<<i<<'\t'<<etaBin[ib]<<'\t'<<etaBin[ib+1]<<'\t'<<secondMaxEta<<'\t'<<secondMaxPt<<'\t'<<maxEta<<'\t'<<maxPt<<'\t'<<avg_pT<<'\t'<<assym<<endl;
	      }
	    }
	  }
	  else if((fabs(secondMaxEta)<1.3 &&fabs(maxEta)<1.3)){
	    float delPhi = TMath::Abs(maxPhi - secondMaxPhi);
	    bool phiCondition = (delPhi >2.7);
	    if(!phiCondition) continue;
	    avg_pT=0.5*(maxPt+secondMaxPt);
	    alpha_Val=thirdMaxPt/avg_pT;
	    bool alphaCondition1103=(alpha_Val<alphaValue);
	    if(!alphaCondition1103) continue;
	    if(avg_pT>=avgPtBin[ipt] &&avg_pT<avgPtBin[ipt+1]){
	      if(maxEta >= etaBin[ib] &&  maxEta < etaBin[ib+1] && secondMaxEta >= etaBin[ib] &&  secondMaxEta < etaBin[ib+1]){
		assym= alternate ? (maxPt-secondMaxPt)/(maxPt+secondMaxPt) : (secondMaxPt-maxPt)/(maxPt+secondMaxPt);
		//assym=(maxPt-secondMaxPt)/(maxPt+secondMaxPt);
		if(fabs(assym)>0.55) continue;	      
		hAssymBins[ipt][ib]->Fill(assym);
		//filled = true;
		//	cout<<"case1: "<<i<<'\t'<<etaBin[ib]<<'\t'<<etaBin[ib+1]<<'\t'<<secondMaxEta<<'\t'<<secondMaxPt<<'\t'<<maxEta<<'\t'<<maxPt<<'\t'<<avg_pT<<'\t'<<assym<<endl;
		alternate = !alternate;
		break;
	      }
	      else if( maxEta>= etaBin[ib] &&  maxEta< etaBin[ib+1]) {
		assym=(maxPt-secondMaxPt)/(maxPt+secondMaxPt);
		if(fabs(assym)>0.55) continue;	      
		hAssymBins[ipt][ib]->Fill(assym);
		//cout<<"case2: "<<i<<'\t'<<etaBin[ib]<<'\t'<<etaBin[ib+1]<<'\t'<<secondMaxEta<<'\t'<<secondMaxPt<<'\t'<<maxEta<<'\t'<<maxPt<<'\t'<<avg_pT<<'\t'<<assym<<endl;
		//break;
	      }
	      else if(secondMaxEta>= etaBin[ib] &&  secondMaxEta< etaBin[ib+1]) {
		assym=(secondMaxPt-maxPt)/(maxPt+secondMaxPt);
		if(fabs(assym)>0.55) continue;	      
		hAssymBins[ipt][ib]->Fill(assym);
		//cout<<"case3: "<<i<<'\t'<<etaBin[ib]<<'\t'<<etaBin[ib+1]<<'\t'<<secondMaxEta<<'\t'<<secondMaxPt<<'\t'<<maxEta<<'\t'<<maxPt<<'\t'<<avg_pT<<'\t'<<assym<<endl;
		//break;
	      }
	    }
	  }	    	
	}
      }
    }
      
#endif
  }
 
#if 1

  // Save Output
  TFile *outputFile = new TFile(Form("/afs/cern.ch/user/u/uacharya/eos/cms_analaysis/JEC_analysis/histo_eff_JEC/dijet_balance_after_JEC_simu_alphaCond0075_final2.root"), "RECREATE");
  for (int ipt = 0; ipt < nptBin; ipt++) {
    for (int ib = 0; ib < nEta; ib++) {
      hAssymBins[ipt][ib]->Write();
      //hAssymBins40[ipt][ib]->Write();
      //hAssymBins80[ipt][ib]->Write();
      //hAssymBins110[ipt][ib]->Write();
    }
  }
  // h_alpha->Write();
  outputFile->Close();
#endif

  
}
