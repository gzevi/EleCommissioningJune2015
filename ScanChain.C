// Usage:
// > root -b doAll.C

// C++
#include <iostream>
#include <vector>

// ROOT
#include "TBenchmark.h"
#include "TChain.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TROOT.h"
#include "TTreeCache.h"

#include "TLorentzVector.h"
#include "Math/Vector4D.h"
#include "Math/LorentzVector.h"

// CMS3
#include "CMS3.cc"
#include "PlotUtilities.h"

using namespace std;
using namespace tas;




void makePlots(std::map<std::string, TH1*> & h_1, TString sel, int iel, float weight = 1) {
  
  float pt = els_p4().at(iel).pt();
  float phi = els_p4().at(iel).phi();
  float eta = els_etaSC().at(iel);
  float seedE = els_eSeed().at(iel);
  float theta = 2.0*TMath::ATan(TMath::Exp(-1.*eta));
  float seedEt = seedE * TMath::Sin(theta);
  float SCrawEt = els_eSCRaw().at(iel) * TMath::Sin(theta);
  
  TString EBEE = "";
  TString EBEEsign = "";
  
  if (fabs(eta) > 1.57) {
    EBEE = "EE";
    if (eta > 0) EBEEsign = "EEpos";
    if (eta < 0) EBEEsign = "EEneg";
  }
  else if (fabs(eta)<1.44){
    EBEE = "EB";
    EBEEsign = "EB";
  }
  else return;
  plot1D(("h"+sel+"_pt"+EBEE).Data(), pt,  weight, h_1, "pt", 50, 0, 100);
  plot1D(("h"+sel+"_seedEt"+EBEE).Data(), seedEt,  weight, h_1, "seed ET", 50, 0, 100);
  plot1D(("h"+sel+"_SCrawEt"+EBEE).Data(), SCrawEt,  weight, h_1, "raw SC ET", 50, 0, 100);
  plot1D(("h"+sel+"_eta").Data(), eta,  weight, h_1, "eta", 50, -2.5, 2.5);
  plot1D(("h"+sel+"_phi").Data(), phi,  weight, h_1, "phi", 50, -3.5, 3.5);
  
  plot1D(("h"+sel+"_relchiso"+EBEE).Data(), els_pfChargedHadronIso().at(iel)/seedEt,  weight, h_1, "PFCh", 100, 0, 1);
  plot1D(("h"+sel+"_relemiso"+EBEE).Data(), els_pfPhotonIso().at(iel)/seedEt,  weight, h_1, "PFEM", 100, 0, 1);
  plot1D(("h"+sel+"_relnhiso"+EBEE).Data(), els_pfNeutralHadronIso().at(iel)/seedEt,  weight, h_1, "PFNh", 100, 0, 1);
  
  plot1D(("h"+sel+"_relECALiso"+EBEE).Data(), els_ecalIso().at(iel)/seedEt,  weight, h_1, "ECAL RelIso", 100, 0, 1);
  plot1D(("h"+sel+"_relHCALiso"+EBEE).Data(), els_hcalIso().at(iel)/seedEt,  weight, h_1, "HCAL RelIso", 100, 0, 1);
  plot1D(("h"+sel+"_relECALHCALiso"+EBEE).Data(), (els_ecalIso().at(iel)+els_hcalIso().at(iel))/seedEt,  weight, h_1, "HCAL RelIso", 100, 0, 1);
  
  plot1D(("h"+sel+"_sieie"+EBEE).Data(), els_sigmaIEtaIEta_full5x5().at(iel),  weight, h_1, "sieie_5x5", 100, 0, EBEE=="EB" ? 0.05 : 0.1);
  plot1D(("h"+sel+"_sipip"+EBEE).Data(), els_sigmaIPhiIPhi_full5x5().at(iel),  weight, h_1, "sipip_5x5", 100, 0, EBEE=="EB" ? 0.05 : 0.1);
  plot1D(("h"+sel+"_deta"+EBEEsign).Data(), els_dEtaIn().at(iel),  weight, h_1, "DeltaEta", 50, -0.1, 0.1);
  if (fabs(eta) > 1.57) plot1D(("h"+sel+"_deta"+EBEE).Data(), els_dEtaIn().at(iel),  weight, h_1, "DeltaEta", 50, -0.1, 0.1);
  plot1D(("h"+sel+"_dphi"+EBEEsign).Data(), els_dPhiIn().at(iel),  weight, h_1, "DeltaPhi", 50, -0.1, 0.1);
  if (fabs(eta) > 1.57) plot1D(("h"+sel+"_dphi"+EBEE).Data(), els_dPhiIn().at(iel),  weight, h_1, "DeltaPhi", 50, -0.1, 0.1);
  plot1D(("h"+sel+"_HoverE"+EBEE).Data(), els_hOverE().at(iel),  weight, h_1, "H/E", 50, 0, 0.5);
  if (fabs(eta) > 1.57) plot1D(("h"+sel+"_detaSeed"+EBEE).Data(), els_dEtaOut().at(iel),  weight, h_1, "DeltaEta", 50, -0.1, 0.1);
  if (fabs(eta) > 1.57) plot1D(("h"+sel+"_dphiSeed"+EBEE).Data(), els_dPhiOut().at(iel),  weight, h_1, "DeltaPhi", 50, -0.1, 0.1);
  plot1D(("h"+sel+"_detaSeed"+EBEEsign).Data(), els_dEtaOut().at(iel),  weight, h_1, "DeltaEta", 50, -0.1, 0.1);
  plot1D(("h"+sel+"_dphiSeed"+EBEEsign).Data(), els_dPhiOut().at(iel),  weight, h_1, "DeltaPhi", 50, -0.1, 0.1);
  
  plot1D(("h"+sel+"_PSEoverRawE"+EBEE).Data(), els_eSCPresh().at(iel) /els_eSCRaw().at(iel) ,  weight, h_1, "EPSoverERawSC", 50, 0, 1);

  return;
  
}

//int ScanChain( TChain* chain, bool fast = true, int nEvents = -1, string skimFilePrefix = "test") {
int ScanChain( TChain* chain, int nEvents = -1, const char* outName = "test", float isData = false) {
  float fast = true;
  // Benchmark
  TBenchmark *bmark = new TBenchmark();
  bmark->Start("benchmark");
  
  // Example Histograms
  TDirectory *rootdir = gDirectory->GetDirectory("Rint:");
  TH1F *samplehisto = new TH1F("samplehisto", "Example histogram", 200,0,200);
  samplehisto->SetDirectory(rootdir);
  
  //  TFile * outfile_ = new TFile("plots.root","RECREATE") ;
  
  // Loop over events to Analyze
  unsigned int nEventsTotal = 0;
  unsigned int nEventsChain = chain->GetEntries();
  if( nEvents >= 0 ) nEventsChain = nEvents;
  TObjArray *listOfFiles = chain->GetListOfFiles();
  TIter fileIter(listOfFiles);
  TFile *currentFile = 0;
  
  int trigEG5 = 0;
  int trigEG20 = 0;
  int trigDiEG = 0;
  int trigDiSCJPsi = 0;
  int trigDiSCUpsilon = 0;
  int trigEleSCJPsi = 0;
  
  std::map<std::string, TH1*> h_1d;
  
  
  float vtxWeight[17];
  vtxWeight[0]=85.7184;
  vtxWeight[1]=35.8861;
  vtxWeight[2]=17.9012;
  vtxWeight[3]=12.8676;
  vtxWeight[4]=6.40865;
  vtxWeight[5]=3.3218;
  vtxWeight[6]=2.42468;
  vtxWeight[7]=1.4026;
  vtxWeight[8]=0.907739;
  vtxWeight[9]=0.559019;
  vtxWeight[10]=0.417363;
  vtxWeight[11]=0.29857;
  vtxWeight[12]=0.225661;
  vtxWeight[13]=0.173138;
  vtxWeight[14]=0.135456;
  vtxWeight[15]=0.106216;
  vtxWeight[16]=0.0842394;
  
  
  // File Loop
  while ( (currentFile = (TFile*)fileIter.Next()) ) {
    TString currfile =currentFile->GetTitle();
    cout<<"Looking at file: "<<currfile<<endl;
    // Get File Content
    TFile *file = new TFile( currentFile->GetTitle() );
    TTree *tree = (TTree*)file->Get("Events");
    if(fast) TTreeCache::SetLearnEntries(10);
    if(fast) tree->SetCacheSize(128*1024*1024);
    cms3.Init(tree);
    
    // Loop over Events in current file
    if( nEventsTotal >= nEventsChain ) continue;
    unsigned int nEventsTree = tree->GetEntriesFast();
    for( unsigned int event = 0; event < nEventsTree; ++event) {
      
      // Get Event Content
      if( nEventsTotal >= nEventsChain ) continue;
      if(fast) tree->LoadTree(event);
      cms3.GetEntry(event);
      ++nEventsTotal;
      
      // Progress
      CMS3::progress( nEventsTotal, nEventsChain );
      
      
      bool passtrigEG5  = false;
      bool passtrigEG20 = false;
      bool passtrigDiEG = false;
      bool passtrigDiSCJPsi = false;
      bool passtrigDiSCUpsilon = false;
      bool passtrigEleSCJPsi = false;
      
      //if (passHLTTrigger("HLT_L1SingleEG5_v1"))                  { passtrigEG5   = true; trigEG5++ ; }
      //if (passHLTTrigger("HLT_L1SingleEG20_v1"))                 { passtrigEG20  = true; trigEG20++; }
      //if (passHLTTrigger("HLT_DiSC30_18_EIso_AND_HE_Mass70_v1")) { passtrigDiEG  = true; trigDiEG++; }
      //if (passHLTTrigger("HLT_DiSC5_JPsi_v1")) { passtrigDiSCJPsi  = true; trigDiSCJPsi++; }
      //if (passHLTTrigger("HLT_DiSC5_Upsilon_v1")) { passtrigDiSCUpsilon  = true; trigDiSCUpsilon++; }
      //if (passHLTTrigger("HLT_Ele5_SC5_JPsi_Mass2to4p5_v2")) { passtrigEleSCJPsi  = true; trigEleSCJPsi++; }
      
      int nvtx =evt_nvtxs();
      plot1D("h_nvtx", nvtx,  1, h_1d, "NVTX", 50, 0, 50);
      plot1D(("h_nvtx_"+currfile).Data(), nvtx,  1, h_1d, "NVTX", 50, 0, 50);
      
      
      // Analysis Code
      for (unsigned int iel = 0; iel < els_p4().size(); ++iel) {
        LorentzVector el_p4 = els_p4().at(iel);
        float pt = el_p4.pt();
        float eta = els_etaSC().at(iel);
        samplehisto->Fill(els_p4().at(iel).pt());
        
        float theta = 2.0*TMath::ATan(TMath::Exp(-1.*eta));
        
        float SCEt = els_eSC().at(iel) * TMath::Sin(theta);
        
        plot1D("h_pt", pt,  1, h_1d, "pT [GeV]", 150, 0, 150);
        if (fabs(eta) > 2.5 || (pt < 5 ) ) continue;
        if ( fabs(eta) < 1.44 && fabs(eta) > 1.57) continue;
        //	if (fabs(eta) > 2.5 || (pt < 10) ) continue;
        float sieie =  els_sigmaIEtaIEta_full5x5().at(iel);
        float sipip =  els_sigmaIPhiIPhi_full5x5().at(iel);
        float deta = els_dEtaIn().at(iel);
        float dphi =  els_dPhiIn().at(iel);
        float detaSeed = els_dEtaOut().at(iel);
        float dphiSeed =  els_dPhiOut().at(iel);
        float hovere = els_hOverE().at(iel);
        float chiso     = els_pfChargedHadronIso().at(iel);
        float nhiso     = els_pfNeutralHadronIso().at(iel);
        float ecaliso = els_ecalIso().at(iel);
        float hcaliso = els_hcalIso().at(iel);
        float emiso     = els_pfPhotonIso().at(iel);
        float deltaBeta = els_pfPUIso().at(iel);
        float iso = chiso + nhiso + emiso - 0.5* deltaBeta;
        float caloiso = nhiso + emiso ;
        float trkshFrac = els_trkshFrac().at(iel);
        if (trkshFrac == -9999.) trkshFrac = -1.;
        float trkdr = els_trkdr().at(iel);
        if (trkdr == -9999.) trkdr = -1.;
        float seedE = els_eSeed().at(iel);
        float seedEt = seedE * TMath::Sin(theta);
        float SCrawE = els_eSCRaw().at(iel);
        float SCrawEt = SCrawE * TMath::Sin(theta);
        
        LorentzVector el_seedCl = el_p4;
        //cout<<"Start with "<<el_seedCl<<endl;
        el_seedCl *= seedE/el_p4.E();
        //cout<<"turn into "<<el_seedCl<<endl;
        LorentzVector el_SC = el_p4;
        //el_SC.SetE(els_eSC().at(iel));
        el_SC *= els_eSC().at(iel) / el_p4.E();
        LorentzVector el_SCraw = el_p4;
        el_SCraw *= els_eSCRaw().at(iel) / el_p4.E();
        
        //	cout<<"seed cluster: E, eta, theta, Et: "<<seedE<<" "<<eta<<" "<<theta<<" "<<seedEt<<endl;
        //	cout<<"supercluster: E, eta, theta, Et: "<<els_eSC().at(iel)<<" "<<eta<<" "<<theta<<" "<<SCEt<<endl;
        
        //if (pt > 150) cout<<"found 150 GeV event/ls/run: "<<evt_event()<<" "<<evt_lumiBlock()<<" "<<evt_run()<<endl;
        
        float invMassForZveto = 0.;
        for (unsigned int jel = 0; jel < els_p4().size(); ++jel) {
          if (iel==jel) continue;
          LorentzVector el2_p4 = els_p4().at(jel);
          LorentzVector el2_seedCl = el2_p4;
          el2_seedCl *= els_eSeed().at(jel)/el2_p4.E();
          LorentzVector massSeedCl = el_seedCl + el2_seedCl;
          invMassForZveto = massSeedCl.M();
        }
        if (invMassForZveto < 70 || invMassForZveto > 110 || invMassForZveto < 12) {
          makePlots(h_1d, "zveto5", iel);
          if (nvtx > 10 && nvtx < 26) makePlots(h_1d, "zveto5nvtxW", iel, isData ? 1 : vtxWeight[nvtx-10]);
        }
        
        
        
        TString sel = "unsel5";
        makePlots(h_1d, sel, iel);
        if (nvtx > 10 && nvtx < 26) makePlots(h_1d, sel+"nvtxW", iel, isData ? 1 : vtxWeight[nvtx-10]);
        
        if (SCrawEt < 10)  continue;
        
        sel = "unsel10";
        makePlots(h_1d, sel, iel);
        
        if (SCrawEt < 20 )  continue;
        
        plot1D("h_trkdr", trkdr,  1, h_1d, "DR (KF,GSF)", 100, -1.1, 0.5);
        plot1D("h_trkshFrac", trkshFrac,  1, h_1d, "Shared Fraction", 100, -1.1, 1.1);
        
        plot1D("h_ElePtMinusSCEt", (pt-SCEt)/SCEt,  1, h_1d, "(pT(ele) - ET(SC)) / ET(SC)", 100, -1.1, 1.1);
        plot1D("h_SCEtMinusSeedEt", (SCEt-seedEt)/seedEt,  1, h_1d, "(ET(SC) - ET(seed)) / ET(seed)", 100, -1.1, 1.1);
        
        
        
        sel = "unsel";
        makePlots(h_1d, sel, iel);
        
        
        bool selected1 = true;
        if ( (ecaliso+hcaliso)/SCrawEt > 0.1) selected1 = false;
        //	if (fabs(eta) < 1.5 && (fabs(deta) > 0.005 || fabs(dphi) > 0.005 || hovere > 0.005 || sieie > 0.012)  ) selected1 = false;
        //	if (fabs(eta) > 1.5 && (fabs(deta) > 0.01  || fabs(dphi) > 0.005 || hovere > 0.005 || sieie > 0.027)  ) selected1 = false;
        if (fabs(eta) < 1.5 && (fabs(deta) > 0.01  || fabs(dphi) > 0.020 || hovere > 0.05 || sieie > 0.012 || sipip > 0.012)  ) selected1 = false;
        if (fabs(eta) > 1.5 && (fabs(deta) > 0.03  || fabs(dphi) > 0.020 || hovere > 0.05 || sieie > 0.030 || sipip > 0.030)  ) selected1 = false;
        
        sel = "sel1";
        if (selected1) makePlots(h_1d, sel, iel);
        
        if (selected1) plot1D("h_ptSelected", pt,  1, h_1d, "pT [GeV]", 150, 0, 150);
        if (selected1 && passtrigEG5   ) plot1D("h_ptSelectedtrigEG5 ", pt,  1, h_1d, "pT [GeV]", 150, 0, 150);
        if (selected1 && passtrigEG20  ) plot1D("h_ptSelectedtrigEG20", pt,  1, h_1d, "pT [GeV]", 150, 0, 150);
        if (selected1 && passtrigDiEG  ) plot1D("h_ptSelectedtrigDiEG", pt,  1, h_1d, "pT [GeV]", 150, 0, 150);
        
        for (unsigned int jel = 0; jel < els_p4().size(); ++jel) {
          if (iel==jel) continue;
          LorentzVector el2_p4 = els_p4().at(jel);
          float pt2 = el2_p4.pt();
          float eta2 = els_etaSC().at(jel);
          float theta2 = 2.0*TMath::ATan(TMath::Exp(-1.*eta2));
          float SCrawEt2 = els_eSCRaw().at(jel) * TMath::Sin(theta2);
          
          if (fabs(eta2) > 2.5 || SCrawEt2 < 20 ) continue;
          //if ( fabs(eta2) < 1.44 && fabs(eta2) > 1.57) continue;
          
          float sieie2 =  els_sigmaIEtaIEta_full5x5().at(jel);
          float sipip2 =  els_sigmaIPhiIPhi_full5x5().at(jel);
          float deta2 = (els_dEtaIn().at(jel));
          float dphi2 =  (els_dPhiIn().at(jel));
          float detaSeed2 = els_dEtaOut().at(jel);
          float dphiSeed2 =  els_dPhiOut().at(jel);
          float hovere2 = els_hOverE().at(jel);
          float chiso2     = els_pfChargedHadronIso().at(jel);
          float nhiso2     = els_pfNeutralHadronIso().at(jel);
          float emiso2     = els_pfPhotonIso().at(jel);
          float ecaliso2 = els_ecalIso().at(jel);
          float hcaliso2 = els_hcalIso().at(jel);
          float deltaBeta2 = els_pfPUIso().at(jel);
          float iso2 = chiso2 + nhiso2 + emiso2 - 0.5* deltaBeta2;
          float caloiso2 = nhiso2 + emiso2 ;

          float seedE2 = els_eSeed().at(jel);
          float seedEt2 = seedE * TMath::Sin(theta2);
 
          LorentzVector el2_seedCl = el2_p4;
          //el2_seedCl.SetE( els_eSeed().at(jel));
          el2_seedCl *= els_eSeed().at(jel)/el2_p4.E();
          
          LorentzVector el2_SC = el2_p4;
          //el2_SC.SetE(els_eSC().at(jel));
          el2_SC *= els_eSC().at(jel) / el2_p4.E();
          
          LorentzVector el2_SCraw = el2_p4;
          el2_SCraw *= els_eSCRaw().at(jel) / el2_p4.E();
          
          
          bool selected2 = true;
          bool selected2loose = true;
          if ((ecaliso2+hcaliso2)/SCrawEt2 > 0.15) selected2 = false;
          if ((ecaliso2+hcaliso2)/SCrawEt2 > 0.3) selected2loose = false;
          if (fabs(eta2) < 1.5 && (fabs(deta2) > 0.01  || fabs(dphi2) > 0.020 || hovere2 > 0.05 || sieie2 > 0.012 || sipip2 > 0.012)  ) selected2 = false;
          if (fabs(eta2) > 1.5 && (fabs(deta2) > 0.03  || fabs(dphi2) > 0.020 || hovere2 > 0.05 || sieie2 > 0.030 || sipip2 > 0.030)  ) selected2 = false;
          //if (fabs(eta2) < 1.5 && (fabs(deta2) > 0.04  || fabs(dphi2) > 0.04 || hovere2 > 0.15 || sieie2 > 0.025)  ) selected2loose = false;
          //if (fabs(eta2) > 1.5 && (fabs(deta2) > 0.04  || fabs(dphi2) > 0.04 || hovere2 > 0.15 || sieie2 > 0.050)  ) selected2loose = false;
          
          sel = "sel2";
          if (selected1 && selected2) makePlots(h_1d, sel, jel);
          
          
          if (selected1 && selected2) plot1D("h_pt2Selected12", pt2,  1, h_1d, "pT [GeV]", 150, 0, 150);
          
          LorentzVector mass = el_p4 + el2_p4;
          LorentzVector massSC = el_SC + el2_SC;
          LorentzVector massSCraw = el_SCraw + el2_SCraw;
          LorentzVector massSeedCl = el_seedCl + el2_seedCl;
          //cout<<"mass, massSC, massSeedCl: "<<mass.M()<<" "<<massSC.M()<<" "<<massSeedCl.M()<<endl;
          
          
          // Plot some distributions for probes inside mass window
          if (selected1 && massSeedCl.M() > 80 && massSeedCl.M() < 100 && selected2loose) {
            
            sel = "tag";
            makePlots(h_1d, sel, iel);
            
            sel = "probe";
            makePlots(h_1d, sel, jel);
            if (nvtx > 10 && nvtx < 26) makePlots(h_1d, sel+"nvtxW", jel, isData ? 1 : vtxWeight[nvtx-10]);
            
            
            
            
          } // if selected1 && mass.M() > 70 && mass.M() < 110
          
          if (iel > jel) continue; // Only make Z's once!
          
          if (selected1 && selected2) {
            plot1D("h_selectedMass", mass.M(),  1, h_1d, "M_ee [GeV]", 200, 0, 200);
            plot1D("h_selectedMassWideBins", mass.M(),  1, h_1d, "M_ee [GeV]", 200, 0, 200);
            if (fabs(eta) > 1.5 && fabs(eta2) > 1.5) plot1D("h_selectedMassEEEE", mass.M(),  1, h_1d, "M_ee [GeV]", 200, 0, 200);
            if (fabs(eta) < 1.5 && fabs(eta2) < 1.5) plot1D("h_selectedMassEBEB", mass.M(),  1, h_1d, "M_ee [GeV]", 200, 0, 200);
            if ( (fabs(eta) > 1.5 && fabs(eta2) < 1.5) || (fabs(eta) < 1.5 && fabs(eta2) > 1.5))
              plot1D("h_selectedMassEBEE", mass.M(),  1, h_1d, "M_ee [GeV]", 200, 0, 200);
            if (fabs(eta) > 1.5 || fabs(eta2) > 1.5) plot1D("h_selectedMassEE", mass.M(),  1, h_1d, "M_ee [GeV]", 200, 0, 200);
            if (fabs(eta) > 1.5 || fabs(eta2) > 1.5) plot1D(("h_selectedMassEE"+currfile).Data(), mass.M(),  1, h_1d, "M_ee [GeV]", 200, 0, 200);
            
            plot1D("h_selectedMassSeed", massSeedCl.M(),  1, h_1d, "M_ee [GeV]", 200, 0, 200);
            if (fabs(eta) > 1.5 && fabs(eta2) > 1.5) plot1D("h_selectedMassSeedEEEE", massSeedCl.M(),  1, h_1d, "M_ee [GeV]", 200, 0, 200);
            if (fabs(eta) < 1.5 && fabs(eta2) < 1.5) plot1D("h_selectedMassSeedEBEB", massSeedCl.M(),  1, h_1d, "M_ee [GeV]", 200, 0, 200);
            if ( (fabs(eta) > 1.5 && fabs(eta2) < 1.5) || (fabs(eta) < 1.5 && fabs(eta2) > 1.5) )
              plot1D("h_selectedMassSeedEBEE", massSeedCl.M(),  1, h_1d, "M_ee [GeV]", 200, 0, 200);
            if (fabs(eta) > 1.5 || fabs(eta2) > 1.5) plot1D("h_selectedMassSeedEE", massSeedCl.M(),  1, h_1d, "M_ee [GeV]", 200, 0, 200);
            
            
            plot1D("h_selectedMassSCraw", massSCraw.M(),  1, h_1d, "M_ee [GeV]", 200, 0, 200);
            plot1D(("h_selectedMassSCraw"+currfile).Data(), massSCraw.M(),  1, h_1d, "M_ee [GeV]", 200, 0, 200);
            if (fabs(eta) < 1.5 && fabs(eta2) < 1.5) plot1D("h_selectedMassSCrawEBEB", massSCraw.M(),  1, h_1d, "M_ee [GeV]", 200, 0, 200);
            if (fabs(eta) > 1.5 || fabs(eta2) > 1.5) plot1D("h_selectedMassSCrawEE", massSCraw.M(),  1, h_1d, "M_ee [GeV]", 200, 0, 200);
            
            
            plot1D("h_selectedMassSC", massSC.M(),  1, h_1d, "M_ee [GeV]", 200, 0, 200);
            if (fabs(eta) < 1.5 && fabs(eta2) < 1.5) plot1D("h_selectedMassSCEBEB", massSC.M(),  1, h_1d, "M_ee [GeV]", 200, 0, 200);
            if (fabs(eta) > 1.5 || fabs(eta2) > 1.5) plot1D("h_selectedMassSCEE", massSC.M(),  1, h_1d, "M_ee [GeV]", 200, 0, 200);
            
            
            //cout<<"found interesting run:ls:event "<<evt_run()<<":"<<evt_lumiBlock()<<":"<<evt_event();
            //cout<<" with invariant mass "<<mass.M()<<endl;
            if (massSeedCl.M() > 70 && massSeedCl.M() < 110 ) {
              plot1D("h_Zrun", evt_run(),  1, h_1d, "Run", 1200, 246900, 248100);
              plot1D(("h_Zrunfile"+currfile).Data(), evt_run(),  1, h_1d, "Run", 1200, 246900, 248100);
              plot1D("h_Znvtx", nvtx,  1, h_1d, "NVTX", 50, 0, 50);
              plot1D(("h_Znvtx"+currfile).Data(), nvtx,  1, h_1d, "NVTX", 50, 0, 50);
              //cout<<"run:ls:event "<<evt_run()<<":"<<evt_lumiBlock()<<":"<<evt_event()<<endl;
              
            }
            
            
          }
          if (selected1 && selected2loose) {
            plot1D("h_halfSelectedMass", massSeedCl.M(),  1, h_1d, "M_ee [GeV]", 200, 0, 200);
            if (fabs(eta) < 1.5 && fabs(eta2) < 1.5) plot1D("h_halfSelectedMassEBEB", massSeedCl.M(),  1, h_1d, "M_ee [GeV]", 200, 0, 200);
            if (fabs(eta) > 1.5 || fabs(eta2) > 1.5) plot1D("h_halfSelectedMassEE", massSeedCl.M(),  1, h_1d, "M_ee [GeV]", 200, 0, 200);
          }
          
          plot1D("h_unselectedMass", massSeedCl.M(),  1, h_1d, "M_ee [GeV]", 200, 0, 200);
          
        } // end of second electron loop
        
        
        
      }// end of electron loop
      
      
    } // end of event loop
    
    cout<<"Trigger report: "<<endl;
    cout<<"trigEG5   :" <<trigEG5  <<endl;
    cout<<"trigEG20  :" <<trigEG20 <<endl;
    cout<<"trigDiEG  :" <<trigDiEG <<endl;
    cout<<"trigDiSCJPsi    " <<trigDiSCJPsi   <<endl;
    cout<<"trigDiSCUpsilon " <<trigDiSCUpsilon<<endl;
    cout<<"trigEleSCJPsi   " <<trigEleSCJPsi  <<endl;
    
    
    
    //outfile_->cd();
    // Clean Up
    delete tree;
    file->Close();
    delete file;
  }
  
  savePlots(h_1d, outName);
  
  if ( nEventsChain != nEventsTotal ) {
    cout << Form( "ERROR: number of events from files (%d) is not equal to total number of events (%d)", nEventsChain, nEventsTotal ) << endl;
  }
  
  // Example Histograms
  //samplehisto->Draw();
  
  // return
  bmark->Stop("benchmark");
  cout << endl;
  cout << nEventsTotal << " Events Processed" << endl;
  cout << "------------------------------" << endl;
  cout << "CPU  Time:	" << Form( "%.01f", bmark->GetCpuTime("benchmark")  ) << endl;
  cout << "Real Time:	" << Form( "%.01f", bmark->GetRealTime("benchmark") ) << endl;
  cout << endl;
  delete bmark;
  return 0;
}
