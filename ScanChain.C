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

// CMS3
#include "CMS3.cc"
#include "PlotUtilities.h"

using namespace std;
using namespace tas;

//int ScanChain( TChain* chain, bool fast = true, int nEvents = -1, string skimFilePrefix = "test") {
int ScanChain( TChain* chain, int nEvents = -1, const char* outName = "test") {
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
  std::map<std::string, TH1*> h_1d;

  int trigEG5 = 0;
  int trigEG20 = 0;
  int trigDiEG = 0;
  int trigDiSCJPsi = 0;
  int trigDiSCUpsilon = 0;
  int trigEleSCJPsi = 0;


  // File Loop
  while ( (currentFile = (TFile*)fileIter.Next()) ) {
    cout<<"Looking at file: "<<currentFile->GetTitle()<<endl;
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


      // Analysis Code
      for (unsigned int iel = 0; iel < els_p4().size(); ++iel) {
	LorentzVector el_p4 = els_p4().at(iel);
	float pt = el_p4.pt();
	float eta = els_etaSC().at(iel);
	samplehisto->Fill(els_p4().at(iel).pt());
	
	float theta = 2.0*TMath::ATan(TMath::Exp(-1.*eta));
	  
	float SCEt = els_eSC().at(iel) * TMath::Sin(theta);

	plot1D("h_pt", pt,  1, h_1d, "pT [GeV]", 150, 0, 150);
	if (fabs(eta) > 2.5 || (pt < 5 && SCEt < 5) ) continue;
	//	if (fabs(eta) > 2.5 || (pt < 10) ) continue;
	float sieie =  els_sigmaIEtaIEta_full5x5().at(iel);
	float deta = els_dEtaIn().at(iel);
	float dphi =  els_dPhiIn().at(iel);
	float detaSeed = els_dEtaOut().at(iel);
	float dphiSeed =  els_dPhiOut().at(iel);
	float hovere = els_hOverE().at(iel);
	float chiso     = els_pfChargedHadronIso().at(iel);
	float nhiso     = els_pfNeutralHadronIso().at(iel);
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
	LorentzVector el_seedCl = el_p4;
	//cout<<"Start with "<<el_seedCl<<endl;
	el_seedCl *= seedE/el_p4.E();
	//cout<<"turn into "<<el_seedCl<<endl;
	LorentzVector el_SC = el_p4;
	//el_SC.SetE(els_eSC().at(iel));
	el_SC *= els_eSC().at(iel) / el_p4.E();

//	cout<<"seed cluster: E, eta, theta, Et: "<<seedE<<" "<<eta<<" "<<theta<<" "<<seedEt<<endl;
//	cout<<"supercluster: E, eta, theta, Et: "<<els_eSC().at(iel)<<" "<<eta<<" "<<theta<<" "<<SCEt<<endl;

	//if (pt > 150) cout<<"found 150 GeV event/ls/run: "<<evt_event()<<" "<<evt_lumiBlock()<<" "<<evt_run()<<endl;
	
	plot1D("h_trkdr", trkdr,  1, h_1d, "DR (KF,GSF)", 100, -1.1, 0.5);
	plot1D("h_trkshFrac", trkshFrac,  1, h_1d, "Shared Fraction", 100, -1.1, 1.1);

	plot1D("h_ElePtMinusSCEt", (pt-SCEt)/SCEt,  1, h_1d, "(pT(ele) - ET(SC)) / ET(SC)", 100, -1.1, 1.1);
	plot1D("h_SCEtMinusSeedEt", (SCEt-seedEt)/seedEt,  1, h_1d, "(ET(SC) - ET(seed)) / ET(seed)", 100, -1.1, 1.1);

	
	if (fabs(eta) > 1.5) {
	  plot1D("h_relchisoEE", chiso/pt,  1, h_1d, "PFCh", 100, 0, 1);
	  plot1D("h_relemisoEE", emiso/pt,  1, h_1d, "PFEM", 100, 0, 1);
	  plot1D("h_relnhisoEE", nhiso/pt,  1, h_1d, "PFNh", 100, 0, 1);
	  
	  plot1D("h_sieieEE", els_sigmaIEtaIEta_full5x5().at(iel),  1, h_1d, "sieie_5x5", 100, 0, 0.05);
	  plot1D("h_sipipEE", els_sigmaIPhiIPhi_full5x5().at(iel),  1, h_1d, "sipip_5x5", 100, 0, 0.1);

	  if (eta > 0) {
	    plot1D("h_detaEEpos", deta,  1, h_1d, "DeltaEta", 50, -0.1, 0.1);
	    plot1D("h_dphiEEpos", dphi,  1, h_1d, "DeltaPhi", 50, -0.1, 0.1);
	    plot1D("h_detaSeedEEpos", detaSeed,  1, h_1d, "DeltaEta", 50, -0.1, 0.1);
	    plot1D("h_dphiSeedEEpos", dphiSeed,  1, h_1d, "DeltaPhi", 50, -0.1, 0.1);
	  }
	  if (eta < 0) {
	    plot1D("h_detaEEneg", deta,  1, h_1d, "DeltaEta", 50, -0.1, 0.1);
	    plot1D("h_dphiEEneg", dphi,  1, h_1d, "DeltaPhi", 50, -0.1, 0.1);
	    plot1D("h_detaSeedEEneg", detaSeed,  1, h_1d, "DeltaEta", 50, -0.1, 0.1);
	    plot1D("h_dphiSeedEEneg", dphiSeed,  1, h_1d, "DeltaPhi", 50, -0.1, 0.1);

	  }
	  
	}
	else {
	  plot1D("h_relchisoEB", chiso/pt,  1, h_1d, "PFCh", 100, 0, 1);
	  plot1D("h_relemisoEB", emiso/pt,  1, h_1d, "PFEM", 100, 0, 1);
	  plot1D("h_relnhisoEB", nhiso/pt,  1, h_1d, "PFNh", 100, 0, 1);

	  plot1D("h_sieieEB", els_sigmaIEtaIEta_full5x5().at(iel),  1, h_1d, "sieie_5x5", 100, 0, 0.05);
	  plot1D("h_sipipEB", els_sigmaIPhiIPhi_full5x5().at(iel),  1, h_1d, "sipip_5x5", 100, 0, 0.1);

	  plot1D("h_detaEB", deta,  1, h_1d, "DeltaEta", 100, -0.1, 0.1);
	  plot1D("h_dphiEB", dphi,  1, h_1d, "DeltaPhi", 100, -0.1, 0.1);
	    

	}


	bool selected1 = true;
	if (caloiso/pt > 0.2) selected1 = false;
//	if (fabs(eta) < 1.5 && (fabs(deta) > 0.005 || fabs(dphi) > 0.005 || hovere > 0.005 || sieie > 0.012)  ) selected1 = false;
//	if (fabs(eta) > 1.5 && (fabs(deta) > 0.01  || fabs(dphi) > 0.005 || hovere > 0.005 || sieie > 0.027)  ) selected1 = false;
	if (fabs(eta) < 1.5 && (fabs(deta) > 0.01  || fabs(dphi) > 0.01 || hovere > 0.05 || sieie > 0.015)  ) selected1 = false;
	if (fabs(eta) > 1.5 && (fabs(deta) > 0.01  || fabs(dphi) > 0.01 || hovere > 0.05 || sieie > 0.030)  ) selected1 = false;

	if (selected1) plot1D("h_ptSelected", pt,  1, h_1d, "pT [GeV]", 150, 0, 150);	
	if (selected1 && passtrigEG5   ) plot1D("h_ptSelectedtrigEG5 ", pt,  1, h_1d, "pT [GeV]", 150, 0, 150);	
	if (selected1 && passtrigEG20  ) plot1D("h_ptSelectedtrigEG20", pt,  1, h_1d, "pT [GeV]", 150, 0, 150);	
	if (selected1 && passtrigDiEG  ) plot1D("h_ptSelectedtrigDiEG", pt,  1, h_1d, "pT [GeV]", 150, 0, 150);	

	for (unsigned int jel = iel+1; jel < els_p4().size(); ++jel) {
	  LorentzVector el2_p4 = els_p4().at(jel);
	  float pt2 = el2_p4.pt();
	  float eta2 = els_etaSC().at(jel);
	  if (fabs(eta2) > 2.5 || pt2 < 5 ) continue;
	  float sieie2 =  els_sigmaIEtaIEta_full5x5().at(jel);
	  float deta2 = fabs(els_dEtaIn().at(jel));
	  float dphi2 =  fabs(els_dPhiIn().at(jel));
	  float hovere2 = els_hOverE().at(jel);
	  float chiso2     = els_pfChargedHadronIso().at(jel);
	  float nhiso2     = els_pfNeutralHadronIso().at(jel);
	  float emiso2     = els_pfPhotonIso().at(jel);
	  float deltaBeta2 = els_pfPUIso().at(jel);
	  float iso2 = chiso2 + nhiso2 + emiso2 - 0.5* deltaBeta2;
	  float caloiso2 = nhiso2 + emiso2 ;
	  LorentzVector el2_seedCl = el2_p4;
	  //el2_seedCl.SetE( els_eSeed().at(jel));
	  el2_seedCl *= els_eSeed().at(jel)/el2_p4.E();

	  LorentzVector el2_SC = el2_p4;
	  //el2_SC.SetE(els_eSC().at(jel));
	  el2_SC *= els_eSC().at(jel) / el2_p4.E();
	  

	  bool selected2 = true;
	  if (caloiso2/pt > 0.2) selected2 = false;
//	  if (fabs(eta2) < 1.5 && (fabs(deta2) > 0.005 || fabs(dphi2) > 0.005 || hovere2 > 0.005 || sieie2 > 0.012)  ) selected2 = false;
//	  if (fabs(eta2) > 1.5 && (fabs(deta2) > 0.01 ||  fabs(dphi2) > 0.005 || hovere2 > 0.005 || sieie2 > 0.027)  )  selected2 = false;
	if (fabs(eta2) < 1.5 && (fabs(deta2) > 0.01  || fabs(dphi2) > 0.01 || hovere2 > 0.05 || sieie2 > 0.015)  ) selected2 = false;
	if (fabs(eta2) > 1.5 && (fabs(deta2) > 0.01  || fabs(dphi2) > 0.01 || hovere2 > 0.05 || sieie2 > 0.030)  ) selected2 = false;


	  if (selected1 && selected2) plot1D("h_pt2Selected12", pt2,  1, h_1d, "pT [GeV]", 150, 0, 150);	

	  LorentzVector mass = el_p4 + el2_p4;
	  LorentzVector massSC = el_SC + el2_SC;
	  LorentzVector massSeedCl = el_seedCl + el2_seedCl;
	  //cout<<"mass, massSC, massSeedCl: "<<mass.M()<<" "<<massSC.M()<<" "<<massSeedCl.M()<<endl;
	  
	  if (selected1 && selected2) {
	    plot1D("h_selectedMass", mass.M(),  1, h_1d, "M_ee [GeV]", 1200, 0, 150);	
	    plot1D("h_selectedMassWideBins", mass.M(),  1, h_1d, "M_ee [GeV]", 150, 0, 150);	
	    plot1D("h_selectedMassSeed", massSeedCl.M(),  1, h_1d, "M_ee [GeV]", 1200, 0, 150);	
	    plot1D("h_selectedMassSC", massSC.M(),  1, h_1d, "M_ee [GeV]", 1200, 0, 150);	
	    if (passtrigEG5    ) plot1D("h_selectedMasstrigEG5 ", mass.M(),  1, h_1d, "M_ee [GeV]", 600, 0, 150);	
	    if (passtrigEG20   ) plot1D("h_selectedMasstrigEG20", mass.M(),  1, h_1d, "M_ee [GeV]", 600, 0, 150);	
	    if (passtrigDiEG   ) plot1D("h_selectedMasstrigDiEG", mass.M(),  1, h_1d, "M_ee [GeV]", 600, 0, 150);      

	    cout<<"found interesting run:ls:event "<<evt_run()<<":"<<evt_lumiBlock()<<":"<<evt_event(); 
	    cout<<" with invariant mass "<<mass.M()<<endl;
	  }
	  if (selected1 || selected2) plot1D("h_halfSelectedMass", mass.M(),  1, h_1d, "M_ee [GeV]", 1200, 0, 150);	
	  plot1D("h_unselectedMass", mass.M(),  1, h_1d, "M_ee [GeV]", 1200, 0, 150);	

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
