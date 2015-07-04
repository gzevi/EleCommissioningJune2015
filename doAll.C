{

  gROOT->ProcessLine(".L PlotUtilities.cc+");
  gROOT->ProcessLine(".L ScanChain.C+");

  TChain *ch7 = new TChain("Events");
//  ch7->Add("ntuple248036SingleAndDiEle.root");
  ch7->Add("ntuple17JuneEGammaDilepton.root");
//  ch7->Add("ntuple18HLTPhysicsDielectron.root");
//  ch7->Add("ntuple17JuneEGMLowPUDilepton.root");
  ScanChain(ch7, -1, "plots.root");

}
