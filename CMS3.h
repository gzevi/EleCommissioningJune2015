// -*- C++ -*-
#ifndef CMS3_H
#define CMS3_H
#include "Math/LorentzVector.h"
#include "Math/Point3D.h"
#include "TMath.h"
#include "TBranch.h"
#include "TTree.h"
#include "TH1F.h"
#include "TFile.h"
#include "TBits.h"
#include <vector> 
#include <unistd.h> 
typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > LorentzVector;

using namespace std; 
class CMS3 {
private: 
protected: 
	unsigned int index;
	TBits hlt_bits_;
	TBranch *hlt_bits_branch;
	bool hlt_bits_isLoaded;
	vector<TString> evt_CMS3tag_;
	TBranch *evt_CMS3tag_branch;
	bool evt_CMS3tag_isLoaded;
	vector<TString> evt_dataset_;
	TBranch *evt_dataset_branch;
	bool evt_dataset_isLoaded;
	vector<TString> hlt_trigNames_;
	TBranch *hlt_trigNames_branch;
	bool hlt_trigNames_isLoaded;
	vector<bool> els_conv_vtx_flag_;
	TBranch *els_conv_vtx_flag_branch;
	bool els_conv_vtx_flag_isLoaded;
	vector<bool> els_isGsfCtfScPixChargeConsistent_;
	TBranch *els_isGsfCtfScPixChargeConsistent_branch;
	bool els_isGsfCtfScPixChargeConsistent_isLoaded;
	vector<bool> els_passingMvaPreselection_;
	TBranch *els_passingMvaPreselection_branch;
	bool els_passingMvaPreselection_isLoaded;
	vector<bool> els_passingPflowPreselection_;
	TBranch *els_passingPflowPreselection_branch;
	bool els_passingPflowPreselection_isLoaded;
	vector<bool> photons_haspixelSeed_;
	TBranch *photons_haspixelSeed_branch;
	bool photons_haspixelSeed_isLoaded;
	float evt_bs_Xwidth_;
	TBranch *evt_bs_Xwidth_branch;
	bool evt_bs_Xwidth_isLoaded;
	float evt_bs_XwidthErr_;
	TBranch *evt_bs_XwidthErr_branch;
	bool evt_bs_XwidthErr_isLoaded;
	float evt_bs_Ywidth_;
	TBranch *evt_bs_Ywidth_branch;
	bool evt_bs_Ywidth_isLoaded;
	float evt_bs_YwidthErr_;
	TBranch *evt_bs_YwidthErr_branch;
	bool evt_bs_YwidthErr_isLoaded;
	float evt_bs_dxdz_;
	TBranch *evt_bs_dxdz_branch;
	bool evt_bs_dxdz_isLoaded;
	float evt_bs_dxdzErr_;
	TBranch *evt_bs_dxdzErr_branch;
	bool evt_bs_dxdzErr_isLoaded;
	float evt_bs_dydz_;
	TBranch *evt_bs_dydz_branch;
	bool evt_bs_dydz_isLoaded;
	float evt_bs_dydzErr_;
	TBranch *evt_bs_dydzErr_branch;
	bool evt_bs_dydzErr_isLoaded;
	float evt_bs_sigmaZ_;
	TBranch *evt_bs_sigmaZ_branch;
	bool evt_bs_sigmaZ_isLoaded;
	float evt_bs_sigmaZErr_;
	TBranch *evt_bs_sigmaZErr_branch;
	bool evt_bs_sigmaZErr_isLoaded;
	float evt_bs_xErr_;
	TBranch *evt_bs_xErr_branch;
	bool evt_bs_xErr_isLoaded;
	float evt_bs_yErr_;
	TBranch *evt_bs_yErr_branch;
	bool evt_bs_yErr_isLoaded;
	float evt_bs_zErr_;
	TBranch *evt_bs_zErr_branch;
	bool evt_bs_zErr_isLoaded;
	float evt_bField_;
	TBranch *evt_bField_branch;
	bool evt_bField_isLoaded;
	float evt_fixgrid_all_rho_;
	TBranch *evt_fixgrid_all_rho_branch;
	bool evt_fixgrid_all_rho_isLoaded;
	float evt_fixgridfastjet_allcalo_rho_;
	TBranch *evt_fixgridfastjet_allcalo_rho_branch;
	bool evt_fixgridfastjet_allcalo_rho_isLoaded;
	float evt_fixgridfastjet_all_rho_;
	TBranch *evt_fixgridfastjet_all_rho_branch;
	bool evt_fixgridfastjet_all_rho_isLoaded;
	float evt_fixgridfastjet_centralcalo_rho_;
	TBranch *evt_fixgridfastjet_centralcalo_rho_branch;
	bool evt_fixgridfastjet_centralcalo_rho_isLoaded;
	float evt_fixgridfastjet_centralchargedpileup_rho_;
	TBranch *evt_fixgridfastjet_centralchargedpileup_rho_branch;
	bool evt_fixgridfastjet_centralchargedpileup_rho_isLoaded;
	float evt_fixgridfastjet_centralneutral_rho_;
	TBranch *evt_fixgridfastjet_centralneutral_rho_branch;
	bool evt_fixgridfastjet_centralneutral_rho_isLoaded;
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >  evt_bsp4_;
	TBranch *evt_bsp4_branch;
	bool evt_bsp4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > els_mc_patMatch_p4_;
	TBranch *els_mc_patMatch_p4_branch;
	bool els_mc_patMatch_p4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > els_p4_;
	TBranch *els_p4_branch;
	bool els_p4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > els_p4In_;
	TBranch *els_p4In_branch;
	bool els_p4In_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > els_p4Out_;
	TBranch *els_p4Out_branch;
	bool els_p4Out_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > els_trk_p4_;
	TBranch *els_trk_p4_branch;
	bool els_trk_p4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > els_trk_vertex_p4_;
	TBranch *els_trk_vertex_p4_branch;
	bool els_trk_vertex_p4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > els_vertex_p4_;
	TBranch *els_vertex_p4_branch;
	bool els_vertex_p4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > photons_p4_;
	TBranch *photons_p4_branch;
	bool photons_p4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > vtxs_position_;
	TBranch *vtxs_position_branch;
	bool vtxs_position_isLoaded;
	vector<vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > > hlt_trigObjs_p4_;
	TBranch *hlt_trigObjs_p4_branch;
	bool hlt_trigObjs_p4_isLoaded;
	vector<float> evt_bs_covMatrix_;
	TBranch *evt_bs_covMatrix_branch;
	bool evt_bs_covMatrix_isLoaded;
	vector<float> els_bs2d_;
	TBranch *els_bs2d_branch;
	bool els_bs2d_isLoaded;
	vector<float> els_bs2derr_;
	TBranch *els_bs2derr_branch;
	bool els_bs2derr_isLoaded;
	vector<float> els_bs3d_;
	TBranch *els_bs3d_branch;
	bool els_bs3d_isLoaded;
	vector<float> els_bs3derr_;
	TBranch *els_bs3derr_branch;
	bool els_bs3derr_isLoaded;
	vector<float> els_chi2_;
	TBranch *els_chi2_branch;
	bool els_chi2_isLoaded;
	vector<float> els_ckf_chi2_;
	TBranch *els_ckf_chi2_branch;
	bool els_ckf_chi2_isLoaded;
	vector<float> els_ckf_ndof_;
	TBranch *els_ckf_ndof_branch;
	bool els_ckf_ndof_isLoaded;
	vector<float> els_d0_;
	TBranch *els_d0_branch;
	bool els_d0_isLoaded;
	vector<float> els_d0Err_;
	TBranch *els_d0Err_branch;
	bool els_d0Err_isLoaded;
	vector<float> els_d0corr_;
	TBranch *els_d0corr_branch;
	bool els_d0corr_isLoaded;
	vector<float> els_d0corrPhi_;
	TBranch *els_d0corrPhi_branch;
	bool els_d0corrPhi_isLoaded;
	vector<float> els_d0phiCov_;
	TBranch *els_d0phiCov_branch;
	bool els_d0phiCov_isLoaded;
	vector<float> els_dEtaIn_;
	TBranch *els_dEtaIn_branch;
	bool els_dEtaIn_isLoaded;
	vector<float> els_dEtaOut_;
	TBranch *els_dEtaOut_branch;
	bool els_dEtaOut_isLoaded;
	vector<float> els_dPhiIn_;
	TBranch *els_dPhiIn_branch;
	bool els_dPhiIn_isLoaded;
	vector<float> els_dPhiInPhiOut_;
	TBranch *els_dPhiInPhiOut_branch;
	bool els_dPhiInPhiOut_isLoaded;
	vector<float> els_dPhiOut_;
	TBranch *els_dPhiOut_branch;
	bool els_dPhiOut_isLoaded;
	vector<float> els_deltaEtaEleClusterTrackAtCalo_;
	TBranch *els_deltaEtaEleClusterTrackAtCalo_branch;
	bool els_deltaEtaEleClusterTrackAtCalo_isLoaded;
	vector<float> els_deltaPhiEleClusterTrackAtCalo_;
	TBranch *els_deltaPhiEleClusterTrackAtCalo_branch;
	bool els_deltaPhiEleClusterTrackAtCalo_isLoaded;
	vector<float> els_dxyPV_;
	TBranch *els_dxyPV_branch;
	bool els_dxyPV_isLoaded;
	vector<float> els_dzPV_;
	TBranch *els_dzPV_branch;
	bool els_dzPV_isLoaded;
	vector<float> els_e1x5_;
	TBranch *els_e1x5_branch;
	bool els_e1x5_isLoaded;
	vector<float> els_e1x5_full5x5_;
	TBranch *els_e1x5_full5x5_branch;
	bool els_e1x5_full5x5_isLoaded;
	vector<float> els_e2x5Max_;
	TBranch *els_e2x5Max_branch;
	bool els_e2x5Max_isLoaded;
	vector<float> els_e2x5Max_full5x5_;
	TBranch *els_e2x5Max_full5x5_branch;
	bool els_e2x5Max_full5x5_isLoaded;
	vector<float> els_e5x5_;
	TBranch *els_e5x5_branch;
	bool els_e5x5_isLoaded;
	vector<float> els_e5x5_full5x5_;
	TBranch *els_e5x5_full5x5_branch;
	bool els_e5x5_full5x5_isLoaded;
	vector<float> els_eOverPIn_;
	TBranch *els_eOverPIn_branch;
	bool els_eOverPIn_isLoaded;
	vector<float> els_eOverPOut_;
	TBranch *els_eOverPOut_branch;
	bool els_eOverPOut_isLoaded;
	vector<float> els_eSC_;
	TBranch *els_eSC_branch;
	bool els_eSC_isLoaded;
	vector<float> els_eSCPresh_;
	TBranch *els_eSCPresh_branch;
	bool els_eSCPresh_isLoaded;
	vector<float> els_eSCRaw_;
	TBranch *els_eSCRaw_branch;
	bool els_eSCRaw_isLoaded;
	vector<float> els_eSeed_;
	TBranch *els_eSeed_branch;
	bool els_eSeed_isLoaded;
	vector<float> els_eSeedOverPIn_;
	TBranch *els_eSeedOverPIn_branch;
	bool els_eSeedOverPIn_isLoaded;
	vector<float> els_eSeedOverPOut_;
	TBranch *els_eSeedOverPOut_branch;
	bool els_eSeedOverPOut_isLoaded;
	vector<float> els_ecalEnergy_;
	TBranch *els_ecalEnergy_branch;
	bool els_ecalEnergy_isLoaded;
	vector<float> els_ecalEnergyError_;
	TBranch *els_ecalEnergyError_branch;
	bool els_ecalEnergyError_isLoaded;
	vector<float> els_ecalIso_;
	TBranch *els_ecalIso_branch;
	bool els_ecalIso_isLoaded;
	vector<float> els_ecalIso04_;
	TBranch *els_ecalIso04_branch;
	bool els_ecalIso04_isLoaded;
	vector<float> els_ecalPFClusterIso_;
	TBranch *els_ecalPFClusterIso_branch;
	bool els_ecalPFClusterIso_isLoaded;
	vector<float> els_etaErr_;
	TBranch *els_etaErr_branch;
	bool els_etaErr_isLoaded;
	vector<float> els_etaSC_;
	TBranch *els_etaSC_branch;
	bool els_etaSC_isLoaded;
	vector<float> els_etaSCwidth_;
	TBranch *els_etaSCwidth_branch;
	bool els_etaSCwidth_isLoaded;
	vector<float> els_fbrem_;
	TBranch *els_fbrem_branch;
	bool els_fbrem_isLoaded;
	vector<float> els_hOverE_;
	TBranch *els_hOverE_branch;
	bool els_hOverE_isLoaded;
	vector<float> els_hOverEBC_;
	TBranch *els_hOverEBC_branch;
	bool els_hOverEBC_isLoaded;
	vector<float> els_hcalDepth1OverEcal_;
	TBranch *els_hcalDepth1OverEcal_branch;
	bool els_hcalDepth1OverEcal_isLoaded;
	vector<float> els_hcalDepth1TowerSumEt_;
	TBranch *els_hcalDepth1TowerSumEt_branch;
	bool els_hcalDepth1TowerSumEt_isLoaded;
	vector<float> els_hcalDepth1TowerSumEt04_;
	TBranch *els_hcalDepth1TowerSumEt04_branch;
	bool els_hcalDepth1TowerSumEt04_isLoaded;
	vector<float> els_hcalDepth2OverEcal_;
	TBranch *els_hcalDepth2OverEcal_branch;
	bool els_hcalDepth2OverEcal_isLoaded;
	vector<float> els_hcalDepth2TowerSumEt_;
	TBranch *els_hcalDepth2TowerSumEt_branch;
	bool els_hcalDepth2TowerSumEt_isLoaded;
	vector<float> els_hcalDepth2TowerSumEt04_;
	TBranch *els_hcalDepth2TowerSumEt04_branch;
	bool els_hcalDepth2TowerSumEt04_isLoaded;
	vector<float> els_hcalIso_;
	TBranch *els_hcalIso_branch;
	bool els_hcalIso_isLoaded;
	vector<float> els_hcalIso04_;
	TBranch *els_hcalIso04_branch;
	bool els_hcalIso04_isLoaded;
	vector<float> els_hcalPFClusterIso_;
	TBranch *els_hcalPFClusterIso_branch;
	bool els_hcalPFClusterIso_isLoaded;
	vector<float> els_ip2d_;
	TBranch *els_ip2d_branch;
	bool els_ip2d_isLoaded;
	vector<float> els_ip2derr_;
	TBranch *els_ip2derr_branch;
	bool els_ip2derr_isLoaded;
	vector<float> els_ip3d_;
	TBranch *els_ip3d_branch;
	bool els_ip3d_isLoaded;
	vector<float> els_ip3derr_;
	TBranch *els_ip3derr_branch;
	bool els_ip3derr_isLoaded;
	vector<float> els_mass_;
	TBranch *els_mass_branch;
	bool els_mass_isLoaded;
	vector<float> els_mc_patMatch_dr_;
	TBranch *els_mc_patMatch_dr_branch;
	bool els_mc_patMatch_dr_isLoaded;
	vector<float> els_miniIso_ch_;
	TBranch *els_miniIso_ch_branch;
	bool els_miniIso_ch_isLoaded;
	vector<float> els_miniIso_db_;
	TBranch *els_miniIso_db_branch;
	bool els_miniIso_db_isLoaded;
	vector<float> els_miniIso_em_;
	TBranch *els_miniIso_em_branch;
	bool els_miniIso_em_isLoaded;
	vector<float> els_miniIso_nh_;
	TBranch *els_miniIso_nh_branch;
	bool els_miniIso_nh_isLoaded;
	vector<float> els_miniIso_uncor_;
	TBranch *els_miniIso_uncor_branch;
	bool els_miniIso_uncor_isLoaded;
	vector<float> els_mva_;
	TBranch *els_mva_branch;
	bool els_mva_isLoaded;
	vector<float> els_ndof_;
	TBranch *els_ndof_branch;
	bool els_ndof_isLoaded;
	vector<float> els_pfChargedHadronIso_;
	TBranch *els_pfChargedHadronIso_branch;
	bool els_pfChargedHadronIso_isLoaded;
	vector<float> els_pfNeutralHadronIso_;
	TBranch *els_pfNeutralHadronIso_branch;
	bool els_pfNeutralHadronIso_isLoaded;
	vector<float> els_pfPUIso_;
	TBranch *els_pfPUIso_branch;
	bool els_pfPUIso_isLoaded;
	vector<float> els_pfPhotonIso_;
	TBranch *els_pfPhotonIso_branch;
	bool els_pfPhotonIso_isLoaded;
	vector<float> els_phiErr_;
	TBranch *els_phiErr_branch;
	bool els_phiErr_isLoaded;
	vector<float> els_phiSC_;
	TBranch *els_phiSC_branch;
	bool els_phiSC_isLoaded;
	vector<float> els_phiSCwidth_;
	TBranch *els_phiSCwidth_branch;
	bool els_phiSCwidth_isLoaded;
	vector<float> els_ptErr_;
	TBranch *els_ptErr_branch;
	bool els_ptErr_isLoaded;
	vector<float> els_ptErrGsf_;
	TBranch *els_ptErrGsf_branch;
	bool els_ptErrGsf_isLoaded;
	vector<float> els_r9_;
	TBranch *els_r9_branch;
	bool els_r9_isLoaded;
	vector<float> els_r9_full5x5_;
	TBranch *els_r9_full5x5_branch;
	bool els_r9_full5x5_isLoaded;
	vector<float> els_sigmaEtaEta_;
	TBranch *els_sigmaEtaEta_branch;
	bool els_sigmaEtaEta_isLoaded;
	vector<float> els_sigmaEtaEta_full5x5_;
	TBranch *els_sigmaEtaEta_full5x5_branch;
	bool els_sigmaEtaEta_full5x5_isLoaded;
	vector<float> els_sigmaIEtaIEta_;
	TBranch *els_sigmaIEtaIEta_branch;
	bool els_sigmaIEtaIEta_isLoaded;
	vector<float> els_sigmaIEtaIEta_full5x5_;
	TBranch *els_sigmaIEtaIEta_full5x5_branch;
	bool els_sigmaIEtaIEta_full5x5_isLoaded;
	vector<float> els_sigmaIPhiIPhi_;
	TBranch *els_sigmaIPhiIPhi_branch;
	bool els_sigmaIPhiIPhi_isLoaded;
	vector<float> els_sigmaIPhiIPhi_full5x5_;
	TBranch *els_sigmaIPhiIPhi_full5x5_branch;
	bool els_sigmaIPhiIPhi_full5x5_isLoaded;
	vector<float> els_sigmaIphiIphi_;
	TBranch *els_sigmaIphiIphi_branch;
	bool els_sigmaIphiIphi_isLoaded;
	vector<float> els_tkIso_;
	TBranch *els_tkIso_branch;
	bool els_tkIso_isLoaded;
	vector<float> els_tkIso04_;
	TBranch *els_tkIso04_branch;
	bool els_tkIso04_isLoaded;
	vector<float> els_trackMomentumError_;
	TBranch *els_trackMomentumError_branch;
	bool els_trackMomentumError_isLoaded;
	vector<float> els_trkdr_;
	TBranch *els_trkdr_branch;
	bool els_trkdr_isLoaded;
	vector<float> els_trkshFrac_;
	TBranch *els_trkshFrac_branch;
	bool els_trkshFrac_isLoaded;
	vector<float> els_z0_;
	TBranch *els_z0_branch;
	bool els_z0_isLoaded;
	vector<float> els_z0Err_;
	TBranch *els_z0Err_branch;
	bool els_z0Err_isLoaded;
	vector<float> els_z0corr_;
	TBranch *els_z0corr_branch;
	bool els_z0corr_isLoaded;
	vector<float> photons_chargedHadronIso_;
	TBranch *photons_chargedHadronIso_branch;
	bool photons_chargedHadronIso_isLoaded;
	vector<float> photons_e1x5_;
	TBranch *photons_e1x5_branch;
	bool photons_e1x5_isLoaded;
	vector<float> photons_e2x5Max_;
	TBranch *photons_e2x5Max_branch;
	bool photons_e2x5Max_isLoaded;
	vector<float> photons_e3x3_;
	TBranch *photons_e3x3_branch;
	bool photons_e3x3_isLoaded;
	vector<float> photons_e5x5_;
	TBranch *photons_e5x5_branch;
	bool photons_e5x5_isLoaded;
	vector<float> photons_eSC_;
	TBranch *photons_eSC_branch;
	bool photons_eSC_isLoaded;
	vector<float> photons_eSCPresh_;
	TBranch *photons_eSCPresh_branch;
	bool photons_eSCPresh_isLoaded;
	vector<float> photons_eSCRaw_;
	TBranch *photons_eSCRaw_branch;
	bool photons_eSCRaw_isLoaded;
	vector<float> photons_ecalIso03_;
	TBranch *photons_ecalIso03_branch;
	bool photons_ecalIso03_isLoaded;
	vector<float> photons_ecalIso04_;
	TBranch *photons_ecalIso04_branch;
	bool photons_ecalIso04_isLoaded;
	vector<float> photons_ecalPFClusterIso_;
	TBranch *photons_ecalPFClusterIso_branch;
	bool photons_ecalPFClusterIso_isLoaded;
	vector<float> photons_etaSC_;
	TBranch *photons_etaSC_branch;
	bool photons_etaSC_isLoaded;
	vector<float> photons_full3x3_e3x3_;
	TBranch *photons_full3x3_e3x3_branch;
	bool photons_full3x3_e3x3_isLoaded;
	vector<float> photons_full5x5_e1x5_;
	TBranch *photons_full5x5_e1x5_branch;
	bool photons_full5x5_e1x5_isLoaded;
	vector<float> photons_full5x5_e2x5Max_;
	TBranch *photons_full5x5_e2x5Max_branch;
	bool photons_full5x5_e2x5Max_isLoaded;
	vector<float> photons_full5x5_e5x5_;
	TBranch *photons_full5x5_e5x5_branch;
	bool photons_full5x5_e5x5_isLoaded;
	vector<float> photons_full5x5_hOverE_;
	TBranch *photons_full5x5_hOverE_branch;
	bool photons_full5x5_hOverE_isLoaded;
	vector<float> photons_full5x5_hOverEtowBC_;
	TBranch *photons_full5x5_hOverEtowBC_branch;
	bool photons_full5x5_hOverEtowBC_isLoaded;
	vector<float> photons_full5x5_r9_;
	TBranch *photons_full5x5_r9_branch;
	bool photons_full5x5_r9_isLoaded;
	vector<float> photons_full5x5_sigmaEtaEta_;
	TBranch *photons_full5x5_sigmaEtaEta_branch;
	bool photons_full5x5_sigmaEtaEta_isLoaded;
	vector<float> photons_full5x5_sigmaIEtaIEta_;
	TBranch *photons_full5x5_sigmaIEtaIEta_branch;
	bool photons_full5x5_sigmaIEtaIEta_isLoaded;
	vector<float> photons_hOverE_;
	TBranch *photons_hOverE_branch;
	bool photons_hOverE_isLoaded;
	vector<float> photons_hOverEtowBC_;
	TBranch *photons_hOverEtowBC_branch;
	bool photons_hOverEtowBC_isLoaded;
	vector<float> photons_hcalDepth1TowerSumEtBcConeDR03_;
	TBranch *photons_hcalDepth1TowerSumEtBcConeDR03_branch;
	bool photons_hcalDepth1TowerSumEtBcConeDR03_isLoaded;
	vector<float> photons_hcalDepth1TowerSumEtBcConeDR04_;
	TBranch *photons_hcalDepth1TowerSumEtBcConeDR04_branch;
	bool photons_hcalDepth1TowerSumEtBcConeDR04_isLoaded;
	vector<float> photons_hcalDepth2TowerSumEtBcConeDR03_;
	TBranch *photons_hcalDepth2TowerSumEtBcConeDR03_branch;
	bool photons_hcalDepth2TowerSumEtBcConeDR03_isLoaded;
	vector<float> photons_hcalDepth2TowerSumEtBcConeDR04_;
	TBranch *photons_hcalDepth2TowerSumEtBcConeDR04_branch;
	bool photons_hcalDepth2TowerSumEtBcConeDR04_isLoaded;
	vector<float> photons_hcalIso03_;
	TBranch *photons_hcalIso03_branch;
	bool photons_hcalIso03_isLoaded;
	vector<float> photons_hcalIso04_;
	TBranch *photons_hcalIso04_branch;
	bool photons_hcalIso04_isLoaded;
	vector<float> photons_hcalPFClusterIso_;
	TBranch *photons_hcalPFClusterIso_branch;
	bool photons_hcalPFClusterIso_isLoaded;
	vector<float> photons_hcalTowerSumEtBcConeDR03_;
	TBranch *photons_hcalTowerSumEtBcConeDR03_branch;
	bool photons_hcalTowerSumEtBcConeDR03_isLoaded;
	vector<float> photons_hcalTowerSumEtBcConeDR04_;
	TBranch *photons_hcalTowerSumEtBcConeDR04_branch;
	bool photons_hcalTowerSumEtBcConeDR04_isLoaded;
	vector<float> photons_mass_;
	TBranch *photons_mass_branch;
	bool photons_mass_isLoaded;
	vector<float> photons_neutralHadronIso_;
	TBranch *photons_neutralHadronIso_branch;
	bool photons_neutralHadronIso_isLoaded;
	vector<float> photons_ntkIsoHollow03_;
	TBranch *photons_ntkIsoHollow03_branch;
	bool photons_ntkIsoHollow03_isLoaded;
	vector<float> photons_ntkIsoHollow04_;
	TBranch *photons_ntkIsoHollow04_branch;
	bool photons_ntkIsoHollow04_isLoaded;
	vector<float> photons_ntkIsoSolid03_;
	TBranch *photons_ntkIsoSolid03_branch;
	bool photons_ntkIsoSolid03_isLoaded;
	vector<float> photons_ntkIsoSolid04_;
	TBranch *photons_ntkIsoSolid04_branch;
	bool photons_ntkIsoSolid04_isLoaded;
	vector<float> photons_phiSC_;
	TBranch *photons_phiSC_branch;
	bool photons_phiSC_isLoaded;
	vector<float> photons_photonIso_;
	TBranch *photons_photonIso_branch;
	bool photons_photonIso_isLoaded;
	vector<float> photons_recoChargedHadronIso_;
	TBranch *photons_recoChargedHadronIso_branch;
	bool photons_recoChargedHadronIso_isLoaded;
	vector<float> photons_recoNeutralHadronIso_;
	TBranch *photons_recoNeutralHadronIso_branch;
	bool photons_recoNeutralHadronIso_isLoaded;
	vector<float> photons_recoPhotonIso_;
	TBranch *photons_recoPhotonIso_branch;
	bool photons_recoPhotonIso_isLoaded;
	vector<float> photons_sigmaEtaEta_;
	TBranch *photons_sigmaEtaEta_branch;
	bool photons_sigmaEtaEta_isLoaded;
	vector<float> photons_sigmaIEtaIEta_;
	TBranch *photons_sigmaIEtaIEta_branch;
	bool photons_sigmaIEtaIEta_isLoaded;
	vector<float> photons_tkIsoHollow03_;
	TBranch *photons_tkIsoHollow03_branch;
	bool photons_tkIsoHollow03_isLoaded;
	vector<float> photons_tkIsoHollow04_;
	TBranch *photons_tkIsoHollow04_branch;
	bool photons_tkIsoHollow04_isLoaded;
	vector<float> photons_tkIsoSolid03_;
	TBranch *photons_tkIsoSolid03_branch;
	bool photons_tkIsoSolid03_isLoaded;
	vector<float> photons_tkIsoSolid04_;
	TBranch *photons_tkIsoSolid04_branch;
	bool photons_tkIsoSolid04_isLoaded;
	vector<float> puInfo_trueNumInteractions_;
	TBranch *puInfo_trueNumInteractions_branch;
	bool puInfo_trueNumInteractions_isLoaded;
	vector<float> vtxs_chi2_;
	TBranch *vtxs_chi2_branch;
	bool vtxs_chi2_isLoaded;
	vector<float> vtxs_ndof_;
	TBranch *vtxs_ndof_branch;
	bool vtxs_ndof_isLoaded;
	vector<float> vtxs_score_;
	TBranch *vtxs_score_branch;
	bool vtxs_score_isLoaded;
	vector<float> vtxs_xError_;
	TBranch *vtxs_xError_branch;
	bool vtxs_xError_isLoaded;
	vector<float> vtxs_yError_;
	TBranch *vtxs_yError_branch;
	bool vtxs_yError_isLoaded;
	vector<float> vtxs_zError_;
	TBranch *vtxs_zError_branch;
	bool vtxs_zError_isLoaded;
	vector<vector<float> > puInfo_instLumi_;
	TBranch *puInfo_instLumi_branch;
	bool puInfo_instLumi_isLoaded;
	vector<vector<float> > vtxs_covMatrix_;
	TBranch *vtxs_covMatrix_branch;
	bool vtxs_covMatrix_isLoaded;
	int evt_bsType_;
	TBranch *evt_bsType_branch;
	bool evt_bsType_isLoaded;
	int evt_bunchCrossing_;
	TBranch *evt_bunchCrossing_branch;
	bool evt_bunchCrossing_isLoaded;
	int evt_experimentType_;
	TBranch *evt_experimentType_branch;
	bool evt_experimentType_isLoaded;
	int evt_isRealData_;
	TBranch *evt_isRealData_branch;
	bool evt_isRealData_isLoaded;
	int evt_orbitNumber_;
	TBranch *evt_orbitNumber_branch;
	bool evt_orbitNumber_isLoaded;
	int evt_storeNumber_;
	TBranch *evt_storeNumber_branch;
	bool evt_storeNumber_isLoaded;
	vector<int> els_category_;
	TBranch *els_category_branch;
	bool els_category_isLoaded;
	vector<int> els_charge_;
	TBranch *els_charge_branch;
	bool els_charge_isLoaded;
	vector<int> els_ckf_charge_;
	TBranch *els_ckf_charge_branch;
	bool els_ckf_charge_isLoaded;
	vector<int> els_ckf_laywithmeas_;
	TBranch *els_ckf_laywithmeas_branch;
	bool els_ckf_laywithmeas_isLoaded;
	vector<int> els_class_;
	TBranch *els_class_branch;
	bool els_class_isLoaded;
	vector<int> els_exp_innerlayers_;
	TBranch *els_exp_innerlayers_branch;
	bool els_exp_innerlayers_isLoaded;
	vector<int> els_exp_outerlayers_;
	TBranch *els_exp_outerlayers_branch;
	bool els_exp_outerlayers_isLoaded;
	vector<int> els_fiduciality_;
	TBranch *els_fiduciality_branch;
	bool els_fiduciality_isLoaded;
	vector<int> els_lostHits_;
	TBranch *els_lostHits_branch;
	bool els_lostHits_isLoaded;
	vector<int> els_lost_pixelhits_;
	TBranch *els_lost_pixelhits_branch;
	bool els_lost_pixelhits_isLoaded;
	vector<int> els_mc_patMatch_id_;
	TBranch *els_mc_patMatch_id_branch;
	bool els_mc_patMatch_id_isLoaded;
	vector<int> els_nSeed_;
	TBranch *els_nSeed_branch;
	bool els_nSeed_isLoaded;
	vector<int> els_nlayers_;
	TBranch *els_nlayers_branch;
	bool els_nlayers_isLoaded;
	vector<int> els_nlayers3D_;
	TBranch *els_nlayers3D_branch;
	bool els_nlayers3D_isLoaded;
	vector<int> els_nlayersLost_;
	TBranch *els_nlayersLost_branch;
	bool els_nlayersLost_isLoaded;
	vector<int> els_sccharge_;
	TBranch *els_sccharge_branch;
	bool els_sccharge_isLoaded;
	vector<int> els_trk_charge_;
	TBranch *els_trk_charge_branch;
	bool els_trk_charge_isLoaded;
	vector<int> els_type_;
	TBranch *els_type_branch;
	bool els_type_isLoaded;
	vector<int> els_validHits_;
	TBranch *els_validHits_branch;
	bool els_validHits_isLoaded;
	vector<int> els_valid_pixelhits_;
	TBranch *els_valid_pixelhits_branch;
	bool els_valid_pixelhits_isLoaded;
	vector<int> els_passLooseId_;
	TBranch *els_passLooseId_branch;
	bool els_passLooseId_isLoaded;
	vector<int> els_passMediumId_;
	TBranch *els_passMediumId_branch;
	bool els_passMediumId_isLoaded;
	vector<int> els_passTightId_;
	TBranch *els_passTightId_branch;
	bool els_passTightId_isLoaded;
	vector<int> els_passVetoId_;
	TBranch *els_passVetoId_branch;
	bool els_passVetoId_isLoaded;
	vector<int> photons_fiduciality_;
	TBranch *photons_fiduciality_branch;
	bool photons_fiduciality_isLoaded;
	vector<int> photons_photonID_loose_;
	TBranch *photons_photonID_loose_branch;
	bool photons_photonID_loose_isLoaded;
	vector<int> photons_photonID_tight_;
	TBranch *photons_photonID_tight_branch;
	bool photons_photonID_tight_isLoaded;
	vector<int> puInfo_bunchCrossing_;
	TBranch *puInfo_bunchCrossing_branch;
	bool puInfo_bunchCrossing_isLoaded;
	vector<int> puInfo_nPUvertices_;
	TBranch *puInfo_nPUvertices_branch;
	bool puInfo_nPUvertices_isLoaded;
	vector<int> vtxs_isFake_;
	TBranch *vtxs_isFake_branch;
	bool vtxs_isFake_isLoaded;
	vector<int> vtxs_isValid_;
	TBranch *vtxs_isValid_branch;
	bool vtxs_isValid_isLoaded;
	vector<vector<int> > els_PFCand_idx_;
	TBranch *els_PFCand_idx_branch;
	bool els_PFCand_idx_isLoaded;
	vector<vector<int> > hlt_trigObjs_id_;
	TBranch *hlt_trigObjs_id_branch;
	bool hlt_trigObjs_id_isLoaded;
	vector<vector<int> > photons_PFCand_idx_;
	TBranch *photons_PFCand_idx_branch;
	bool photons_PFCand_idx_isLoaded;
	unsigned int evt_nels_;
	TBranch *evt_nels_branch;
	bool evt_nels_isLoaded;
	unsigned int evt_detectorStatus_;
	TBranch *evt_detectorStatus_branch;
	bool evt_detectorStatus_isLoaded;
	unsigned int evt_lumiBlock_;
	TBranch *evt_lumiBlock_branch;
	bool evt_lumiBlock_isLoaded;
  unsigned long long evt_event_;
  TBranch *evt_event_branch;
  bool evt_event_isLoaded;
	unsigned int evt_run_;
	TBranch *evt_run_branch;
	bool evt_run_isLoaded;
	unsigned int evt_nphotons_;
	TBranch *evt_nphotons_branch;
	bool evt_nphotons_isLoaded;
	unsigned int evt_nvtxs_;
	TBranch *evt_nvtxs_branch;
	bool evt_nvtxs_isLoaded;
	vector<unsigned int> hlt_prescales_;
	TBranch *hlt_prescales_branch;
	bool hlt_prescales_isLoaded;
public: 
void Init(TTree *tree) {
	hlt_bits_branch = 0;
	if (tree->GetAlias("hlt_bits") != 0) {
		hlt_bits_branch = tree->GetBranch(tree->GetAlias("hlt_bits"));
		if (hlt_bits_branch) {hlt_bits_branch->SetAddress(&hlt_bits_);}
	}
	evt_bsp4_branch = 0;
	if (tree->GetAlias("evt_bsp4") != 0) {
		evt_bsp4_branch = tree->GetBranch(tree->GetAlias("evt_bsp4"));
		if (evt_bsp4_branch) {evt_bsp4_branch->SetAddress(&evt_bsp4_);}
	}
	els_mc_patMatch_p4_branch = 0;
	if (tree->GetAlias("els_mc_patMatch_p4") != 0) {
		els_mc_patMatch_p4_branch = tree->GetBranch(tree->GetAlias("els_mc_patMatch_p4"));
		if (els_mc_patMatch_p4_branch) {els_mc_patMatch_p4_branch->SetAddress(&els_mc_patMatch_p4_);}
	}
	els_p4_branch = 0;
	if (tree->GetAlias("els_p4") != 0) {
		els_p4_branch = tree->GetBranch(tree->GetAlias("els_p4"));
		if (els_p4_branch) {els_p4_branch->SetAddress(&els_p4_);}
	}
	els_p4In_branch = 0;
	if (tree->GetAlias("els_p4In") != 0) {
		els_p4In_branch = tree->GetBranch(tree->GetAlias("els_p4In"));
		if (els_p4In_branch) {els_p4In_branch->SetAddress(&els_p4In_);}
	}
	els_p4Out_branch = 0;
	if (tree->GetAlias("els_p4Out") != 0) {
		els_p4Out_branch = tree->GetBranch(tree->GetAlias("els_p4Out"));
		if (els_p4Out_branch) {els_p4Out_branch->SetAddress(&els_p4Out_);}
	}
	els_trk_p4_branch = 0;
	if (tree->GetAlias("els_trk_p4") != 0) {
		els_trk_p4_branch = tree->GetBranch(tree->GetAlias("els_trk_p4"));
		if (els_trk_p4_branch) {els_trk_p4_branch->SetAddress(&els_trk_p4_);}
	}
	els_trk_vertex_p4_branch = 0;
	if (tree->GetAlias("els_trk_vertex_p4") != 0) {
		els_trk_vertex_p4_branch = tree->GetBranch(tree->GetAlias("els_trk_vertex_p4"));
		if (els_trk_vertex_p4_branch) {els_trk_vertex_p4_branch->SetAddress(&els_trk_vertex_p4_);}
	}
	els_vertex_p4_branch = 0;
	if (tree->GetAlias("els_vertex_p4") != 0) {
		els_vertex_p4_branch = tree->GetBranch(tree->GetAlias("els_vertex_p4"));
		if (els_vertex_p4_branch) {els_vertex_p4_branch->SetAddress(&els_vertex_p4_);}
	}
	photons_p4_branch = 0;
	if (tree->GetAlias("photons_p4") != 0) {
		photons_p4_branch = tree->GetBranch(tree->GetAlias("photons_p4"));
		if (photons_p4_branch) {photons_p4_branch->SetAddress(&photons_p4_);}
	}
	vtxs_position_branch = 0;
	if (tree->GetAlias("vtxs_position") != 0) {
		vtxs_position_branch = tree->GetBranch(tree->GetAlias("vtxs_position"));
		if (vtxs_position_branch) {vtxs_position_branch->SetAddress(&vtxs_position_);}
	}
  tree->SetMakeClass(1);
	evt_CMS3tag_branch = 0;
	if (tree->GetAlias("evt_CMS3tag") != 0) {
		evt_CMS3tag_branch = tree->GetBranch(tree->GetAlias("evt_CMS3tag"));
		if (evt_CMS3tag_branch) {evt_CMS3tag_branch->SetAddress(&evt_CMS3tag_);}
	}
	evt_dataset_branch = 0;
	if (tree->GetAlias("evt_dataset") != 0) {
		evt_dataset_branch = tree->GetBranch(tree->GetAlias("evt_dataset"));
		if (evt_dataset_branch) {evt_dataset_branch->SetAddress(&evt_dataset_);}
	}
	hlt_trigNames_branch = 0;
	if (tree->GetAlias("hlt_trigNames") != 0) {
		hlt_trigNames_branch = tree->GetBranch(tree->GetAlias("hlt_trigNames"));
		if (hlt_trigNames_branch) {hlt_trigNames_branch->SetAddress(&hlt_trigNames_);}
	}
	els_conv_vtx_flag_branch = 0;
	if (tree->GetAlias("els_conv_vtx_flag") != 0) {
		els_conv_vtx_flag_branch = tree->GetBranch(tree->GetAlias("els_conv_vtx_flag"));
		if (els_conv_vtx_flag_branch) {els_conv_vtx_flag_branch->SetAddress(&els_conv_vtx_flag_);}
	}
	els_isGsfCtfScPixChargeConsistent_branch = 0;
	if (tree->GetAlias("els_isGsfCtfScPixChargeConsistent") != 0) {
		els_isGsfCtfScPixChargeConsistent_branch = tree->GetBranch(tree->GetAlias("els_isGsfCtfScPixChargeConsistent"));
		if (els_isGsfCtfScPixChargeConsistent_branch) {els_isGsfCtfScPixChargeConsistent_branch->SetAddress(&els_isGsfCtfScPixChargeConsistent_);}
	}
	els_passingMvaPreselection_branch = 0;
	if (tree->GetAlias("els_passingMvaPreselection") != 0) {
		els_passingMvaPreselection_branch = tree->GetBranch(tree->GetAlias("els_passingMvaPreselection"));
		if (els_passingMvaPreselection_branch) {els_passingMvaPreselection_branch->SetAddress(&els_passingMvaPreselection_);}
	}
	els_passingPflowPreselection_branch = 0;
	if (tree->GetAlias("els_passingPflowPreselection") != 0) {
		els_passingPflowPreselection_branch = tree->GetBranch(tree->GetAlias("els_passingPflowPreselection"));
		if (els_passingPflowPreselection_branch) {els_passingPflowPreselection_branch->SetAddress(&els_passingPflowPreselection_);}
	}
	photons_haspixelSeed_branch = 0;
	if (tree->GetAlias("photons_haspixelSeed") != 0) {
		photons_haspixelSeed_branch = tree->GetBranch(tree->GetAlias("photons_haspixelSeed"));
		if (photons_haspixelSeed_branch) {photons_haspixelSeed_branch->SetAddress(&photons_haspixelSeed_);}
	}
	evt_bs_Xwidth_branch = 0;
	if (tree->GetAlias("evt_bs_Xwidth") != 0) {
		evt_bs_Xwidth_branch = tree->GetBranch(tree->GetAlias("evt_bs_Xwidth"));
		if (evt_bs_Xwidth_branch) {evt_bs_Xwidth_branch->SetAddress(&evt_bs_Xwidth_);}
	}
	evt_bs_XwidthErr_branch = 0;
	if (tree->GetAlias("evt_bs_XwidthErr") != 0) {
		evt_bs_XwidthErr_branch = tree->GetBranch(tree->GetAlias("evt_bs_XwidthErr"));
		if (evt_bs_XwidthErr_branch) {evt_bs_XwidthErr_branch->SetAddress(&evt_bs_XwidthErr_);}
	}
	evt_bs_Ywidth_branch = 0;
	if (tree->GetAlias("evt_bs_Ywidth") != 0) {
		evt_bs_Ywidth_branch = tree->GetBranch(tree->GetAlias("evt_bs_Ywidth"));
		if (evt_bs_Ywidth_branch) {evt_bs_Ywidth_branch->SetAddress(&evt_bs_Ywidth_);}
	}
	evt_bs_YwidthErr_branch = 0;
	if (tree->GetAlias("evt_bs_YwidthErr") != 0) {
		evt_bs_YwidthErr_branch = tree->GetBranch(tree->GetAlias("evt_bs_YwidthErr"));
		if (evt_bs_YwidthErr_branch) {evt_bs_YwidthErr_branch->SetAddress(&evt_bs_YwidthErr_);}
	}
	evt_bs_dxdz_branch = 0;
	if (tree->GetAlias("evt_bs_dxdz") != 0) {
		evt_bs_dxdz_branch = tree->GetBranch(tree->GetAlias("evt_bs_dxdz"));
		if (evt_bs_dxdz_branch) {evt_bs_dxdz_branch->SetAddress(&evt_bs_dxdz_);}
	}
	evt_bs_dxdzErr_branch = 0;
	if (tree->GetAlias("evt_bs_dxdzErr") != 0) {
		evt_bs_dxdzErr_branch = tree->GetBranch(tree->GetAlias("evt_bs_dxdzErr"));
		if (evt_bs_dxdzErr_branch) {evt_bs_dxdzErr_branch->SetAddress(&evt_bs_dxdzErr_);}
	}
	evt_bs_dydz_branch = 0;
	if (tree->GetAlias("evt_bs_dydz") != 0) {
		evt_bs_dydz_branch = tree->GetBranch(tree->GetAlias("evt_bs_dydz"));
		if (evt_bs_dydz_branch) {evt_bs_dydz_branch->SetAddress(&evt_bs_dydz_);}
	}
	evt_bs_dydzErr_branch = 0;
	if (tree->GetAlias("evt_bs_dydzErr") != 0) {
		evt_bs_dydzErr_branch = tree->GetBranch(tree->GetAlias("evt_bs_dydzErr"));
		if (evt_bs_dydzErr_branch) {evt_bs_dydzErr_branch->SetAddress(&evt_bs_dydzErr_);}
	}
	evt_bs_sigmaZ_branch = 0;
	if (tree->GetAlias("evt_bs_sigmaZ") != 0) {
		evt_bs_sigmaZ_branch = tree->GetBranch(tree->GetAlias("evt_bs_sigmaZ"));
		if (evt_bs_sigmaZ_branch) {evt_bs_sigmaZ_branch->SetAddress(&evt_bs_sigmaZ_);}
	}
	evt_bs_sigmaZErr_branch = 0;
	if (tree->GetAlias("evt_bs_sigmaZErr") != 0) {
		evt_bs_sigmaZErr_branch = tree->GetBranch(tree->GetAlias("evt_bs_sigmaZErr"));
		if (evt_bs_sigmaZErr_branch) {evt_bs_sigmaZErr_branch->SetAddress(&evt_bs_sigmaZErr_);}
	}
	evt_bs_xErr_branch = 0;
	if (tree->GetAlias("evt_bs_xErr") != 0) {
		evt_bs_xErr_branch = tree->GetBranch(tree->GetAlias("evt_bs_xErr"));
		if (evt_bs_xErr_branch) {evt_bs_xErr_branch->SetAddress(&evt_bs_xErr_);}
	}
	evt_bs_yErr_branch = 0;
	if (tree->GetAlias("evt_bs_yErr") != 0) {
		evt_bs_yErr_branch = tree->GetBranch(tree->GetAlias("evt_bs_yErr"));
		if (evt_bs_yErr_branch) {evt_bs_yErr_branch->SetAddress(&evt_bs_yErr_);}
	}
	evt_bs_zErr_branch = 0;
	if (tree->GetAlias("evt_bs_zErr") != 0) {
		evt_bs_zErr_branch = tree->GetBranch(tree->GetAlias("evt_bs_zErr"));
		if (evt_bs_zErr_branch) {evt_bs_zErr_branch->SetAddress(&evt_bs_zErr_);}
	}
	evt_bField_branch = 0;
	if (tree->GetAlias("evt_bField") != 0) {
		evt_bField_branch = tree->GetBranch(tree->GetAlias("evt_bField"));
		if (evt_bField_branch) {evt_bField_branch->SetAddress(&evt_bField_);}
	}
	evt_fixgrid_all_rho_branch = 0;
	if (tree->GetAlias("evt_fixgrid_all_rho") != 0) {
		evt_fixgrid_all_rho_branch = tree->GetBranch(tree->GetAlias("evt_fixgrid_all_rho"));
		if (evt_fixgrid_all_rho_branch) {evt_fixgrid_all_rho_branch->SetAddress(&evt_fixgrid_all_rho_);}
	}
	evt_fixgridfastjet_allcalo_rho_branch = 0;
	if (tree->GetAlias("evt_fixgridfastjet_allcalo_rho") != 0) {
		evt_fixgridfastjet_allcalo_rho_branch = tree->GetBranch(tree->GetAlias("evt_fixgridfastjet_allcalo_rho"));
		if (evt_fixgridfastjet_allcalo_rho_branch) {evt_fixgridfastjet_allcalo_rho_branch->SetAddress(&evt_fixgridfastjet_allcalo_rho_);}
	}
	evt_fixgridfastjet_all_rho_branch = 0;
	if (tree->GetAlias("evt_fixgridfastjet_all_rho") != 0) {
		evt_fixgridfastjet_all_rho_branch = tree->GetBranch(tree->GetAlias("evt_fixgridfastjet_all_rho"));
		if (evt_fixgridfastjet_all_rho_branch) {evt_fixgridfastjet_all_rho_branch->SetAddress(&evt_fixgridfastjet_all_rho_);}
	}
	evt_fixgridfastjet_centralcalo_rho_branch = 0;
	if (tree->GetAlias("evt_fixgridfastjet_centralcalo_rho") != 0) {
		evt_fixgridfastjet_centralcalo_rho_branch = tree->GetBranch(tree->GetAlias("evt_fixgridfastjet_centralcalo_rho"));
		if (evt_fixgridfastjet_centralcalo_rho_branch) {evt_fixgridfastjet_centralcalo_rho_branch->SetAddress(&evt_fixgridfastjet_centralcalo_rho_);}
	}
	evt_fixgridfastjet_centralchargedpileup_rho_branch = 0;
	if (tree->GetAlias("evt_fixgridfastjet_centralchargedpileup_rho") != 0) {
		evt_fixgridfastjet_centralchargedpileup_rho_branch = tree->GetBranch(tree->GetAlias("evt_fixgridfastjet_centralchargedpileup_rho"));
		if (evt_fixgridfastjet_centralchargedpileup_rho_branch) {evt_fixgridfastjet_centralchargedpileup_rho_branch->SetAddress(&evt_fixgridfastjet_centralchargedpileup_rho_);}
	}
	evt_fixgridfastjet_centralneutral_rho_branch = 0;
	if (tree->GetAlias("evt_fixgridfastjet_centralneutral_rho") != 0) {
		evt_fixgridfastjet_centralneutral_rho_branch = tree->GetBranch(tree->GetAlias("evt_fixgridfastjet_centralneutral_rho"));
		if (evt_fixgridfastjet_centralneutral_rho_branch) {evt_fixgridfastjet_centralneutral_rho_branch->SetAddress(&evt_fixgridfastjet_centralneutral_rho_);}
	}
	hlt_trigObjs_p4_branch = 0;
	if (tree->GetAlias("hlt_trigObjs_p4") != 0) {
		hlt_trigObjs_p4_branch = tree->GetBranch(tree->GetAlias("hlt_trigObjs_p4"));
		if (hlt_trigObjs_p4_branch) {hlt_trigObjs_p4_branch->SetAddress(&hlt_trigObjs_p4_);}
	}
	evt_bs_covMatrix_branch = 0;
	if (tree->GetAlias("evt_bs_covMatrix") != 0) {
		evt_bs_covMatrix_branch = tree->GetBranch(tree->GetAlias("evt_bs_covMatrix"));
		if (evt_bs_covMatrix_branch) {evt_bs_covMatrix_branch->SetAddress(&evt_bs_covMatrix_);}
	}
	els_bs2d_branch = 0;
	if (tree->GetAlias("els_bs2d") != 0) {
		els_bs2d_branch = tree->GetBranch(tree->GetAlias("els_bs2d"));
		if (els_bs2d_branch) {els_bs2d_branch->SetAddress(&els_bs2d_);}
	}
	els_bs2derr_branch = 0;
	if (tree->GetAlias("els_bs2derr") != 0) {
		els_bs2derr_branch = tree->GetBranch(tree->GetAlias("els_bs2derr"));
		if (els_bs2derr_branch) {els_bs2derr_branch->SetAddress(&els_bs2derr_);}
	}
	els_bs3d_branch = 0;
	if (tree->GetAlias("els_bs3d") != 0) {
		els_bs3d_branch = tree->GetBranch(tree->GetAlias("els_bs3d"));
		if (els_bs3d_branch) {els_bs3d_branch->SetAddress(&els_bs3d_);}
	}
	els_bs3derr_branch = 0;
	if (tree->GetAlias("els_bs3derr") != 0) {
		els_bs3derr_branch = tree->GetBranch(tree->GetAlias("els_bs3derr"));
		if (els_bs3derr_branch) {els_bs3derr_branch->SetAddress(&els_bs3derr_);}
	}
	els_chi2_branch = 0;
	if (tree->GetAlias("els_chi2") != 0) {
		els_chi2_branch = tree->GetBranch(tree->GetAlias("els_chi2"));
		if (els_chi2_branch) {els_chi2_branch->SetAddress(&els_chi2_);}
	}
	els_ckf_chi2_branch = 0;
	if (tree->GetAlias("els_ckf_chi2") != 0) {
		els_ckf_chi2_branch = tree->GetBranch(tree->GetAlias("els_ckf_chi2"));
		if (els_ckf_chi2_branch) {els_ckf_chi2_branch->SetAddress(&els_ckf_chi2_);}
	}
	els_ckf_ndof_branch = 0;
	if (tree->GetAlias("els_ckf_ndof") != 0) {
		els_ckf_ndof_branch = tree->GetBranch(tree->GetAlias("els_ckf_ndof"));
		if (els_ckf_ndof_branch) {els_ckf_ndof_branch->SetAddress(&els_ckf_ndof_);}
	}
	els_d0_branch = 0;
	if (tree->GetAlias("els_d0") != 0) {
		els_d0_branch = tree->GetBranch(tree->GetAlias("els_d0"));
		if (els_d0_branch) {els_d0_branch->SetAddress(&els_d0_);}
	}
	els_d0Err_branch = 0;
	if (tree->GetAlias("els_d0Err") != 0) {
		els_d0Err_branch = tree->GetBranch(tree->GetAlias("els_d0Err"));
		if (els_d0Err_branch) {els_d0Err_branch->SetAddress(&els_d0Err_);}
	}
	els_d0corr_branch = 0;
	if (tree->GetAlias("els_d0corr") != 0) {
		els_d0corr_branch = tree->GetBranch(tree->GetAlias("els_d0corr"));
		if (els_d0corr_branch) {els_d0corr_branch->SetAddress(&els_d0corr_);}
	}
	els_d0corrPhi_branch = 0;
	if (tree->GetAlias("els_d0corrPhi") != 0) {
		els_d0corrPhi_branch = tree->GetBranch(tree->GetAlias("els_d0corrPhi"));
		if (els_d0corrPhi_branch) {els_d0corrPhi_branch->SetAddress(&els_d0corrPhi_);}
	}
	els_d0phiCov_branch = 0;
	if (tree->GetAlias("els_d0phiCov") != 0) {
		els_d0phiCov_branch = tree->GetBranch(tree->GetAlias("els_d0phiCov"));
		if (els_d0phiCov_branch) {els_d0phiCov_branch->SetAddress(&els_d0phiCov_);}
	}
	els_dEtaIn_branch = 0;
	if (tree->GetAlias("els_dEtaIn") != 0) {
		els_dEtaIn_branch = tree->GetBranch(tree->GetAlias("els_dEtaIn"));
		if (els_dEtaIn_branch) {els_dEtaIn_branch->SetAddress(&els_dEtaIn_);}
	}
	els_dEtaOut_branch = 0;
	if (tree->GetAlias("els_dEtaOut") != 0) {
		els_dEtaOut_branch = tree->GetBranch(tree->GetAlias("els_dEtaOut"));
		if (els_dEtaOut_branch) {els_dEtaOut_branch->SetAddress(&els_dEtaOut_);}
	}
	els_dPhiIn_branch = 0;
	if (tree->GetAlias("els_dPhiIn") != 0) {
		els_dPhiIn_branch = tree->GetBranch(tree->GetAlias("els_dPhiIn"));
		if (els_dPhiIn_branch) {els_dPhiIn_branch->SetAddress(&els_dPhiIn_);}
	}
	els_dPhiInPhiOut_branch = 0;
	if (tree->GetAlias("els_dPhiInPhiOut") != 0) {
		els_dPhiInPhiOut_branch = tree->GetBranch(tree->GetAlias("els_dPhiInPhiOut"));
		if (els_dPhiInPhiOut_branch) {els_dPhiInPhiOut_branch->SetAddress(&els_dPhiInPhiOut_);}
	}
	els_dPhiOut_branch = 0;
	if (tree->GetAlias("els_dPhiOut") != 0) {
		els_dPhiOut_branch = tree->GetBranch(tree->GetAlias("els_dPhiOut"));
		if (els_dPhiOut_branch) {els_dPhiOut_branch->SetAddress(&els_dPhiOut_);}
	}
	els_deltaEtaEleClusterTrackAtCalo_branch = 0;
	if (tree->GetAlias("els_deltaEtaEleClusterTrackAtCalo") != 0) {
		els_deltaEtaEleClusterTrackAtCalo_branch = tree->GetBranch(tree->GetAlias("els_deltaEtaEleClusterTrackAtCalo"));
		if (els_deltaEtaEleClusterTrackAtCalo_branch) {els_deltaEtaEleClusterTrackAtCalo_branch->SetAddress(&els_deltaEtaEleClusterTrackAtCalo_);}
	}
	els_deltaPhiEleClusterTrackAtCalo_branch = 0;
	if (tree->GetAlias("els_deltaPhiEleClusterTrackAtCalo") != 0) {
		els_deltaPhiEleClusterTrackAtCalo_branch = tree->GetBranch(tree->GetAlias("els_deltaPhiEleClusterTrackAtCalo"));
		if (els_deltaPhiEleClusterTrackAtCalo_branch) {els_deltaPhiEleClusterTrackAtCalo_branch->SetAddress(&els_deltaPhiEleClusterTrackAtCalo_);}
	}
	els_dxyPV_branch = 0;
	if (tree->GetAlias("els_dxyPV") != 0) {
		els_dxyPV_branch = tree->GetBranch(tree->GetAlias("els_dxyPV"));
		if (els_dxyPV_branch) {els_dxyPV_branch->SetAddress(&els_dxyPV_);}
	}
	els_dzPV_branch = 0;
	if (tree->GetAlias("els_dzPV") != 0) {
		els_dzPV_branch = tree->GetBranch(tree->GetAlias("els_dzPV"));
		if (els_dzPV_branch) {els_dzPV_branch->SetAddress(&els_dzPV_);}
	}
	els_e1x5_branch = 0;
	if (tree->GetAlias("els_e1x5") != 0) {
		els_e1x5_branch = tree->GetBranch(tree->GetAlias("els_e1x5"));
		if (els_e1x5_branch) {els_e1x5_branch->SetAddress(&els_e1x5_);}
	}
	els_e1x5_full5x5_branch = 0;
	if (tree->GetAlias("els_e1x5_full5x5") != 0) {
		els_e1x5_full5x5_branch = tree->GetBranch(tree->GetAlias("els_e1x5_full5x5"));
		if (els_e1x5_full5x5_branch) {els_e1x5_full5x5_branch->SetAddress(&els_e1x5_full5x5_);}
	}
	els_e2x5Max_branch = 0;
	if (tree->GetAlias("els_e2x5Max") != 0) {
		els_e2x5Max_branch = tree->GetBranch(tree->GetAlias("els_e2x5Max"));
		if (els_e2x5Max_branch) {els_e2x5Max_branch->SetAddress(&els_e2x5Max_);}
	}
	els_e2x5Max_full5x5_branch = 0;
	if (tree->GetAlias("els_e2x5Max_full5x5") != 0) {
		els_e2x5Max_full5x5_branch = tree->GetBranch(tree->GetAlias("els_e2x5Max_full5x5"));
		if (els_e2x5Max_full5x5_branch) {els_e2x5Max_full5x5_branch->SetAddress(&els_e2x5Max_full5x5_);}
	}
	els_e5x5_branch = 0;
	if (tree->GetAlias("els_e5x5") != 0) {
		els_e5x5_branch = tree->GetBranch(tree->GetAlias("els_e5x5"));
		if (els_e5x5_branch) {els_e5x5_branch->SetAddress(&els_e5x5_);}
	}
	els_e5x5_full5x5_branch = 0;
	if (tree->GetAlias("els_e5x5_full5x5") != 0) {
		els_e5x5_full5x5_branch = tree->GetBranch(tree->GetAlias("els_e5x5_full5x5"));
		if (els_e5x5_full5x5_branch) {els_e5x5_full5x5_branch->SetAddress(&els_e5x5_full5x5_);}
	}
	els_eOverPIn_branch = 0;
	if (tree->GetAlias("els_eOverPIn") != 0) {
		els_eOverPIn_branch = tree->GetBranch(tree->GetAlias("els_eOverPIn"));
		if (els_eOverPIn_branch) {els_eOverPIn_branch->SetAddress(&els_eOverPIn_);}
	}
	els_eOverPOut_branch = 0;
	if (tree->GetAlias("els_eOverPOut") != 0) {
		els_eOverPOut_branch = tree->GetBranch(tree->GetAlias("els_eOverPOut"));
		if (els_eOverPOut_branch) {els_eOverPOut_branch->SetAddress(&els_eOverPOut_);}
	}
	els_eSC_branch = 0;
	if (tree->GetAlias("els_eSC") != 0) {
		els_eSC_branch = tree->GetBranch(tree->GetAlias("els_eSC"));
		if (els_eSC_branch) {els_eSC_branch->SetAddress(&els_eSC_);}
	}
	els_eSCPresh_branch = 0;
	if (tree->GetAlias("els_eSCPresh") != 0) {
		els_eSCPresh_branch = tree->GetBranch(tree->GetAlias("els_eSCPresh"));
		if (els_eSCPresh_branch) {els_eSCPresh_branch->SetAddress(&els_eSCPresh_);}
	}
	els_eSCRaw_branch = 0;
	if (tree->GetAlias("els_eSCRaw") != 0) {
		els_eSCRaw_branch = tree->GetBranch(tree->GetAlias("els_eSCRaw"));
		if (els_eSCRaw_branch) {els_eSCRaw_branch->SetAddress(&els_eSCRaw_);}
	}
	els_eSeed_branch = 0;
	if (tree->GetAlias("els_eSeed") != 0) {
		els_eSeed_branch = tree->GetBranch(tree->GetAlias("els_eSeed"));
		if (els_eSeed_branch) {els_eSeed_branch->SetAddress(&els_eSeed_);}
	}
	els_eSeedOverPIn_branch = 0;
	if (tree->GetAlias("els_eSeedOverPIn") != 0) {
		els_eSeedOverPIn_branch = tree->GetBranch(tree->GetAlias("els_eSeedOverPIn"));
		if (els_eSeedOverPIn_branch) {els_eSeedOverPIn_branch->SetAddress(&els_eSeedOverPIn_);}
	}
	els_eSeedOverPOut_branch = 0;
	if (tree->GetAlias("els_eSeedOverPOut") != 0) {
		els_eSeedOverPOut_branch = tree->GetBranch(tree->GetAlias("els_eSeedOverPOut"));
		if (els_eSeedOverPOut_branch) {els_eSeedOverPOut_branch->SetAddress(&els_eSeedOverPOut_);}
	}
	els_ecalEnergy_branch = 0;
	if (tree->GetAlias("els_ecalEnergy") != 0) {
		els_ecalEnergy_branch = tree->GetBranch(tree->GetAlias("els_ecalEnergy"));
		if (els_ecalEnergy_branch) {els_ecalEnergy_branch->SetAddress(&els_ecalEnergy_);}
	}
	els_ecalEnergyError_branch = 0;
	if (tree->GetAlias("els_ecalEnergyError") != 0) {
		els_ecalEnergyError_branch = tree->GetBranch(tree->GetAlias("els_ecalEnergyError"));
		if (els_ecalEnergyError_branch) {els_ecalEnergyError_branch->SetAddress(&els_ecalEnergyError_);}
	}
	els_ecalIso_branch = 0;
	if (tree->GetAlias("els_ecalIso") != 0) {
		els_ecalIso_branch = tree->GetBranch(tree->GetAlias("els_ecalIso"));
		if (els_ecalIso_branch) {els_ecalIso_branch->SetAddress(&els_ecalIso_);}
	}
	els_ecalIso04_branch = 0;
	if (tree->GetAlias("els_ecalIso04") != 0) {
		els_ecalIso04_branch = tree->GetBranch(tree->GetAlias("els_ecalIso04"));
		if (els_ecalIso04_branch) {els_ecalIso04_branch->SetAddress(&els_ecalIso04_);}
	}
	els_ecalPFClusterIso_branch = 0;
	if (tree->GetAlias("els_ecalPFClusterIso") != 0) {
		els_ecalPFClusterIso_branch = tree->GetBranch(tree->GetAlias("els_ecalPFClusterIso"));
		if (els_ecalPFClusterIso_branch) {els_ecalPFClusterIso_branch->SetAddress(&els_ecalPFClusterIso_);}
	}
	els_etaErr_branch = 0;
	if (tree->GetAlias("els_etaErr") != 0) {
		els_etaErr_branch = tree->GetBranch(tree->GetAlias("els_etaErr"));
		if (els_etaErr_branch) {els_etaErr_branch->SetAddress(&els_etaErr_);}
	}
	els_etaSC_branch = 0;
	if (tree->GetAlias("els_etaSC") != 0) {
		els_etaSC_branch = tree->GetBranch(tree->GetAlias("els_etaSC"));
		if (els_etaSC_branch) {els_etaSC_branch->SetAddress(&els_etaSC_);}
	}
	els_etaSCwidth_branch = 0;
	if (tree->GetAlias("els_etaSCwidth") != 0) {
		els_etaSCwidth_branch = tree->GetBranch(tree->GetAlias("els_etaSCwidth"));
		if (els_etaSCwidth_branch) {els_etaSCwidth_branch->SetAddress(&els_etaSCwidth_);}
	}
	els_fbrem_branch = 0;
	if (tree->GetAlias("els_fbrem") != 0) {
		els_fbrem_branch = tree->GetBranch(tree->GetAlias("els_fbrem"));
		if (els_fbrem_branch) {els_fbrem_branch->SetAddress(&els_fbrem_);}
	}
	els_hOverE_branch = 0;
	if (tree->GetAlias("els_hOverE") != 0) {
		els_hOverE_branch = tree->GetBranch(tree->GetAlias("els_hOverE"));
		if (els_hOverE_branch) {els_hOverE_branch->SetAddress(&els_hOverE_);}
	}
	els_hOverEBC_branch = 0;
	if (tree->GetAlias("els_hOverEBC") != 0) {
		els_hOverEBC_branch = tree->GetBranch(tree->GetAlias("els_hOverEBC"));
		if (els_hOverEBC_branch) {els_hOverEBC_branch->SetAddress(&els_hOverEBC_);}
	}
	els_hcalDepth1OverEcal_branch = 0;
	if (tree->GetAlias("els_hcalDepth1OverEcal") != 0) {
		els_hcalDepth1OverEcal_branch = tree->GetBranch(tree->GetAlias("els_hcalDepth1OverEcal"));
		if (els_hcalDepth1OverEcal_branch) {els_hcalDepth1OverEcal_branch->SetAddress(&els_hcalDepth1OverEcal_);}
	}
	els_hcalDepth1TowerSumEt_branch = 0;
	if (tree->GetAlias("els_hcalDepth1TowerSumEt") != 0) {
		els_hcalDepth1TowerSumEt_branch = tree->GetBranch(tree->GetAlias("els_hcalDepth1TowerSumEt"));
		if (els_hcalDepth1TowerSumEt_branch) {els_hcalDepth1TowerSumEt_branch->SetAddress(&els_hcalDepth1TowerSumEt_);}
	}
	els_hcalDepth1TowerSumEt04_branch = 0;
	if (tree->GetAlias("els_hcalDepth1TowerSumEt04") != 0) {
		els_hcalDepth1TowerSumEt04_branch = tree->GetBranch(tree->GetAlias("els_hcalDepth1TowerSumEt04"));
		if (els_hcalDepth1TowerSumEt04_branch) {els_hcalDepth1TowerSumEt04_branch->SetAddress(&els_hcalDepth1TowerSumEt04_);}
	}
	els_hcalDepth2OverEcal_branch = 0;
	if (tree->GetAlias("els_hcalDepth2OverEcal") != 0) {
		els_hcalDepth2OverEcal_branch = tree->GetBranch(tree->GetAlias("els_hcalDepth2OverEcal"));
		if (els_hcalDepth2OverEcal_branch) {els_hcalDepth2OverEcal_branch->SetAddress(&els_hcalDepth2OverEcal_);}
	}
	els_hcalDepth2TowerSumEt_branch = 0;
	if (tree->GetAlias("els_hcalDepth2TowerSumEt") != 0) {
		els_hcalDepth2TowerSumEt_branch = tree->GetBranch(tree->GetAlias("els_hcalDepth2TowerSumEt"));
		if (els_hcalDepth2TowerSumEt_branch) {els_hcalDepth2TowerSumEt_branch->SetAddress(&els_hcalDepth2TowerSumEt_);}
	}
	els_hcalDepth2TowerSumEt04_branch = 0;
	if (tree->GetAlias("els_hcalDepth2TowerSumEt04") != 0) {
		els_hcalDepth2TowerSumEt04_branch = tree->GetBranch(tree->GetAlias("els_hcalDepth2TowerSumEt04"));
		if (els_hcalDepth2TowerSumEt04_branch) {els_hcalDepth2TowerSumEt04_branch->SetAddress(&els_hcalDepth2TowerSumEt04_);}
	}
	els_hcalIso_branch = 0;
	if (tree->GetAlias("els_hcalIso") != 0) {
		els_hcalIso_branch = tree->GetBranch(tree->GetAlias("els_hcalIso"));
		if (els_hcalIso_branch) {els_hcalIso_branch->SetAddress(&els_hcalIso_);}
	}
	els_hcalIso04_branch = 0;
	if (tree->GetAlias("els_hcalIso04") != 0) {
		els_hcalIso04_branch = tree->GetBranch(tree->GetAlias("els_hcalIso04"));
		if (els_hcalIso04_branch) {els_hcalIso04_branch->SetAddress(&els_hcalIso04_);}
	}
	els_hcalPFClusterIso_branch = 0;
	if (tree->GetAlias("els_hcalPFClusterIso") != 0) {
		els_hcalPFClusterIso_branch = tree->GetBranch(tree->GetAlias("els_hcalPFClusterIso"));
		if (els_hcalPFClusterIso_branch) {els_hcalPFClusterIso_branch->SetAddress(&els_hcalPFClusterIso_);}
	}
	els_ip2d_branch = 0;
	if (tree->GetAlias("els_ip2d") != 0) {
		els_ip2d_branch = tree->GetBranch(tree->GetAlias("els_ip2d"));
		if (els_ip2d_branch) {els_ip2d_branch->SetAddress(&els_ip2d_);}
	}
	els_ip2derr_branch = 0;
	if (tree->GetAlias("els_ip2derr") != 0) {
		els_ip2derr_branch = tree->GetBranch(tree->GetAlias("els_ip2derr"));
		if (els_ip2derr_branch) {els_ip2derr_branch->SetAddress(&els_ip2derr_);}
	}
	els_ip3d_branch = 0;
	if (tree->GetAlias("els_ip3d") != 0) {
		els_ip3d_branch = tree->GetBranch(tree->GetAlias("els_ip3d"));
		if (els_ip3d_branch) {els_ip3d_branch->SetAddress(&els_ip3d_);}
	}
	els_ip3derr_branch = 0;
	if (tree->GetAlias("els_ip3derr") != 0) {
		els_ip3derr_branch = tree->GetBranch(tree->GetAlias("els_ip3derr"));
		if (els_ip3derr_branch) {els_ip3derr_branch->SetAddress(&els_ip3derr_);}
	}
	els_mass_branch = 0;
	if (tree->GetAlias("els_mass") != 0) {
		els_mass_branch = tree->GetBranch(tree->GetAlias("els_mass"));
		if (els_mass_branch) {els_mass_branch->SetAddress(&els_mass_);}
	}
	els_mc_patMatch_dr_branch = 0;
	if (tree->GetAlias("els_mc_patMatch_dr") != 0) {
		els_mc_patMatch_dr_branch = tree->GetBranch(tree->GetAlias("els_mc_patMatch_dr"));
		if (els_mc_patMatch_dr_branch) {els_mc_patMatch_dr_branch->SetAddress(&els_mc_patMatch_dr_);}
	}
	els_miniIso_ch_branch = 0;
	if (tree->GetAlias("els_miniIso_ch") != 0) {
		els_miniIso_ch_branch = tree->GetBranch(tree->GetAlias("els_miniIso_ch"));
		if (els_miniIso_ch_branch) {els_miniIso_ch_branch->SetAddress(&els_miniIso_ch_);}
	}
	els_miniIso_db_branch = 0;
	if (tree->GetAlias("els_miniIso_db") != 0) {
		els_miniIso_db_branch = tree->GetBranch(tree->GetAlias("els_miniIso_db"));
		if (els_miniIso_db_branch) {els_miniIso_db_branch->SetAddress(&els_miniIso_db_);}
	}
	els_miniIso_em_branch = 0;
	if (tree->GetAlias("els_miniIso_em") != 0) {
		els_miniIso_em_branch = tree->GetBranch(tree->GetAlias("els_miniIso_em"));
		if (els_miniIso_em_branch) {els_miniIso_em_branch->SetAddress(&els_miniIso_em_);}
	}
	els_miniIso_nh_branch = 0;
	if (tree->GetAlias("els_miniIso_nh") != 0) {
		els_miniIso_nh_branch = tree->GetBranch(tree->GetAlias("els_miniIso_nh"));
		if (els_miniIso_nh_branch) {els_miniIso_nh_branch->SetAddress(&els_miniIso_nh_);}
	}
	els_miniIso_uncor_branch = 0;
	if (tree->GetAlias("els_miniIso_uncor") != 0) {
		els_miniIso_uncor_branch = tree->GetBranch(tree->GetAlias("els_miniIso_uncor"));
		if (els_miniIso_uncor_branch) {els_miniIso_uncor_branch->SetAddress(&els_miniIso_uncor_);}
	}
	els_mva_branch = 0;
	if (tree->GetAlias("els_mva") != 0) {
		els_mva_branch = tree->GetBranch(tree->GetAlias("els_mva"));
		if (els_mva_branch) {els_mva_branch->SetAddress(&els_mva_);}
	}
	els_ndof_branch = 0;
	if (tree->GetAlias("els_ndof") != 0) {
		els_ndof_branch = tree->GetBranch(tree->GetAlias("els_ndof"));
		if (els_ndof_branch) {els_ndof_branch->SetAddress(&els_ndof_);}
	}
	els_pfChargedHadronIso_branch = 0;
	if (tree->GetAlias("els_pfChargedHadronIso") != 0) {
		els_pfChargedHadronIso_branch = tree->GetBranch(tree->GetAlias("els_pfChargedHadronIso"));
		if (els_pfChargedHadronIso_branch) {els_pfChargedHadronIso_branch->SetAddress(&els_pfChargedHadronIso_);}
	}
	els_pfNeutralHadronIso_branch = 0;
	if (tree->GetAlias("els_pfNeutralHadronIso") != 0) {
		els_pfNeutralHadronIso_branch = tree->GetBranch(tree->GetAlias("els_pfNeutralHadronIso"));
		if (els_pfNeutralHadronIso_branch) {els_pfNeutralHadronIso_branch->SetAddress(&els_pfNeutralHadronIso_);}
	}
	els_pfPUIso_branch = 0;
	if (tree->GetAlias("els_pfPUIso") != 0) {
		els_pfPUIso_branch = tree->GetBranch(tree->GetAlias("els_pfPUIso"));
		if (els_pfPUIso_branch) {els_pfPUIso_branch->SetAddress(&els_pfPUIso_);}
	}
	els_pfPhotonIso_branch = 0;
	if (tree->GetAlias("els_pfPhotonIso") != 0) {
		els_pfPhotonIso_branch = tree->GetBranch(tree->GetAlias("els_pfPhotonIso"));
		if (els_pfPhotonIso_branch) {els_pfPhotonIso_branch->SetAddress(&els_pfPhotonIso_);}
	}
	els_phiErr_branch = 0;
	if (tree->GetAlias("els_phiErr") != 0) {
		els_phiErr_branch = tree->GetBranch(tree->GetAlias("els_phiErr"));
		if (els_phiErr_branch) {els_phiErr_branch->SetAddress(&els_phiErr_);}
	}
	els_phiSC_branch = 0;
	if (tree->GetAlias("els_phiSC") != 0) {
		els_phiSC_branch = tree->GetBranch(tree->GetAlias("els_phiSC"));
		if (els_phiSC_branch) {els_phiSC_branch->SetAddress(&els_phiSC_);}
	}
	els_phiSCwidth_branch = 0;
	if (tree->GetAlias("els_phiSCwidth") != 0) {
		els_phiSCwidth_branch = tree->GetBranch(tree->GetAlias("els_phiSCwidth"));
		if (els_phiSCwidth_branch) {els_phiSCwidth_branch->SetAddress(&els_phiSCwidth_);}
	}
	els_ptErr_branch = 0;
	if (tree->GetAlias("els_ptErr") != 0) {
		els_ptErr_branch = tree->GetBranch(tree->GetAlias("els_ptErr"));
		if (els_ptErr_branch) {els_ptErr_branch->SetAddress(&els_ptErr_);}
	}
	els_ptErrGsf_branch = 0;
	if (tree->GetAlias("els_ptErrGsf") != 0) {
		els_ptErrGsf_branch = tree->GetBranch(tree->GetAlias("els_ptErrGsf"));
		if (els_ptErrGsf_branch) {els_ptErrGsf_branch->SetAddress(&els_ptErrGsf_);}
	}
	els_r9_branch = 0;
	if (tree->GetAlias("els_r9") != 0) {
		els_r9_branch = tree->GetBranch(tree->GetAlias("els_r9"));
		if (els_r9_branch) {els_r9_branch->SetAddress(&els_r9_);}
	}
	els_r9_full5x5_branch = 0;
	if (tree->GetAlias("els_r9_full5x5") != 0) {
		els_r9_full5x5_branch = tree->GetBranch(tree->GetAlias("els_r9_full5x5"));
		if (els_r9_full5x5_branch) {els_r9_full5x5_branch->SetAddress(&els_r9_full5x5_);}
	}
	els_sigmaEtaEta_branch = 0;
	if (tree->GetAlias("els_sigmaEtaEta") != 0) {
		els_sigmaEtaEta_branch = tree->GetBranch(tree->GetAlias("els_sigmaEtaEta"));
		if (els_sigmaEtaEta_branch) {els_sigmaEtaEta_branch->SetAddress(&els_sigmaEtaEta_);}
	}
	els_sigmaEtaEta_full5x5_branch = 0;
	if (tree->GetAlias("els_sigmaEtaEta_full5x5") != 0) {
		els_sigmaEtaEta_full5x5_branch = tree->GetBranch(tree->GetAlias("els_sigmaEtaEta_full5x5"));
		if (els_sigmaEtaEta_full5x5_branch) {els_sigmaEtaEta_full5x5_branch->SetAddress(&els_sigmaEtaEta_full5x5_);}
	}
	els_sigmaIEtaIEta_branch = 0;
	if (tree->GetAlias("els_sigmaIEtaIEta") != 0) {
		els_sigmaIEtaIEta_branch = tree->GetBranch(tree->GetAlias("els_sigmaIEtaIEta"));
		if (els_sigmaIEtaIEta_branch) {els_sigmaIEtaIEta_branch->SetAddress(&els_sigmaIEtaIEta_);}
	}
	els_sigmaIEtaIEta_full5x5_branch = 0;
	if (tree->GetAlias("els_sigmaIEtaIEta_full5x5") != 0) {
		els_sigmaIEtaIEta_full5x5_branch = tree->GetBranch(tree->GetAlias("els_sigmaIEtaIEta_full5x5"));
		if (els_sigmaIEtaIEta_full5x5_branch) {els_sigmaIEtaIEta_full5x5_branch->SetAddress(&els_sigmaIEtaIEta_full5x5_);}
	}
	els_sigmaIPhiIPhi_branch = 0;
	if (tree->GetAlias("els_sigmaIPhiIPhi") != 0) {
		els_sigmaIPhiIPhi_branch = tree->GetBranch(tree->GetAlias("els_sigmaIPhiIPhi"));
		if (els_sigmaIPhiIPhi_branch) {els_sigmaIPhiIPhi_branch->SetAddress(&els_sigmaIPhiIPhi_);}
	}
	els_sigmaIPhiIPhi_full5x5_branch = 0;
	if (tree->GetAlias("els_sigmaIPhiIPhi_full5x5") != 0) {
		els_sigmaIPhiIPhi_full5x5_branch = tree->GetBranch(tree->GetAlias("els_sigmaIPhiIPhi_full5x5"));
		if (els_sigmaIPhiIPhi_full5x5_branch) {els_sigmaIPhiIPhi_full5x5_branch->SetAddress(&els_sigmaIPhiIPhi_full5x5_);}
	}
	els_sigmaIphiIphi_branch = 0;
	if (tree->GetAlias("els_sigmaIphiIphi") != 0) {
		els_sigmaIphiIphi_branch = tree->GetBranch(tree->GetAlias("els_sigmaIphiIphi"));
		if (els_sigmaIphiIphi_branch) {els_sigmaIphiIphi_branch->SetAddress(&els_sigmaIphiIphi_);}
	}
	els_tkIso_branch = 0;
	if (tree->GetAlias("els_tkIso") != 0) {
		els_tkIso_branch = tree->GetBranch(tree->GetAlias("els_tkIso"));
		if (els_tkIso_branch) {els_tkIso_branch->SetAddress(&els_tkIso_);}
	}
	els_tkIso04_branch = 0;
	if (tree->GetAlias("els_tkIso04") != 0) {
		els_tkIso04_branch = tree->GetBranch(tree->GetAlias("els_tkIso04"));
		if (els_tkIso04_branch) {els_tkIso04_branch->SetAddress(&els_tkIso04_);}
	}
	els_trackMomentumError_branch = 0;
	if (tree->GetAlias("els_trackMomentumError") != 0) {
		els_trackMomentumError_branch = tree->GetBranch(tree->GetAlias("els_trackMomentumError"));
		if (els_trackMomentumError_branch) {els_trackMomentumError_branch->SetAddress(&els_trackMomentumError_);}
	}
	els_trkdr_branch = 0;
	if (tree->GetAlias("els_trkdr") != 0) {
		els_trkdr_branch = tree->GetBranch(tree->GetAlias("els_trkdr"));
		if (els_trkdr_branch) {els_trkdr_branch->SetAddress(&els_trkdr_);}
	}
	els_trkshFrac_branch = 0;
	if (tree->GetAlias("els_trkshFrac") != 0) {
		els_trkshFrac_branch = tree->GetBranch(tree->GetAlias("els_trkshFrac"));
		if (els_trkshFrac_branch) {els_trkshFrac_branch->SetAddress(&els_trkshFrac_);}
	}
	els_z0_branch = 0;
	if (tree->GetAlias("els_z0") != 0) {
		els_z0_branch = tree->GetBranch(tree->GetAlias("els_z0"));
		if (els_z0_branch) {els_z0_branch->SetAddress(&els_z0_);}
	}
	els_z0Err_branch = 0;
	if (tree->GetAlias("els_z0Err") != 0) {
		els_z0Err_branch = tree->GetBranch(tree->GetAlias("els_z0Err"));
		if (els_z0Err_branch) {els_z0Err_branch->SetAddress(&els_z0Err_);}
	}
	els_z0corr_branch = 0;
	if (tree->GetAlias("els_z0corr") != 0) {
		els_z0corr_branch = tree->GetBranch(tree->GetAlias("els_z0corr"));
		if (els_z0corr_branch) {els_z0corr_branch->SetAddress(&els_z0corr_);}
	}
	photons_chargedHadronIso_branch = 0;
	if (tree->GetAlias("photons_chargedHadronIso") != 0) {
		photons_chargedHadronIso_branch = tree->GetBranch(tree->GetAlias("photons_chargedHadronIso"));
		if (photons_chargedHadronIso_branch) {photons_chargedHadronIso_branch->SetAddress(&photons_chargedHadronIso_);}
	}
	photons_e1x5_branch = 0;
	if (tree->GetAlias("photons_e1x5") != 0) {
		photons_e1x5_branch = tree->GetBranch(tree->GetAlias("photons_e1x5"));
		if (photons_e1x5_branch) {photons_e1x5_branch->SetAddress(&photons_e1x5_);}
	}
	photons_e2x5Max_branch = 0;
	if (tree->GetAlias("photons_e2x5Max") != 0) {
		photons_e2x5Max_branch = tree->GetBranch(tree->GetAlias("photons_e2x5Max"));
		if (photons_e2x5Max_branch) {photons_e2x5Max_branch->SetAddress(&photons_e2x5Max_);}
	}
	photons_e3x3_branch = 0;
	if (tree->GetAlias("photons_e3x3") != 0) {
		photons_e3x3_branch = tree->GetBranch(tree->GetAlias("photons_e3x3"));
		if (photons_e3x3_branch) {photons_e3x3_branch->SetAddress(&photons_e3x3_);}
	}
	photons_e5x5_branch = 0;
	if (tree->GetAlias("photons_e5x5") != 0) {
		photons_e5x5_branch = tree->GetBranch(tree->GetAlias("photons_e5x5"));
		if (photons_e5x5_branch) {photons_e5x5_branch->SetAddress(&photons_e5x5_);}
	}
	photons_eSC_branch = 0;
	if (tree->GetAlias("photons_eSC") != 0) {
		photons_eSC_branch = tree->GetBranch(tree->GetAlias("photons_eSC"));
		if (photons_eSC_branch) {photons_eSC_branch->SetAddress(&photons_eSC_);}
	}
	photons_eSCPresh_branch = 0;
	if (tree->GetAlias("photons_eSCPresh") != 0) {
		photons_eSCPresh_branch = tree->GetBranch(tree->GetAlias("photons_eSCPresh"));
		if (photons_eSCPresh_branch) {photons_eSCPresh_branch->SetAddress(&photons_eSCPresh_);}
	}
	photons_eSCRaw_branch = 0;
	if (tree->GetAlias("photons_eSCRaw") != 0) {
		photons_eSCRaw_branch = tree->GetBranch(tree->GetAlias("photons_eSCRaw"));
		if (photons_eSCRaw_branch) {photons_eSCRaw_branch->SetAddress(&photons_eSCRaw_);}
	}
	photons_ecalIso03_branch = 0;
	if (tree->GetAlias("photons_ecalIso03") != 0) {
		photons_ecalIso03_branch = tree->GetBranch(tree->GetAlias("photons_ecalIso03"));
		if (photons_ecalIso03_branch) {photons_ecalIso03_branch->SetAddress(&photons_ecalIso03_);}
	}
	photons_ecalIso04_branch = 0;
	if (tree->GetAlias("photons_ecalIso04") != 0) {
		photons_ecalIso04_branch = tree->GetBranch(tree->GetAlias("photons_ecalIso04"));
		if (photons_ecalIso04_branch) {photons_ecalIso04_branch->SetAddress(&photons_ecalIso04_);}
	}
	photons_ecalPFClusterIso_branch = 0;
	if (tree->GetAlias("photons_ecalPFClusterIso") != 0) {
		photons_ecalPFClusterIso_branch = tree->GetBranch(tree->GetAlias("photons_ecalPFClusterIso"));
		if (photons_ecalPFClusterIso_branch) {photons_ecalPFClusterIso_branch->SetAddress(&photons_ecalPFClusterIso_);}
	}
	photons_etaSC_branch = 0;
	if (tree->GetAlias("photons_etaSC") != 0) {
		photons_etaSC_branch = tree->GetBranch(tree->GetAlias("photons_etaSC"));
		if (photons_etaSC_branch) {photons_etaSC_branch->SetAddress(&photons_etaSC_);}
	}
	photons_full3x3_e3x3_branch = 0;
	if (tree->GetAlias("photons_full3x3_e3x3") != 0) {
		photons_full3x3_e3x3_branch = tree->GetBranch(tree->GetAlias("photons_full3x3_e3x3"));
		if (photons_full3x3_e3x3_branch) {photons_full3x3_e3x3_branch->SetAddress(&photons_full3x3_e3x3_);}
	}
	photons_full5x5_e1x5_branch = 0;
	if (tree->GetAlias("photons_full5x5_e1x5") != 0) {
		photons_full5x5_e1x5_branch = tree->GetBranch(tree->GetAlias("photons_full5x5_e1x5"));
		if (photons_full5x5_e1x5_branch) {photons_full5x5_e1x5_branch->SetAddress(&photons_full5x5_e1x5_);}
	}
	photons_full5x5_e2x5Max_branch = 0;
	if (tree->GetAlias("photons_full5x5_e2x5Max") != 0) {
		photons_full5x5_e2x5Max_branch = tree->GetBranch(tree->GetAlias("photons_full5x5_e2x5Max"));
		if (photons_full5x5_e2x5Max_branch) {photons_full5x5_e2x5Max_branch->SetAddress(&photons_full5x5_e2x5Max_);}
	}
	photons_full5x5_e5x5_branch = 0;
	if (tree->GetAlias("photons_full5x5_e5x5") != 0) {
		photons_full5x5_e5x5_branch = tree->GetBranch(tree->GetAlias("photons_full5x5_e5x5"));
		if (photons_full5x5_e5x5_branch) {photons_full5x5_e5x5_branch->SetAddress(&photons_full5x5_e5x5_);}
	}
	photons_full5x5_hOverE_branch = 0;
	if (tree->GetAlias("photons_full5x5_hOverE") != 0) {
		photons_full5x5_hOverE_branch = tree->GetBranch(tree->GetAlias("photons_full5x5_hOverE"));
		if (photons_full5x5_hOverE_branch) {photons_full5x5_hOverE_branch->SetAddress(&photons_full5x5_hOverE_);}
	}
	photons_full5x5_hOverEtowBC_branch = 0;
	if (tree->GetAlias("photons_full5x5_hOverEtowBC") != 0) {
		photons_full5x5_hOverEtowBC_branch = tree->GetBranch(tree->GetAlias("photons_full5x5_hOverEtowBC"));
		if (photons_full5x5_hOverEtowBC_branch) {photons_full5x5_hOverEtowBC_branch->SetAddress(&photons_full5x5_hOverEtowBC_);}
	}
	photons_full5x5_r9_branch = 0;
	if (tree->GetAlias("photons_full5x5_r9") != 0) {
		photons_full5x5_r9_branch = tree->GetBranch(tree->GetAlias("photons_full5x5_r9"));
		if (photons_full5x5_r9_branch) {photons_full5x5_r9_branch->SetAddress(&photons_full5x5_r9_);}
	}
	photons_full5x5_sigmaEtaEta_branch = 0;
	if (tree->GetAlias("photons_full5x5_sigmaEtaEta") != 0) {
		photons_full5x5_sigmaEtaEta_branch = tree->GetBranch(tree->GetAlias("photons_full5x5_sigmaEtaEta"));
		if (photons_full5x5_sigmaEtaEta_branch) {photons_full5x5_sigmaEtaEta_branch->SetAddress(&photons_full5x5_sigmaEtaEta_);}
	}
	photons_full5x5_sigmaIEtaIEta_branch = 0;
	if (tree->GetAlias("photons_full5x5_sigmaIEtaIEta") != 0) {
		photons_full5x5_sigmaIEtaIEta_branch = tree->GetBranch(tree->GetAlias("photons_full5x5_sigmaIEtaIEta"));
		if (photons_full5x5_sigmaIEtaIEta_branch) {photons_full5x5_sigmaIEtaIEta_branch->SetAddress(&photons_full5x5_sigmaIEtaIEta_);}
	}
	photons_hOverE_branch = 0;
	if (tree->GetAlias("photons_hOverE") != 0) {
		photons_hOverE_branch = tree->GetBranch(tree->GetAlias("photons_hOverE"));
		if (photons_hOverE_branch) {photons_hOverE_branch->SetAddress(&photons_hOverE_);}
	}
	photons_hOverEtowBC_branch = 0;
	if (tree->GetAlias("photons_hOverEtowBC") != 0) {
		photons_hOverEtowBC_branch = tree->GetBranch(tree->GetAlias("photons_hOverEtowBC"));
		if (photons_hOverEtowBC_branch) {photons_hOverEtowBC_branch->SetAddress(&photons_hOverEtowBC_);}
	}
	photons_hcalDepth1TowerSumEtBcConeDR03_branch = 0;
	if (tree->GetAlias("photons_hcalDepth1TowerSumEtBcConeDR03") != 0) {
		photons_hcalDepth1TowerSumEtBcConeDR03_branch = tree->GetBranch(tree->GetAlias("photons_hcalDepth1TowerSumEtBcConeDR03"));
		if (photons_hcalDepth1TowerSumEtBcConeDR03_branch) {photons_hcalDepth1TowerSumEtBcConeDR03_branch->SetAddress(&photons_hcalDepth1TowerSumEtBcConeDR03_);}
	}
	photons_hcalDepth1TowerSumEtBcConeDR04_branch = 0;
	if (tree->GetAlias("photons_hcalDepth1TowerSumEtBcConeDR04") != 0) {
		photons_hcalDepth1TowerSumEtBcConeDR04_branch = tree->GetBranch(tree->GetAlias("photons_hcalDepth1TowerSumEtBcConeDR04"));
		if (photons_hcalDepth1TowerSumEtBcConeDR04_branch) {photons_hcalDepth1TowerSumEtBcConeDR04_branch->SetAddress(&photons_hcalDepth1TowerSumEtBcConeDR04_);}
	}
	photons_hcalDepth2TowerSumEtBcConeDR03_branch = 0;
	if (tree->GetAlias("photons_hcalDepth2TowerSumEtBcConeDR03") != 0) {
		photons_hcalDepth2TowerSumEtBcConeDR03_branch = tree->GetBranch(tree->GetAlias("photons_hcalDepth2TowerSumEtBcConeDR03"));
		if (photons_hcalDepth2TowerSumEtBcConeDR03_branch) {photons_hcalDepth2TowerSumEtBcConeDR03_branch->SetAddress(&photons_hcalDepth2TowerSumEtBcConeDR03_);}
	}
	photons_hcalDepth2TowerSumEtBcConeDR04_branch = 0;
	if (tree->GetAlias("photons_hcalDepth2TowerSumEtBcConeDR04") != 0) {
		photons_hcalDepth2TowerSumEtBcConeDR04_branch = tree->GetBranch(tree->GetAlias("photons_hcalDepth2TowerSumEtBcConeDR04"));
		if (photons_hcalDepth2TowerSumEtBcConeDR04_branch) {photons_hcalDepth2TowerSumEtBcConeDR04_branch->SetAddress(&photons_hcalDepth2TowerSumEtBcConeDR04_);}
	}
	photons_hcalIso03_branch = 0;
	if (tree->GetAlias("photons_hcalIso03") != 0) {
		photons_hcalIso03_branch = tree->GetBranch(tree->GetAlias("photons_hcalIso03"));
		if (photons_hcalIso03_branch) {photons_hcalIso03_branch->SetAddress(&photons_hcalIso03_);}
	}
	photons_hcalIso04_branch = 0;
	if (tree->GetAlias("photons_hcalIso04") != 0) {
		photons_hcalIso04_branch = tree->GetBranch(tree->GetAlias("photons_hcalIso04"));
		if (photons_hcalIso04_branch) {photons_hcalIso04_branch->SetAddress(&photons_hcalIso04_);}
	}
	photons_hcalPFClusterIso_branch = 0;
	if (tree->GetAlias("photons_hcalPFClusterIso") != 0) {
		photons_hcalPFClusterIso_branch = tree->GetBranch(tree->GetAlias("photons_hcalPFClusterIso"));
		if (photons_hcalPFClusterIso_branch) {photons_hcalPFClusterIso_branch->SetAddress(&photons_hcalPFClusterIso_);}
	}
	photons_hcalTowerSumEtBcConeDR03_branch = 0;
	if (tree->GetAlias("photons_hcalTowerSumEtBcConeDR03") != 0) {
		photons_hcalTowerSumEtBcConeDR03_branch = tree->GetBranch(tree->GetAlias("photons_hcalTowerSumEtBcConeDR03"));
		if (photons_hcalTowerSumEtBcConeDR03_branch) {photons_hcalTowerSumEtBcConeDR03_branch->SetAddress(&photons_hcalTowerSumEtBcConeDR03_);}
	}
	photons_hcalTowerSumEtBcConeDR04_branch = 0;
	if (tree->GetAlias("photons_hcalTowerSumEtBcConeDR04") != 0) {
		photons_hcalTowerSumEtBcConeDR04_branch = tree->GetBranch(tree->GetAlias("photons_hcalTowerSumEtBcConeDR04"));
		if (photons_hcalTowerSumEtBcConeDR04_branch) {photons_hcalTowerSumEtBcConeDR04_branch->SetAddress(&photons_hcalTowerSumEtBcConeDR04_);}
	}
	photons_mass_branch = 0;
	if (tree->GetAlias("photons_mass") != 0) {
		photons_mass_branch = tree->GetBranch(tree->GetAlias("photons_mass"));
		if (photons_mass_branch) {photons_mass_branch->SetAddress(&photons_mass_);}
	}
	photons_neutralHadronIso_branch = 0;
	if (tree->GetAlias("photons_neutralHadronIso") != 0) {
		photons_neutralHadronIso_branch = tree->GetBranch(tree->GetAlias("photons_neutralHadronIso"));
		if (photons_neutralHadronIso_branch) {photons_neutralHadronIso_branch->SetAddress(&photons_neutralHadronIso_);}
	}
	photons_ntkIsoHollow03_branch = 0;
	if (tree->GetAlias("photons_ntkIsoHollow03") != 0) {
		photons_ntkIsoHollow03_branch = tree->GetBranch(tree->GetAlias("photons_ntkIsoHollow03"));
		if (photons_ntkIsoHollow03_branch) {photons_ntkIsoHollow03_branch->SetAddress(&photons_ntkIsoHollow03_);}
	}
	photons_ntkIsoHollow04_branch = 0;
	if (tree->GetAlias("photons_ntkIsoHollow04") != 0) {
		photons_ntkIsoHollow04_branch = tree->GetBranch(tree->GetAlias("photons_ntkIsoHollow04"));
		if (photons_ntkIsoHollow04_branch) {photons_ntkIsoHollow04_branch->SetAddress(&photons_ntkIsoHollow04_);}
	}
	photons_ntkIsoSolid03_branch = 0;
	if (tree->GetAlias("photons_ntkIsoSolid03") != 0) {
		photons_ntkIsoSolid03_branch = tree->GetBranch(tree->GetAlias("photons_ntkIsoSolid03"));
		if (photons_ntkIsoSolid03_branch) {photons_ntkIsoSolid03_branch->SetAddress(&photons_ntkIsoSolid03_);}
	}
	photons_ntkIsoSolid04_branch = 0;
	if (tree->GetAlias("photons_ntkIsoSolid04") != 0) {
		photons_ntkIsoSolid04_branch = tree->GetBranch(tree->GetAlias("photons_ntkIsoSolid04"));
		if (photons_ntkIsoSolid04_branch) {photons_ntkIsoSolid04_branch->SetAddress(&photons_ntkIsoSolid04_);}
	}
	photons_phiSC_branch = 0;
	if (tree->GetAlias("photons_phiSC") != 0) {
		photons_phiSC_branch = tree->GetBranch(tree->GetAlias("photons_phiSC"));
		if (photons_phiSC_branch) {photons_phiSC_branch->SetAddress(&photons_phiSC_);}
	}
	photons_photonIso_branch = 0;
	if (tree->GetAlias("photons_photonIso") != 0) {
		photons_photonIso_branch = tree->GetBranch(tree->GetAlias("photons_photonIso"));
		if (photons_photonIso_branch) {photons_photonIso_branch->SetAddress(&photons_photonIso_);}
	}
	photons_recoChargedHadronIso_branch = 0;
	if (tree->GetAlias("photons_recoChargedHadronIso") != 0) {
		photons_recoChargedHadronIso_branch = tree->GetBranch(tree->GetAlias("photons_recoChargedHadronIso"));
		if (photons_recoChargedHadronIso_branch) {photons_recoChargedHadronIso_branch->SetAddress(&photons_recoChargedHadronIso_);}
	}
	photons_recoNeutralHadronIso_branch = 0;
	if (tree->GetAlias("photons_recoNeutralHadronIso") != 0) {
		photons_recoNeutralHadronIso_branch = tree->GetBranch(tree->GetAlias("photons_recoNeutralHadronIso"));
		if (photons_recoNeutralHadronIso_branch) {photons_recoNeutralHadronIso_branch->SetAddress(&photons_recoNeutralHadronIso_);}
	}
	photons_recoPhotonIso_branch = 0;
	if (tree->GetAlias("photons_recoPhotonIso") != 0) {
		photons_recoPhotonIso_branch = tree->GetBranch(tree->GetAlias("photons_recoPhotonIso"));
		if (photons_recoPhotonIso_branch) {photons_recoPhotonIso_branch->SetAddress(&photons_recoPhotonIso_);}
	}
	photons_sigmaEtaEta_branch = 0;
	if (tree->GetAlias("photons_sigmaEtaEta") != 0) {
		photons_sigmaEtaEta_branch = tree->GetBranch(tree->GetAlias("photons_sigmaEtaEta"));
		if (photons_sigmaEtaEta_branch) {photons_sigmaEtaEta_branch->SetAddress(&photons_sigmaEtaEta_);}
	}
	photons_sigmaIEtaIEta_branch = 0;
	if (tree->GetAlias("photons_sigmaIEtaIEta") != 0) {
		photons_sigmaIEtaIEta_branch = tree->GetBranch(tree->GetAlias("photons_sigmaIEtaIEta"));
		if (photons_sigmaIEtaIEta_branch) {photons_sigmaIEtaIEta_branch->SetAddress(&photons_sigmaIEtaIEta_);}
	}
	photons_tkIsoHollow03_branch = 0;
	if (tree->GetAlias("photons_tkIsoHollow03") != 0) {
		photons_tkIsoHollow03_branch = tree->GetBranch(tree->GetAlias("photons_tkIsoHollow03"));
		if (photons_tkIsoHollow03_branch) {photons_tkIsoHollow03_branch->SetAddress(&photons_tkIsoHollow03_);}
	}
	photons_tkIsoHollow04_branch = 0;
	if (tree->GetAlias("photons_tkIsoHollow04") != 0) {
		photons_tkIsoHollow04_branch = tree->GetBranch(tree->GetAlias("photons_tkIsoHollow04"));
		if (photons_tkIsoHollow04_branch) {photons_tkIsoHollow04_branch->SetAddress(&photons_tkIsoHollow04_);}
	}
	photons_tkIsoSolid03_branch = 0;
	if (tree->GetAlias("photons_tkIsoSolid03") != 0) {
		photons_tkIsoSolid03_branch = tree->GetBranch(tree->GetAlias("photons_tkIsoSolid03"));
		if (photons_tkIsoSolid03_branch) {photons_tkIsoSolid03_branch->SetAddress(&photons_tkIsoSolid03_);}
	}
	photons_tkIsoSolid04_branch = 0;
	if (tree->GetAlias("photons_tkIsoSolid04") != 0) {
		photons_tkIsoSolid04_branch = tree->GetBranch(tree->GetAlias("photons_tkIsoSolid04"));
		if (photons_tkIsoSolid04_branch) {photons_tkIsoSolid04_branch->SetAddress(&photons_tkIsoSolid04_);}
	}
	puInfo_trueNumInteractions_branch = 0;
	if (tree->GetAlias("puInfo_trueNumInteractions") != 0) {
		puInfo_trueNumInteractions_branch = tree->GetBranch(tree->GetAlias("puInfo_trueNumInteractions"));
		if (puInfo_trueNumInteractions_branch) {puInfo_trueNumInteractions_branch->SetAddress(&puInfo_trueNumInteractions_);}
	}
	vtxs_chi2_branch = 0;
	if (tree->GetAlias("vtxs_chi2") != 0) {
		vtxs_chi2_branch = tree->GetBranch(tree->GetAlias("vtxs_chi2"));
		if (vtxs_chi2_branch) {vtxs_chi2_branch->SetAddress(&vtxs_chi2_);}
	}
	vtxs_ndof_branch = 0;
	if (tree->GetAlias("vtxs_ndof") != 0) {
		vtxs_ndof_branch = tree->GetBranch(tree->GetAlias("vtxs_ndof"));
		if (vtxs_ndof_branch) {vtxs_ndof_branch->SetAddress(&vtxs_ndof_);}
	}
	vtxs_score_branch = 0;
	if (tree->GetAlias("vtxs_score") != 0) {
		vtxs_score_branch = tree->GetBranch(tree->GetAlias("vtxs_score"));
		if (vtxs_score_branch) {vtxs_score_branch->SetAddress(&vtxs_score_);}
	}
	vtxs_xError_branch = 0;
	if (tree->GetAlias("vtxs_xError") != 0) {
		vtxs_xError_branch = tree->GetBranch(tree->GetAlias("vtxs_xError"));
		if (vtxs_xError_branch) {vtxs_xError_branch->SetAddress(&vtxs_xError_);}
	}
	vtxs_yError_branch = 0;
	if (tree->GetAlias("vtxs_yError") != 0) {
		vtxs_yError_branch = tree->GetBranch(tree->GetAlias("vtxs_yError"));
		if (vtxs_yError_branch) {vtxs_yError_branch->SetAddress(&vtxs_yError_);}
	}
	vtxs_zError_branch = 0;
	if (tree->GetAlias("vtxs_zError") != 0) {
		vtxs_zError_branch = tree->GetBranch(tree->GetAlias("vtxs_zError"));
		if (vtxs_zError_branch) {vtxs_zError_branch->SetAddress(&vtxs_zError_);}
	}
	puInfo_instLumi_branch = 0;
	if (tree->GetAlias("puInfo_instLumi") != 0) {
		puInfo_instLumi_branch = tree->GetBranch(tree->GetAlias("puInfo_instLumi"));
		if (puInfo_instLumi_branch) {puInfo_instLumi_branch->SetAddress(&puInfo_instLumi_);}
	}
	vtxs_covMatrix_branch = 0;
	if (tree->GetAlias("vtxs_covMatrix") != 0) {
		vtxs_covMatrix_branch = tree->GetBranch(tree->GetAlias("vtxs_covMatrix"));
		if (vtxs_covMatrix_branch) {vtxs_covMatrix_branch->SetAddress(&vtxs_covMatrix_);}
	}
	evt_bsType_branch = 0;
	if (tree->GetAlias("evt_bsType") != 0) {
		evt_bsType_branch = tree->GetBranch(tree->GetAlias("evt_bsType"));
		if (evt_bsType_branch) {evt_bsType_branch->SetAddress(&evt_bsType_);}
	}
	evt_bunchCrossing_branch = 0;
	if (tree->GetAlias("evt_bunchCrossing") != 0) {
		evt_bunchCrossing_branch = tree->GetBranch(tree->GetAlias("evt_bunchCrossing"));
		if (evt_bunchCrossing_branch) {evt_bunchCrossing_branch->SetAddress(&evt_bunchCrossing_);}
	}
	evt_experimentType_branch = 0;
	if (tree->GetAlias("evt_experimentType") != 0) {
		evt_experimentType_branch = tree->GetBranch(tree->GetAlias("evt_experimentType"));
		if (evt_experimentType_branch) {evt_experimentType_branch->SetAddress(&evt_experimentType_);}
	}
	evt_isRealData_branch = 0;
	if (tree->GetAlias("evt_isRealData") != 0) {
		evt_isRealData_branch = tree->GetBranch(tree->GetAlias("evt_isRealData"));
		if (evt_isRealData_branch) {evt_isRealData_branch->SetAddress(&evt_isRealData_);}
	}
	evt_orbitNumber_branch = 0;
	if (tree->GetAlias("evt_orbitNumber") != 0) {
		evt_orbitNumber_branch = tree->GetBranch(tree->GetAlias("evt_orbitNumber"));
		if (evt_orbitNumber_branch) {evt_orbitNumber_branch->SetAddress(&evt_orbitNumber_);}
	}
	evt_storeNumber_branch = 0;
	if (tree->GetAlias("evt_storeNumber") != 0) {
		evt_storeNumber_branch = tree->GetBranch(tree->GetAlias("evt_storeNumber"));
		if (evt_storeNumber_branch) {evt_storeNumber_branch->SetAddress(&evt_storeNumber_);}
	}
	els_category_branch = 0;
	if (tree->GetAlias("els_category") != 0) {
		els_category_branch = tree->GetBranch(tree->GetAlias("els_category"));
		if (els_category_branch) {els_category_branch->SetAddress(&els_category_);}
	}
	els_charge_branch = 0;
	if (tree->GetAlias("els_charge") != 0) {
		els_charge_branch = tree->GetBranch(tree->GetAlias("els_charge"));
		if (els_charge_branch) {els_charge_branch->SetAddress(&els_charge_);}
	}
	els_ckf_charge_branch = 0;
	if (tree->GetAlias("els_ckf_charge") != 0) {
		els_ckf_charge_branch = tree->GetBranch(tree->GetAlias("els_ckf_charge"));
		if (els_ckf_charge_branch) {els_ckf_charge_branch->SetAddress(&els_ckf_charge_);}
	}
	els_ckf_laywithmeas_branch = 0;
	if (tree->GetAlias("els_ckf_laywithmeas") != 0) {
		els_ckf_laywithmeas_branch = tree->GetBranch(tree->GetAlias("els_ckf_laywithmeas"));
		if (els_ckf_laywithmeas_branch) {els_ckf_laywithmeas_branch->SetAddress(&els_ckf_laywithmeas_);}
	}
	els_class_branch = 0;
	if (tree->GetAlias("els_class") != 0) {
		els_class_branch = tree->GetBranch(tree->GetAlias("els_class"));
		if (els_class_branch) {els_class_branch->SetAddress(&els_class_);}
	}
	els_exp_innerlayers_branch = 0;
	if (tree->GetAlias("els_exp_innerlayers") != 0) {
		els_exp_innerlayers_branch = tree->GetBranch(tree->GetAlias("els_exp_innerlayers"));
		if (els_exp_innerlayers_branch) {els_exp_innerlayers_branch->SetAddress(&els_exp_innerlayers_);}
	}
	els_exp_outerlayers_branch = 0;
	if (tree->GetAlias("els_exp_outerlayers") != 0) {
		els_exp_outerlayers_branch = tree->GetBranch(tree->GetAlias("els_exp_outerlayers"));
		if (els_exp_outerlayers_branch) {els_exp_outerlayers_branch->SetAddress(&els_exp_outerlayers_);}
	}
	els_fiduciality_branch = 0;
	if (tree->GetAlias("els_fiduciality") != 0) {
		els_fiduciality_branch = tree->GetBranch(tree->GetAlias("els_fiduciality"));
		if (els_fiduciality_branch) {els_fiduciality_branch->SetAddress(&els_fiduciality_);}
	}
	els_lostHits_branch = 0;
	if (tree->GetAlias("els_lostHits") != 0) {
		els_lostHits_branch = tree->GetBranch(tree->GetAlias("els_lostHits"));
		if (els_lostHits_branch) {els_lostHits_branch->SetAddress(&els_lostHits_);}
	}
	els_lost_pixelhits_branch = 0;
	if (tree->GetAlias("els_lost_pixelhits") != 0) {
		els_lost_pixelhits_branch = tree->GetBranch(tree->GetAlias("els_lost_pixelhits"));
		if (els_lost_pixelhits_branch) {els_lost_pixelhits_branch->SetAddress(&els_lost_pixelhits_);}
	}
	els_mc_patMatch_id_branch = 0;
	if (tree->GetAlias("els_mc_patMatch_id") != 0) {
		els_mc_patMatch_id_branch = tree->GetBranch(tree->GetAlias("els_mc_patMatch_id"));
		if (els_mc_patMatch_id_branch) {els_mc_patMatch_id_branch->SetAddress(&els_mc_patMatch_id_);}
	}
	els_nSeed_branch = 0;
	if (tree->GetAlias("els_nSeed") != 0) {
		els_nSeed_branch = tree->GetBranch(tree->GetAlias("els_nSeed"));
		if (els_nSeed_branch) {els_nSeed_branch->SetAddress(&els_nSeed_);}
	}
	els_nlayers_branch = 0;
	if (tree->GetAlias("els_nlayers") != 0) {
		els_nlayers_branch = tree->GetBranch(tree->GetAlias("els_nlayers"));
		if (els_nlayers_branch) {els_nlayers_branch->SetAddress(&els_nlayers_);}
	}
	els_nlayers3D_branch = 0;
	if (tree->GetAlias("els_nlayers3D") != 0) {
		els_nlayers3D_branch = tree->GetBranch(tree->GetAlias("els_nlayers3D"));
		if (els_nlayers3D_branch) {els_nlayers3D_branch->SetAddress(&els_nlayers3D_);}
	}
	els_nlayersLost_branch = 0;
	if (tree->GetAlias("els_nlayersLost") != 0) {
		els_nlayersLost_branch = tree->GetBranch(tree->GetAlias("els_nlayersLost"));
		if (els_nlayersLost_branch) {els_nlayersLost_branch->SetAddress(&els_nlayersLost_);}
	}
	els_sccharge_branch = 0;
	if (tree->GetAlias("els_sccharge") != 0) {
		els_sccharge_branch = tree->GetBranch(tree->GetAlias("els_sccharge"));
		if (els_sccharge_branch) {els_sccharge_branch->SetAddress(&els_sccharge_);}
	}
	els_trk_charge_branch = 0;
	if (tree->GetAlias("els_trk_charge") != 0) {
		els_trk_charge_branch = tree->GetBranch(tree->GetAlias("els_trk_charge"));
		if (els_trk_charge_branch) {els_trk_charge_branch->SetAddress(&els_trk_charge_);}
	}
	els_type_branch = 0;
	if (tree->GetAlias("els_type") != 0) {
		els_type_branch = tree->GetBranch(tree->GetAlias("els_type"));
		if (els_type_branch) {els_type_branch->SetAddress(&els_type_);}
	}
	els_validHits_branch = 0;
	if (tree->GetAlias("els_validHits") != 0) {
		els_validHits_branch = tree->GetBranch(tree->GetAlias("els_validHits"));
		if (els_validHits_branch) {els_validHits_branch->SetAddress(&els_validHits_);}
	}
	els_valid_pixelhits_branch = 0;
	if (tree->GetAlias("els_valid_pixelhits") != 0) {
		els_valid_pixelhits_branch = tree->GetBranch(tree->GetAlias("els_valid_pixelhits"));
		if (els_valid_pixelhits_branch) {els_valid_pixelhits_branch->SetAddress(&els_valid_pixelhits_);}
	}
	els_passLooseId_branch = 0;
	if (tree->GetAlias("els_passLooseId") != 0) {
		els_passLooseId_branch = tree->GetBranch(tree->GetAlias("els_passLooseId"));
		if (els_passLooseId_branch) {els_passLooseId_branch->SetAddress(&els_passLooseId_);}
	}
	els_passMediumId_branch = 0;
	if (tree->GetAlias("els_passMediumId") != 0) {
		els_passMediumId_branch = tree->GetBranch(tree->GetAlias("els_passMediumId"));
		if (els_passMediumId_branch) {els_passMediumId_branch->SetAddress(&els_passMediumId_);}
	}
	els_passTightId_branch = 0;
	if (tree->GetAlias("els_passTightId") != 0) {
		els_passTightId_branch = tree->GetBranch(tree->GetAlias("els_passTightId"));
		if (els_passTightId_branch) {els_passTightId_branch->SetAddress(&els_passTightId_);}
	}
	els_passVetoId_branch = 0;
	if (tree->GetAlias("els_passVetoId") != 0) {
		els_passVetoId_branch = tree->GetBranch(tree->GetAlias("els_passVetoId"));
		if (els_passVetoId_branch) {els_passVetoId_branch->SetAddress(&els_passVetoId_);}
	}
	photons_fiduciality_branch = 0;
	if (tree->GetAlias("photons_fiduciality") != 0) {
		photons_fiduciality_branch = tree->GetBranch(tree->GetAlias("photons_fiduciality"));
		if (photons_fiduciality_branch) {photons_fiduciality_branch->SetAddress(&photons_fiduciality_);}
	}
	photons_photonID_loose_branch = 0;
	if (tree->GetAlias("photons_photonID_loose") != 0) {
		photons_photonID_loose_branch = tree->GetBranch(tree->GetAlias("photons_photonID_loose"));
		if (photons_photonID_loose_branch) {photons_photonID_loose_branch->SetAddress(&photons_photonID_loose_);}
	}
	photons_photonID_tight_branch = 0;
	if (tree->GetAlias("photons_photonID_tight") != 0) {
		photons_photonID_tight_branch = tree->GetBranch(tree->GetAlias("photons_photonID_tight"));
		if (photons_photonID_tight_branch) {photons_photonID_tight_branch->SetAddress(&photons_photonID_tight_);}
	}
	puInfo_bunchCrossing_branch = 0;
	if (tree->GetAlias("puInfo_bunchCrossing") != 0) {
		puInfo_bunchCrossing_branch = tree->GetBranch(tree->GetAlias("puInfo_bunchCrossing"));
		if (puInfo_bunchCrossing_branch) {puInfo_bunchCrossing_branch->SetAddress(&puInfo_bunchCrossing_);}
	}
	puInfo_nPUvertices_branch = 0;
	if (tree->GetAlias("puInfo_nPUvertices") != 0) {
		puInfo_nPUvertices_branch = tree->GetBranch(tree->GetAlias("puInfo_nPUvertices"));
		if (puInfo_nPUvertices_branch) {puInfo_nPUvertices_branch->SetAddress(&puInfo_nPUvertices_);}
	}
	vtxs_isFake_branch = 0;
	if (tree->GetAlias("vtxs_isFake") != 0) {
		vtxs_isFake_branch = tree->GetBranch(tree->GetAlias("vtxs_isFake"));
		if (vtxs_isFake_branch) {vtxs_isFake_branch->SetAddress(&vtxs_isFake_);}
	}
	vtxs_isValid_branch = 0;
	if (tree->GetAlias("vtxs_isValid") != 0) {
		vtxs_isValid_branch = tree->GetBranch(tree->GetAlias("vtxs_isValid"));
		if (vtxs_isValid_branch) {vtxs_isValid_branch->SetAddress(&vtxs_isValid_);}
	}
	els_PFCand_idx_branch = 0;
	if (tree->GetAlias("els_PFCand_idx") != 0) {
		els_PFCand_idx_branch = tree->GetBranch(tree->GetAlias("els_PFCand_idx"));
		if (els_PFCand_idx_branch) {els_PFCand_idx_branch->SetAddress(&els_PFCand_idx_);}
	}
	hlt_trigObjs_id_branch = 0;
	if (tree->GetAlias("hlt_trigObjs_id") != 0) {
		hlt_trigObjs_id_branch = tree->GetBranch(tree->GetAlias("hlt_trigObjs_id"));
		if (hlt_trigObjs_id_branch) {hlt_trigObjs_id_branch->SetAddress(&hlt_trigObjs_id_);}
	}
	photons_PFCand_idx_branch = 0;
	if (tree->GetAlias("photons_PFCand_idx") != 0) {
		photons_PFCand_idx_branch = tree->GetBranch(tree->GetAlias("photons_PFCand_idx"));
		if (photons_PFCand_idx_branch) {photons_PFCand_idx_branch->SetAddress(&photons_PFCand_idx_);}
	}
	evt_nels_branch = 0;
	if (tree->GetAlias("evt_nels") != 0) {
		evt_nels_branch = tree->GetBranch(tree->GetAlias("evt_nels"));
		if (evt_nels_branch) {evt_nels_branch->SetAddress(&evt_nels_);}
	}
	evt_detectorStatus_branch = 0;
	if (tree->GetAlias("evt_detectorStatus") != 0) {
		evt_detectorStatus_branch = tree->GetBranch(tree->GetAlias("evt_detectorStatus"));
		if (evt_detectorStatus_branch) {evt_detectorStatus_branch->SetAddress(&evt_detectorStatus_);}
	}
	evt_lumiBlock_branch = 0;
	if (tree->GetAlias("evt_lumiBlock") != 0) {
		evt_lumiBlock_branch = tree->GetBranch(tree->GetAlias("evt_lumiBlock"));
		if (evt_lumiBlock_branch) {evt_lumiBlock_branch->SetAddress(&evt_lumiBlock_);}
	}
	evt_event_branch = 0;
	if (tree->GetAlias("evt_event") != 0) {
	  evt_event_branch = tree->GetBranch(tree->GetAlias("evt_event"));
	  if (evt_event_branch) {evt_event_branch->SetAddress(&evt_event_);}
	}
	evt_run_branch = 0;
	if (tree->GetAlias("evt_run") != 0) {
		evt_run_branch = tree->GetBranch(tree->GetAlias("evt_run"));
		if (evt_run_branch) {evt_run_branch->SetAddress(&evt_run_);}
	}
	evt_nphotons_branch = 0;
	if (tree->GetAlias("evt_nphotons") != 0) {
		evt_nphotons_branch = tree->GetBranch(tree->GetAlias("evt_nphotons"));
		if (evt_nphotons_branch) {evt_nphotons_branch->SetAddress(&evt_nphotons_);}
	}
	evt_nvtxs_branch = 0;
	if (tree->GetAlias("evt_nvtxs") != 0) {
		evt_nvtxs_branch = tree->GetBranch(tree->GetAlias("evt_nvtxs"));
		if (evt_nvtxs_branch) {evt_nvtxs_branch->SetAddress(&evt_nvtxs_);}
	}
	hlt_prescales_branch = 0;
	if (tree->GetAlias("hlt_prescales") != 0) {
		hlt_prescales_branch = tree->GetBranch(tree->GetAlias("hlt_prescales"));
		if (hlt_prescales_branch) {hlt_prescales_branch->SetAddress(&hlt_prescales_);}
	}
  tree->SetMakeClass(0);
}
void GetEntry(unsigned int idx) 
	// this only marks branches as not loaded, saving a lot of time
	{
		index = idx;
		hlt_bits_isLoaded = false;
		evt_CMS3tag_isLoaded = false;
		evt_dataset_isLoaded = false;
		hlt_trigNames_isLoaded = false;
		els_conv_vtx_flag_isLoaded = false;
		els_isGsfCtfScPixChargeConsistent_isLoaded = false;
		els_passingMvaPreselection_isLoaded = false;
		els_passingPflowPreselection_isLoaded = false;
		photons_haspixelSeed_isLoaded = false;
		evt_bs_Xwidth_isLoaded = false;
		evt_bs_XwidthErr_isLoaded = false;
		evt_bs_Ywidth_isLoaded = false;
		evt_bs_YwidthErr_isLoaded = false;
		evt_bs_dxdz_isLoaded = false;
		evt_bs_dxdzErr_isLoaded = false;
		evt_bs_dydz_isLoaded = false;
		evt_bs_dydzErr_isLoaded = false;
		evt_bs_sigmaZ_isLoaded = false;
		evt_bs_sigmaZErr_isLoaded = false;
		evt_bs_xErr_isLoaded = false;
		evt_bs_yErr_isLoaded = false;
		evt_bs_zErr_isLoaded = false;
		evt_bField_isLoaded = false;
		evt_fixgrid_all_rho_isLoaded = false;
		evt_fixgridfastjet_allcalo_rho_isLoaded = false;
		evt_fixgridfastjet_all_rho_isLoaded = false;
		evt_fixgridfastjet_centralcalo_rho_isLoaded = false;
		evt_fixgridfastjet_centralchargedpileup_rho_isLoaded = false;
		evt_fixgridfastjet_centralneutral_rho_isLoaded = false;
		evt_bsp4_isLoaded = false;
		els_mc_patMatch_p4_isLoaded = false;
		els_p4_isLoaded = false;
		els_p4In_isLoaded = false;
		els_p4Out_isLoaded = false;
		els_trk_p4_isLoaded = false;
		els_trk_vertex_p4_isLoaded = false;
		els_vertex_p4_isLoaded = false;
		photons_p4_isLoaded = false;
		vtxs_position_isLoaded = false;
		hlt_trigObjs_p4_isLoaded = false;
		evt_bs_covMatrix_isLoaded = false;
		els_bs2d_isLoaded = false;
		els_bs2derr_isLoaded = false;
		els_bs3d_isLoaded = false;
		els_bs3derr_isLoaded = false;
		els_chi2_isLoaded = false;
		els_ckf_chi2_isLoaded = false;
		els_ckf_ndof_isLoaded = false;
		els_d0_isLoaded = false;
		els_d0Err_isLoaded = false;
		els_d0corr_isLoaded = false;
		els_d0corrPhi_isLoaded = false;
		els_d0phiCov_isLoaded = false;
		els_dEtaIn_isLoaded = false;
		els_dEtaOut_isLoaded = false;
		els_dPhiIn_isLoaded = false;
		els_dPhiInPhiOut_isLoaded = false;
		els_dPhiOut_isLoaded = false;
		els_deltaEtaEleClusterTrackAtCalo_isLoaded = false;
		els_deltaPhiEleClusterTrackAtCalo_isLoaded = false;
		els_dxyPV_isLoaded = false;
		els_dzPV_isLoaded = false;
		els_e1x5_isLoaded = false;
		els_e1x5_full5x5_isLoaded = false;
		els_e2x5Max_isLoaded = false;
		els_e2x5Max_full5x5_isLoaded = false;
		els_e5x5_isLoaded = false;
		els_e5x5_full5x5_isLoaded = false;
		els_eOverPIn_isLoaded = false;
		els_eOverPOut_isLoaded = false;
		els_eSC_isLoaded = false;
		els_eSCPresh_isLoaded = false;
		els_eSCRaw_isLoaded = false;
		els_eSeed_isLoaded = false;
		els_eSeedOverPIn_isLoaded = false;
		els_eSeedOverPOut_isLoaded = false;
		els_ecalEnergy_isLoaded = false;
		els_ecalEnergyError_isLoaded = false;
		els_ecalIso_isLoaded = false;
		els_ecalIso04_isLoaded = false;
		els_ecalPFClusterIso_isLoaded = false;
		els_etaErr_isLoaded = false;
		els_etaSC_isLoaded = false;
		els_etaSCwidth_isLoaded = false;
		els_fbrem_isLoaded = false;
		els_hOverE_isLoaded = false;
		els_hOverEBC_isLoaded = false;
		els_hcalDepth1OverEcal_isLoaded = false;
		els_hcalDepth1TowerSumEt_isLoaded = false;
		els_hcalDepth1TowerSumEt04_isLoaded = false;
		els_hcalDepth2OverEcal_isLoaded = false;
		els_hcalDepth2TowerSumEt_isLoaded = false;
		els_hcalDepth2TowerSumEt04_isLoaded = false;
		els_hcalIso_isLoaded = false;
		els_hcalIso04_isLoaded = false;
		els_hcalPFClusterIso_isLoaded = false;
		els_ip2d_isLoaded = false;
		els_ip2derr_isLoaded = false;
		els_ip3d_isLoaded = false;
		els_ip3derr_isLoaded = false;
		els_mass_isLoaded = false;
		els_mc_patMatch_dr_isLoaded = false;
		els_miniIso_ch_isLoaded = false;
		els_miniIso_db_isLoaded = false;
		els_miniIso_em_isLoaded = false;
		els_miniIso_nh_isLoaded = false;
		els_miniIso_uncor_isLoaded = false;
		els_mva_isLoaded = false;
		els_ndof_isLoaded = false;
		els_pfChargedHadronIso_isLoaded = false;
		els_pfNeutralHadronIso_isLoaded = false;
		els_pfPUIso_isLoaded = false;
		els_pfPhotonIso_isLoaded = false;
		els_phiErr_isLoaded = false;
		els_phiSC_isLoaded = false;
		els_phiSCwidth_isLoaded = false;
		els_ptErr_isLoaded = false;
		els_ptErrGsf_isLoaded = false;
		els_r9_isLoaded = false;
		els_r9_full5x5_isLoaded = false;
		els_sigmaEtaEta_isLoaded = false;
		els_sigmaEtaEta_full5x5_isLoaded = false;
		els_sigmaIEtaIEta_isLoaded = false;
		els_sigmaIEtaIEta_full5x5_isLoaded = false;
		els_sigmaIPhiIPhi_isLoaded = false;
		els_sigmaIPhiIPhi_full5x5_isLoaded = false;
		els_sigmaIphiIphi_isLoaded = false;
		els_tkIso_isLoaded = false;
		els_tkIso04_isLoaded = false;
		els_trackMomentumError_isLoaded = false;
		els_trkdr_isLoaded = false;
		els_trkshFrac_isLoaded = false;
		els_z0_isLoaded = false;
		els_z0Err_isLoaded = false;
		els_z0corr_isLoaded = false;
		photons_chargedHadronIso_isLoaded = false;
		photons_e1x5_isLoaded = false;
		photons_e2x5Max_isLoaded = false;
		photons_e3x3_isLoaded = false;
		photons_e5x5_isLoaded = false;
		photons_eSC_isLoaded = false;
		photons_eSCPresh_isLoaded = false;
		photons_eSCRaw_isLoaded = false;
		photons_ecalIso03_isLoaded = false;
		photons_ecalIso04_isLoaded = false;
		photons_ecalPFClusterIso_isLoaded = false;
		photons_etaSC_isLoaded = false;
		photons_full3x3_e3x3_isLoaded = false;
		photons_full5x5_e1x5_isLoaded = false;
		photons_full5x5_e2x5Max_isLoaded = false;
		photons_full5x5_e5x5_isLoaded = false;
		photons_full5x5_hOverE_isLoaded = false;
		photons_full5x5_hOverEtowBC_isLoaded = false;
		photons_full5x5_r9_isLoaded = false;
		photons_full5x5_sigmaEtaEta_isLoaded = false;
		photons_full5x5_sigmaIEtaIEta_isLoaded = false;
		photons_hOverE_isLoaded = false;
		photons_hOverEtowBC_isLoaded = false;
		photons_hcalDepth1TowerSumEtBcConeDR03_isLoaded = false;
		photons_hcalDepth1TowerSumEtBcConeDR04_isLoaded = false;
		photons_hcalDepth2TowerSumEtBcConeDR03_isLoaded = false;
		photons_hcalDepth2TowerSumEtBcConeDR04_isLoaded = false;
		photons_hcalIso03_isLoaded = false;
		photons_hcalIso04_isLoaded = false;
		photons_hcalPFClusterIso_isLoaded = false;
		photons_hcalTowerSumEtBcConeDR03_isLoaded = false;
		photons_hcalTowerSumEtBcConeDR04_isLoaded = false;
		photons_mass_isLoaded = false;
		photons_neutralHadronIso_isLoaded = false;
		photons_ntkIsoHollow03_isLoaded = false;
		photons_ntkIsoHollow04_isLoaded = false;
		photons_ntkIsoSolid03_isLoaded = false;
		photons_ntkIsoSolid04_isLoaded = false;
		photons_phiSC_isLoaded = false;
		photons_photonIso_isLoaded = false;
		photons_recoChargedHadronIso_isLoaded = false;
		photons_recoNeutralHadronIso_isLoaded = false;
		photons_recoPhotonIso_isLoaded = false;
		photons_sigmaEtaEta_isLoaded = false;
		photons_sigmaIEtaIEta_isLoaded = false;
		photons_tkIsoHollow03_isLoaded = false;
		photons_tkIsoHollow04_isLoaded = false;
		photons_tkIsoSolid03_isLoaded = false;
		photons_tkIsoSolid04_isLoaded = false;
		puInfo_trueNumInteractions_isLoaded = false;
		vtxs_chi2_isLoaded = false;
		vtxs_ndof_isLoaded = false;
		vtxs_score_isLoaded = false;
		vtxs_xError_isLoaded = false;
		vtxs_yError_isLoaded = false;
		vtxs_zError_isLoaded = false;
		puInfo_instLumi_isLoaded = false;
		vtxs_covMatrix_isLoaded = false;
		evt_bsType_isLoaded = false;
		evt_bunchCrossing_isLoaded = false;
		evt_experimentType_isLoaded = false;
		evt_isRealData_isLoaded = false;
		evt_orbitNumber_isLoaded = false;
		evt_storeNumber_isLoaded = false;
		els_category_isLoaded = false;
		els_charge_isLoaded = false;
		els_ckf_charge_isLoaded = false;
		els_ckf_laywithmeas_isLoaded = false;
		els_class_isLoaded = false;
		els_exp_innerlayers_isLoaded = false;
		els_exp_outerlayers_isLoaded = false;
		els_fiduciality_isLoaded = false;
		els_lostHits_isLoaded = false;
		els_lost_pixelhits_isLoaded = false;
		els_mc_patMatch_id_isLoaded = false;
		els_nSeed_isLoaded = false;
		els_nlayers_isLoaded = false;
		els_nlayers3D_isLoaded = false;
		els_nlayersLost_isLoaded = false;
		els_sccharge_isLoaded = false;
		els_trk_charge_isLoaded = false;
		els_type_isLoaded = false;
		els_validHits_isLoaded = false;
		els_valid_pixelhits_isLoaded = false;
		els_passLooseId_isLoaded = false;
		els_passMediumId_isLoaded = false;
		els_passTightId_isLoaded = false;
		els_passVetoId_isLoaded = false;
		photons_fiduciality_isLoaded = false;
		photons_photonID_loose_isLoaded = false;
		photons_photonID_tight_isLoaded = false;
		puInfo_bunchCrossing_isLoaded = false;
		puInfo_nPUvertices_isLoaded = false;
		vtxs_isFake_isLoaded = false;
		vtxs_isValid_isLoaded = false;
		els_PFCand_idx_isLoaded = false;
		hlt_trigObjs_id_isLoaded = false;
		photons_PFCand_idx_isLoaded = false;
		evt_nels_isLoaded = false;
		evt_detectorStatus_isLoaded = false;
		evt_lumiBlock_isLoaded = false;
		evt_event_isLoaded = false;
		evt_run_isLoaded = false;
		evt_nphotons_isLoaded = false;
		evt_nvtxs_isLoaded = false;
		hlt_prescales_isLoaded = false;
	}

void LoadAllBranches() 
	// load all branches
{
	if (hlt_bits_branch != 0) hlt_bits();
	if (evt_CMS3tag_branch != 0) evt_CMS3tag();
	if (evt_dataset_branch != 0) evt_dataset();
	if (hlt_trigNames_branch != 0) hlt_trigNames();
	if (els_conv_vtx_flag_branch != 0) els_conv_vtx_flag();
	if (els_isGsfCtfScPixChargeConsistent_branch != 0) els_isGsfCtfScPixChargeConsistent();
	if (els_passingMvaPreselection_branch != 0) els_passingMvaPreselection();
	if (els_passingPflowPreselection_branch != 0) els_passingPflowPreselection();
	if (photons_haspixelSeed_branch != 0) photons_haspixelSeed();
	if (evt_bs_Xwidth_branch != 0) evt_bs_Xwidth();
	if (evt_bs_XwidthErr_branch != 0) evt_bs_XwidthErr();
	if (evt_bs_Ywidth_branch != 0) evt_bs_Ywidth();
	if (evt_bs_YwidthErr_branch != 0) evt_bs_YwidthErr();
	if (evt_bs_dxdz_branch != 0) evt_bs_dxdz();
	if (evt_bs_dxdzErr_branch != 0) evt_bs_dxdzErr();
	if (evt_bs_dydz_branch != 0) evt_bs_dydz();
	if (evt_bs_dydzErr_branch != 0) evt_bs_dydzErr();
	if (evt_bs_sigmaZ_branch != 0) evt_bs_sigmaZ();
	if (evt_bs_sigmaZErr_branch != 0) evt_bs_sigmaZErr();
	if (evt_bs_xErr_branch != 0) evt_bs_xErr();
	if (evt_bs_yErr_branch != 0) evt_bs_yErr();
	if (evt_bs_zErr_branch != 0) evt_bs_zErr();
	if (evt_bField_branch != 0) evt_bField();
	if (evt_fixgrid_all_rho_branch != 0) evt_fixgrid_all_rho();
	if (evt_fixgridfastjet_allcalo_rho_branch != 0) evt_fixgridfastjet_allcalo_rho();
	if (evt_fixgridfastjet_all_rho_branch != 0) evt_fixgridfastjet_all_rho();
	if (evt_fixgridfastjet_centralcalo_rho_branch != 0) evt_fixgridfastjet_centralcalo_rho();
	if (evt_fixgridfastjet_centralchargedpileup_rho_branch != 0) evt_fixgridfastjet_centralchargedpileup_rho();
	if (evt_fixgridfastjet_centralneutral_rho_branch != 0) evt_fixgridfastjet_centralneutral_rho();
	if (evt_bsp4_branch != 0) evt_bsp4();
	if (els_mc_patMatch_p4_branch != 0) els_mc_patMatch_p4();
	if (els_p4_branch != 0) els_p4();
	if (els_p4In_branch != 0) els_p4In();
	if (els_p4Out_branch != 0) els_p4Out();
	if (els_trk_p4_branch != 0) els_trk_p4();
	if (els_trk_vertex_p4_branch != 0) els_trk_vertex_p4();
	if (els_vertex_p4_branch != 0) els_vertex_p4();
	if (photons_p4_branch != 0) photons_p4();
	if (vtxs_position_branch != 0) vtxs_position();
	if (hlt_trigObjs_p4_branch != 0) hlt_trigObjs_p4();
	if (evt_bs_covMatrix_branch != 0) evt_bs_covMatrix();
	if (els_bs2d_branch != 0) els_bs2d();
	if (els_bs2derr_branch != 0) els_bs2derr();
	if (els_bs3d_branch != 0) els_bs3d();
	if (els_bs3derr_branch != 0) els_bs3derr();
	if (els_chi2_branch != 0) els_chi2();
	if (els_ckf_chi2_branch != 0) els_ckf_chi2();
	if (els_ckf_ndof_branch != 0) els_ckf_ndof();
	if (els_d0_branch != 0) els_d0();
	if (els_d0Err_branch != 0) els_d0Err();
	if (els_d0corr_branch != 0) els_d0corr();
	if (els_d0corrPhi_branch != 0) els_d0corrPhi();
	if (els_d0phiCov_branch != 0) els_d0phiCov();
	if (els_dEtaIn_branch != 0) els_dEtaIn();
	if (els_dEtaOut_branch != 0) els_dEtaOut();
	if (els_dPhiIn_branch != 0) els_dPhiIn();
	if (els_dPhiInPhiOut_branch != 0) els_dPhiInPhiOut();
	if (els_dPhiOut_branch != 0) els_dPhiOut();
	if (els_deltaEtaEleClusterTrackAtCalo_branch != 0) els_deltaEtaEleClusterTrackAtCalo();
	if (els_deltaPhiEleClusterTrackAtCalo_branch != 0) els_deltaPhiEleClusterTrackAtCalo();
	if (els_dxyPV_branch != 0) els_dxyPV();
	if (els_dzPV_branch != 0) els_dzPV();
	if (els_e1x5_branch != 0) els_e1x5();
	if (els_e1x5_full5x5_branch != 0) els_e1x5_full5x5();
	if (els_e2x5Max_branch != 0) els_e2x5Max();
	if (els_e2x5Max_full5x5_branch != 0) els_e2x5Max_full5x5();
	if (els_e5x5_branch != 0) els_e5x5();
	if (els_e5x5_full5x5_branch != 0) els_e5x5_full5x5();
	if (els_eOverPIn_branch != 0) els_eOverPIn();
	if (els_eOverPOut_branch != 0) els_eOverPOut();
	if (els_eSC_branch != 0) els_eSC();
	if (els_eSCPresh_branch != 0) els_eSCPresh();
	if (els_eSCRaw_branch != 0) els_eSCRaw();
	if (els_eSeed_branch != 0) els_eSeed();
	if (els_eSeedOverPIn_branch != 0) els_eSeedOverPIn();
	if (els_eSeedOverPOut_branch != 0) els_eSeedOverPOut();
	if (els_ecalEnergy_branch != 0) els_ecalEnergy();
	if (els_ecalEnergyError_branch != 0) els_ecalEnergyError();
	if (els_ecalIso_branch != 0) els_ecalIso();
	if (els_ecalIso04_branch != 0) els_ecalIso04();
	if (els_ecalPFClusterIso_branch != 0) els_ecalPFClusterIso();
	if (els_etaErr_branch != 0) els_etaErr();
	if (els_etaSC_branch != 0) els_etaSC();
	if (els_etaSCwidth_branch != 0) els_etaSCwidth();
	if (els_fbrem_branch != 0) els_fbrem();
	if (els_hOverE_branch != 0) els_hOverE();
	if (els_hOverEBC_branch != 0) els_hOverEBC();
	if (els_hcalDepth1OverEcal_branch != 0) els_hcalDepth1OverEcal();
	if (els_hcalDepth1TowerSumEt_branch != 0) els_hcalDepth1TowerSumEt();
	if (els_hcalDepth1TowerSumEt04_branch != 0) els_hcalDepth1TowerSumEt04();
	if (els_hcalDepth2OverEcal_branch != 0) els_hcalDepth2OverEcal();
	if (els_hcalDepth2TowerSumEt_branch != 0) els_hcalDepth2TowerSumEt();
	if (els_hcalDepth2TowerSumEt04_branch != 0) els_hcalDepth2TowerSumEt04();
	if (els_hcalIso_branch != 0) els_hcalIso();
	if (els_hcalIso04_branch != 0) els_hcalIso04();
	if (els_hcalPFClusterIso_branch != 0) els_hcalPFClusterIso();
	if (els_ip2d_branch != 0) els_ip2d();
	if (els_ip2derr_branch != 0) els_ip2derr();
	if (els_ip3d_branch != 0) els_ip3d();
	if (els_ip3derr_branch != 0) els_ip3derr();
	if (els_mass_branch != 0) els_mass();
	if (els_mc_patMatch_dr_branch != 0) els_mc_patMatch_dr();
	if (els_miniIso_ch_branch != 0) els_miniIso_ch();
	if (els_miniIso_db_branch != 0) els_miniIso_db();
	if (els_miniIso_em_branch != 0) els_miniIso_em();
	if (els_miniIso_nh_branch != 0) els_miniIso_nh();
	if (els_miniIso_uncor_branch != 0) els_miniIso_uncor();
	if (els_mva_branch != 0) els_mva();
	if (els_ndof_branch != 0) els_ndof();
	if (els_pfChargedHadronIso_branch != 0) els_pfChargedHadronIso();
	if (els_pfNeutralHadronIso_branch != 0) els_pfNeutralHadronIso();
	if (els_pfPUIso_branch != 0) els_pfPUIso();
	if (els_pfPhotonIso_branch != 0) els_pfPhotonIso();
	if (els_phiErr_branch != 0) els_phiErr();
	if (els_phiSC_branch != 0) els_phiSC();
	if (els_phiSCwidth_branch != 0) els_phiSCwidth();
	if (els_ptErr_branch != 0) els_ptErr();
	if (els_ptErrGsf_branch != 0) els_ptErrGsf();
	if (els_r9_branch != 0) els_r9();
	if (els_r9_full5x5_branch != 0) els_r9_full5x5();
	if (els_sigmaEtaEta_branch != 0) els_sigmaEtaEta();
	if (els_sigmaEtaEta_full5x5_branch != 0) els_sigmaEtaEta_full5x5();
	if (els_sigmaIEtaIEta_branch != 0) els_sigmaIEtaIEta();
	if (els_sigmaIEtaIEta_full5x5_branch != 0) els_sigmaIEtaIEta_full5x5();
	if (els_sigmaIPhiIPhi_branch != 0) els_sigmaIPhiIPhi();
	if (els_sigmaIPhiIPhi_full5x5_branch != 0) els_sigmaIPhiIPhi_full5x5();
	if (els_sigmaIphiIphi_branch != 0) els_sigmaIphiIphi();
	if (els_tkIso_branch != 0) els_tkIso();
	if (els_tkIso04_branch != 0) els_tkIso04();
	if (els_trackMomentumError_branch != 0) els_trackMomentumError();
	if (els_trkdr_branch != 0) els_trkdr();
	if (els_trkshFrac_branch != 0) els_trkshFrac();
	if (els_z0_branch != 0) els_z0();
	if (els_z0Err_branch != 0) els_z0Err();
	if (els_z0corr_branch != 0) els_z0corr();
	if (photons_chargedHadronIso_branch != 0) photons_chargedHadronIso();
	if (photons_e1x5_branch != 0) photons_e1x5();
	if (photons_e2x5Max_branch != 0) photons_e2x5Max();
	if (photons_e3x3_branch != 0) photons_e3x3();
	if (photons_e5x5_branch != 0) photons_e5x5();
	if (photons_eSC_branch != 0) photons_eSC();
	if (photons_eSCPresh_branch != 0) photons_eSCPresh();
	if (photons_eSCRaw_branch != 0) photons_eSCRaw();
	if (photons_ecalIso03_branch != 0) photons_ecalIso03();
	if (photons_ecalIso04_branch != 0) photons_ecalIso04();
	if (photons_ecalPFClusterIso_branch != 0) photons_ecalPFClusterIso();
	if (photons_etaSC_branch != 0) photons_etaSC();
	if (photons_full3x3_e3x3_branch != 0) photons_full3x3_e3x3();
	if (photons_full5x5_e1x5_branch != 0) photons_full5x5_e1x5();
	if (photons_full5x5_e2x5Max_branch != 0) photons_full5x5_e2x5Max();
	if (photons_full5x5_e5x5_branch != 0) photons_full5x5_e5x5();
	if (photons_full5x5_hOverE_branch != 0) photons_full5x5_hOverE();
	if (photons_full5x5_hOverEtowBC_branch != 0) photons_full5x5_hOverEtowBC();
	if (photons_full5x5_r9_branch != 0) photons_full5x5_r9();
	if (photons_full5x5_sigmaEtaEta_branch != 0) photons_full5x5_sigmaEtaEta();
	if (photons_full5x5_sigmaIEtaIEta_branch != 0) photons_full5x5_sigmaIEtaIEta();
	if (photons_hOverE_branch != 0) photons_hOverE();
	if (photons_hOverEtowBC_branch != 0) photons_hOverEtowBC();
	if (photons_hcalDepth1TowerSumEtBcConeDR03_branch != 0) photons_hcalDepth1TowerSumEtBcConeDR03();
	if (photons_hcalDepth1TowerSumEtBcConeDR04_branch != 0) photons_hcalDepth1TowerSumEtBcConeDR04();
	if (photons_hcalDepth2TowerSumEtBcConeDR03_branch != 0) photons_hcalDepth2TowerSumEtBcConeDR03();
	if (photons_hcalDepth2TowerSumEtBcConeDR04_branch != 0) photons_hcalDepth2TowerSumEtBcConeDR04();
	if (photons_hcalIso03_branch != 0) photons_hcalIso03();
	if (photons_hcalIso04_branch != 0) photons_hcalIso04();
	if (photons_hcalPFClusterIso_branch != 0) photons_hcalPFClusterIso();
	if (photons_hcalTowerSumEtBcConeDR03_branch != 0) photons_hcalTowerSumEtBcConeDR03();
	if (photons_hcalTowerSumEtBcConeDR04_branch != 0) photons_hcalTowerSumEtBcConeDR04();
	if (photons_mass_branch != 0) photons_mass();
	if (photons_neutralHadronIso_branch != 0) photons_neutralHadronIso();
	if (photons_ntkIsoHollow03_branch != 0) photons_ntkIsoHollow03();
	if (photons_ntkIsoHollow04_branch != 0) photons_ntkIsoHollow04();
	if (photons_ntkIsoSolid03_branch != 0) photons_ntkIsoSolid03();
	if (photons_ntkIsoSolid04_branch != 0) photons_ntkIsoSolid04();
	if (photons_phiSC_branch != 0) photons_phiSC();
	if (photons_photonIso_branch != 0) photons_photonIso();
	if (photons_recoChargedHadronIso_branch != 0) photons_recoChargedHadronIso();
	if (photons_recoNeutralHadronIso_branch != 0) photons_recoNeutralHadronIso();
	if (photons_recoPhotonIso_branch != 0) photons_recoPhotonIso();
	if (photons_sigmaEtaEta_branch != 0) photons_sigmaEtaEta();
	if (photons_sigmaIEtaIEta_branch != 0) photons_sigmaIEtaIEta();
	if (photons_tkIsoHollow03_branch != 0) photons_tkIsoHollow03();
	if (photons_tkIsoHollow04_branch != 0) photons_tkIsoHollow04();
	if (photons_tkIsoSolid03_branch != 0) photons_tkIsoSolid03();
	if (photons_tkIsoSolid04_branch != 0) photons_tkIsoSolid04();
	if (puInfo_trueNumInteractions_branch != 0) puInfo_trueNumInteractions();
	if (vtxs_chi2_branch != 0) vtxs_chi2();
	if (vtxs_ndof_branch != 0) vtxs_ndof();
	if (vtxs_score_branch != 0) vtxs_score();
	if (vtxs_xError_branch != 0) vtxs_xError();
	if (vtxs_yError_branch != 0) vtxs_yError();
	if (vtxs_zError_branch != 0) vtxs_zError();
	if (puInfo_instLumi_branch != 0) puInfo_instLumi();
	if (vtxs_covMatrix_branch != 0) vtxs_covMatrix();
	if (evt_bsType_branch != 0) evt_bsType();
	if (evt_bunchCrossing_branch != 0) evt_bunchCrossing();
	if (evt_experimentType_branch != 0) evt_experimentType();
	if (evt_isRealData_branch != 0) evt_isRealData();
	if (evt_orbitNumber_branch != 0) evt_orbitNumber();
	if (evt_storeNumber_branch != 0) evt_storeNumber();
	if (els_category_branch != 0) els_category();
	if (els_charge_branch != 0) els_charge();
	if (els_ckf_charge_branch != 0) els_ckf_charge();
	if (els_ckf_laywithmeas_branch != 0) els_ckf_laywithmeas();
	if (els_class_branch != 0) els_class();
	if (els_exp_innerlayers_branch != 0) els_exp_innerlayers();
	if (els_exp_outerlayers_branch != 0) els_exp_outerlayers();
	if (els_fiduciality_branch != 0) els_fiduciality();
	if (els_lostHits_branch != 0) els_lostHits();
	if (els_lost_pixelhits_branch != 0) els_lost_pixelhits();
	if (els_mc_patMatch_id_branch != 0) els_mc_patMatch_id();
	if (els_nSeed_branch != 0) els_nSeed();
	if (els_nlayers_branch != 0) els_nlayers();
	if (els_nlayers3D_branch != 0) els_nlayers3D();
	if (els_nlayersLost_branch != 0) els_nlayersLost();
	if (els_sccharge_branch != 0) els_sccharge();
	if (els_trk_charge_branch != 0) els_trk_charge();
	if (els_type_branch != 0) els_type();
	if (els_validHits_branch != 0) els_validHits();
	if (els_valid_pixelhits_branch != 0) els_valid_pixelhits();
	if (els_passLooseId_branch != 0) els_passLooseId();
	if (els_passMediumId_branch != 0) els_passMediumId();
	if (els_passTightId_branch != 0) els_passTightId();
	if (els_passVetoId_branch != 0) els_passVetoId();
	if (photons_fiduciality_branch != 0) photons_fiduciality();
	if (photons_photonID_loose_branch != 0) photons_photonID_loose();
	if (photons_photonID_tight_branch != 0) photons_photonID_tight();
	if (puInfo_bunchCrossing_branch != 0) puInfo_bunchCrossing();
	if (puInfo_nPUvertices_branch != 0) puInfo_nPUvertices();
	if (vtxs_isFake_branch != 0) vtxs_isFake();
	if (vtxs_isValid_branch != 0) vtxs_isValid();
	if (els_PFCand_idx_branch != 0) els_PFCand_idx();
	if (hlt_trigObjs_id_branch != 0) hlt_trigObjs_id();
	if (photons_PFCand_idx_branch != 0) photons_PFCand_idx();
	if (evt_nels_branch != 0) evt_nels();
	if (evt_detectorStatus_branch != 0) evt_detectorStatus();
	if (evt_lumiBlock_branch != 0) evt_lumiBlock();
	if (evt_event_branch != 0) evt_event();
	if (evt_run_branch != 0) evt_run();
	if (evt_nphotons_branch != 0) evt_nphotons();
	if (evt_nvtxs_branch != 0) evt_nvtxs();
	if (hlt_prescales_branch != 0) hlt_prescales();
}

	TBits &hlt_bits()
	{
		if (not hlt_bits_isLoaded) {
			if (hlt_bits_branch != 0) {
				hlt_bits_branch->GetEntry(index);
			} else { 
				printf("branch hlt_bits_branch does not exist!\n");
				exit(1);
			}
			hlt_bits_isLoaded = true;
		}
		return hlt_bits_;
	}
	const vector<TString> &evt_CMS3tag()
	{
		if (not evt_CMS3tag_isLoaded) {
			if (evt_CMS3tag_branch != 0) {
				evt_CMS3tag_branch->GetEntry(index);
			} else { 
				printf("branch evt_CMS3tag_branch does not exist!\n");
				exit(1);
			}
			evt_CMS3tag_isLoaded = true;
		}
		return evt_CMS3tag_;
	}
	const vector<TString> &evt_dataset()
	{
		if (not evt_dataset_isLoaded) {
			if (evt_dataset_branch != 0) {
				evt_dataset_branch->GetEntry(index);
			} else { 
				printf("branch evt_dataset_branch does not exist!\n");
				exit(1);
			}
			evt_dataset_isLoaded = true;
		}
		return evt_dataset_;
	}
	const vector<TString> &hlt_trigNames()
	{
		if (not hlt_trigNames_isLoaded) {
			if (hlt_trigNames_branch != 0) {
				hlt_trigNames_branch->GetEntry(index);
			} else { 
				printf("branch hlt_trigNames_branch does not exist!\n");
				exit(1);
			}
			hlt_trigNames_isLoaded = true;
		}
		return hlt_trigNames_;
	}
	const vector<bool> &els_conv_vtx_flag()
	{
		if (not els_conv_vtx_flag_isLoaded) {
			if (els_conv_vtx_flag_branch != 0) {
				els_conv_vtx_flag_branch->GetEntry(index);
			} else { 
				printf("branch els_conv_vtx_flag_branch does not exist!\n");
				exit(1);
			}
			els_conv_vtx_flag_isLoaded = true;
		}
		return els_conv_vtx_flag_;
	}
	const vector<bool> &els_isGsfCtfScPixChargeConsistent()
	{
		if (not els_isGsfCtfScPixChargeConsistent_isLoaded) {
			if (els_isGsfCtfScPixChargeConsistent_branch != 0) {
				els_isGsfCtfScPixChargeConsistent_branch->GetEntry(index);
			} else { 
				printf("branch els_isGsfCtfScPixChargeConsistent_branch does not exist!\n");
				exit(1);
			}
			els_isGsfCtfScPixChargeConsistent_isLoaded = true;
		}
		return els_isGsfCtfScPixChargeConsistent_;
	}
	const vector<bool> &els_passingMvaPreselection()
	{
		if (not els_passingMvaPreselection_isLoaded) {
			if (els_passingMvaPreselection_branch != 0) {
				els_passingMvaPreselection_branch->GetEntry(index);
			} else { 
				printf("branch els_passingMvaPreselection_branch does not exist!\n");
				exit(1);
			}
			els_passingMvaPreselection_isLoaded = true;
		}
		return els_passingMvaPreselection_;
	}
	const vector<bool> &els_passingPflowPreselection()
	{
		if (not els_passingPflowPreselection_isLoaded) {
			if (els_passingPflowPreselection_branch != 0) {
				els_passingPflowPreselection_branch->GetEntry(index);
			} else { 
				printf("branch els_passingPflowPreselection_branch does not exist!\n");
				exit(1);
			}
			els_passingPflowPreselection_isLoaded = true;
		}
		return els_passingPflowPreselection_;
	}
	const vector<bool> &photons_haspixelSeed()
	{
		if (not photons_haspixelSeed_isLoaded) {
			if (photons_haspixelSeed_branch != 0) {
				photons_haspixelSeed_branch->GetEntry(index);
			} else { 
				printf("branch photons_haspixelSeed_branch does not exist!\n");
				exit(1);
			}
			photons_haspixelSeed_isLoaded = true;
		}
		return photons_haspixelSeed_;
	}
	float &evt_bs_Xwidth()
	{
		if (not evt_bs_Xwidth_isLoaded) {
			if (evt_bs_Xwidth_branch != 0) {
				evt_bs_Xwidth_branch->GetEntry(index);
			} else { 
				printf("branch evt_bs_Xwidth_branch does not exist!\n");
				exit(1);
			}
			evt_bs_Xwidth_isLoaded = true;
		}
		return evt_bs_Xwidth_;
	}
	float &evt_bs_XwidthErr()
	{
		if (not evt_bs_XwidthErr_isLoaded) {
			if (evt_bs_XwidthErr_branch != 0) {
				evt_bs_XwidthErr_branch->GetEntry(index);
			} else { 
				printf("branch evt_bs_XwidthErr_branch does not exist!\n");
				exit(1);
			}
			evt_bs_XwidthErr_isLoaded = true;
		}
		return evt_bs_XwidthErr_;
	}
	float &evt_bs_Ywidth()
	{
		if (not evt_bs_Ywidth_isLoaded) {
			if (evt_bs_Ywidth_branch != 0) {
				evt_bs_Ywidth_branch->GetEntry(index);
			} else { 
				printf("branch evt_bs_Ywidth_branch does not exist!\n");
				exit(1);
			}
			evt_bs_Ywidth_isLoaded = true;
		}
		return evt_bs_Ywidth_;
	}
	float &evt_bs_YwidthErr()
	{
		if (not evt_bs_YwidthErr_isLoaded) {
			if (evt_bs_YwidthErr_branch != 0) {
				evt_bs_YwidthErr_branch->GetEntry(index);
			} else { 
				printf("branch evt_bs_YwidthErr_branch does not exist!\n");
				exit(1);
			}
			evt_bs_YwidthErr_isLoaded = true;
		}
		return evt_bs_YwidthErr_;
	}
	float &evt_bs_dxdz()
	{
		if (not evt_bs_dxdz_isLoaded) {
			if (evt_bs_dxdz_branch != 0) {
				evt_bs_dxdz_branch->GetEntry(index);
			} else { 
				printf("branch evt_bs_dxdz_branch does not exist!\n");
				exit(1);
			}
			evt_bs_dxdz_isLoaded = true;
		}
		return evt_bs_dxdz_;
	}
	float &evt_bs_dxdzErr()
	{
		if (not evt_bs_dxdzErr_isLoaded) {
			if (evt_bs_dxdzErr_branch != 0) {
				evt_bs_dxdzErr_branch->GetEntry(index);
			} else { 
				printf("branch evt_bs_dxdzErr_branch does not exist!\n");
				exit(1);
			}
			evt_bs_dxdzErr_isLoaded = true;
		}
		return evt_bs_dxdzErr_;
	}
	float &evt_bs_dydz()
	{
		if (not evt_bs_dydz_isLoaded) {
			if (evt_bs_dydz_branch != 0) {
				evt_bs_dydz_branch->GetEntry(index);
			} else { 
				printf("branch evt_bs_dydz_branch does not exist!\n");
				exit(1);
			}
			evt_bs_dydz_isLoaded = true;
		}
		return evt_bs_dydz_;
	}
	float &evt_bs_dydzErr()
	{
		if (not evt_bs_dydzErr_isLoaded) {
			if (evt_bs_dydzErr_branch != 0) {
				evt_bs_dydzErr_branch->GetEntry(index);
			} else { 
				printf("branch evt_bs_dydzErr_branch does not exist!\n");
				exit(1);
			}
			evt_bs_dydzErr_isLoaded = true;
		}
		return evt_bs_dydzErr_;
	}
	float &evt_bs_sigmaZ()
	{
		if (not evt_bs_sigmaZ_isLoaded) {
			if (evt_bs_sigmaZ_branch != 0) {
				evt_bs_sigmaZ_branch->GetEntry(index);
			} else { 
				printf("branch evt_bs_sigmaZ_branch does not exist!\n");
				exit(1);
			}
			evt_bs_sigmaZ_isLoaded = true;
		}
		return evt_bs_sigmaZ_;
	}
	float &evt_bs_sigmaZErr()
	{
		if (not evt_bs_sigmaZErr_isLoaded) {
			if (evt_bs_sigmaZErr_branch != 0) {
				evt_bs_sigmaZErr_branch->GetEntry(index);
			} else { 
				printf("branch evt_bs_sigmaZErr_branch does not exist!\n");
				exit(1);
			}
			evt_bs_sigmaZErr_isLoaded = true;
		}
		return evt_bs_sigmaZErr_;
	}
	float &evt_bs_xErr()
	{
		if (not evt_bs_xErr_isLoaded) {
			if (evt_bs_xErr_branch != 0) {
				evt_bs_xErr_branch->GetEntry(index);
			} else { 
				printf("branch evt_bs_xErr_branch does not exist!\n");
				exit(1);
			}
			evt_bs_xErr_isLoaded = true;
		}
		return evt_bs_xErr_;
	}
	float &evt_bs_yErr()
	{
		if (not evt_bs_yErr_isLoaded) {
			if (evt_bs_yErr_branch != 0) {
				evt_bs_yErr_branch->GetEntry(index);
			} else { 
				printf("branch evt_bs_yErr_branch does not exist!\n");
				exit(1);
			}
			evt_bs_yErr_isLoaded = true;
		}
		return evt_bs_yErr_;
	}
	float &evt_bs_zErr()
	{
		if (not evt_bs_zErr_isLoaded) {
			if (evt_bs_zErr_branch != 0) {
				evt_bs_zErr_branch->GetEntry(index);
			} else { 
				printf("branch evt_bs_zErr_branch does not exist!\n");
				exit(1);
			}
			evt_bs_zErr_isLoaded = true;
		}
		return evt_bs_zErr_;
	}
	float &evt_bField()
	{
		if (not evt_bField_isLoaded) {
			if (evt_bField_branch != 0) {
				evt_bField_branch->GetEntry(index);
			} else { 
				printf("branch evt_bField_branch does not exist!\n");
				exit(1);
			}
			evt_bField_isLoaded = true;
		}
		return evt_bField_;
	}
	float &evt_fixgrid_all_rho()
	{
		if (not evt_fixgrid_all_rho_isLoaded) {
			if (evt_fixgrid_all_rho_branch != 0) {
				evt_fixgrid_all_rho_branch->GetEntry(index);
			} else { 
				printf("branch evt_fixgrid_all_rho_branch does not exist!\n");
				exit(1);
			}
			evt_fixgrid_all_rho_isLoaded = true;
		}
		return evt_fixgrid_all_rho_;
	}
	float &evt_fixgridfastjet_allcalo_rho()
	{
		if (not evt_fixgridfastjet_allcalo_rho_isLoaded) {
			if (evt_fixgridfastjet_allcalo_rho_branch != 0) {
				evt_fixgridfastjet_allcalo_rho_branch->GetEntry(index);
			} else { 
				printf("branch evt_fixgridfastjet_allcalo_rho_branch does not exist!\n");
				exit(1);
			}
			evt_fixgridfastjet_allcalo_rho_isLoaded = true;
		}
		return evt_fixgridfastjet_allcalo_rho_;
	}
	float &evt_fixgridfastjet_all_rho()
	{
		if (not evt_fixgridfastjet_all_rho_isLoaded) {
			if (evt_fixgridfastjet_all_rho_branch != 0) {
				evt_fixgridfastjet_all_rho_branch->GetEntry(index);
			} else { 
				printf("branch evt_fixgridfastjet_all_rho_branch does not exist!\n");
				exit(1);
			}
			evt_fixgridfastjet_all_rho_isLoaded = true;
		}
		return evt_fixgridfastjet_all_rho_;
	}
	float &evt_fixgridfastjet_centralcalo_rho()
	{
		if (not evt_fixgridfastjet_centralcalo_rho_isLoaded) {
			if (evt_fixgridfastjet_centralcalo_rho_branch != 0) {
				evt_fixgridfastjet_centralcalo_rho_branch->GetEntry(index);
			} else { 
				printf("branch evt_fixgridfastjet_centralcalo_rho_branch does not exist!\n");
				exit(1);
			}
			evt_fixgridfastjet_centralcalo_rho_isLoaded = true;
		}
		return evt_fixgridfastjet_centralcalo_rho_;
	}
	float &evt_fixgridfastjet_centralchargedpileup_rho()
	{
		if (not evt_fixgridfastjet_centralchargedpileup_rho_isLoaded) {
			if (evt_fixgridfastjet_centralchargedpileup_rho_branch != 0) {
				evt_fixgridfastjet_centralchargedpileup_rho_branch->GetEntry(index);
			} else { 
				printf("branch evt_fixgridfastjet_centralchargedpileup_rho_branch does not exist!\n");
				exit(1);
			}
			evt_fixgridfastjet_centralchargedpileup_rho_isLoaded = true;
		}
		return evt_fixgridfastjet_centralchargedpileup_rho_;
	}
	float &evt_fixgridfastjet_centralneutral_rho()
	{
		if (not evt_fixgridfastjet_centralneutral_rho_isLoaded) {
			if (evt_fixgridfastjet_centralneutral_rho_branch != 0) {
				evt_fixgridfastjet_centralneutral_rho_branch->GetEntry(index);
			} else { 
				printf("branch evt_fixgridfastjet_centralneutral_rho_branch does not exist!\n");
				exit(1);
			}
			evt_fixgridfastjet_centralneutral_rho_isLoaded = true;
		}
		return evt_fixgridfastjet_centralneutral_rho_;
	}
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >  &evt_bsp4()
	{
		if (not evt_bsp4_isLoaded) {
			if (evt_bsp4_branch != 0) {
				evt_bsp4_branch->GetEntry(index);
			} else { 
				printf("branch evt_bsp4_branch does not exist!\n");
				exit(1);
			}
			evt_bsp4_isLoaded = true;
		}
		return evt_bsp4_;
	}
	const vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &els_mc_patMatch_p4()
	{
		if (not els_mc_patMatch_p4_isLoaded) {
			if (els_mc_patMatch_p4_branch != 0) {
				els_mc_patMatch_p4_branch->GetEntry(index);
			} else { 
				printf("branch els_mc_patMatch_p4_branch does not exist!\n");
				exit(1);
			}
			els_mc_patMatch_p4_isLoaded = true;
		}
		return els_mc_patMatch_p4_;
	}
	const vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &els_p4()
	{
		if (not els_p4_isLoaded) {
			if (els_p4_branch != 0) {
				els_p4_branch->GetEntry(index);
			} else { 
				printf("branch els_p4_branch does not exist!\n");
				exit(1);
			}
			els_p4_isLoaded = true;
		}
		return els_p4_;
	}
	const vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &els_p4In()
	{
		if (not els_p4In_isLoaded) {
			if (els_p4In_branch != 0) {
				els_p4In_branch->GetEntry(index);
			} else { 
				printf("branch els_p4In_branch does not exist!\n");
				exit(1);
			}
			els_p4In_isLoaded = true;
		}
		return els_p4In_;
	}
	const vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &els_p4Out()
	{
		if (not els_p4Out_isLoaded) {
			if (els_p4Out_branch != 0) {
				els_p4Out_branch->GetEntry(index);
			} else { 
				printf("branch els_p4Out_branch does not exist!\n");
				exit(1);
			}
			els_p4Out_isLoaded = true;
		}
		return els_p4Out_;
	}
	const vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &els_trk_p4()
	{
		if (not els_trk_p4_isLoaded) {
			if (els_trk_p4_branch != 0) {
				els_trk_p4_branch->GetEntry(index);
			} else { 
				printf("branch els_trk_p4_branch does not exist!\n");
				exit(1);
			}
			els_trk_p4_isLoaded = true;
		}
		return els_trk_p4_;
	}
	const vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &els_trk_vertex_p4()
	{
		if (not els_trk_vertex_p4_isLoaded) {
			if (els_trk_vertex_p4_branch != 0) {
				els_trk_vertex_p4_branch->GetEntry(index);
			} else { 
				printf("branch els_trk_vertex_p4_branch does not exist!\n");
				exit(1);
			}
			els_trk_vertex_p4_isLoaded = true;
		}
		return els_trk_vertex_p4_;
	}
	const vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &els_vertex_p4()
	{
		if (not els_vertex_p4_isLoaded) {
			if (els_vertex_p4_branch != 0) {
				els_vertex_p4_branch->GetEntry(index);
			} else { 
				printf("branch els_vertex_p4_branch does not exist!\n");
				exit(1);
			}
			els_vertex_p4_isLoaded = true;
		}
		return els_vertex_p4_;
	}
	const vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &photons_p4()
	{
		if (not photons_p4_isLoaded) {
			if (photons_p4_branch != 0) {
				photons_p4_branch->GetEntry(index);
			} else { 
				printf("branch photons_p4_branch does not exist!\n");
				exit(1);
			}
			photons_p4_isLoaded = true;
		}
		return photons_p4_;
	}
	const vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &vtxs_position()
	{
		if (not vtxs_position_isLoaded) {
			if (vtxs_position_branch != 0) {
				vtxs_position_branch->GetEntry(index);
			} else { 
				printf("branch vtxs_position_branch does not exist!\n");
				exit(1);
			}
			vtxs_position_isLoaded = true;
		}
		return vtxs_position_;
	}
	const vector<vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > > &hlt_trigObjs_p4()
	{
		if (not hlt_trigObjs_p4_isLoaded) {
			if (hlt_trigObjs_p4_branch != 0) {
				hlt_trigObjs_p4_branch->GetEntry(index);
			} else { 
				printf("branch hlt_trigObjs_p4_branch does not exist!\n");
				exit(1);
			}
			hlt_trigObjs_p4_isLoaded = true;
		}
		return hlt_trigObjs_p4_;
	}
	const vector<float> &evt_bs_covMatrix()
	{
		if (not evt_bs_covMatrix_isLoaded) {
			if (evt_bs_covMatrix_branch != 0) {
				evt_bs_covMatrix_branch->GetEntry(index);
			} else { 
				printf("branch evt_bs_covMatrix_branch does not exist!\n");
				exit(1);
			}
			evt_bs_covMatrix_isLoaded = true;
		}
		return evt_bs_covMatrix_;
	}
	const vector<float> &els_bs2d()
	{
		if (not els_bs2d_isLoaded) {
			if (els_bs2d_branch != 0) {
				els_bs2d_branch->GetEntry(index);
			} else { 
				printf("branch els_bs2d_branch does not exist!\n");
				exit(1);
			}
			els_bs2d_isLoaded = true;
		}
		return els_bs2d_;
	}
	const vector<float> &els_bs2derr()
	{
		if (not els_bs2derr_isLoaded) {
			if (els_bs2derr_branch != 0) {
				els_bs2derr_branch->GetEntry(index);
			} else { 
				printf("branch els_bs2derr_branch does not exist!\n");
				exit(1);
			}
			els_bs2derr_isLoaded = true;
		}
		return els_bs2derr_;
	}
	const vector<float> &els_bs3d()
	{
		if (not els_bs3d_isLoaded) {
			if (els_bs3d_branch != 0) {
				els_bs3d_branch->GetEntry(index);
			} else { 
				printf("branch els_bs3d_branch does not exist!\n");
				exit(1);
			}
			els_bs3d_isLoaded = true;
		}
		return els_bs3d_;
	}
	const vector<float> &els_bs3derr()
	{
		if (not els_bs3derr_isLoaded) {
			if (els_bs3derr_branch != 0) {
				els_bs3derr_branch->GetEntry(index);
			} else { 
				printf("branch els_bs3derr_branch does not exist!\n");
				exit(1);
			}
			els_bs3derr_isLoaded = true;
		}
		return els_bs3derr_;
	}
	const vector<float> &els_chi2()
	{
		if (not els_chi2_isLoaded) {
			if (els_chi2_branch != 0) {
				els_chi2_branch->GetEntry(index);
			} else { 
				printf("branch els_chi2_branch does not exist!\n");
				exit(1);
			}
			els_chi2_isLoaded = true;
		}
		return els_chi2_;
	}
	const vector<float> &els_ckf_chi2()
	{
		if (not els_ckf_chi2_isLoaded) {
			if (els_ckf_chi2_branch != 0) {
				els_ckf_chi2_branch->GetEntry(index);
			} else { 
				printf("branch els_ckf_chi2_branch does not exist!\n");
				exit(1);
			}
			els_ckf_chi2_isLoaded = true;
		}
		return els_ckf_chi2_;
	}
	const vector<float> &els_ckf_ndof()
	{
		if (not els_ckf_ndof_isLoaded) {
			if (els_ckf_ndof_branch != 0) {
				els_ckf_ndof_branch->GetEntry(index);
			} else { 
				printf("branch els_ckf_ndof_branch does not exist!\n");
				exit(1);
			}
			els_ckf_ndof_isLoaded = true;
		}
		return els_ckf_ndof_;
	}
	const vector<float> &els_d0()
	{
		if (not els_d0_isLoaded) {
			if (els_d0_branch != 0) {
				els_d0_branch->GetEntry(index);
			} else { 
				printf("branch els_d0_branch does not exist!\n");
				exit(1);
			}
			els_d0_isLoaded = true;
		}
		return els_d0_;
	}
	const vector<float> &els_d0Err()
	{
		if (not els_d0Err_isLoaded) {
			if (els_d0Err_branch != 0) {
				els_d0Err_branch->GetEntry(index);
			} else { 
				printf("branch els_d0Err_branch does not exist!\n");
				exit(1);
			}
			els_d0Err_isLoaded = true;
		}
		return els_d0Err_;
	}
	const vector<float> &els_d0corr()
	{
		if (not els_d0corr_isLoaded) {
			if (els_d0corr_branch != 0) {
				els_d0corr_branch->GetEntry(index);
			} else { 
				printf("branch els_d0corr_branch does not exist!\n");
				exit(1);
			}
			els_d0corr_isLoaded = true;
		}
		return els_d0corr_;
	}
	const vector<float> &els_d0corrPhi()
	{
		if (not els_d0corrPhi_isLoaded) {
			if (els_d0corrPhi_branch != 0) {
				els_d0corrPhi_branch->GetEntry(index);
			} else { 
				printf("branch els_d0corrPhi_branch does not exist!\n");
				exit(1);
			}
			els_d0corrPhi_isLoaded = true;
		}
		return els_d0corrPhi_;
	}
	const vector<float> &els_d0phiCov()
	{
		if (not els_d0phiCov_isLoaded) {
			if (els_d0phiCov_branch != 0) {
				els_d0phiCov_branch->GetEntry(index);
			} else { 
				printf("branch els_d0phiCov_branch does not exist!\n");
				exit(1);
			}
			els_d0phiCov_isLoaded = true;
		}
		return els_d0phiCov_;
	}
	const vector<float> &els_dEtaIn()
	{
		if (not els_dEtaIn_isLoaded) {
			if (els_dEtaIn_branch != 0) {
				els_dEtaIn_branch->GetEntry(index);
			} else { 
				printf("branch els_dEtaIn_branch does not exist!\n");
				exit(1);
			}
			els_dEtaIn_isLoaded = true;
		}
		return els_dEtaIn_;
	}
	const vector<float> &els_dEtaOut()
	{
		if (not els_dEtaOut_isLoaded) {
			if (els_dEtaOut_branch != 0) {
				els_dEtaOut_branch->GetEntry(index);
			} else { 
				printf("branch els_dEtaOut_branch does not exist!\n");
				exit(1);
			}
			els_dEtaOut_isLoaded = true;
		}
		return els_dEtaOut_;
	}
	const vector<float> &els_dPhiIn()
	{
		if (not els_dPhiIn_isLoaded) {
			if (els_dPhiIn_branch != 0) {
				els_dPhiIn_branch->GetEntry(index);
			} else { 
				printf("branch els_dPhiIn_branch does not exist!\n");
				exit(1);
			}
			els_dPhiIn_isLoaded = true;
		}
		return els_dPhiIn_;
	}
	const vector<float> &els_dPhiInPhiOut()
	{
		if (not els_dPhiInPhiOut_isLoaded) {
			if (els_dPhiInPhiOut_branch != 0) {
				els_dPhiInPhiOut_branch->GetEntry(index);
			} else { 
				printf("branch els_dPhiInPhiOut_branch does not exist!\n");
				exit(1);
			}
			els_dPhiInPhiOut_isLoaded = true;
		}
		return els_dPhiInPhiOut_;
	}
	const vector<float> &els_dPhiOut()
	{
		if (not els_dPhiOut_isLoaded) {
			if (els_dPhiOut_branch != 0) {
				els_dPhiOut_branch->GetEntry(index);
			} else { 
				printf("branch els_dPhiOut_branch does not exist!\n");
				exit(1);
			}
			els_dPhiOut_isLoaded = true;
		}
		return els_dPhiOut_;
	}
	const vector<float> &els_deltaEtaEleClusterTrackAtCalo()
	{
		if (not els_deltaEtaEleClusterTrackAtCalo_isLoaded) {
			if (els_deltaEtaEleClusterTrackAtCalo_branch != 0) {
				els_deltaEtaEleClusterTrackAtCalo_branch->GetEntry(index);
			} else { 
				printf("branch els_deltaEtaEleClusterTrackAtCalo_branch does not exist!\n");
				exit(1);
			}
			els_deltaEtaEleClusterTrackAtCalo_isLoaded = true;
		}
		return els_deltaEtaEleClusterTrackAtCalo_;
	}
	const vector<float> &els_deltaPhiEleClusterTrackAtCalo()
	{
		if (not els_deltaPhiEleClusterTrackAtCalo_isLoaded) {
			if (els_deltaPhiEleClusterTrackAtCalo_branch != 0) {
				els_deltaPhiEleClusterTrackAtCalo_branch->GetEntry(index);
			} else { 
				printf("branch els_deltaPhiEleClusterTrackAtCalo_branch does not exist!\n");
				exit(1);
			}
			els_deltaPhiEleClusterTrackAtCalo_isLoaded = true;
		}
		return els_deltaPhiEleClusterTrackAtCalo_;
	}
	const vector<float> &els_dxyPV()
	{
		if (not els_dxyPV_isLoaded) {
			if (els_dxyPV_branch != 0) {
				els_dxyPV_branch->GetEntry(index);
			} else { 
				printf("branch els_dxyPV_branch does not exist!\n");
				exit(1);
			}
			els_dxyPV_isLoaded = true;
		}
		return els_dxyPV_;
	}
	const vector<float> &els_dzPV()
	{
		if (not els_dzPV_isLoaded) {
			if (els_dzPV_branch != 0) {
				els_dzPV_branch->GetEntry(index);
			} else { 
				printf("branch els_dzPV_branch does not exist!\n");
				exit(1);
			}
			els_dzPV_isLoaded = true;
		}
		return els_dzPV_;
	}
	const vector<float> &els_e1x5()
	{
		if (not els_e1x5_isLoaded) {
			if (els_e1x5_branch != 0) {
				els_e1x5_branch->GetEntry(index);
			} else { 
				printf("branch els_e1x5_branch does not exist!\n");
				exit(1);
			}
			els_e1x5_isLoaded = true;
		}
		return els_e1x5_;
	}
	const vector<float> &els_e1x5_full5x5()
	{
		if (not els_e1x5_full5x5_isLoaded) {
			if (els_e1x5_full5x5_branch != 0) {
				els_e1x5_full5x5_branch->GetEntry(index);
			} else { 
				printf("branch els_e1x5_full5x5_branch does not exist!\n");
				exit(1);
			}
			els_e1x5_full5x5_isLoaded = true;
		}
		return els_e1x5_full5x5_;
	}
	const vector<float> &els_e2x5Max()
	{
		if (not els_e2x5Max_isLoaded) {
			if (els_e2x5Max_branch != 0) {
				els_e2x5Max_branch->GetEntry(index);
			} else { 
				printf("branch els_e2x5Max_branch does not exist!\n");
				exit(1);
			}
			els_e2x5Max_isLoaded = true;
		}
		return els_e2x5Max_;
	}
	const vector<float> &els_e2x5Max_full5x5()
	{
		if (not els_e2x5Max_full5x5_isLoaded) {
			if (els_e2x5Max_full5x5_branch != 0) {
				els_e2x5Max_full5x5_branch->GetEntry(index);
			} else { 
				printf("branch els_e2x5Max_full5x5_branch does not exist!\n");
				exit(1);
			}
			els_e2x5Max_full5x5_isLoaded = true;
		}
		return els_e2x5Max_full5x5_;
	}
	const vector<float> &els_e5x5()
	{
		if (not els_e5x5_isLoaded) {
			if (els_e5x5_branch != 0) {
				els_e5x5_branch->GetEntry(index);
			} else { 
				printf("branch els_e5x5_branch does not exist!\n");
				exit(1);
			}
			els_e5x5_isLoaded = true;
		}
		return els_e5x5_;
	}
	const vector<float> &els_e5x5_full5x5()
	{
		if (not els_e5x5_full5x5_isLoaded) {
			if (els_e5x5_full5x5_branch != 0) {
				els_e5x5_full5x5_branch->GetEntry(index);
			} else { 
				printf("branch els_e5x5_full5x5_branch does not exist!\n");
				exit(1);
			}
			els_e5x5_full5x5_isLoaded = true;
		}
		return els_e5x5_full5x5_;
	}
	const vector<float> &els_eOverPIn()
	{
		if (not els_eOverPIn_isLoaded) {
			if (els_eOverPIn_branch != 0) {
				els_eOverPIn_branch->GetEntry(index);
			} else { 
				printf("branch els_eOverPIn_branch does not exist!\n");
				exit(1);
			}
			els_eOverPIn_isLoaded = true;
		}
		return els_eOverPIn_;
	}
	const vector<float> &els_eOverPOut()
	{
		if (not els_eOverPOut_isLoaded) {
			if (els_eOverPOut_branch != 0) {
				els_eOverPOut_branch->GetEntry(index);
			} else { 
				printf("branch els_eOverPOut_branch does not exist!\n");
				exit(1);
			}
			els_eOverPOut_isLoaded = true;
		}
		return els_eOverPOut_;
	}
	const vector<float> &els_eSC()
	{
		if (not els_eSC_isLoaded) {
			if (els_eSC_branch != 0) {
				els_eSC_branch->GetEntry(index);
			} else { 
				printf("branch els_eSC_branch does not exist!\n");
				exit(1);
			}
			els_eSC_isLoaded = true;
		}
		return els_eSC_;
	}
	const vector<float> &els_eSCPresh()
	{
		if (not els_eSCPresh_isLoaded) {
			if (els_eSCPresh_branch != 0) {
				els_eSCPresh_branch->GetEntry(index);
			} else { 
				printf("branch els_eSCPresh_branch does not exist!\n");
				exit(1);
			}
			els_eSCPresh_isLoaded = true;
		}
		return els_eSCPresh_;
	}
	const vector<float> &els_eSCRaw()
	{
		if (not els_eSCRaw_isLoaded) {
			if (els_eSCRaw_branch != 0) {
				els_eSCRaw_branch->GetEntry(index);
			} else { 
				printf("branch els_eSCRaw_branch does not exist!\n");
				exit(1);
			}
			els_eSCRaw_isLoaded = true;
		}
		return els_eSCRaw_;
	}
	const vector<float> &els_eSeed()
	{
		if (not els_eSeed_isLoaded) {
			if (els_eSeed_branch != 0) {
				els_eSeed_branch->GetEntry(index);
			} else { 
				printf("branch els_eSeed_branch does not exist!\n");
				exit(1);
			}
			els_eSeed_isLoaded = true;
		}
		return els_eSeed_;
	}
	const vector<float> &els_eSeedOverPIn()
	{
		if (not els_eSeedOverPIn_isLoaded) {
			if (els_eSeedOverPIn_branch != 0) {
				els_eSeedOverPIn_branch->GetEntry(index);
			} else { 
				printf("branch els_eSeedOverPIn_branch does not exist!\n");
				exit(1);
			}
			els_eSeedOverPIn_isLoaded = true;
		}
		return els_eSeedOverPIn_;
	}
	const vector<float> &els_eSeedOverPOut()
	{
		if (not els_eSeedOverPOut_isLoaded) {
			if (els_eSeedOverPOut_branch != 0) {
				els_eSeedOverPOut_branch->GetEntry(index);
			} else { 
				printf("branch els_eSeedOverPOut_branch does not exist!\n");
				exit(1);
			}
			els_eSeedOverPOut_isLoaded = true;
		}
		return els_eSeedOverPOut_;
	}
	const vector<float> &els_ecalEnergy()
	{
		if (not els_ecalEnergy_isLoaded) {
			if (els_ecalEnergy_branch != 0) {
				els_ecalEnergy_branch->GetEntry(index);
			} else { 
				printf("branch els_ecalEnergy_branch does not exist!\n");
				exit(1);
			}
			els_ecalEnergy_isLoaded = true;
		}
		return els_ecalEnergy_;
	}
	const vector<float> &els_ecalEnergyError()
	{
		if (not els_ecalEnergyError_isLoaded) {
			if (els_ecalEnergyError_branch != 0) {
				els_ecalEnergyError_branch->GetEntry(index);
			} else { 
				printf("branch els_ecalEnergyError_branch does not exist!\n");
				exit(1);
			}
			els_ecalEnergyError_isLoaded = true;
		}
		return els_ecalEnergyError_;
	}
	const vector<float> &els_ecalIso()
	{
		if (not els_ecalIso_isLoaded) {
			if (els_ecalIso_branch != 0) {
				els_ecalIso_branch->GetEntry(index);
			} else { 
				printf("branch els_ecalIso_branch does not exist!\n");
				exit(1);
			}
			els_ecalIso_isLoaded = true;
		}
		return els_ecalIso_;
	}
	const vector<float> &els_ecalIso04()
	{
		if (not els_ecalIso04_isLoaded) {
			if (els_ecalIso04_branch != 0) {
				els_ecalIso04_branch->GetEntry(index);
			} else { 
				printf("branch els_ecalIso04_branch does not exist!\n");
				exit(1);
			}
			els_ecalIso04_isLoaded = true;
		}
		return els_ecalIso04_;
	}
	const vector<float> &els_ecalPFClusterIso()
	{
		if (not els_ecalPFClusterIso_isLoaded) {
			if (els_ecalPFClusterIso_branch != 0) {
				els_ecalPFClusterIso_branch->GetEntry(index);
			} else { 
				printf("branch els_ecalPFClusterIso_branch does not exist!\n");
				exit(1);
			}
			els_ecalPFClusterIso_isLoaded = true;
		}
		return els_ecalPFClusterIso_;
	}
	const vector<float> &els_etaErr()
	{
		if (not els_etaErr_isLoaded) {
			if (els_etaErr_branch != 0) {
				els_etaErr_branch->GetEntry(index);
			} else { 
				printf("branch els_etaErr_branch does not exist!\n");
				exit(1);
			}
			els_etaErr_isLoaded = true;
		}
		return els_etaErr_;
	}
	const vector<float> &els_etaSC()
	{
		if (not els_etaSC_isLoaded) {
			if (els_etaSC_branch != 0) {
				els_etaSC_branch->GetEntry(index);
			} else { 
				printf("branch els_etaSC_branch does not exist!\n");
				exit(1);
			}
			els_etaSC_isLoaded = true;
		}
		return els_etaSC_;
	}
	const vector<float> &els_etaSCwidth()
	{
		if (not els_etaSCwidth_isLoaded) {
			if (els_etaSCwidth_branch != 0) {
				els_etaSCwidth_branch->GetEntry(index);
			} else { 
				printf("branch els_etaSCwidth_branch does not exist!\n");
				exit(1);
			}
			els_etaSCwidth_isLoaded = true;
		}
		return els_etaSCwidth_;
	}
	const vector<float> &els_fbrem()
	{
		if (not els_fbrem_isLoaded) {
			if (els_fbrem_branch != 0) {
				els_fbrem_branch->GetEntry(index);
			} else { 
				printf("branch els_fbrem_branch does not exist!\n");
				exit(1);
			}
			els_fbrem_isLoaded = true;
		}
		return els_fbrem_;
	}
	const vector<float> &els_hOverE()
	{
		if (not els_hOverE_isLoaded) {
			if (els_hOverE_branch != 0) {
				els_hOverE_branch->GetEntry(index);
			} else { 
				printf("branch els_hOverE_branch does not exist!\n");
				exit(1);
			}
			els_hOverE_isLoaded = true;
		}
		return els_hOverE_;
	}
	const vector<float> &els_hOverEBC()
	{
		if (not els_hOverEBC_isLoaded) {
			if (els_hOverEBC_branch != 0) {
				els_hOverEBC_branch->GetEntry(index);
			} else { 
				printf("branch els_hOverEBC_branch does not exist!\n");
				exit(1);
			}
			els_hOverEBC_isLoaded = true;
		}
		return els_hOverEBC_;
	}
	const vector<float> &els_hcalDepth1OverEcal()
	{
		if (not els_hcalDepth1OverEcal_isLoaded) {
			if (els_hcalDepth1OverEcal_branch != 0) {
				els_hcalDepth1OverEcal_branch->GetEntry(index);
			} else { 
				printf("branch els_hcalDepth1OverEcal_branch does not exist!\n");
				exit(1);
			}
			els_hcalDepth1OverEcal_isLoaded = true;
		}
		return els_hcalDepth1OverEcal_;
	}
	const vector<float> &els_hcalDepth1TowerSumEt()
	{
		if (not els_hcalDepth1TowerSumEt_isLoaded) {
			if (els_hcalDepth1TowerSumEt_branch != 0) {
				els_hcalDepth1TowerSumEt_branch->GetEntry(index);
			} else { 
				printf("branch els_hcalDepth1TowerSumEt_branch does not exist!\n");
				exit(1);
			}
			els_hcalDepth1TowerSumEt_isLoaded = true;
		}
		return els_hcalDepth1TowerSumEt_;
	}
	const vector<float> &els_hcalDepth1TowerSumEt04()
	{
		if (not els_hcalDepth1TowerSumEt04_isLoaded) {
			if (els_hcalDepth1TowerSumEt04_branch != 0) {
				els_hcalDepth1TowerSumEt04_branch->GetEntry(index);
			} else { 
				printf("branch els_hcalDepth1TowerSumEt04_branch does not exist!\n");
				exit(1);
			}
			els_hcalDepth1TowerSumEt04_isLoaded = true;
		}
		return els_hcalDepth1TowerSumEt04_;
	}
	const vector<float> &els_hcalDepth2OverEcal()
	{
		if (not els_hcalDepth2OverEcal_isLoaded) {
			if (els_hcalDepth2OverEcal_branch != 0) {
				els_hcalDepth2OverEcal_branch->GetEntry(index);
			} else { 
				printf("branch els_hcalDepth2OverEcal_branch does not exist!\n");
				exit(1);
			}
			els_hcalDepth2OverEcal_isLoaded = true;
		}
		return els_hcalDepth2OverEcal_;
	}
	const vector<float> &els_hcalDepth2TowerSumEt()
	{
		if (not els_hcalDepth2TowerSumEt_isLoaded) {
			if (els_hcalDepth2TowerSumEt_branch != 0) {
				els_hcalDepth2TowerSumEt_branch->GetEntry(index);
			} else { 
				printf("branch els_hcalDepth2TowerSumEt_branch does not exist!\n");
				exit(1);
			}
			els_hcalDepth2TowerSumEt_isLoaded = true;
		}
		return els_hcalDepth2TowerSumEt_;
	}
	const vector<float> &els_hcalDepth2TowerSumEt04()
	{
		if (not els_hcalDepth2TowerSumEt04_isLoaded) {
			if (els_hcalDepth2TowerSumEt04_branch != 0) {
				els_hcalDepth2TowerSumEt04_branch->GetEntry(index);
			} else { 
				printf("branch els_hcalDepth2TowerSumEt04_branch does not exist!\n");
				exit(1);
			}
			els_hcalDepth2TowerSumEt04_isLoaded = true;
		}
		return els_hcalDepth2TowerSumEt04_;
	}
	const vector<float> &els_hcalIso()
	{
		if (not els_hcalIso_isLoaded) {
			if (els_hcalIso_branch != 0) {
				els_hcalIso_branch->GetEntry(index);
			} else { 
				printf("branch els_hcalIso_branch does not exist!\n");
				exit(1);
			}
			els_hcalIso_isLoaded = true;
		}
		return els_hcalIso_;
	}
	const vector<float> &els_hcalIso04()
	{
		if (not els_hcalIso04_isLoaded) {
			if (els_hcalIso04_branch != 0) {
				els_hcalIso04_branch->GetEntry(index);
			} else { 
				printf("branch els_hcalIso04_branch does not exist!\n");
				exit(1);
			}
			els_hcalIso04_isLoaded = true;
		}
		return els_hcalIso04_;
	}
	const vector<float> &els_hcalPFClusterIso()
	{
		if (not els_hcalPFClusterIso_isLoaded) {
			if (els_hcalPFClusterIso_branch != 0) {
				els_hcalPFClusterIso_branch->GetEntry(index);
			} else { 
				printf("branch els_hcalPFClusterIso_branch does not exist!\n");
				exit(1);
			}
			els_hcalPFClusterIso_isLoaded = true;
		}
		return els_hcalPFClusterIso_;
	}
	const vector<float> &els_ip2d()
	{
		if (not els_ip2d_isLoaded) {
			if (els_ip2d_branch != 0) {
				els_ip2d_branch->GetEntry(index);
			} else { 
				printf("branch els_ip2d_branch does not exist!\n");
				exit(1);
			}
			els_ip2d_isLoaded = true;
		}
		return els_ip2d_;
	}
	const vector<float> &els_ip2derr()
	{
		if (not els_ip2derr_isLoaded) {
			if (els_ip2derr_branch != 0) {
				els_ip2derr_branch->GetEntry(index);
			} else { 
				printf("branch els_ip2derr_branch does not exist!\n");
				exit(1);
			}
			els_ip2derr_isLoaded = true;
		}
		return els_ip2derr_;
	}
	const vector<float> &els_ip3d()
	{
		if (not els_ip3d_isLoaded) {
			if (els_ip3d_branch != 0) {
				els_ip3d_branch->GetEntry(index);
			} else { 
				printf("branch els_ip3d_branch does not exist!\n");
				exit(1);
			}
			els_ip3d_isLoaded = true;
		}
		return els_ip3d_;
	}
	const vector<float> &els_ip3derr()
	{
		if (not els_ip3derr_isLoaded) {
			if (els_ip3derr_branch != 0) {
				els_ip3derr_branch->GetEntry(index);
			} else { 
				printf("branch els_ip3derr_branch does not exist!\n");
				exit(1);
			}
			els_ip3derr_isLoaded = true;
		}
		return els_ip3derr_;
	}
	const vector<float> &els_mass()
	{
		if (not els_mass_isLoaded) {
			if (els_mass_branch != 0) {
				els_mass_branch->GetEntry(index);
			} else { 
				printf("branch els_mass_branch does not exist!\n");
				exit(1);
			}
			els_mass_isLoaded = true;
		}
		return els_mass_;
	}
	const vector<float> &els_mc_patMatch_dr()
	{
		if (not els_mc_patMatch_dr_isLoaded) {
			if (els_mc_patMatch_dr_branch != 0) {
				els_mc_patMatch_dr_branch->GetEntry(index);
			} else { 
				printf("branch els_mc_patMatch_dr_branch does not exist!\n");
				exit(1);
			}
			els_mc_patMatch_dr_isLoaded = true;
		}
		return els_mc_patMatch_dr_;
	}
	const vector<float> &els_miniIso_ch()
	{
		if (not els_miniIso_ch_isLoaded) {
			if (els_miniIso_ch_branch != 0) {
				els_miniIso_ch_branch->GetEntry(index);
			} else { 
				printf("branch els_miniIso_ch_branch does not exist!\n");
				exit(1);
			}
			els_miniIso_ch_isLoaded = true;
		}
		return els_miniIso_ch_;
	}
	const vector<float> &els_miniIso_db()
	{
		if (not els_miniIso_db_isLoaded) {
			if (els_miniIso_db_branch != 0) {
				els_miniIso_db_branch->GetEntry(index);
			} else { 
				printf("branch els_miniIso_db_branch does not exist!\n");
				exit(1);
			}
			els_miniIso_db_isLoaded = true;
		}
		return els_miniIso_db_;
	}
	const vector<float> &els_miniIso_em()
	{
		if (not els_miniIso_em_isLoaded) {
			if (els_miniIso_em_branch != 0) {
				els_miniIso_em_branch->GetEntry(index);
			} else { 
				printf("branch els_miniIso_em_branch does not exist!\n");
				exit(1);
			}
			els_miniIso_em_isLoaded = true;
		}
		return els_miniIso_em_;
	}
	const vector<float> &els_miniIso_nh()
	{
		if (not els_miniIso_nh_isLoaded) {
			if (els_miniIso_nh_branch != 0) {
				els_miniIso_nh_branch->GetEntry(index);
			} else { 
				printf("branch els_miniIso_nh_branch does not exist!\n");
				exit(1);
			}
			els_miniIso_nh_isLoaded = true;
		}
		return els_miniIso_nh_;
	}
	const vector<float> &els_miniIso_uncor()
	{
		if (not els_miniIso_uncor_isLoaded) {
			if (els_miniIso_uncor_branch != 0) {
				els_miniIso_uncor_branch->GetEntry(index);
			} else { 
				printf("branch els_miniIso_uncor_branch does not exist!\n");
				exit(1);
			}
			els_miniIso_uncor_isLoaded = true;
		}
		return els_miniIso_uncor_;
	}
	const vector<float> &els_mva()
	{
		if (not els_mva_isLoaded) {
			if (els_mva_branch != 0) {
				els_mva_branch->GetEntry(index);
			} else { 
				printf("branch els_mva_branch does not exist!\n");
				exit(1);
			}
			els_mva_isLoaded = true;
		}
		return els_mva_;
	}
	const vector<float> &els_ndof()
	{
		if (not els_ndof_isLoaded) {
			if (els_ndof_branch != 0) {
				els_ndof_branch->GetEntry(index);
			} else { 
				printf("branch els_ndof_branch does not exist!\n");
				exit(1);
			}
			els_ndof_isLoaded = true;
		}
		return els_ndof_;
	}
	const vector<float> &els_pfChargedHadronIso()
	{
		if (not els_pfChargedHadronIso_isLoaded) {
			if (els_pfChargedHadronIso_branch != 0) {
				els_pfChargedHadronIso_branch->GetEntry(index);
			} else { 
				printf("branch els_pfChargedHadronIso_branch does not exist!\n");
				exit(1);
			}
			els_pfChargedHadronIso_isLoaded = true;
		}
		return els_pfChargedHadronIso_;
	}
	const vector<float> &els_pfNeutralHadronIso()
	{
		if (not els_pfNeutralHadronIso_isLoaded) {
			if (els_pfNeutralHadronIso_branch != 0) {
				els_pfNeutralHadronIso_branch->GetEntry(index);
			} else { 
				printf("branch els_pfNeutralHadronIso_branch does not exist!\n");
				exit(1);
			}
			els_pfNeutralHadronIso_isLoaded = true;
		}
		return els_pfNeutralHadronIso_;
	}
	const vector<float> &els_pfPUIso()
	{
		if (not els_pfPUIso_isLoaded) {
			if (els_pfPUIso_branch != 0) {
				els_pfPUIso_branch->GetEntry(index);
			} else { 
				printf("branch els_pfPUIso_branch does not exist!\n");
				exit(1);
			}
			els_pfPUIso_isLoaded = true;
		}
		return els_pfPUIso_;
	}
	const vector<float> &els_pfPhotonIso()
	{
		if (not els_pfPhotonIso_isLoaded) {
			if (els_pfPhotonIso_branch != 0) {
				els_pfPhotonIso_branch->GetEntry(index);
			} else { 
				printf("branch els_pfPhotonIso_branch does not exist!\n");
				exit(1);
			}
			els_pfPhotonIso_isLoaded = true;
		}
		return els_pfPhotonIso_;
	}
	const vector<float> &els_phiErr()
	{
		if (not els_phiErr_isLoaded) {
			if (els_phiErr_branch != 0) {
				els_phiErr_branch->GetEntry(index);
			} else { 
				printf("branch els_phiErr_branch does not exist!\n");
				exit(1);
			}
			els_phiErr_isLoaded = true;
		}
		return els_phiErr_;
	}
	const vector<float> &els_phiSC()
	{
		if (not els_phiSC_isLoaded) {
			if (els_phiSC_branch != 0) {
				els_phiSC_branch->GetEntry(index);
			} else { 
				printf("branch els_phiSC_branch does not exist!\n");
				exit(1);
			}
			els_phiSC_isLoaded = true;
		}
		return els_phiSC_;
	}
	const vector<float> &els_phiSCwidth()
	{
		if (not els_phiSCwidth_isLoaded) {
			if (els_phiSCwidth_branch != 0) {
				els_phiSCwidth_branch->GetEntry(index);
			} else { 
				printf("branch els_phiSCwidth_branch does not exist!\n");
				exit(1);
			}
			els_phiSCwidth_isLoaded = true;
		}
		return els_phiSCwidth_;
	}
	const vector<float> &els_ptErr()
	{
		if (not els_ptErr_isLoaded) {
			if (els_ptErr_branch != 0) {
				els_ptErr_branch->GetEntry(index);
			} else { 
				printf("branch els_ptErr_branch does not exist!\n");
				exit(1);
			}
			els_ptErr_isLoaded = true;
		}
		return els_ptErr_;
	}
	const vector<float> &els_ptErrGsf()
	{
		if (not els_ptErrGsf_isLoaded) {
			if (els_ptErrGsf_branch != 0) {
				els_ptErrGsf_branch->GetEntry(index);
			} else { 
				printf("branch els_ptErrGsf_branch does not exist!\n");
				exit(1);
			}
			els_ptErrGsf_isLoaded = true;
		}
		return els_ptErrGsf_;
	}
	const vector<float> &els_r9()
	{
		if (not els_r9_isLoaded) {
			if (els_r9_branch != 0) {
				els_r9_branch->GetEntry(index);
			} else { 
				printf("branch els_r9_branch does not exist!\n");
				exit(1);
			}
			els_r9_isLoaded = true;
		}
		return els_r9_;
	}
	const vector<float> &els_r9_full5x5()
	{
		if (not els_r9_full5x5_isLoaded) {
			if (els_r9_full5x5_branch != 0) {
				els_r9_full5x5_branch->GetEntry(index);
			} else { 
				printf("branch els_r9_full5x5_branch does not exist!\n");
				exit(1);
			}
			els_r9_full5x5_isLoaded = true;
		}
		return els_r9_full5x5_;
	}
	const vector<float> &els_sigmaEtaEta()
	{
		if (not els_sigmaEtaEta_isLoaded) {
			if (els_sigmaEtaEta_branch != 0) {
				els_sigmaEtaEta_branch->GetEntry(index);
			} else { 
				printf("branch els_sigmaEtaEta_branch does not exist!\n");
				exit(1);
			}
			els_sigmaEtaEta_isLoaded = true;
		}
		return els_sigmaEtaEta_;
	}
	const vector<float> &els_sigmaEtaEta_full5x5()
	{
		if (not els_sigmaEtaEta_full5x5_isLoaded) {
			if (els_sigmaEtaEta_full5x5_branch != 0) {
				els_sigmaEtaEta_full5x5_branch->GetEntry(index);
			} else { 
				printf("branch els_sigmaEtaEta_full5x5_branch does not exist!\n");
				exit(1);
			}
			els_sigmaEtaEta_full5x5_isLoaded = true;
		}
		return els_sigmaEtaEta_full5x5_;
	}
	const vector<float> &els_sigmaIEtaIEta()
	{
		if (not els_sigmaIEtaIEta_isLoaded) {
			if (els_sigmaIEtaIEta_branch != 0) {
				els_sigmaIEtaIEta_branch->GetEntry(index);
			} else { 
				printf("branch els_sigmaIEtaIEta_branch does not exist!\n");
				exit(1);
			}
			els_sigmaIEtaIEta_isLoaded = true;
		}
		return els_sigmaIEtaIEta_;
	}
	const vector<float> &els_sigmaIEtaIEta_full5x5()
	{
		if (not els_sigmaIEtaIEta_full5x5_isLoaded) {
			if (els_sigmaIEtaIEta_full5x5_branch != 0) {
				els_sigmaIEtaIEta_full5x5_branch->GetEntry(index);
			} else { 
				printf("branch els_sigmaIEtaIEta_full5x5_branch does not exist!\n");
				exit(1);
			}
			els_sigmaIEtaIEta_full5x5_isLoaded = true;
		}
		return els_sigmaIEtaIEta_full5x5_;
	}
	const vector<float> &els_sigmaIPhiIPhi()
	{
		if (not els_sigmaIPhiIPhi_isLoaded) {
			if (els_sigmaIPhiIPhi_branch != 0) {
				els_sigmaIPhiIPhi_branch->GetEntry(index);
			} else { 
				printf("branch els_sigmaIPhiIPhi_branch does not exist!\n");
				exit(1);
			}
			els_sigmaIPhiIPhi_isLoaded = true;
		}
		return els_sigmaIPhiIPhi_;
	}
	const vector<float> &els_sigmaIPhiIPhi_full5x5()
	{
		if (not els_sigmaIPhiIPhi_full5x5_isLoaded) {
			if (els_sigmaIPhiIPhi_full5x5_branch != 0) {
				els_sigmaIPhiIPhi_full5x5_branch->GetEntry(index);
			} else { 
				printf("branch els_sigmaIPhiIPhi_full5x5_branch does not exist!\n");
				exit(1);
			}
			els_sigmaIPhiIPhi_full5x5_isLoaded = true;
		}
		return els_sigmaIPhiIPhi_full5x5_;
	}
	const vector<float> &els_sigmaIphiIphi()
	{
		if (not els_sigmaIphiIphi_isLoaded) {
			if (els_sigmaIphiIphi_branch != 0) {
				els_sigmaIphiIphi_branch->GetEntry(index);
			} else { 
				printf("branch els_sigmaIphiIphi_branch does not exist!\n");
				exit(1);
			}
			els_sigmaIphiIphi_isLoaded = true;
		}
		return els_sigmaIphiIphi_;
	}
	const vector<float> &els_tkIso()
	{
		if (not els_tkIso_isLoaded) {
			if (els_tkIso_branch != 0) {
				els_tkIso_branch->GetEntry(index);
			} else { 
				printf("branch els_tkIso_branch does not exist!\n");
				exit(1);
			}
			els_tkIso_isLoaded = true;
		}
		return els_tkIso_;
	}
	const vector<float> &els_tkIso04()
	{
		if (not els_tkIso04_isLoaded) {
			if (els_tkIso04_branch != 0) {
				els_tkIso04_branch->GetEntry(index);
			} else { 
				printf("branch els_tkIso04_branch does not exist!\n");
				exit(1);
			}
			els_tkIso04_isLoaded = true;
		}
		return els_tkIso04_;
	}
	const vector<float> &els_trackMomentumError()
	{
		if (not els_trackMomentumError_isLoaded) {
			if (els_trackMomentumError_branch != 0) {
				els_trackMomentumError_branch->GetEntry(index);
			} else { 
				printf("branch els_trackMomentumError_branch does not exist!\n");
				exit(1);
			}
			els_trackMomentumError_isLoaded = true;
		}
		return els_trackMomentumError_;
	}
	const vector<float> &els_trkdr()
	{
		if (not els_trkdr_isLoaded) {
			if (els_trkdr_branch != 0) {
				els_trkdr_branch->GetEntry(index);
			} else { 
				printf("branch els_trkdr_branch does not exist!\n");
				exit(1);
			}
			els_trkdr_isLoaded = true;
		}
		return els_trkdr_;
	}
	const vector<float> &els_trkshFrac()
	{
		if (not els_trkshFrac_isLoaded) {
			if (els_trkshFrac_branch != 0) {
				els_trkshFrac_branch->GetEntry(index);
			} else { 
				printf("branch els_trkshFrac_branch does not exist!\n");
				exit(1);
			}
			els_trkshFrac_isLoaded = true;
		}
		return els_trkshFrac_;
	}
	const vector<float> &els_z0()
	{
		if (not els_z0_isLoaded) {
			if (els_z0_branch != 0) {
				els_z0_branch->GetEntry(index);
			} else { 
				printf("branch els_z0_branch does not exist!\n");
				exit(1);
			}
			els_z0_isLoaded = true;
		}
		return els_z0_;
	}
	const vector<float> &els_z0Err()
	{
		if (not els_z0Err_isLoaded) {
			if (els_z0Err_branch != 0) {
				els_z0Err_branch->GetEntry(index);
			} else { 
				printf("branch els_z0Err_branch does not exist!\n");
				exit(1);
			}
			els_z0Err_isLoaded = true;
		}
		return els_z0Err_;
	}
	const vector<float> &els_z0corr()
	{
		if (not els_z0corr_isLoaded) {
			if (els_z0corr_branch != 0) {
				els_z0corr_branch->GetEntry(index);
			} else { 
				printf("branch els_z0corr_branch does not exist!\n");
				exit(1);
			}
			els_z0corr_isLoaded = true;
		}
		return els_z0corr_;
	}
	const vector<float> &photons_chargedHadronIso()
	{
		if (not photons_chargedHadronIso_isLoaded) {
			if (photons_chargedHadronIso_branch != 0) {
				photons_chargedHadronIso_branch->GetEntry(index);
			} else { 
				printf("branch photons_chargedHadronIso_branch does not exist!\n");
				exit(1);
			}
			photons_chargedHadronIso_isLoaded = true;
		}
		return photons_chargedHadronIso_;
	}
	const vector<float> &photons_e1x5()
	{
		if (not photons_e1x5_isLoaded) {
			if (photons_e1x5_branch != 0) {
				photons_e1x5_branch->GetEntry(index);
			} else { 
				printf("branch photons_e1x5_branch does not exist!\n");
				exit(1);
			}
			photons_e1x5_isLoaded = true;
		}
		return photons_e1x5_;
	}
	const vector<float> &photons_e2x5Max()
	{
		if (not photons_e2x5Max_isLoaded) {
			if (photons_e2x5Max_branch != 0) {
				photons_e2x5Max_branch->GetEntry(index);
			} else { 
				printf("branch photons_e2x5Max_branch does not exist!\n");
				exit(1);
			}
			photons_e2x5Max_isLoaded = true;
		}
		return photons_e2x5Max_;
	}
	const vector<float> &photons_e3x3()
	{
		if (not photons_e3x3_isLoaded) {
			if (photons_e3x3_branch != 0) {
				photons_e3x3_branch->GetEntry(index);
			} else { 
				printf("branch photons_e3x3_branch does not exist!\n");
				exit(1);
			}
			photons_e3x3_isLoaded = true;
		}
		return photons_e3x3_;
	}
	const vector<float> &photons_e5x5()
	{
		if (not photons_e5x5_isLoaded) {
			if (photons_e5x5_branch != 0) {
				photons_e5x5_branch->GetEntry(index);
			} else { 
				printf("branch photons_e5x5_branch does not exist!\n");
				exit(1);
			}
			photons_e5x5_isLoaded = true;
		}
		return photons_e5x5_;
	}
	const vector<float> &photons_eSC()
	{
		if (not photons_eSC_isLoaded) {
			if (photons_eSC_branch != 0) {
				photons_eSC_branch->GetEntry(index);
			} else { 
				printf("branch photons_eSC_branch does not exist!\n");
				exit(1);
			}
			photons_eSC_isLoaded = true;
		}
		return photons_eSC_;
	}
	const vector<float> &photons_eSCPresh()
	{
		if (not photons_eSCPresh_isLoaded) {
			if (photons_eSCPresh_branch != 0) {
				photons_eSCPresh_branch->GetEntry(index);
			} else { 
				printf("branch photons_eSCPresh_branch does not exist!\n");
				exit(1);
			}
			photons_eSCPresh_isLoaded = true;
		}
		return photons_eSCPresh_;
	}
	const vector<float> &photons_eSCRaw()
	{
		if (not photons_eSCRaw_isLoaded) {
			if (photons_eSCRaw_branch != 0) {
				photons_eSCRaw_branch->GetEntry(index);
			} else { 
				printf("branch photons_eSCRaw_branch does not exist!\n");
				exit(1);
			}
			photons_eSCRaw_isLoaded = true;
		}
		return photons_eSCRaw_;
	}
	const vector<float> &photons_ecalIso03()
	{
		if (not photons_ecalIso03_isLoaded) {
			if (photons_ecalIso03_branch != 0) {
				photons_ecalIso03_branch->GetEntry(index);
			} else { 
				printf("branch photons_ecalIso03_branch does not exist!\n");
				exit(1);
			}
			photons_ecalIso03_isLoaded = true;
		}
		return photons_ecalIso03_;
	}
	const vector<float> &photons_ecalIso04()
	{
		if (not photons_ecalIso04_isLoaded) {
			if (photons_ecalIso04_branch != 0) {
				photons_ecalIso04_branch->GetEntry(index);
			} else { 
				printf("branch photons_ecalIso04_branch does not exist!\n");
				exit(1);
			}
			photons_ecalIso04_isLoaded = true;
		}
		return photons_ecalIso04_;
	}
	const vector<float> &photons_ecalPFClusterIso()
	{
		if (not photons_ecalPFClusterIso_isLoaded) {
			if (photons_ecalPFClusterIso_branch != 0) {
				photons_ecalPFClusterIso_branch->GetEntry(index);
			} else { 
				printf("branch photons_ecalPFClusterIso_branch does not exist!\n");
				exit(1);
			}
			photons_ecalPFClusterIso_isLoaded = true;
		}
		return photons_ecalPFClusterIso_;
	}
	const vector<float> &photons_etaSC()
	{
		if (not photons_etaSC_isLoaded) {
			if (photons_etaSC_branch != 0) {
				photons_etaSC_branch->GetEntry(index);
			} else { 
				printf("branch photons_etaSC_branch does not exist!\n");
				exit(1);
			}
			photons_etaSC_isLoaded = true;
		}
		return photons_etaSC_;
	}
	const vector<float> &photons_full3x3_e3x3()
	{
		if (not photons_full3x3_e3x3_isLoaded) {
			if (photons_full3x3_e3x3_branch != 0) {
				photons_full3x3_e3x3_branch->GetEntry(index);
			} else { 
				printf("branch photons_full3x3_e3x3_branch does not exist!\n");
				exit(1);
			}
			photons_full3x3_e3x3_isLoaded = true;
		}
		return photons_full3x3_e3x3_;
	}
	const vector<float> &photons_full5x5_e1x5()
	{
		if (not photons_full5x5_e1x5_isLoaded) {
			if (photons_full5x5_e1x5_branch != 0) {
				photons_full5x5_e1x5_branch->GetEntry(index);
			} else { 
				printf("branch photons_full5x5_e1x5_branch does not exist!\n");
				exit(1);
			}
			photons_full5x5_e1x5_isLoaded = true;
		}
		return photons_full5x5_e1x5_;
	}
	const vector<float> &photons_full5x5_e2x5Max()
	{
		if (not photons_full5x5_e2x5Max_isLoaded) {
			if (photons_full5x5_e2x5Max_branch != 0) {
				photons_full5x5_e2x5Max_branch->GetEntry(index);
			} else { 
				printf("branch photons_full5x5_e2x5Max_branch does not exist!\n");
				exit(1);
			}
			photons_full5x5_e2x5Max_isLoaded = true;
		}
		return photons_full5x5_e2x5Max_;
	}
	const vector<float> &photons_full5x5_e5x5()
	{
		if (not photons_full5x5_e5x5_isLoaded) {
			if (photons_full5x5_e5x5_branch != 0) {
				photons_full5x5_e5x5_branch->GetEntry(index);
			} else { 
				printf("branch photons_full5x5_e5x5_branch does not exist!\n");
				exit(1);
			}
			photons_full5x5_e5x5_isLoaded = true;
		}
		return photons_full5x5_e5x5_;
	}
	const vector<float> &photons_full5x5_hOverE()
	{
		if (not photons_full5x5_hOverE_isLoaded) {
			if (photons_full5x5_hOverE_branch != 0) {
				photons_full5x5_hOverE_branch->GetEntry(index);
			} else { 
				printf("branch photons_full5x5_hOverE_branch does not exist!\n");
				exit(1);
			}
			photons_full5x5_hOverE_isLoaded = true;
		}
		return photons_full5x5_hOverE_;
	}
	const vector<float> &photons_full5x5_hOverEtowBC()
	{
		if (not photons_full5x5_hOverEtowBC_isLoaded) {
			if (photons_full5x5_hOverEtowBC_branch != 0) {
				photons_full5x5_hOverEtowBC_branch->GetEntry(index);
			} else { 
				printf("branch photons_full5x5_hOverEtowBC_branch does not exist!\n");
				exit(1);
			}
			photons_full5x5_hOverEtowBC_isLoaded = true;
		}
		return photons_full5x5_hOverEtowBC_;
	}
	const vector<float> &photons_full5x5_r9()
	{
		if (not photons_full5x5_r9_isLoaded) {
			if (photons_full5x5_r9_branch != 0) {
				photons_full5x5_r9_branch->GetEntry(index);
			} else { 
				printf("branch photons_full5x5_r9_branch does not exist!\n");
				exit(1);
			}
			photons_full5x5_r9_isLoaded = true;
		}
		return photons_full5x5_r9_;
	}
	const vector<float> &photons_full5x5_sigmaEtaEta()
	{
		if (not photons_full5x5_sigmaEtaEta_isLoaded) {
			if (photons_full5x5_sigmaEtaEta_branch != 0) {
				photons_full5x5_sigmaEtaEta_branch->GetEntry(index);
			} else { 
				printf("branch photons_full5x5_sigmaEtaEta_branch does not exist!\n");
				exit(1);
			}
			photons_full5x5_sigmaEtaEta_isLoaded = true;
		}
		return photons_full5x5_sigmaEtaEta_;
	}
	const vector<float> &photons_full5x5_sigmaIEtaIEta()
	{
		if (not photons_full5x5_sigmaIEtaIEta_isLoaded) {
			if (photons_full5x5_sigmaIEtaIEta_branch != 0) {
				photons_full5x5_sigmaIEtaIEta_branch->GetEntry(index);
			} else { 
				printf("branch photons_full5x5_sigmaIEtaIEta_branch does not exist!\n");
				exit(1);
			}
			photons_full5x5_sigmaIEtaIEta_isLoaded = true;
		}
		return photons_full5x5_sigmaIEtaIEta_;
	}
	const vector<float> &photons_hOverE()
	{
		if (not photons_hOverE_isLoaded) {
			if (photons_hOverE_branch != 0) {
				photons_hOverE_branch->GetEntry(index);
			} else { 
				printf("branch photons_hOverE_branch does not exist!\n");
				exit(1);
			}
			photons_hOverE_isLoaded = true;
		}
		return photons_hOverE_;
	}
	const vector<float> &photons_hOverEtowBC()
	{
		if (not photons_hOverEtowBC_isLoaded) {
			if (photons_hOverEtowBC_branch != 0) {
				photons_hOverEtowBC_branch->GetEntry(index);
			} else { 
				printf("branch photons_hOverEtowBC_branch does not exist!\n");
				exit(1);
			}
			photons_hOverEtowBC_isLoaded = true;
		}
		return photons_hOverEtowBC_;
	}
	const vector<float> &photons_hcalDepth1TowerSumEtBcConeDR03()
	{
		if (not photons_hcalDepth1TowerSumEtBcConeDR03_isLoaded) {
			if (photons_hcalDepth1TowerSumEtBcConeDR03_branch != 0) {
				photons_hcalDepth1TowerSumEtBcConeDR03_branch->GetEntry(index);
			} else { 
				printf("branch photons_hcalDepth1TowerSumEtBcConeDR03_branch does not exist!\n");
				exit(1);
			}
			photons_hcalDepth1TowerSumEtBcConeDR03_isLoaded = true;
		}
		return photons_hcalDepth1TowerSumEtBcConeDR03_;
	}
	const vector<float> &photons_hcalDepth1TowerSumEtBcConeDR04()
	{
		if (not photons_hcalDepth1TowerSumEtBcConeDR04_isLoaded) {
			if (photons_hcalDepth1TowerSumEtBcConeDR04_branch != 0) {
				photons_hcalDepth1TowerSumEtBcConeDR04_branch->GetEntry(index);
			} else { 
				printf("branch photons_hcalDepth1TowerSumEtBcConeDR04_branch does not exist!\n");
				exit(1);
			}
			photons_hcalDepth1TowerSumEtBcConeDR04_isLoaded = true;
		}
		return photons_hcalDepth1TowerSumEtBcConeDR04_;
	}
	const vector<float> &photons_hcalDepth2TowerSumEtBcConeDR03()
	{
		if (not photons_hcalDepth2TowerSumEtBcConeDR03_isLoaded) {
			if (photons_hcalDepth2TowerSumEtBcConeDR03_branch != 0) {
				photons_hcalDepth2TowerSumEtBcConeDR03_branch->GetEntry(index);
			} else { 
				printf("branch photons_hcalDepth2TowerSumEtBcConeDR03_branch does not exist!\n");
				exit(1);
			}
			photons_hcalDepth2TowerSumEtBcConeDR03_isLoaded = true;
		}
		return photons_hcalDepth2TowerSumEtBcConeDR03_;
	}
	const vector<float> &photons_hcalDepth2TowerSumEtBcConeDR04()
	{
		if (not photons_hcalDepth2TowerSumEtBcConeDR04_isLoaded) {
			if (photons_hcalDepth2TowerSumEtBcConeDR04_branch != 0) {
				photons_hcalDepth2TowerSumEtBcConeDR04_branch->GetEntry(index);
			} else { 
				printf("branch photons_hcalDepth2TowerSumEtBcConeDR04_branch does not exist!\n");
				exit(1);
			}
			photons_hcalDepth2TowerSumEtBcConeDR04_isLoaded = true;
		}
		return photons_hcalDepth2TowerSumEtBcConeDR04_;
	}
	const vector<float> &photons_hcalIso03()
	{
		if (not photons_hcalIso03_isLoaded) {
			if (photons_hcalIso03_branch != 0) {
				photons_hcalIso03_branch->GetEntry(index);
			} else { 
				printf("branch photons_hcalIso03_branch does not exist!\n");
				exit(1);
			}
			photons_hcalIso03_isLoaded = true;
		}
		return photons_hcalIso03_;
	}
	const vector<float> &photons_hcalIso04()
	{
		if (not photons_hcalIso04_isLoaded) {
			if (photons_hcalIso04_branch != 0) {
				photons_hcalIso04_branch->GetEntry(index);
			} else { 
				printf("branch photons_hcalIso04_branch does not exist!\n");
				exit(1);
			}
			photons_hcalIso04_isLoaded = true;
		}
		return photons_hcalIso04_;
	}
	const vector<float> &photons_hcalPFClusterIso()
	{
		if (not photons_hcalPFClusterIso_isLoaded) {
			if (photons_hcalPFClusterIso_branch != 0) {
				photons_hcalPFClusterIso_branch->GetEntry(index);
			} else { 
				printf("branch photons_hcalPFClusterIso_branch does not exist!\n");
				exit(1);
			}
			photons_hcalPFClusterIso_isLoaded = true;
		}
		return photons_hcalPFClusterIso_;
	}
	const vector<float> &photons_hcalTowerSumEtBcConeDR03()
	{
		if (not photons_hcalTowerSumEtBcConeDR03_isLoaded) {
			if (photons_hcalTowerSumEtBcConeDR03_branch != 0) {
				photons_hcalTowerSumEtBcConeDR03_branch->GetEntry(index);
			} else { 
				printf("branch photons_hcalTowerSumEtBcConeDR03_branch does not exist!\n");
				exit(1);
			}
			photons_hcalTowerSumEtBcConeDR03_isLoaded = true;
		}
		return photons_hcalTowerSumEtBcConeDR03_;
	}
	const vector<float> &photons_hcalTowerSumEtBcConeDR04()
	{
		if (not photons_hcalTowerSumEtBcConeDR04_isLoaded) {
			if (photons_hcalTowerSumEtBcConeDR04_branch != 0) {
				photons_hcalTowerSumEtBcConeDR04_branch->GetEntry(index);
			} else { 
				printf("branch photons_hcalTowerSumEtBcConeDR04_branch does not exist!\n");
				exit(1);
			}
			photons_hcalTowerSumEtBcConeDR04_isLoaded = true;
		}
		return photons_hcalTowerSumEtBcConeDR04_;
	}
	const vector<float> &photons_mass()
	{
		if (not photons_mass_isLoaded) {
			if (photons_mass_branch != 0) {
				photons_mass_branch->GetEntry(index);
			} else { 
				printf("branch photons_mass_branch does not exist!\n");
				exit(1);
			}
			photons_mass_isLoaded = true;
		}
		return photons_mass_;
	}
	const vector<float> &photons_neutralHadronIso()
	{
		if (not photons_neutralHadronIso_isLoaded) {
			if (photons_neutralHadronIso_branch != 0) {
				photons_neutralHadronIso_branch->GetEntry(index);
			} else { 
				printf("branch photons_neutralHadronIso_branch does not exist!\n");
				exit(1);
			}
			photons_neutralHadronIso_isLoaded = true;
		}
		return photons_neutralHadronIso_;
	}
	const vector<float> &photons_ntkIsoHollow03()
	{
		if (not photons_ntkIsoHollow03_isLoaded) {
			if (photons_ntkIsoHollow03_branch != 0) {
				photons_ntkIsoHollow03_branch->GetEntry(index);
			} else { 
				printf("branch photons_ntkIsoHollow03_branch does not exist!\n");
				exit(1);
			}
			photons_ntkIsoHollow03_isLoaded = true;
		}
		return photons_ntkIsoHollow03_;
	}
	const vector<float> &photons_ntkIsoHollow04()
	{
		if (not photons_ntkIsoHollow04_isLoaded) {
			if (photons_ntkIsoHollow04_branch != 0) {
				photons_ntkIsoHollow04_branch->GetEntry(index);
			} else { 
				printf("branch photons_ntkIsoHollow04_branch does not exist!\n");
				exit(1);
			}
			photons_ntkIsoHollow04_isLoaded = true;
		}
		return photons_ntkIsoHollow04_;
	}
	const vector<float> &photons_ntkIsoSolid03()
	{
		if (not photons_ntkIsoSolid03_isLoaded) {
			if (photons_ntkIsoSolid03_branch != 0) {
				photons_ntkIsoSolid03_branch->GetEntry(index);
			} else { 
				printf("branch photons_ntkIsoSolid03_branch does not exist!\n");
				exit(1);
			}
			photons_ntkIsoSolid03_isLoaded = true;
		}
		return photons_ntkIsoSolid03_;
	}
	const vector<float> &photons_ntkIsoSolid04()
	{
		if (not photons_ntkIsoSolid04_isLoaded) {
			if (photons_ntkIsoSolid04_branch != 0) {
				photons_ntkIsoSolid04_branch->GetEntry(index);
			} else { 
				printf("branch photons_ntkIsoSolid04_branch does not exist!\n");
				exit(1);
			}
			photons_ntkIsoSolid04_isLoaded = true;
		}
		return photons_ntkIsoSolid04_;
	}
	const vector<float> &photons_phiSC()
	{
		if (not photons_phiSC_isLoaded) {
			if (photons_phiSC_branch != 0) {
				photons_phiSC_branch->GetEntry(index);
			} else { 
				printf("branch photons_phiSC_branch does not exist!\n");
				exit(1);
			}
			photons_phiSC_isLoaded = true;
		}
		return photons_phiSC_;
	}
	const vector<float> &photons_photonIso()
	{
		if (not photons_photonIso_isLoaded) {
			if (photons_photonIso_branch != 0) {
				photons_photonIso_branch->GetEntry(index);
			} else { 
				printf("branch photons_photonIso_branch does not exist!\n");
				exit(1);
			}
			photons_photonIso_isLoaded = true;
		}
		return photons_photonIso_;
	}
	const vector<float> &photons_recoChargedHadronIso()
	{
		if (not photons_recoChargedHadronIso_isLoaded) {
			if (photons_recoChargedHadronIso_branch != 0) {
				photons_recoChargedHadronIso_branch->GetEntry(index);
			} else { 
				printf("branch photons_recoChargedHadronIso_branch does not exist!\n");
				exit(1);
			}
			photons_recoChargedHadronIso_isLoaded = true;
		}
		return photons_recoChargedHadronIso_;
	}
	const vector<float> &photons_recoNeutralHadronIso()
	{
		if (not photons_recoNeutralHadronIso_isLoaded) {
			if (photons_recoNeutralHadronIso_branch != 0) {
				photons_recoNeutralHadronIso_branch->GetEntry(index);
			} else { 
				printf("branch photons_recoNeutralHadronIso_branch does not exist!\n");
				exit(1);
			}
			photons_recoNeutralHadronIso_isLoaded = true;
		}
		return photons_recoNeutralHadronIso_;
	}
	const vector<float> &photons_recoPhotonIso()
	{
		if (not photons_recoPhotonIso_isLoaded) {
			if (photons_recoPhotonIso_branch != 0) {
				photons_recoPhotonIso_branch->GetEntry(index);
			} else { 
				printf("branch photons_recoPhotonIso_branch does not exist!\n");
				exit(1);
			}
			photons_recoPhotonIso_isLoaded = true;
		}
		return photons_recoPhotonIso_;
	}
	const vector<float> &photons_sigmaEtaEta()
	{
		if (not photons_sigmaEtaEta_isLoaded) {
			if (photons_sigmaEtaEta_branch != 0) {
				photons_sigmaEtaEta_branch->GetEntry(index);
			} else { 
				printf("branch photons_sigmaEtaEta_branch does not exist!\n");
				exit(1);
			}
			photons_sigmaEtaEta_isLoaded = true;
		}
		return photons_sigmaEtaEta_;
	}
	const vector<float> &photons_sigmaIEtaIEta()
	{
		if (not photons_sigmaIEtaIEta_isLoaded) {
			if (photons_sigmaIEtaIEta_branch != 0) {
				photons_sigmaIEtaIEta_branch->GetEntry(index);
			} else { 
				printf("branch photons_sigmaIEtaIEta_branch does not exist!\n");
				exit(1);
			}
			photons_sigmaIEtaIEta_isLoaded = true;
		}
		return photons_sigmaIEtaIEta_;
	}
	const vector<float> &photons_tkIsoHollow03()
	{
		if (not photons_tkIsoHollow03_isLoaded) {
			if (photons_tkIsoHollow03_branch != 0) {
				photons_tkIsoHollow03_branch->GetEntry(index);
			} else { 
				printf("branch photons_tkIsoHollow03_branch does not exist!\n");
				exit(1);
			}
			photons_tkIsoHollow03_isLoaded = true;
		}
		return photons_tkIsoHollow03_;
	}
	const vector<float> &photons_tkIsoHollow04()
	{
		if (not photons_tkIsoHollow04_isLoaded) {
			if (photons_tkIsoHollow04_branch != 0) {
				photons_tkIsoHollow04_branch->GetEntry(index);
			} else { 
				printf("branch photons_tkIsoHollow04_branch does not exist!\n");
				exit(1);
			}
			photons_tkIsoHollow04_isLoaded = true;
		}
		return photons_tkIsoHollow04_;
	}
	const vector<float> &photons_tkIsoSolid03()
	{
		if (not photons_tkIsoSolid03_isLoaded) {
			if (photons_tkIsoSolid03_branch != 0) {
				photons_tkIsoSolid03_branch->GetEntry(index);
			} else { 
				printf("branch photons_tkIsoSolid03_branch does not exist!\n");
				exit(1);
			}
			photons_tkIsoSolid03_isLoaded = true;
		}
		return photons_tkIsoSolid03_;
	}
	const vector<float> &photons_tkIsoSolid04()
	{
		if (not photons_tkIsoSolid04_isLoaded) {
			if (photons_tkIsoSolid04_branch != 0) {
				photons_tkIsoSolid04_branch->GetEntry(index);
			} else { 
				printf("branch photons_tkIsoSolid04_branch does not exist!\n");
				exit(1);
			}
			photons_tkIsoSolid04_isLoaded = true;
		}
		return photons_tkIsoSolid04_;
	}
	const vector<float> &puInfo_trueNumInteractions()
	{
		if (not puInfo_trueNumInteractions_isLoaded) {
			if (puInfo_trueNumInteractions_branch != 0) {
				puInfo_trueNumInteractions_branch->GetEntry(index);
			} else { 
				printf("branch puInfo_trueNumInteractions_branch does not exist!\n");
				exit(1);
			}
			puInfo_trueNumInteractions_isLoaded = true;
		}
		return puInfo_trueNumInteractions_;
	}
	const vector<float> &vtxs_chi2()
	{
		if (not vtxs_chi2_isLoaded) {
			if (vtxs_chi2_branch != 0) {
				vtxs_chi2_branch->GetEntry(index);
			} else { 
				printf("branch vtxs_chi2_branch does not exist!\n");
				exit(1);
			}
			vtxs_chi2_isLoaded = true;
		}
		return vtxs_chi2_;
	}
	const vector<float> &vtxs_ndof()
	{
		if (not vtxs_ndof_isLoaded) {
			if (vtxs_ndof_branch != 0) {
				vtxs_ndof_branch->GetEntry(index);
			} else { 
				printf("branch vtxs_ndof_branch does not exist!\n");
				exit(1);
			}
			vtxs_ndof_isLoaded = true;
		}
		return vtxs_ndof_;
	}
	const vector<float> &vtxs_score()
	{
		if (not vtxs_score_isLoaded) {
			if (vtxs_score_branch != 0) {
				vtxs_score_branch->GetEntry(index);
			} else { 
				printf("branch vtxs_score_branch does not exist!\n");
				exit(1);
			}
			vtxs_score_isLoaded = true;
		}
		return vtxs_score_;
	}
	const vector<float> &vtxs_xError()
	{
		if (not vtxs_xError_isLoaded) {
			if (vtxs_xError_branch != 0) {
				vtxs_xError_branch->GetEntry(index);
			} else { 
				printf("branch vtxs_xError_branch does not exist!\n");
				exit(1);
			}
			vtxs_xError_isLoaded = true;
		}
		return vtxs_xError_;
	}
	const vector<float> &vtxs_yError()
	{
		if (not vtxs_yError_isLoaded) {
			if (vtxs_yError_branch != 0) {
				vtxs_yError_branch->GetEntry(index);
			} else { 
				printf("branch vtxs_yError_branch does not exist!\n");
				exit(1);
			}
			vtxs_yError_isLoaded = true;
		}
		return vtxs_yError_;
	}
	const vector<float> &vtxs_zError()
	{
		if (not vtxs_zError_isLoaded) {
			if (vtxs_zError_branch != 0) {
				vtxs_zError_branch->GetEntry(index);
			} else { 
				printf("branch vtxs_zError_branch does not exist!\n");
				exit(1);
			}
			vtxs_zError_isLoaded = true;
		}
		return vtxs_zError_;
	}
	const vector<vector<float> > &puInfo_instLumi()
	{
		if (not puInfo_instLumi_isLoaded) {
			if (puInfo_instLumi_branch != 0) {
				puInfo_instLumi_branch->GetEntry(index);
			} else { 
				printf("branch puInfo_instLumi_branch does not exist!\n");
				exit(1);
			}
			puInfo_instLumi_isLoaded = true;
		}
		return puInfo_instLumi_;
	}
	const vector<vector<float> > &vtxs_covMatrix()
	{
		if (not vtxs_covMatrix_isLoaded) {
			if (vtxs_covMatrix_branch != 0) {
				vtxs_covMatrix_branch->GetEntry(index);
			} else { 
				printf("branch vtxs_covMatrix_branch does not exist!\n");
				exit(1);
			}
			vtxs_covMatrix_isLoaded = true;
		}
		return vtxs_covMatrix_;
	}
	int &evt_bsType()
	{
		if (not evt_bsType_isLoaded) {
			if (evt_bsType_branch != 0) {
				evt_bsType_branch->GetEntry(index);
			} else { 
				printf("branch evt_bsType_branch does not exist!\n");
				exit(1);
			}
			evt_bsType_isLoaded = true;
		}
		return evt_bsType_;
	}
	int &evt_bunchCrossing()
	{
		if (not evt_bunchCrossing_isLoaded) {
			if (evt_bunchCrossing_branch != 0) {
				evt_bunchCrossing_branch->GetEntry(index);
			} else { 
				printf("branch evt_bunchCrossing_branch does not exist!\n");
				exit(1);
			}
			evt_bunchCrossing_isLoaded = true;
		}
		return evt_bunchCrossing_;
	}
	int &evt_experimentType()
	{
		if (not evt_experimentType_isLoaded) {
			if (evt_experimentType_branch != 0) {
				evt_experimentType_branch->GetEntry(index);
			} else { 
				printf("branch evt_experimentType_branch does not exist!\n");
				exit(1);
			}
			evt_experimentType_isLoaded = true;
		}
		return evt_experimentType_;
	}
	int &evt_isRealData()
	{
		if (not evt_isRealData_isLoaded) {
			if (evt_isRealData_branch != 0) {
				evt_isRealData_branch->GetEntry(index);
			} else { 
				printf("branch evt_isRealData_branch does not exist!\n");
				exit(1);
			}
			evt_isRealData_isLoaded = true;
		}
		return evt_isRealData_;
	}
	int &evt_orbitNumber()
	{
		if (not evt_orbitNumber_isLoaded) {
			if (evt_orbitNumber_branch != 0) {
				evt_orbitNumber_branch->GetEntry(index);
			} else { 
				printf("branch evt_orbitNumber_branch does not exist!\n");
				exit(1);
			}
			evt_orbitNumber_isLoaded = true;
		}
		return evt_orbitNumber_;
	}
	int &evt_storeNumber()
	{
		if (not evt_storeNumber_isLoaded) {
			if (evt_storeNumber_branch != 0) {
				evt_storeNumber_branch->GetEntry(index);
			} else { 
				printf("branch evt_storeNumber_branch does not exist!\n");
				exit(1);
			}
			evt_storeNumber_isLoaded = true;
		}
		return evt_storeNumber_;
	}
	const vector<int> &els_category()
	{
		if (not els_category_isLoaded) {
			if (els_category_branch != 0) {
				els_category_branch->GetEntry(index);
			} else { 
				printf("branch els_category_branch does not exist!\n");
				exit(1);
			}
			els_category_isLoaded = true;
		}
		return els_category_;
	}
	const vector<int> &els_charge()
	{
		if (not els_charge_isLoaded) {
			if (els_charge_branch != 0) {
				els_charge_branch->GetEntry(index);
			} else { 
				printf("branch els_charge_branch does not exist!\n");
				exit(1);
			}
			els_charge_isLoaded = true;
		}
		return els_charge_;
	}
	const vector<int> &els_ckf_charge()
	{
		if (not els_ckf_charge_isLoaded) {
			if (els_ckf_charge_branch != 0) {
				els_ckf_charge_branch->GetEntry(index);
			} else { 
				printf("branch els_ckf_charge_branch does not exist!\n");
				exit(1);
			}
			els_ckf_charge_isLoaded = true;
		}
		return els_ckf_charge_;
	}
	const vector<int> &els_ckf_laywithmeas()
	{
		if (not els_ckf_laywithmeas_isLoaded) {
			if (els_ckf_laywithmeas_branch != 0) {
				els_ckf_laywithmeas_branch->GetEntry(index);
			} else { 
				printf("branch els_ckf_laywithmeas_branch does not exist!\n");
				exit(1);
			}
			els_ckf_laywithmeas_isLoaded = true;
		}
		return els_ckf_laywithmeas_;
	}
	const vector<int> &els_class()
	{
		if (not els_class_isLoaded) {
			if (els_class_branch != 0) {
				els_class_branch->GetEntry(index);
			} else { 
				printf("branch els_class_branch does not exist!\n");
				exit(1);
			}
			els_class_isLoaded = true;
		}
		return els_class_;
	}
	const vector<int> &els_exp_innerlayers()
	{
		if (not els_exp_innerlayers_isLoaded) {
			if (els_exp_innerlayers_branch != 0) {
				els_exp_innerlayers_branch->GetEntry(index);
			} else { 
				printf("branch els_exp_innerlayers_branch does not exist!\n");
				exit(1);
			}
			els_exp_innerlayers_isLoaded = true;
		}
		return els_exp_innerlayers_;
	}
	const vector<int> &els_exp_outerlayers()
	{
		if (not els_exp_outerlayers_isLoaded) {
			if (els_exp_outerlayers_branch != 0) {
				els_exp_outerlayers_branch->GetEntry(index);
			} else { 
				printf("branch els_exp_outerlayers_branch does not exist!\n");
				exit(1);
			}
			els_exp_outerlayers_isLoaded = true;
		}
		return els_exp_outerlayers_;
	}
	const vector<int> &els_fiduciality()
	{
		if (not els_fiduciality_isLoaded) {
			if (els_fiduciality_branch != 0) {
				els_fiduciality_branch->GetEntry(index);
			} else { 
				printf("branch els_fiduciality_branch does not exist!\n");
				exit(1);
			}
			els_fiduciality_isLoaded = true;
		}
		return els_fiduciality_;
	}
	const vector<int> &els_lostHits()
	{
		if (not els_lostHits_isLoaded) {
			if (els_lostHits_branch != 0) {
				els_lostHits_branch->GetEntry(index);
			} else { 
				printf("branch els_lostHits_branch does not exist!\n");
				exit(1);
			}
			els_lostHits_isLoaded = true;
		}
		return els_lostHits_;
	}
	const vector<int> &els_lost_pixelhits()
	{
		if (not els_lost_pixelhits_isLoaded) {
			if (els_lost_pixelhits_branch != 0) {
				els_lost_pixelhits_branch->GetEntry(index);
			} else { 
				printf("branch els_lost_pixelhits_branch does not exist!\n");
				exit(1);
			}
			els_lost_pixelhits_isLoaded = true;
		}
		return els_lost_pixelhits_;
	}
	const vector<int> &els_mc_patMatch_id()
	{
		if (not els_mc_patMatch_id_isLoaded) {
			if (els_mc_patMatch_id_branch != 0) {
				els_mc_patMatch_id_branch->GetEntry(index);
			} else { 
				printf("branch els_mc_patMatch_id_branch does not exist!\n");
				exit(1);
			}
			els_mc_patMatch_id_isLoaded = true;
		}
		return els_mc_patMatch_id_;
	}
	const vector<int> &els_nSeed()
	{
		if (not els_nSeed_isLoaded) {
			if (els_nSeed_branch != 0) {
				els_nSeed_branch->GetEntry(index);
			} else { 
				printf("branch els_nSeed_branch does not exist!\n");
				exit(1);
			}
			els_nSeed_isLoaded = true;
		}
		return els_nSeed_;
	}
	const vector<int> &els_nlayers()
	{
		if (not els_nlayers_isLoaded) {
			if (els_nlayers_branch != 0) {
				els_nlayers_branch->GetEntry(index);
			} else { 
				printf("branch els_nlayers_branch does not exist!\n");
				exit(1);
			}
			els_nlayers_isLoaded = true;
		}
		return els_nlayers_;
	}
	const vector<int> &els_nlayers3D()
	{
		if (not els_nlayers3D_isLoaded) {
			if (els_nlayers3D_branch != 0) {
				els_nlayers3D_branch->GetEntry(index);
			} else { 
				printf("branch els_nlayers3D_branch does not exist!\n");
				exit(1);
			}
			els_nlayers3D_isLoaded = true;
		}
		return els_nlayers3D_;
	}
	const vector<int> &els_nlayersLost()
	{
		if (not els_nlayersLost_isLoaded) {
			if (els_nlayersLost_branch != 0) {
				els_nlayersLost_branch->GetEntry(index);
			} else { 
				printf("branch els_nlayersLost_branch does not exist!\n");
				exit(1);
			}
			els_nlayersLost_isLoaded = true;
		}
		return els_nlayersLost_;
	}
	const vector<int> &els_sccharge()
	{
		if (not els_sccharge_isLoaded) {
			if (els_sccharge_branch != 0) {
				els_sccharge_branch->GetEntry(index);
			} else { 
				printf("branch els_sccharge_branch does not exist!\n");
				exit(1);
			}
			els_sccharge_isLoaded = true;
		}
		return els_sccharge_;
	}
	const vector<int> &els_trk_charge()
	{
		if (not els_trk_charge_isLoaded) {
			if (els_trk_charge_branch != 0) {
				els_trk_charge_branch->GetEntry(index);
			} else { 
				printf("branch els_trk_charge_branch does not exist!\n");
				exit(1);
			}
			els_trk_charge_isLoaded = true;
		}
		return els_trk_charge_;
	}
	const vector<int> &els_type()
	{
		if (not els_type_isLoaded) {
			if (els_type_branch != 0) {
				els_type_branch->GetEntry(index);
			} else { 
				printf("branch els_type_branch does not exist!\n");
				exit(1);
			}
			els_type_isLoaded = true;
		}
		return els_type_;
	}
	const vector<int> &els_validHits()
	{
		if (not els_validHits_isLoaded) {
			if (els_validHits_branch != 0) {
				els_validHits_branch->GetEntry(index);
			} else { 
				printf("branch els_validHits_branch does not exist!\n");
				exit(1);
			}
			els_validHits_isLoaded = true;
		}
		return els_validHits_;
	}
	const vector<int> &els_valid_pixelhits()
	{
		if (not els_valid_pixelhits_isLoaded) {
			if (els_valid_pixelhits_branch != 0) {
				els_valid_pixelhits_branch->GetEntry(index);
			} else { 
				printf("branch els_valid_pixelhits_branch does not exist!\n");
				exit(1);
			}
			els_valid_pixelhits_isLoaded = true;
		}
		return els_valid_pixelhits_;
	}
	const vector<int> &els_passLooseId()
	{
		if (not els_passLooseId_isLoaded) {
			if (els_passLooseId_branch != 0) {
				els_passLooseId_branch->GetEntry(index);
			} else { 
				printf("branch els_passLooseId_branch does not exist!\n");
				exit(1);
			}
			els_passLooseId_isLoaded = true;
		}
		return els_passLooseId_;
	}
	const vector<int> &els_passMediumId()
	{
		if (not els_passMediumId_isLoaded) {
			if (els_passMediumId_branch != 0) {
				els_passMediumId_branch->GetEntry(index);
			} else { 
				printf("branch els_passMediumId_branch does not exist!\n");
				exit(1);
			}
			els_passMediumId_isLoaded = true;
		}
		return els_passMediumId_;
	}
	const vector<int> &els_passTightId()
	{
		if (not els_passTightId_isLoaded) {
			if (els_passTightId_branch != 0) {
				els_passTightId_branch->GetEntry(index);
			} else { 
				printf("branch els_passTightId_branch does not exist!\n");
				exit(1);
			}
			els_passTightId_isLoaded = true;
		}
		return els_passTightId_;
	}
	const vector<int> &els_passVetoId()
	{
		if (not els_passVetoId_isLoaded) {
			if (els_passVetoId_branch != 0) {
				els_passVetoId_branch->GetEntry(index);
			} else { 
				printf("branch els_passVetoId_branch does not exist!\n");
				exit(1);
			}
			els_passVetoId_isLoaded = true;
		}
		return els_passVetoId_;
	}
	const vector<int> &photons_fiduciality()
	{
		if (not photons_fiduciality_isLoaded) {
			if (photons_fiduciality_branch != 0) {
				photons_fiduciality_branch->GetEntry(index);
			} else { 
				printf("branch photons_fiduciality_branch does not exist!\n");
				exit(1);
			}
			photons_fiduciality_isLoaded = true;
		}
		return photons_fiduciality_;
	}
	const vector<int> &photons_photonID_loose()
	{
		if (not photons_photonID_loose_isLoaded) {
			if (photons_photonID_loose_branch != 0) {
				photons_photonID_loose_branch->GetEntry(index);
			} else { 
				printf("branch photons_photonID_loose_branch does not exist!\n");
				exit(1);
			}
			photons_photonID_loose_isLoaded = true;
		}
		return photons_photonID_loose_;
	}
	const vector<int> &photons_photonID_tight()
	{
		if (not photons_photonID_tight_isLoaded) {
			if (photons_photonID_tight_branch != 0) {
				photons_photonID_tight_branch->GetEntry(index);
			} else { 
				printf("branch photons_photonID_tight_branch does not exist!\n");
				exit(1);
			}
			photons_photonID_tight_isLoaded = true;
		}
		return photons_photonID_tight_;
	}
	const vector<int> &puInfo_bunchCrossing()
	{
		if (not puInfo_bunchCrossing_isLoaded) {
			if (puInfo_bunchCrossing_branch != 0) {
				puInfo_bunchCrossing_branch->GetEntry(index);
			} else { 
				printf("branch puInfo_bunchCrossing_branch does not exist!\n");
				exit(1);
			}
			puInfo_bunchCrossing_isLoaded = true;
		}
		return puInfo_bunchCrossing_;
	}
	const vector<int> &puInfo_nPUvertices()
	{
		if (not puInfo_nPUvertices_isLoaded) {
			if (puInfo_nPUvertices_branch != 0) {
				puInfo_nPUvertices_branch->GetEntry(index);
			} else { 
				printf("branch puInfo_nPUvertices_branch does not exist!\n");
				exit(1);
			}
			puInfo_nPUvertices_isLoaded = true;
		}
		return puInfo_nPUvertices_;
	}
	const vector<int> &vtxs_isFake()
	{
		if (not vtxs_isFake_isLoaded) {
			if (vtxs_isFake_branch != 0) {
				vtxs_isFake_branch->GetEntry(index);
			} else { 
				printf("branch vtxs_isFake_branch does not exist!\n");
				exit(1);
			}
			vtxs_isFake_isLoaded = true;
		}
		return vtxs_isFake_;
	}
	const vector<int> &vtxs_isValid()
	{
		if (not vtxs_isValid_isLoaded) {
			if (vtxs_isValid_branch != 0) {
				vtxs_isValid_branch->GetEntry(index);
			} else { 
				printf("branch vtxs_isValid_branch does not exist!\n");
				exit(1);
			}
			vtxs_isValid_isLoaded = true;
		}
		return vtxs_isValid_;
	}
	const vector<vector<int> > &els_PFCand_idx()
	{
		if (not els_PFCand_idx_isLoaded) {
			if (els_PFCand_idx_branch != 0) {
				els_PFCand_idx_branch->GetEntry(index);
			} else { 
				printf("branch els_PFCand_idx_branch does not exist!\n");
				exit(1);
			}
			els_PFCand_idx_isLoaded = true;
		}
		return els_PFCand_idx_;
	}
	const vector<vector<int> > &hlt_trigObjs_id()
	{
		if (not hlt_trigObjs_id_isLoaded) {
			if (hlt_trigObjs_id_branch != 0) {
				hlt_trigObjs_id_branch->GetEntry(index);
			} else { 
				printf("branch hlt_trigObjs_id_branch does not exist!\n");
				exit(1);
			}
			hlt_trigObjs_id_isLoaded = true;
		}
		return hlt_trigObjs_id_;
	}
	const vector<vector<int> > &photons_PFCand_idx()
	{
		if (not photons_PFCand_idx_isLoaded) {
			if (photons_PFCand_idx_branch != 0) {
				photons_PFCand_idx_branch->GetEntry(index);
			} else { 
				printf("branch photons_PFCand_idx_branch does not exist!\n");
				exit(1);
			}
			photons_PFCand_idx_isLoaded = true;
		}
		return photons_PFCand_idx_;
	}
	unsigned int &evt_nels()
	{
		if (not evt_nels_isLoaded) {
			if (evt_nels_branch != 0) {
				evt_nels_branch->GetEntry(index);
			} else { 
				printf("branch evt_nels_branch does not exist!\n");
				exit(1);
			}
			evt_nels_isLoaded = true;
		}
		return evt_nels_;
	}
	unsigned int &evt_detectorStatus()
	{
		if (not evt_detectorStatus_isLoaded) {
			if (evt_detectorStatus_branch != 0) {
				evt_detectorStatus_branch->GetEntry(index);
			} else { 
				printf("branch evt_detectorStatus_branch does not exist!\n");
				exit(1);
			}
			evt_detectorStatus_isLoaded = true;
		}
		return evt_detectorStatus_;
	}
	unsigned int &evt_lumiBlock()
	{
		if (not evt_lumiBlock_isLoaded) {
			if (evt_lumiBlock_branch != 0) {
				evt_lumiBlock_branch->GetEntry(index);
			} else { 
				printf("branch evt_lumiBlock_branch does not exist!\n");
				exit(1);
			}
			evt_lumiBlock_isLoaded = true;
		}
		return evt_lumiBlock_;
	}
  unsigned long long &evt_event()
  {
    if (not evt_event_isLoaded) {
      if (evt_event_branch != 0) {
	evt_event_branch->GetEntry(index);
      } else { 
	printf("branch evt_event_branch does not exist!\n");
	exit(1);
      }
      evt_event_isLoaded = true;
    }
    return evt_event_;
  }
	unsigned int &evt_run()
	{
		if (not evt_run_isLoaded) {
			if (evt_run_branch != 0) {
				evt_run_branch->GetEntry(index);
			} else { 
				printf("branch evt_run_branch does not exist!\n");
				exit(1);
			}
			evt_run_isLoaded = true;
		}
		return evt_run_;
	}
	unsigned int &evt_nphotons()
	{
		if (not evt_nphotons_isLoaded) {
			if (evt_nphotons_branch != 0) {
				evt_nphotons_branch->GetEntry(index);
			} else { 
				printf("branch evt_nphotons_branch does not exist!\n");
				exit(1);
			}
			evt_nphotons_isLoaded = true;
		}
		return evt_nphotons_;
	}
	unsigned int &evt_nvtxs()
	{
		if (not evt_nvtxs_isLoaded) {
			if (evt_nvtxs_branch != 0) {
				evt_nvtxs_branch->GetEntry(index);
			} else { 
				printf("branch evt_nvtxs_branch does not exist!\n");
				exit(1);
			}
			evt_nvtxs_isLoaded = true;
		}
		return evt_nvtxs_;
	}
	const vector<unsigned int> &hlt_prescales()
	{
		if (not hlt_prescales_isLoaded) {
			if (hlt_prescales_branch != 0) {
				hlt_prescales_branch->GetEntry(index);
			} else { 
				printf("branch hlt_prescales_branch does not exist!\n");
				exit(1);
			}
			hlt_prescales_isLoaded = true;
		}
		return hlt_prescales_;
	}
	bool passHLTTrigger(TString trigName) {
		int trigIndx;
		vector<TString>::const_iterator begin_it = hlt_trigNames().begin();
		vector<TString>::const_iterator end_it = hlt_trigNames().end();
		vector<TString>::const_iterator found_it = find(begin_it, end_it, trigName);
		if(found_it != end_it)
			trigIndx = found_it - begin_it;
		else {
			cout << "Cannot find Trigger " << trigName << endl; 
			return 0;
		}

	return hlt_bits().TestBitNumber(trigIndx);
	}

  static void progress( int nEventsTotal, int nEventsChain ){
    int period = 1000;
    if(nEventsTotal%1000 == 0) {
      // xterm magic from L. Vacavant and A. Cerri
      if (isatty(1)) {
        if( ( nEventsChain - nEventsTotal ) > period ){
          float frac = (float)nEventsTotal/(nEventsChain*0.01);
          printf("\015\033[32m ---> \033[1m\033[31m%4.1f%%"
               "\033[0m\033[32m <---\033[0m\015", frac);
          fflush(stdout);
        }
        else {
          printf("\015\033[32m ---> \033[1m\033[31m%4.1f%%"
                 "\033[0m\033[32m <---\033[0m\015", 100.);
          cout << endl;
        }
      }
    }
  }
  
};

#ifndef __CINT__
extern CMS3 cms3;
#endif

namespace tas {
	const TBits &hlt_bits();
	const vector<TString> &evt_CMS3tag();
	const vector<TString> &evt_dataset();
	const vector<TString> &hlt_trigNames();
	const vector<bool> &els_conv_vtx_flag();
	const vector<bool> &els_isGsfCtfScPixChargeConsistent();
	const vector<bool> &els_passingMvaPreselection();
	const vector<bool> &els_passingPflowPreselection();
	const vector<bool> &photons_haspixelSeed();
	const float &evt_bs_Xwidth();
	const float &evt_bs_XwidthErr();
	const float &evt_bs_Ywidth();
	const float &evt_bs_YwidthErr();
	const float &evt_bs_dxdz();
	const float &evt_bs_dxdzErr();
	const float &evt_bs_dydz();
	const float &evt_bs_dydzErr();
	const float &evt_bs_sigmaZ();
	const float &evt_bs_sigmaZErr();
	const float &evt_bs_xErr();
	const float &evt_bs_yErr();
	const float &evt_bs_zErr();
	const float &evt_bField();
	const float &evt_fixgrid_all_rho();
	const float &evt_fixgridfastjet_allcalo_rho();
	const float &evt_fixgridfastjet_all_rho();
	const float &evt_fixgridfastjet_centralcalo_rho();
	const float &evt_fixgridfastjet_centralchargedpileup_rho();
	const float &evt_fixgridfastjet_centralneutral_rho();
	const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >  &evt_bsp4();
	const vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &els_mc_patMatch_p4();
	const vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &els_p4();
	const vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &els_p4In();
	const vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &els_p4Out();
	const vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &els_trk_p4();
	const vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &els_trk_vertex_p4();
	const vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &els_vertex_p4();
	const vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &photons_p4();
	const vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &vtxs_position();
	const vector<vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > > &hlt_trigObjs_p4();
	const vector<float> &evt_bs_covMatrix();
	const vector<float> &els_bs2d();
	const vector<float> &els_bs2derr();
	const vector<float> &els_bs3d();
	const vector<float> &els_bs3derr();
	const vector<float> &els_chi2();
	const vector<float> &els_ckf_chi2();
	const vector<float> &els_ckf_ndof();
	const vector<float> &els_d0();
	const vector<float> &els_d0Err();
	const vector<float> &els_d0corr();
	const vector<float> &els_d0corrPhi();
	const vector<float> &els_d0phiCov();
	const vector<float> &els_dEtaIn();
	const vector<float> &els_dEtaOut();
	const vector<float> &els_dPhiIn();
	const vector<float> &els_dPhiInPhiOut();
	const vector<float> &els_dPhiOut();
	const vector<float> &els_deltaEtaEleClusterTrackAtCalo();
	const vector<float> &els_deltaPhiEleClusterTrackAtCalo();
	const vector<float> &els_dxyPV();
	const vector<float> &els_dzPV();
	const vector<float> &els_e1x5();
	const vector<float> &els_e1x5_full5x5();
	const vector<float> &els_e2x5Max();
	const vector<float> &els_e2x5Max_full5x5();
	const vector<float> &els_e5x5();
	const vector<float> &els_e5x5_full5x5();
	const vector<float> &els_eOverPIn();
	const vector<float> &els_eOverPOut();
	const vector<float> &els_eSC();
	const vector<float> &els_eSCPresh();
	const vector<float> &els_eSCRaw();
	const vector<float> &els_eSeed();
	const vector<float> &els_eSeedOverPIn();
	const vector<float> &els_eSeedOverPOut();
	const vector<float> &els_ecalEnergy();
	const vector<float> &els_ecalEnergyError();
	const vector<float> &els_ecalIso();
	const vector<float> &els_ecalIso04();
	const vector<float> &els_ecalPFClusterIso();
	const vector<float> &els_etaErr();
	const vector<float> &els_etaSC();
	const vector<float> &els_etaSCwidth();
	const vector<float> &els_fbrem();
	const vector<float> &els_hOverE();
	const vector<float> &els_hOverEBC();
	const vector<float> &els_hcalDepth1OverEcal();
	const vector<float> &els_hcalDepth1TowerSumEt();
	const vector<float> &els_hcalDepth1TowerSumEt04();
	const vector<float> &els_hcalDepth2OverEcal();
	const vector<float> &els_hcalDepth2TowerSumEt();
	const vector<float> &els_hcalDepth2TowerSumEt04();
	const vector<float> &els_hcalIso();
	const vector<float> &els_hcalIso04();
	const vector<float> &els_hcalPFClusterIso();
	const vector<float> &els_ip2d();
	const vector<float> &els_ip2derr();
	const vector<float> &els_ip3d();
	const vector<float> &els_ip3derr();
	const vector<float> &els_mass();
	const vector<float> &els_mc_patMatch_dr();
	const vector<float> &els_miniIso_ch();
	const vector<float> &els_miniIso_db();
	const vector<float> &els_miniIso_em();
	const vector<float> &els_miniIso_nh();
	const vector<float> &els_miniIso_uncor();
	const vector<float> &els_mva();
	const vector<float> &els_ndof();
	const vector<float> &els_pfChargedHadronIso();
	const vector<float> &els_pfNeutralHadronIso();
	const vector<float> &els_pfPUIso();
	const vector<float> &els_pfPhotonIso();
	const vector<float> &els_phiErr();
	const vector<float> &els_phiSC();
	const vector<float> &els_phiSCwidth();
	const vector<float> &els_ptErr();
	const vector<float> &els_ptErrGsf();
	const vector<float> &els_r9();
	const vector<float> &els_r9_full5x5();
	const vector<float> &els_sigmaEtaEta();
	const vector<float> &els_sigmaEtaEta_full5x5();
	const vector<float> &els_sigmaIEtaIEta();
	const vector<float> &els_sigmaIEtaIEta_full5x5();
	const vector<float> &els_sigmaIPhiIPhi();
	const vector<float> &els_sigmaIPhiIPhi_full5x5();
	const vector<float> &els_sigmaIphiIphi();
	const vector<float> &els_tkIso();
	const vector<float> &els_tkIso04();
	const vector<float> &els_trackMomentumError();
	const vector<float> &els_trkdr();
	const vector<float> &els_trkshFrac();
	const vector<float> &els_z0();
	const vector<float> &els_z0Err();
	const vector<float> &els_z0corr();
	const vector<float> &photons_chargedHadronIso();
	const vector<float> &photons_e1x5();
	const vector<float> &photons_e2x5Max();
	const vector<float> &photons_e3x3();
	const vector<float> &photons_e5x5();
	const vector<float> &photons_eSC();
	const vector<float> &photons_eSCPresh();
	const vector<float> &photons_eSCRaw();
	const vector<float> &photons_ecalIso03();
	const vector<float> &photons_ecalIso04();
	const vector<float> &photons_ecalPFClusterIso();
	const vector<float> &photons_etaSC();
	const vector<float> &photons_full3x3_e3x3();
	const vector<float> &photons_full5x5_e1x5();
	const vector<float> &photons_full5x5_e2x5Max();
	const vector<float> &photons_full5x5_e5x5();
	const vector<float> &photons_full5x5_hOverE();
	const vector<float> &photons_full5x5_hOverEtowBC();
	const vector<float> &photons_full5x5_r9();
	const vector<float> &photons_full5x5_sigmaEtaEta();
	const vector<float> &photons_full5x5_sigmaIEtaIEta();
	const vector<float> &photons_hOverE();
	const vector<float> &photons_hOverEtowBC();
	const vector<float> &photons_hcalDepth1TowerSumEtBcConeDR03();
	const vector<float> &photons_hcalDepth1TowerSumEtBcConeDR04();
	const vector<float> &photons_hcalDepth2TowerSumEtBcConeDR03();
	const vector<float> &photons_hcalDepth2TowerSumEtBcConeDR04();
	const vector<float> &photons_hcalIso03();
	const vector<float> &photons_hcalIso04();
	const vector<float> &photons_hcalPFClusterIso();
	const vector<float> &photons_hcalTowerSumEtBcConeDR03();
	const vector<float> &photons_hcalTowerSumEtBcConeDR04();
	const vector<float> &photons_mass();
	const vector<float> &photons_neutralHadronIso();
	const vector<float> &photons_ntkIsoHollow03();
	const vector<float> &photons_ntkIsoHollow04();
	const vector<float> &photons_ntkIsoSolid03();
	const vector<float> &photons_ntkIsoSolid04();
	const vector<float> &photons_phiSC();
	const vector<float> &photons_photonIso();
	const vector<float> &photons_recoChargedHadronIso();
	const vector<float> &photons_recoNeutralHadronIso();
	const vector<float> &photons_recoPhotonIso();
	const vector<float> &photons_sigmaEtaEta();
	const vector<float> &photons_sigmaIEtaIEta();
	const vector<float> &photons_tkIsoHollow03();
	const vector<float> &photons_tkIsoHollow04();
	const vector<float> &photons_tkIsoSolid03();
	const vector<float> &photons_tkIsoSolid04();
	const vector<float> &puInfo_trueNumInteractions();
	const vector<float> &vtxs_chi2();
	const vector<float> &vtxs_ndof();
	const vector<float> &vtxs_score();
	const vector<float> &vtxs_xError();
	const vector<float> &vtxs_yError();
	const vector<float> &vtxs_zError();
	const vector<vector<float> > &puInfo_instLumi();
	const vector<vector<float> > &vtxs_covMatrix();
	const int &evt_bsType();
	const int &evt_bunchCrossing();
	const int &evt_experimentType();
	const int &evt_isRealData();
	const int &evt_orbitNumber();
	const int &evt_storeNumber();
	const vector<int> &els_category();
	const vector<int> &els_charge();
	const vector<int> &els_ckf_charge();
	const vector<int> &els_ckf_laywithmeas();
	const vector<int> &els_class();
	const vector<int> &els_exp_innerlayers();
	const vector<int> &els_exp_outerlayers();
	const vector<int> &els_fiduciality();
	const vector<int> &els_lostHits();
	const vector<int> &els_lost_pixelhits();
	const vector<int> &els_mc_patMatch_id();
	const vector<int> &els_nSeed();
	const vector<int> &els_nlayers();
	const vector<int> &els_nlayers3D();
	const vector<int> &els_nlayersLost();
	const vector<int> &els_sccharge();
	const vector<int> &els_trk_charge();
	const vector<int> &els_type();
	const vector<int> &els_validHits();
	const vector<int> &els_valid_pixelhits();
	const vector<int> &els_passLooseId();
	const vector<int> &els_passMediumId();
	const vector<int> &els_passTightId();
	const vector<int> &els_passVetoId();
	const vector<int> &photons_fiduciality();
	const vector<int> &photons_photonID_loose();
	const vector<int> &photons_photonID_tight();
	const vector<int> &puInfo_bunchCrossing();
	const vector<int> &puInfo_nPUvertices();
	const vector<int> &vtxs_isFake();
	const vector<int> &vtxs_isValid();
	const vector<vector<int> > &els_PFCand_idx();
	const vector<vector<int> > &hlt_trigObjs_id();
	const vector<vector<int> > &photons_PFCand_idx();
	const unsigned int &evt_nels();
	const unsigned int &evt_detectorStatus();
	const unsigned int &evt_lumiBlock();
  const unsigned long long &evt_event();
	const unsigned int &evt_run();
	const unsigned int &evt_nphotons();
	const unsigned int &evt_nvtxs();
	const vector<unsigned int> &hlt_prescales();
	bool passHLTTrigger(TString trigName);
}
#endif
