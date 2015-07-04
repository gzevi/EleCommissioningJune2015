#include "CMS3.h"
CMS3 cms3;
namespace tas {
	const TBits &hlt_bits() { return cms3.hlt_bits(); }
	const vector<TString> &evt_CMS3tag() { return cms3.evt_CMS3tag(); }
	const vector<TString> &evt_dataset() { return cms3.evt_dataset(); }
	const vector<TString> &hlt_trigNames() { return cms3.hlt_trigNames(); }
	const vector<bool> &els_conv_vtx_flag() { return cms3.els_conv_vtx_flag(); }
	const vector<bool> &els_isGsfCtfScPixChargeConsistent() { return cms3.els_isGsfCtfScPixChargeConsistent(); }
	const vector<bool> &els_passingMvaPreselection() { return cms3.els_passingMvaPreselection(); }
	const vector<bool> &els_passingPflowPreselection() { return cms3.els_passingPflowPreselection(); }
	const vector<bool> &photons_haspixelSeed() { return cms3.photons_haspixelSeed(); }
	const float &evt_bs_Xwidth() { return cms3.evt_bs_Xwidth(); }
	const float &evt_bs_XwidthErr() { return cms3.evt_bs_XwidthErr(); }
	const float &evt_bs_Ywidth() { return cms3.evt_bs_Ywidth(); }
	const float &evt_bs_YwidthErr() { return cms3.evt_bs_YwidthErr(); }
	const float &evt_bs_dxdz() { return cms3.evt_bs_dxdz(); }
	const float &evt_bs_dxdzErr() { return cms3.evt_bs_dxdzErr(); }
	const float &evt_bs_dydz() { return cms3.evt_bs_dydz(); }
	const float &evt_bs_dydzErr() { return cms3.evt_bs_dydzErr(); }
	const float &evt_bs_sigmaZ() { return cms3.evt_bs_sigmaZ(); }
	const float &evt_bs_sigmaZErr() { return cms3.evt_bs_sigmaZErr(); }
	const float &evt_bs_xErr() { return cms3.evt_bs_xErr(); }
	const float &evt_bs_yErr() { return cms3.evt_bs_yErr(); }
	const float &evt_bs_zErr() { return cms3.evt_bs_zErr(); }
	const float &evt_bField() { return cms3.evt_bField(); }
	const float &evt_fixgrid_all_rho() { return cms3.evt_fixgrid_all_rho(); }
	const float &evt_fixgridfastjet_allcalo_rho() { return cms3.evt_fixgridfastjet_allcalo_rho(); }
	const float &evt_fixgridfastjet_all_rho() { return cms3.evt_fixgridfastjet_all_rho(); }
	const float &evt_fixgridfastjet_centralcalo_rho() { return cms3.evt_fixgridfastjet_centralcalo_rho(); }
	const float &evt_fixgridfastjet_centralchargedpileup_rho() { return cms3.evt_fixgridfastjet_centralchargedpileup_rho(); }
	const float &evt_fixgridfastjet_centralneutral_rho() { return cms3.evt_fixgridfastjet_centralneutral_rho(); }
	const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >  &evt_bsp4() { return cms3.evt_bsp4(); }
	const vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &els_mc_patMatch_p4() { return cms3.els_mc_patMatch_p4(); }
	const vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &els_p4() { return cms3.els_p4(); }
	const vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &els_p4In() { return cms3.els_p4In(); }
	const vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &els_p4Out() { return cms3.els_p4Out(); }
	const vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &els_trk_p4() { return cms3.els_trk_p4(); }
	const vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &els_trk_vertex_p4() { return cms3.els_trk_vertex_p4(); }
	const vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &els_vertex_p4() { return cms3.els_vertex_p4(); }
	const vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &photons_p4() { return cms3.photons_p4(); }
	const vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &vtxs_position() { return cms3.vtxs_position(); }
	const vector<vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > > &hlt_trigObjs_p4() { return cms3.hlt_trigObjs_p4(); }
	const vector<float> &evt_bs_covMatrix() { return cms3.evt_bs_covMatrix(); }
	const vector<float> &els_bs2d() { return cms3.els_bs2d(); }
	const vector<float> &els_bs2derr() { return cms3.els_bs2derr(); }
	const vector<float> &els_bs3d() { return cms3.els_bs3d(); }
	const vector<float> &els_bs3derr() { return cms3.els_bs3derr(); }
	const vector<float> &els_chi2() { return cms3.els_chi2(); }
	const vector<float> &els_ckf_chi2() { return cms3.els_ckf_chi2(); }
	const vector<float> &els_ckf_ndof() { return cms3.els_ckf_ndof(); }
	const vector<float> &els_d0() { return cms3.els_d0(); }
	const vector<float> &els_d0Err() { return cms3.els_d0Err(); }
	const vector<float> &els_d0corr() { return cms3.els_d0corr(); }
	const vector<float> &els_d0corrPhi() { return cms3.els_d0corrPhi(); }
	const vector<float> &els_d0phiCov() { return cms3.els_d0phiCov(); }
	const vector<float> &els_dEtaIn() { return cms3.els_dEtaIn(); }
	const vector<float> &els_dEtaOut() { return cms3.els_dEtaOut(); }
	const vector<float> &els_dPhiIn() { return cms3.els_dPhiIn(); }
	const vector<float> &els_dPhiInPhiOut() { return cms3.els_dPhiInPhiOut(); }
	const vector<float> &els_dPhiOut() { return cms3.els_dPhiOut(); }
	const vector<float> &els_deltaEtaEleClusterTrackAtCalo() { return cms3.els_deltaEtaEleClusterTrackAtCalo(); }
	const vector<float> &els_deltaPhiEleClusterTrackAtCalo() { return cms3.els_deltaPhiEleClusterTrackAtCalo(); }
	const vector<float> &els_dxyPV() { return cms3.els_dxyPV(); }
	const vector<float> &els_dzPV() { return cms3.els_dzPV(); }
	const vector<float> &els_e1x5() { return cms3.els_e1x5(); }
	const vector<float> &els_e1x5_full5x5() { return cms3.els_e1x5_full5x5(); }
	const vector<float> &els_e2x5Max() { return cms3.els_e2x5Max(); }
	const vector<float> &els_e2x5Max_full5x5() { return cms3.els_e2x5Max_full5x5(); }
	const vector<float> &els_e5x5() { return cms3.els_e5x5(); }
	const vector<float> &els_e5x5_full5x5() { return cms3.els_e5x5_full5x5(); }
	const vector<float> &els_eOverPIn() { return cms3.els_eOverPIn(); }
	const vector<float> &els_eOverPOut() { return cms3.els_eOverPOut(); }
	const vector<float> &els_eSC() { return cms3.els_eSC(); }
	const vector<float> &els_eSCPresh() { return cms3.els_eSCPresh(); }
	const vector<float> &els_eSCRaw() { return cms3.els_eSCRaw(); }
	const vector<float> &els_eSeed() { return cms3.els_eSeed(); }
	const vector<float> &els_eSeedOverPIn() { return cms3.els_eSeedOverPIn(); }
	const vector<float> &els_eSeedOverPOut() { return cms3.els_eSeedOverPOut(); }
	const vector<float> &els_ecalEnergy() { return cms3.els_ecalEnergy(); }
	const vector<float> &els_ecalEnergyError() { return cms3.els_ecalEnergyError(); }
	const vector<float> &els_ecalIso() { return cms3.els_ecalIso(); }
	const vector<float> &els_ecalIso04() { return cms3.els_ecalIso04(); }
	const vector<float> &els_ecalPFClusterIso() { return cms3.els_ecalPFClusterIso(); }
	const vector<float> &els_etaErr() { return cms3.els_etaErr(); }
	const vector<float> &els_etaSC() { return cms3.els_etaSC(); }
	const vector<float> &els_etaSCwidth() { return cms3.els_etaSCwidth(); }
	const vector<float> &els_fbrem() { return cms3.els_fbrem(); }
	const vector<float> &els_hOverE() { return cms3.els_hOverE(); }
	const vector<float> &els_hOverEBC() { return cms3.els_hOverEBC(); }
	const vector<float> &els_hcalDepth1OverEcal() { return cms3.els_hcalDepth1OverEcal(); }
	const vector<float> &els_hcalDepth1TowerSumEt() { return cms3.els_hcalDepth1TowerSumEt(); }
	const vector<float> &els_hcalDepth1TowerSumEt04() { return cms3.els_hcalDepth1TowerSumEt04(); }
	const vector<float> &els_hcalDepth2OverEcal() { return cms3.els_hcalDepth2OverEcal(); }
	const vector<float> &els_hcalDepth2TowerSumEt() { return cms3.els_hcalDepth2TowerSumEt(); }
	const vector<float> &els_hcalDepth2TowerSumEt04() { return cms3.els_hcalDepth2TowerSumEt04(); }
	const vector<float> &els_hcalIso() { return cms3.els_hcalIso(); }
	const vector<float> &els_hcalIso04() { return cms3.els_hcalIso04(); }
	const vector<float> &els_hcalPFClusterIso() { return cms3.els_hcalPFClusterIso(); }
	const vector<float> &els_ip2d() { return cms3.els_ip2d(); }
	const vector<float> &els_ip2derr() { return cms3.els_ip2derr(); }
	const vector<float> &els_ip3d() { return cms3.els_ip3d(); }
	const vector<float> &els_ip3derr() { return cms3.els_ip3derr(); }
	const vector<float> &els_mass() { return cms3.els_mass(); }
	const vector<float> &els_mc_patMatch_dr() { return cms3.els_mc_patMatch_dr(); }
	const vector<float> &els_miniIso_ch() { return cms3.els_miniIso_ch(); }
	const vector<float> &els_miniIso_db() { return cms3.els_miniIso_db(); }
	const vector<float> &els_miniIso_em() { return cms3.els_miniIso_em(); }
	const vector<float> &els_miniIso_nh() { return cms3.els_miniIso_nh(); }
	const vector<float> &els_miniIso_uncor() { return cms3.els_miniIso_uncor(); }
	const vector<float> &els_mva() { return cms3.els_mva(); }
	const vector<float> &els_ndof() { return cms3.els_ndof(); }
	const vector<float> &els_pfChargedHadronIso() { return cms3.els_pfChargedHadronIso(); }
	const vector<float> &els_pfNeutralHadronIso() { return cms3.els_pfNeutralHadronIso(); }
	const vector<float> &els_pfPUIso() { return cms3.els_pfPUIso(); }
	const vector<float> &els_pfPhotonIso() { return cms3.els_pfPhotonIso(); }
	const vector<float> &els_phiErr() { return cms3.els_phiErr(); }
	const vector<float> &els_phiSC() { return cms3.els_phiSC(); }
	const vector<float> &els_phiSCwidth() { return cms3.els_phiSCwidth(); }
	const vector<float> &els_ptErr() { return cms3.els_ptErr(); }
	const vector<float> &els_ptErrGsf() { return cms3.els_ptErrGsf(); }
	const vector<float> &els_r9() { return cms3.els_r9(); }
	const vector<float> &els_r9_full5x5() { return cms3.els_r9_full5x5(); }
	const vector<float> &els_sigmaEtaEta() { return cms3.els_sigmaEtaEta(); }
	const vector<float> &els_sigmaEtaEta_full5x5() { return cms3.els_sigmaEtaEta_full5x5(); }
	const vector<float> &els_sigmaIEtaIEta() { return cms3.els_sigmaIEtaIEta(); }
	const vector<float> &els_sigmaIEtaIEta_full5x5() { return cms3.els_sigmaIEtaIEta_full5x5(); }
	const vector<float> &els_sigmaIPhiIPhi() { return cms3.els_sigmaIPhiIPhi(); }
	const vector<float> &els_sigmaIPhiIPhi_full5x5() { return cms3.els_sigmaIPhiIPhi_full5x5(); }
	const vector<float> &els_sigmaIphiIphi() { return cms3.els_sigmaIphiIphi(); }
	const vector<float> &els_tkIso() { return cms3.els_tkIso(); }
	const vector<float> &els_tkIso04() { return cms3.els_tkIso04(); }
	const vector<float> &els_trackMomentumError() { return cms3.els_trackMomentumError(); }
	const vector<float> &els_trkdr() { return cms3.els_trkdr(); }
	const vector<float> &els_trkshFrac() { return cms3.els_trkshFrac(); }
	const vector<float> &els_z0() { return cms3.els_z0(); }
	const vector<float> &els_z0Err() { return cms3.els_z0Err(); }
	const vector<float> &els_z0corr() { return cms3.els_z0corr(); }
	const vector<float> &photons_chargedHadronIso() { return cms3.photons_chargedHadronIso(); }
	const vector<float> &photons_e1x5() { return cms3.photons_e1x5(); }
	const vector<float> &photons_e2x5Max() { return cms3.photons_e2x5Max(); }
	const vector<float> &photons_e3x3() { return cms3.photons_e3x3(); }
	const vector<float> &photons_e5x5() { return cms3.photons_e5x5(); }
	const vector<float> &photons_eSC() { return cms3.photons_eSC(); }
	const vector<float> &photons_eSCPresh() { return cms3.photons_eSCPresh(); }
	const vector<float> &photons_eSCRaw() { return cms3.photons_eSCRaw(); }
	const vector<float> &photons_ecalIso03() { return cms3.photons_ecalIso03(); }
	const vector<float> &photons_ecalIso04() { return cms3.photons_ecalIso04(); }
	const vector<float> &photons_ecalPFClusterIso() { return cms3.photons_ecalPFClusterIso(); }
	const vector<float> &photons_etaSC() { return cms3.photons_etaSC(); }
	const vector<float> &photons_full3x3_e3x3() { return cms3.photons_full3x3_e3x3(); }
	const vector<float> &photons_full5x5_e1x5() { return cms3.photons_full5x5_e1x5(); }
	const vector<float> &photons_full5x5_e2x5Max() { return cms3.photons_full5x5_e2x5Max(); }
	const vector<float> &photons_full5x5_e5x5() { return cms3.photons_full5x5_e5x5(); }
	const vector<float> &photons_full5x5_hOverE() { return cms3.photons_full5x5_hOverE(); }
	const vector<float> &photons_full5x5_hOverEtowBC() { return cms3.photons_full5x5_hOverEtowBC(); }
	const vector<float> &photons_full5x5_r9() { return cms3.photons_full5x5_r9(); }
	const vector<float> &photons_full5x5_sigmaEtaEta() { return cms3.photons_full5x5_sigmaEtaEta(); }
	const vector<float> &photons_full5x5_sigmaIEtaIEta() { return cms3.photons_full5x5_sigmaIEtaIEta(); }
	const vector<float> &photons_hOverE() { return cms3.photons_hOverE(); }
	const vector<float> &photons_hOverEtowBC() { return cms3.photons_hOverEtowBC(); }
	const vector<float> &photons_hcalDepth1TowerSumEtBcConeDR03() { return cms3.photons_hcalDepth1TowerSumEtBcConeDR03(); }
	const vector<float> &photons_hcalDepth1TowerSumEtBcConeDR04() { return cms3.photons_hcalDepth1TowerSumEtBcConeDR04(); }
	const vector<float> &photons_hcalDepth2TowerSumEtBcConeDR03() { return cms3.photons_hcalDepth2TowerSumEtBcConeDR03(); }
	const vector<float> &photons_hcalDepth2TowerSumEtBcConeDR04() { return cms3.photons_hcalDepth2TowerSumEtBcConeDR04(); }
	const vector<float> &photons_hcalIso03() { return cms3.photons_hcalIso03(); }
	const vector<float> &photons_hcalIso04() { return cms3.photons_hcalIso04(); }
	const vector<float> &photons_hcalPFClusterIso() { return cms3.photons_hcalPFClusterIso(); }
	const vector<float> &photons_hcalTowerSumEtBcConeDR03() { return cms3.photons_hcalTowerSumEtBcConeDR03(); }
	const vector<float> &photons_hcalTowerSumEtBcConeDR04() { return cms3.photons_hcalTowerSumEtBcConeDR04(); }
	const vector<float> &photons_mass() { return cms3.photons_mass(); }
	const vector<float> &photons_neutralHadronIso() { return cms3.photons_neutralHadronIso(); }
	const vector<float> &photons_ntkIsoHollow03() { return cms3.photons_ntkIsoHollow03(); }
	const vector<float> &photons_ntkIsoHollow04() { return cms3.photons_ntkIsoHollow04(); }
	const vector<float> &photons_ntkIsoSolid03() { return cms3.photons_ntkIsoSolid03(); }
	const vector<float> &photons_ntkIsoSolid04() { return cms3.photons_ntkIsoSolid04(); }
	const vector<float> &photons_phiSC() { return cms3.photons_phiSC(); }
	const vector<float> &photons_photonIso() { return cms3.photons_photonIso(); }
	const vector<float> &photons_recoChargedHadronIso() { return cms3.photons_recoChargedHadronIso(); }
	const vector<float> &photons_recoNeutralHadronIso() { return cms3.photons_recoNeutralHadronIso(); }
	const vector<float> &photons_recoPhotonIso() { return cms3.photons_recoPhotonIso(); }
	const vector<float> &photons_sigmaEtaEta() { return cms3.photons_sigmaEtaEta(); }
	const vector<float> &photons_sigmaIEtaIEta() { return cms3.photons_sigmaIEtaIEta(); }
	const vector<float> &photons_tkIsoHollow03() { return cms3.photons_tkIsoHollow03(); }
	const vector<float> &photons_tkIsoHollow04() { return cms3.photons_tkIsoHollow04(); }
	const vector<float> &photons_tkIsoSolid03() { return cms3.photons_tkIsoSolid03(); }
	const vector<float> &photons_tkIsoSolid04() { return cms3.photons_tkIsoSolid04(); }
	const vector<float> &puInfo_trueNumInteractions() { return cms3.puInfo_trueNumInteractions(); }
	const vector<float> &vtxs_chi2() { return cms3.vtxs_chi2(); }
	const vector<float> &vtxs_ndof() { return cms3.vtxs_ndof(); }
	const vector<float> &vtxs_score() { return cms3.vtxs_score(); }
	const vector<float> &vtxs_xError() { return cms3.vtxs_xError(); }
	const vector<float> &vtxs_yError() { return cms3.vtxs_yError(); }
	const vector<float> &vtxs_zError() { return cms3.vtxs_zError(); }
	const vector<vector<float> > &puInfo_instLumi() { return cms3.puInfo_instLumi(); }
	const vector<vector<float> > &vtxs_covMatrix() { return cms3.vtxs_covMatrix(); }
	const int &evt_bsType() { return cms3.evt_bsType(); }
	const int &evt_bunchCrossing() { return cms3.evt_bunchCrossing(); }
	const int &evt_experimentType() { return cms3.evt_experimentType(); }
	const int &evt_isRealData() { return cms3.evt_isRealData(); }
	const int &evt_orbitNumber() { return cms3.evt_orbitNumber(); }
	const int &evt_storeNumber() { return cms3.evt_storeNumber(); }
	const vector<int> &els_category() { return cms3.els_category(); }
	const vector<int> &els_charge() { return cms3.els_charge(); }
	const vector<int> &els_ckf_charge() { return cms3.els_ckf_charge(); }
	const vector<int> &els_ckf_laywithmeas() { return cms3.els_ckf_laywithmeas(); }
	const vector<int> &els_class() { return cms3.els_class(); }
	const vector<int> &els_exp_innerlayers() { return cms3.els_exp_innerlayers(); }
	const vector<int> &els_exp_outerlayers() { return cms3.els_exp_outerlayers(); }
	const vector<int> &els_fiduciality() { return cms3.els_fiduciality(); }
	const vector<int> &els_lostHits() { return cms3.els_lostHits(); }
	const vector<int> &els_lost_pixelhits() { return cms3.els_lost_pixelhits(); }
	const vector<int> &els_mc_patMatch_id() { return cms3.els_mc_patMatch_id(); }
	const vector<int> &els_nSeed() { return cms3.els_nSeed(); }
	const vector<int> &els_nlayers() { return cms3.els_nlayers(); }
	const vector<int> &els_nlayers3D() { return cms3.els_nlayers3D(); }
	const vector<int> &els_nlayersLost() { return cms3.els_nlayersLost(); }
	const vector<int> &els_sccharge() { return cms3.els_sccharge(); }
	const vector<int> &els_trk_charge() { return cms3.els_trk_charge(); }
	const vector<int> &els_type() { return cms3.els_type(); }
	const vector<int> &els_validHits() { return cms3.els_validHits(); }
	const vector<int> &els_valid_pixelhits() { return cms3.els_valid_pixelhits(); }
	const vector<int> &els_passLooseId() { return cms3.els_passLooseId(); }
	const vector<int> &els_passMediumId() { return cms3.els_passMediumId(); }
	const vector<int> &els_passTightId() { return cms3.els_passTightId(); }
	const vector<int> &els_passVetoId() { return cms3.els_passVetoId(); }
	const vector<int> &photons_fiduciality() { return cms3.photons_fiduciality(); }
	const vector<int> &photons_photonID_loose() { return cms3.photons_photonID_loose(); }
	const vector<int> &photons_photonID_tight() { return cms3.photons_photonID_tight(); }
	const vector<int> &puInfo_bunchCrossing() { return cms3.puInfo_bunchCrossing(); }
	const vector<int> &puInfo_nPUvertices() { return cms3.puInfo_nPUvertices(); }
	const vector<int> &vtxs_isFake() { return cms3.vtxs_isFake(); }
	const vector<int> &vtxs_isValid() { return cms3.vtxs_isValid(); }
	const vector<vector<int> > &els_PFCand_idx() { return cms3.els_PFCand_idx(); }
	const vector<vector<int> > &hlt_trigObjs_id() { return cms3.hlt_trigObjs_id(); }
	const vector<vector<int> > &photons_PFCand_idx() { return cms3.photons_PFCand_idx(); }
	const unsigned int &evt_nels() { return cms3.evt_nels(); }
	const unsigned int &evt_detectorStatus() { return cms3.evt_detectorStatus(); }
	const unsigned int &evt_lumiBlock() { return cms3.evt_lumiBlock(); }
  const unsigned long long &evt_event() { return cms3.evt_event(); }

	const unsigned int &evt_run() { return cms3.evt_run(); }
	const unsigned int &evt_nphotons() { return cms3.evt_nphotons(); }
	const unsigned int &evt_nvtxs() { return cms3.evt_nvtxs(); }
	const vector<unsigned int> &hlt_prescales() { return cms3.hlt_prescales(); }
	bool passHLTTrigger(TString trigName) { return cms3.passHLTTrigger(trigName); }
}
