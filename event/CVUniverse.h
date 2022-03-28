// =============================================================================
// Base class for an un-systematically shifted (i.e. CV) universe. Implement
// "Get" functions for all the quantities that you need for your analysis.
//
// This class inherits from PU::MinervaUniverse, which in turn inherits from
// PU::BaseUniverse. PU::BU defines the interface with anatuples.
// 
// Within the class, "WeightFunctions" and "MuonFunctions" are included to gain
// access to standardized weight and muon variable getters. See:
// https://cdcvs.fnal.gov/redmine/projects/minerva-sw/wiki/MinervaUniverse_Structure_
// for a full list of standardized functions you can use. In general, if a
// standard version of a function is available, you should be using it.
// =============================================================================
#ifndef CVUNIVERSE_H
#define CVUNIVERSE_H

#include <iostream>

#include "PlotUtils/MinervaUniverse.h"

class CVUniverse : public PlotUtils::MinervaUniverse {

public:
#include "PlotUtils/MuonFunctions.h" // GetMinosEfficiencyWeight
#include "PlotUtils/TruthFunctions.h" //Getq3True
    // ========================================================================
    // Constructor/Destructor
    // ========================================================================
    CVUniverse(PlotUtils::ChainWrapper* chw, double nsigma = 0)
	: PlotUtils::MinervaUniverse(chw, nsigma) {}

    virtual ~CVUniverse() {}

    // ========================================================================
    // Quantities defined here as constants for the sake of below. Definition
    // matched to Dan's CCQENuInclusiveME variables from:
    // `/minerva/app/users/drut1186/cmtuser/Minerva_v22r1p1_OrigCCQENuInc/Ana/CCQENu/ana_common/include/CCQENuUtils.h`
    // ========================================================================
    // C++ note: `constexpr` must be evaluated at compile time.
    //           `const` evaluation can be deferred until runtime.
    // ========================================================================
    static constexpr double M_n = 939.56536;
    static constexpr double M_p = 938.272013;
    static constexpr double M_nucleon = (1.5*M_n+M_p)/2.5;

    static constexpr int PDG_n = 2112;
    static constexpr int PDG_p = 2212;

    static constexpr double rad2deg = 180. / TMath::Pi();
    static constexpr double deg2rad = 1. / rad2deg;

    // ========================================================================
    // Write a "Get" function for all quantities access by your analysis.
    // For composite quantities (e.g. Enu) use a calculator function.
    //
    // In order to properly calculate muon variables and systematics use the
    // various functions defined in MinervaUniverse.
    // E.g. GetPmu, GetEmu, etc.
    // ========================================================================

    
    //============================================
    // My accessor functions - John
    //============================================

    double GetThetamuDeg() const {
	return GetThetamu() * rad2deg;
    }

    double GetThetalepTrueDeg() const {
	return GetThetalepTrue() * rad2deg;
    }

    double GetEnuTrueGeV() const {
	return GetEnuTrue() * 1e-3;
    }

    // Q2 GeV^2
 
    double GetQ2GeV() const
    {
	return GetQ2() / 1e+6;
    }

    double GetQ2TrueGeV() const
    {
	return GetTrueExperimentersQ2()/(1000. * 1000.);
    }

    // must restructure pion variables to query as vector element!
    // iterate over hadron candidates first and get leading hadron
    // i.e. pion energy

    int GetLeadingHadronIndex() const
    {
	const int len = ( GetVecDouble( "MasterAnaDev_pion_E" ) ).size( );
	double E = -1.0; int idx = -1;

	for( Int_t i = 0; i < len; i++ ){
	    if( GetVecElem( "MasterAnaDev_pion_E", i ) > E ){
		idx = i; E = GetVecElem( "MasterAnaDev_pion_E", i );
	    }
	}

	return idx;
    }

    // alternative to leading hadron index is the best LLR candidate
    int GetBestPionIndex() const
    {
	const int len = ( GetVecDouble( "MasterAnaDev_hadron_piFit_scoreLLR" ) ).size( );
	double LLR = -999.0; int idx = -1;

	for( Int_t i = 0; i < len; i++ ){
	    if( GetVecElem( "MasterAnaDev_hadron_piFit_scoreLLR", i ) > LLR ){
		idx = i; LLR = GetVecElem( "MasterAnaDev_hadron_piFit_scoreLLR", i );
	    }
	}

	return idx;
    }

    // where is primary hadron going

    double HadronIsExitingID() const
    {
	//const int leadIdx = GetLeadingHadronIndex();
	const int leadIdx = GetBestPionIndex();
	return GetVecElem( (GetAnaToolName()+"_hadron_isExiting").c_str(), leadIdx );
    }

    // Topology

    double GetPionScore() const
    {
	return GetDouble((GetAnaToolName()+"_pion_score").c_str());
    }

    double GetProtonScore() const
    {
	return GetDouble((GetAnaToolName()+"_proton_score").c_str());
    }

    double GetHadronLLR() const
    {
	//const int leadIdx = GetLeadingHadronIndex();
	const int leadIdx = GetBestPionIndex();
	return GetVecElem((GetAnaToolName()+"_hadron_piFit_scoreLLR").c_str(), leadIdx);
    }

    int GetNSecondaryHadrons() const
    {
	return GetInt((GetAnaToolName()+"_sec_protons_pion_scores_sz").c_str());
    }

    int GetNHadrons() const
    {
	return GetInt( (GetAnaToolName()+"_hadron_number").c_str() );
    }

    int GetNPartTrue(const int PDGCode) const
    {
	int npart = 0;

	std::string pdgBranch = "mc_FSPartPDG";
	for(int i = 0; i < GetNFSPart(); i++){
	    if(GetVecElem(pdgBranch.c_str(), i) == PDGCode){npart++;}
	}
	return npart;
    }

    int GetNPiPlusTrue() const
    {
	const int npi = GetInt( "truth_N_pip" );
	return npi;
    }

    int GetNPiMinusTrue() const
    {
	const int npi = GetInt( "truth_N_pim" );
	return npi;
    }

    int GetNPiZeroTrue() const
    {
	const int npi = GetInt( "truth_N_pi0" );
	return npi;
    }

    int GetNProtonsTrue() const
    {
	const int kProtonPDG = 2212;
	int np = 0;
	
	std::string pdgBranch = "mc_FSPartPDG";
	for(int i = 0; i < GetNFSPart(); i++){
	    if(GetVecElem(pdgBranch.c_str(), i) == kProtonPDG) np++;
	}
	return np;
    }

    int GetNNeutral() const
    {
	const int kNeutronPDG = 2112;
	const int kP0PDG      = 111;
	const int kK0PDG      = 311;
	const int kK0LPDG     = 130;
	const int kK0SPDG     = 310;
	const int kGAMPDG     = 22;
	int nneut = 0;
	
	std::string pdgBranch = "mc_FSPartPDG";
	for(int i = 0; i < GetNFSPart(); i++){
	    int tmpPDG = GetVecElem(pdgBranch.c_str(), i);
	    switch(tmpPDG){
		case kNeutronPDG: case kP0PDG: case kK0PDG: case kK0LPDG: case kK0SPDG: case kGAMPDG: nneut++; break;
		default: if(false){bool dd = true;}
	    }
	}
	return nneut;
    }

    double GetNNegativeHadrons() const
    {
	const int kPDG_pim = 211;
	const int kPDG_Km  = 321;

	int nmin = 0;
	std::string pdgBranch = "mc_FSPartPDG";
	for(int i = 0; i < GetNFSPart(); i++){
	    int tmpPDG = GetVecElem(pdgBranch.c_str(), i);
	    switch(tmpPDG){
		case kPDG_pim: case kPDG_Km: nmin++; break;
		default: if(false){bool dd = true;}
	    }
	}
	return nmin;
    }

    double GetNNeutronClusters() const
    {
	return GetInt("NneutronClusters3"); // least aggressive cut?
	//return GetInt((GetAnaToolName()+"_EvtHasNBlobIncTracks").c_str());
    }

    int GetNNeutronsTrue() const
    {
	const int kNeutronPDG = 2112;
	
	int nneut = 0;
	std::string pdgBranch = "mc_FSPartPDG";
	for(int i = 0; i < GetNFSPart(); i++){
	    if(GetVecElem(pdgBranch.c_str(), i) == kNeutronPDG) nneut++;
	}
	return nneut;
    }

    int GetNFSPart() const
    {
	return GetInt("mc_nFSPart");
    }

    int GetPDGFromRecord(const int pos) const
    {
	if(pos >= GetNFSPart()) return -1;

	std::string pdgBranch = "mc_FSPartPDG";
	return GetVecElem(pdgBranch.c_str(), pos);
    }

    // reco

    double GetEVtx() const // MeV
    {
	double EVTX = GetDouble("vtx_blobs_energy");
	double EISO = GetDouble("vtx_iso_blobs_energy_outside_radius"); // 150 mm by default
	return EVTX - EISO;
    }

    double GetERecoilVtx0mm() const // MeV
    {
	double EVNM = GetDouble("recoil_energy_nonmuon_vtx0mm");
	return EVNM;
    }

    double GetERecoilVtx50mm() const // MeV
    {
	double EVNM = GetDouble("recoil_energy_nonmuon_vtx50mm");
	return EVNM;
    }

    double GetERecoilVtx100mm() const // MeV
    {
	double EVNM = GetDouble("recoil_energy_nonmuon_vtx100mm");
	return EVNM;
    }

    double GetERecoilVtx150mm() const // MeV
    {
	double EVNM = GetDouble("recoil_energy_nonmuon_vtx150mm");
	return EVNM;
    }
    
    double GetERecoilVtx200mm() const // MeV
    {
	double EVNM = GetDouble("recoil_energy_nonmuon_vtx200mm");
	return EVNM;
    }

    double GetERecoilVtx250mm() const // MeV
    {
	double EVNM = GetDouble("recoil_energy_nonmuon_vtx250mm");
	return EVNM;
    }

    double GetERecoilVtx300mm() const // MeV
    {
	double EVNM = GetDouble("recoil_energy_nonmuon_vtx300mm");
	return EVNM;
    }

    double GetEhad() const //MeV
    {
	//double hadE = GetVecElem((GetAnaToolName()+"_HadronE").c_str(), 3);
	// VITAL: THIS IS NEUTRON ENERGY. NOT PROTON!!!!!
	double hadE = GetDouble((GetAnaToolName()+"_proton_E_fromdEdx").c_str());
	//double hadE = GetDouble((GetAnaToolName()+"_proton_calib_energy").c_str());
        return hadE;
    }

    double GetEhadGeV() const
    {
	return GetEhad() / 1e+3;
    }

    double GetPhad(const double Mhad) const //MeV
    {
	double Ehad = GetEhad();
	return sqrt(Ehad*Ehad - Mhad*Mhad); 
    }

    double GetPhadGeV(const double Mhad) const // Mhad in MeV
    {
	return GetPhad(Mhad) / 1e+3;
    }

    double GetHadronThetaX() const {
	std::string thetaBranch = GetAnaToolName() + "_proton_thetaX"; //means "hadron".
	return GetDouble(thetaBranch.c_str()) + GetBeamAngleOffsetX();
    }
    
    double GetHadronThetaY() const {
	std::string thetaBranch = GetAnaToolName() + "_proton_thetaY"; //means "hadron".
	return GetDouble(thetaBranch.c_str()) + GetBeamAngleOffsetY();
    }

    double GetHadronTheta() const {
	//return theta3D(GetHadronThetaX(), GetHadronThetaY());
	return GetDouble((GetAnaToolName() + "_proton_theta").c_str());
    }

    double GetHadronPhi() const {
	return GetDouble((GetAnaToolName() + "_proton_phi").c_str());
    }

    double GetHadronThetaDeg() const {
	return GetHadronTheta() * rad2deg;
    }
    
    ROOT::Math::PxPyPzEVector GetHadron4V(const double Mhad) const {
	double Phad = GetPhad(Mhad);
	double theta = GetHadronTheta();
	double phi = GetHadronPhi();
	double px = Phad * std::sin(theta) * std::cos(phi);
	double py = Phad * std::sin(theta) * std::sin(phi);
	double pz = Phad * std::cos(theta);
	double Ehad = GetEhad();
	return ROOT::Math::PxPyPzEVector(px, py, pz, Ehad);
    }

    // pion branches here

    double GetEpiCorr( ) const // MeV, passive-material-corrected
    {
	//const int leadIdx = GetLeadingHadronIndex();
	const int leadIdx = GetBestPionIndex();
	return GetVecElem( ( GetAnaToolName()+"_hadron_pion_E_corr" ).c_str(), leadIdx );
    }

    double GetEpiCorrGeV( ) const
    {
	return GetEpiCorr( ) * 1.0e-3;
    }

    double GetPxpiCorr( ) const // wrt nu beam
    {
	//const int leadIdx = GetLeadingHadronIndex();
	/* Broken branches, gotta fix
	const int leadIdx = GetBestPionIndex();
	return GetVecElem( ( GetAnaToolName()+"_hadron_pion_px_corr" ).c_str(), leadIdx );
	*/

	return GetPpiCorr() * GetPxpiWrtBeam() / GetPpi();
    }
    
    double GetPypiCorr( ) const // wrt nu beam
    {
	//const int leadIdx = GetLeadingHadronIndex();
	/* broken branches
	const int leadIdx = GetBestPionIndex();
	return GetVecElem( ( GetAnaToolName()+"_hadron_pion_py_corr" ).c_str(), leadIdx );
	*/

	return GetPpiCorr() * GetPypiWrtBeam() / GetPpi();
    }

    double GetPzpiCorr( ) const // wrt nu beam
    {
	//const int leadIdx = GetLeadingHadronIndex();
	/* broken branches
	const int leadIdx = GetBestPionIndex();
	return GetVecElem( ( GetAnaToolName()+"_hadron_pion_pz_corr" ).c_str(), leadIdx );
	*/
	return GetPpiCorr() * GetPzpiWrtBeam() / GetPpi();
    }

    double GetPpiCorr( ) const // MeV
    {
	//const int leadIdx = GetLeadingHadronIndex();
	/* broken branches
	const int leadIdx = GetBestPionIndex();
	return GetVecElem( ( GetAnaToolName()+"_hadron_pion_p_corr" ).c_str(), leadIdx );
	*/

	return std::sqrt( GetEpiCorr() * GetEpiCorr() - MinervaUnits::M_pion * MinervaUnits::M_pion );
    }

    double GetThetapi( ) const // rad, wrt nu beam
    {
	//const int leadIdx = GetLeadingHadronIndex();
	const int leadIdx = GetBestPionIndex();
	return GetVecElem( ( GetAnaToolName()+"_pion_theta" ).c_str(), leadIdx );
    }

    double GetPhipi( ) const // rad, wrt nu beam
    {
	//const int leadIdx = GetLeadingHadronIndex();
	const int leadIdx = GetBestPionIndex();
	return GetVecElem( ( GetAnaToolName()+"_pion_phi" ).c_str(), leadIdx );
    }

    // let's do some muon branches for inspection
    double GetEmuMAD( ) const
    {
	return GetDouble( "MasterAnaDev_muon_E" );
    }

    double GetPxmu( ) const
    {
	return GetDouble( "MasterAnaDev_muon_Px" );
    }

    double GetPymu( ) const
    {
	return GetDouble( "MasterAnaDev_muon_Py" );
    }

    double GetPzmu( ) const
    {
	return GetDouble( "MasterAnaDev_muon_Pz" );
    }

    double GetPxmuWrtNuBeam() const
    {
	return GetPxmu( );
    }

    double GetPymuWrtNuBeam() const
    {
	return GetPymu( ) * std::cos( MinervaUnits::numi_beam_angle_rad )
	    - GetPzmu( ) * std::sin( MinervaUnits::numi_beam_angle_rad );
    }

    double GetPzmuWrtNuBeam() const
    {
	return GetPymu( ) * std::sin( MinervaUnits::numi_beam_angle_rad )
	    + GetPzmu( ) * std::cos( MinervaUnits::numi_beam_angle_rad );
    }

    // uncorrected branches here... safe with MAD 1.53

    // change this branch for different definitions of Epi
    double GetEpi( ) const // MeV
    {
	return GetEpiRaw();
    }
    
    double GetEpiRaw( ) const // MeV
    {
	// enforce selection on ~*leading* hadron~ best pion candidate!
	//const int leadIdx = GetLeadingHadronIndex();
	const int leadIdx = GetBestPionIndex();
	return GetVecElem( ( GetAnaToolName() + "_pion_E" ).c_str(), leadIdx ); // dEdx, bad!

	//return GetDouble( ( GetAnaToolName() + "_recoil_E" ).c_str() ); // better estimator!
	//return GetDouble( ( GetAnaToolName() + "_hadron_recoil_two_track" ).c_str() );

	//return GetEpiMatchedTrue( ); // to see what happens!
    }

    double GetEpiRawScaled( ) const // MeV, calorimetric scaling factor applied here
    {
	const double alpha = 2.11609; // from getScaleFactor
	return GetEpiRaw() * alpha;
    }

    double GetEpiCorrScaled( ) const // MeV
    {
	const double alpha = 2.07734;
	return GetEpiCorr() * alpha;
    }

    double GetEpiGeV( ) const
    {
	return GetEpi( ) * 1.0e-3;
    }

    double GetPxpi( ) const // wrt MINERvA -- rotate by numi_beam_angle_rad for nu!
    {
	//const int leadIdx = GetLeadingHadronIndex();
	const int leadIdx = GetBestPionIndex();
	//return GetVecElem( ( GetAnaToolName() + "_pion_Px" ).c_str(), leadIdx );

	const double px = GetVecElem( ( GetAnaToolName() + "_pion_Px" ).c_str(), leadIdx );
	const double p3 = GetVecElem( ( GetAnaToolName() + "_pion_P" ).c_str(), leadIdx );

	return GetPpi() * px/p3;
    }
    
    double GetPypi( ) const // wrt MINERvA -- rotate by numi_beam_angle_rad for nu!
    {
	//const int leadIdx = GetLeadingHadronIndex();
	const int leadIdx = GetBestPionIndex();
	//return GetVecElem( ( GetAnaToolName() + "_pion_Py" ).c_str(), leadIdx );


	const double py = GetVecElem( ( GetAnaToolName() + "_pion_Py" ).c_str(), leadIdx );
	const double p3 = GetVecElem( ( GetAnaToolName() + "_pion_P" ).c_str(), leadIdx );

	return GetPpi() * py/p3;
    }

    double GetPzpi( ) const // wrt MINERvA -- rotate by numi_beam_angle_rad for nu!
    {
	//const int leadIdx = GetLeadingHadronIndex();
	const int leadIdx = GetBestPionIndex();
	//return GetVecElem( ( GetAnaToolName() + "_pion_Pz" ).c_str(), leadIdx );

	const double pz = GetVecElem( ( GetAnaToolName() + "_pion_Pz" ).c_str(), leadIdx );
	const double p3 = GetVecElem( ( GetAnaToolName() + "_pion_P" ).c_str(), leadIdx );

	return GetPpi() * pz/p3;
    }

    double GetPpi( ) const // MeV
    {
	//const int leadIdx = GetLeadingHadronIndex();
	//return GetVecElem( ( GetAnaToolName() + "_pion_P" ).c_str(), leadIdx );

	const double Epi = GetEpi( );
	const double mpi = MinervaUnits::M_pion;
	
	return std::sqrt( Epi * Epi - mpi * mpi );
    }

    // let's explicitly rotate this on the yz plane
    double GetPxpiWrtBeam( ) const
    {
	return GetPxpi( ); //same component
    }

    double GetPypiWrtBeam( ) const
    {
	return GetPypi( ) * TMath::Cos( MinervaUnits::numi_beam_angle_rad )
	    - GetPzpi( ) * TMath::Sin( MinervaUnits::numi_beam_angle_rad );
	    
    }

    double GetPzpiWrtBeam( ) const
    {
	return GetPypi( ) * TMath::Sin( MinervaUnits::numi_beam_angle_rad )
	    + GetPzpi( ) * TMath::Cos( MinervaUnits::numi_beam_angle_rad );
    }

    ROOT::Math::PxPyPzEVector GetPion4V( ) const { // wrt nu beam!
	double Px = GetPxpiWrtBeam( );
	double Py = GetPypiWrtBeam( );
	double Pz = GetPzpiWrtBeam( );
	double E = GetEpi( );

	return ROOT::Math::PxPyPzEVector( Px, Py, Pz, E );
    }

    // -- let's add muon truth, shall we?

    double GetEmuMADTrue( ) const
    {
	return GetDouble( "truth_muon_E" );
    }

    double GetPxmuTrue( ) const
    {
	return GetDouble( "truth_muon_px" );
    }

    double GetPymuTrue( ) const
    {
	return GetDouble( "truth_muon_py" );
    }

    double GetPzmuTrue( ) const
    {
	return GetDouble( "truth_muon_pz" );
    }

    double GetPxmuTrueWrtNuBeam( ) const
    {
	return GetPxmuTrue( );
    }

    double GetPymuTrueWrtNuBeam( ) const
    {
	return GetPymuTrue( ) * std::cos( MinervaUnits::numi_beam_angle_rad )
	    - GetPzmuTrue( ) * std::sin( MinervaUnits::numi_beam_angle_rad );
    }

    double GetPzmuTrueWrtNuBeam( ) const
    {
	return GetPymuTrue( ) * std::sin( MinervaUnits::numi_beam_angle_rad )
	    + GetPzmuTrue( ) * std::cos( MinervaUnits::numi_beam_angle_rad );
    }

    // -- let's add pion truth, shall we?

    double GetTpiMatchedTrue( ) const
    {
	//const int leadIdx = GetLeadingHadronIndex();
	const int leadIdx = GetBestPionIndex();
	const int truePDG = GetVecElem( "MasterAnaDev_hadron_tm_PDGCode", leadIdx );
	if( truePDG != 211 ) return -1.0; // not pion!

	return GetVecElem( "MasterAnaDev_hadron_tm_beginKE", leadIdx );
    }

    double GetEpiMatchedTrue( ) const
    {
	const double trueTpi = GetTpiMatchedTrue( );
	if( trueTpi < 0.0 ) return -1.0;

	return trueTpi + MinervaUnits::M_pion;
    }

    double GetPpiMatchedTrue( ) const
    {
	const double trueEpi = GetEpiMatchedTrue( );
	if( trueEpi < 0.0 ) return -1.0;

	const double mpi = MinervaUnits::M_pion;
	return std::sqrt( trueEpi*trueEpi - mpi*mpi );
    }
	
    double GetEpiTrue( ) const
    {
        const char * PDGBranch = "mc_FSPartPDG";
	const char * piEBranch = "mc_FSPartE";
	const char * varBranch = "mc_FSPartE";
	const int pipPDG = 211;
	const int len = (GetVecDouble( PDGBranch ) ).size();
	double maxE = -1.0e+6; int iBest = -1;
	for( Int_t i = 0; i < len; i++ ){ // return leading pion
	    double thisE = GetVecElem( piEBranch, i );
	    int thisPDG = GetVecElem( PDGBranch, i );
	    if( thisE > maxE && thisPDG == pipPDG ){ maxE = thisE; iBest = i; }
	}
	return GetVecElem( varBranch, iBest );
    }

    double GetEpiTrueGeV( ) const
    {
	return GetEpiTrue( ) * 1.0e-3;
    }

    double GetPxpiTrue( ) const // wrt MINERvA -- rotate by numi_beam_angle_rad for nu!
    {
	const char * PDGBranch = "mc_FSPartPDG";
	const char * piEBranch = "mc_FSPartE";
	const char * varBranch = "mc_FSPartPx";
	const int pipPDG = 211;
	const int len = (GetVecDouble( PDGBranch ) ).size();
	double maxE = -1.0e+6; int iBest = -1;
	for( Int_t i = 0; i < len; i++ ){ // return leading pion
	    double thisE = GetVecElem( piEBranch, i );
	    int thisPDG = GetVecElem( PDGBranch, i );
	    if( thisE > maxE && thisPDG == pipPDG ){ maxE = thisE; iBest = i; }
	}
	return GetVecElem( varBranch, iBest );
    }

    double GetPypiTrue( ) const // wrt MINERvA -- rotate by numi_beam_angle_rad for nu!
    {
        const char * PDGBranch = "mc_FSPartPDG";
	const char * piEBranch = "mc_FSPartE";
	const char * varBranch = "mc_FSPartPy";
	const int pipPDG = 211;
	const int len = (GetVecDouble( PDGBranch ) ).size();
	double maxE = -1.0e+6; int iBest = -1;
	for( Int_t i = 0; i < len; i++ ){ // return leading pion
	    double thisE = GetVecElem( piEBranch, i );
	    int thisPDG = GetVecElem( PDGBranch, i );
	    if( thisE > maxE && thisPDG == pipPDG ){ maxE = thisE; iBest = i; }
	}
	return GetVecElem( varBranch, iBest );
    }

    double GetPzpiTrue( ) const // wrt MINERvA -- rotate by numi_beam_angle_rad for nu!
    {
        const char * PDGBranch = "mc_FSPartPDG";
	const char * piEBranch = "mc_FSPartE";
	const char * varBranch = "mc_FSPartPz";
	const int pipPDG = 211;
	const int len = (GetVecDouble( PDGBranch ) ).size();
	double maxE = -1.0e+6; int iBest = -1;
	for( Int_t i = 0; i < len; i++ ){ // return leading pion
	    double thisE = GetVecElem( piEBranch, i );
	    int thisPDG = GetVecElem( PDGBranch, i );
	    if( thisE > maxE && thisPDG == pipPDG ){ maxE = thisE; iBest = i; }
	}
	return GetVecElem( varBranch, iBest );
    }

    double GetPpiTrue( ) const
    {
	return TMath::Sqrt(
	    std::pow( GetPxpiTrue( ), 2.0)
	    + std::pow( GetPypiTrue( ), 2.0)
	    + std::pow( GetPzpiTrue( ), 2.0)
	    );
    }

    // and explicitly rotate on yz plane

    double GetPxpiTrueWrtBeam( ) const
    {
	return GetPxpiTrue( ); // same component
    }

    double GetPypiTrueWrtBeam( ) const
    {
	return GetPypiTrue( ) * TMath::Cos( MinervaUnits::numi_beam_angle_rad )
	    - GetPzpiTrue( ) * TMath::Sin( MinervaUnits::numi_beam_angle_rad );
    }

    double GetPzpiTrueWrtBeam( ) const
    {
	return GetPypiTrue( ) * TMath::Sin( MinervaUnits::numi_beam_angle_rad )
	    + GetPzpiTrue( ) * TMath::Cos( MinervaUnits::numi_beam_angle_rad );
    }

    ROOT::Math::PxPyPzEVector GetHadron4VTrue(const int PDG) const {
	int idx = GetIndexOfHad(PDG);
	double Px = GetVecElem("mc_FSPartPx",idx);
	double Py = GetVecElem("mc_FSPartPy",idx);
	double Pz = GetVecElem("mc_FSPartPz",idx);
	double E = GetVecElem("mc_FSPartE",idx);

	return ROOT::Math::PxPyPzEVector(Px,Py,Pz,E);
    }

    ROOT::Math::PxPyPzEVector GetPion4VTrue( const int PDG ) const {
	const int pip_PDG = 211;
	return GetHadron4VTrue( pip_PDG );
    }

    double GetMuHadAngle(const double Mhad) const{
	ROOT::Math::PxPyPzEVector p4mu  = GetMuon4V();
	ROOT::Math::PxPyPzEVector p4had = GetHadron4V(Mhad);

	double mupx = p4mu.Px(); double hadpx = p4had.Px();
	double mupy = p4mu.Py(); double hadpy = p4had.Py();
	double mupz = p4mu.Pz(); double hadpz = p4had.Pz();
	double mup  = p4mu.P();  double hadp  = p4had.P();

	double num = mupx*hadpx + mupy*hadpy + mupz*hadpz;
	double den = mup*hadp;

	return std::acos(num/den);
    }

    // and the same func *without arguments* for Variable declaration
    double GetThetaMuHadVar() const {
	const double Mhad = ( GetHadronLLR() >= 0.0 ) ? MinervaUnits::M_pion : MinervaUnits::M_p;
	return GetMuHadAngle(Mhad) * 180. / TMath::Pi();
    }

    // simpler |t| calc using A. Mislivec's Eq. (8.6)
    // |t| = |(p_\nu - p_\mu - p_\pi)^2|

    double GetAbsT() const { //MeV^2
	ROOT::Math::PxPyPzEVector mu4V = GetMuon4V(); // wrt nu direction, i.e. along z = z_nu
	//const double Mhad = ( GetHadronLLR() >= 0.0 ) ? MinervaUnits::M_pion : MinervaUnits::M_p ;
	//ROOT::Math::PxPyPzEVector had4V = GetHadron4V(Mhad);
	ROOT::Math::PxPyPzEVector pi4V = GetPion4V();
	const double Enu = GetEnu();
	const double py = 0.0; const double pz = Enu;
	ROOT::Math::PxPyPzEVector nu4V(0.,py,pz,Enu);

	ROOT::Math::PxPyPzEVector sys4V = nu4V - mu4V - pi4V;
	return std::abs(sys4V.M2());
    }

    // and a resolution

    double GetResAbsT() const {
	const double treco = GetAbsT(); // MeV^2
	const double ttrue = GetAbsTTrue(); // MeV^2

	return ( treco - ttrue ) / ttrue;
    }

    // some diagnostics for sys4V
    double GetSys4VE() const { // MeV
	const double muE = GetEmuMAD();
	const double piE = GetEpi();
	const double nuE = GetEnu();
	return nuE - muE - piE;
    }

    double GetSys4VPx() const { // MeV
	const double muPx = GetPxmuWrtNuBeam();
	const double piPx = GetPxpiWrtBeam();
	const double nuPx = 0.0;
	return nuPx - muPx - piPx;
    }

    double GetSys4VPy() const { // MeV
	const double muPy = GetPymuWrtNuBeam();
	const double piPy = GetPypiWrtBeam();
	const double nuPy = 0.0;
	return nuPy - muPy - piPy;
    }

    double GetSys4VPz() const { // MeV
	const double muPz = GetPzmuWrtNuBeam();
	const double piPz = GetPzpiWrtBeam();
	const double nuPz = GetEnu();
	return nuPz - muPz - piPz;
    }

    double GetAbsTGeV() const {
	return GetAbsT() / 1e+6;
    }
    
    // now construct |t| using Alex's Eq. (6.1)
    double GetAlexAbsT() const { //MeV^2
	//const double Mhad = ( GetHadronLLR() >= 0.0 ) ?  MinervaUnits::M_pion : MinervaUnits::M_proton; //need to declare in func
	    
	const double Emu      = GetEmuMAD();
	//double Ehad     = GetEhad();
	const double Epi      = GetEpi();
	/*
	double Pmu      = GetPmu();
	double Phad     = GetPhad(Mhad);
	double thetamu  = GetThetamu();
	double thetahad = GetHadronTheta();
	double thetaMuHad = GetMuHadAngle(Mhad);
	*/

	const double Pxmu = GetPxmuWrtNuBeam();
	const double Pymu = GetPymuWrtNuBeam();
	const double Pzmu = GetPzmuWrtNuBeam();
	const double Pmu = std::sqrt( Pxmu*Pxmu + Pymu*Pymu + Pzmu*Pzmu );

	const double Pxpi = GetPxpiWrtBeam();
	const double Pypi = GetPypiWrtBeam();
	const double Pzpi = GetPzpiWrtBeam();
	const double Ppi = std::sqrt( Pxpi*Pxpi + Pypi*Pypi + Pzpi*Pzpi );

	const double inprodNum = Pxmu*Pxpi + Pymu*Pypi + Pzmu*Pzpi;
	const double inprodDen = Pmu * Ppi;
	const double thetaMuPi = TMath::ACos( inprodNum / inprodDen );
	
	const double Mmu = MinervaUnits::M_mu;
	const double Mpi = MinervaUnits::M_pion;

	double sumE = Emu + Epi;
	double term1 = -2. * sumE * ( Emu - Pzmu );
	double term2 = Mmu * Mmu;
	double term3 = -2. * ( Epi*Epi - sumE*Pzpi + inprodNum );
	double term4 = Mpi * Mpi;

	return std::abs(term1 + term2 + term3 + term4);
    }

    double GetAlexAbsTGeV() const {
	return 1.0e-6 * GetAlexAbsT();
    }

    // calculator function for Enu under COH
    // Equation (3.7) in Alex's thesis

    double GetCOHEnu() const {
	return GetEmu() + GetEpi();
    }
    
    double GetCOHEnuGeV() const {
	return GetCOHEnu() * 1e-3;
    }

    double GetCOHQ2() const { // MeV^2 for consistency
	return calcCOHQ2() * 1e+6;
    }

    double GetCOHQ2GeV() const { // GeV^2
	return calcCOHQ2();
    }

    // calc COH Q2
    double calcCOHQ2() const { // GeV^2
	const double mmu = MinervaUnits::M_mu * 1e-3; //GeV
	double Emu     = GetEmuGeV(); // GeV
	double Pmu     = std::sqrt( Emu*Emu - mmu*mmu ); // GeV
	double Pzmu    = GetPzmuWrtNuBeam() * 1.0e-3; // GeV
	double thetamu = TMath::ACos( Pzmu / Pmu );
	double Enu     = GetCOHEnuGeV();

	return 2. * Enu * (Emu - Pmu * std::cos(thetamu)) - mmu * mmu;
    }

    // get true COH Q2
    double GetCOHQ2TrueGeV() const { // GeV^2
	const double Enu  = GetEnuTrue() * 1.0e-3; //GeV
	const double Emu  = GetEmuMADTrue() * 1.0e-3; //GeV
	const double mmu  = MinervaUnits::M_mu * 1.0e-3; //GeV
	
	const double pzmu = GetPzmuTrueWrtNuBeam( ) * 1.0e-3; //GeV
	const double pmu  = std::sqrt( Emu*Emu - mmu*mmu ); // GeV
	const double theta = TMath::ACos( pzmu / pmu );

	return 2. * Enu * ( Emu - pzmu ) - mmu * mmu;
    }

    //better definition of Q2 for higher |t|
    //see docDB 16281 slide 42, Anne Norrick
    //NB Anne is missing a prefactor 4, see GetTrueExperimentersQ2
    double GetQ2() const { //MeV^2
	const double Emu = GetEmu();
	const double Enu = GetEnu();
	const double thetamu = GetThetamu();
	return CalcTrueExperimentersQ2(Enu, Emu, thetamu);
    }

    double GetW() const { //needed for sidebands, MeV
	std::string Wbranch = GetAnaToolName() + "_W";
	return GetDouble(Wbranch.c_str());
    }

    double GetWGeV() const {
	return GetW() * 1e-3;
    }

    double GetEnu() const {
	return GetEmu() + GetEpi(); //good assuming measurement of Ehad is good
    }

    double GetEnuGeV() const {
	return GetEnu() * 1e-3;
    }

    // truth

    int GetIndexOfHad(const int PDG) const
    { // mc_FSPart{} stuff is a vector of all FS particles. Search for
	// correct PDG and return location of 1st entry.
	for(int idx = 0; idx < GetInt("mc_nFSPart"); idx++){
	    if(GetVecElem("mc_FSPartPDG",idx) == PDG) return idx;
	}
	return -1;
    }

    double GetEhadTrue() const // MeV
    {
	const int PDG_pi = 211;
	const int PDG_p  = 2212;
	const int PDG = (GetNPiPlusTrue() > GetNProtonsTrue()) ? PDG_pi : PDG_p;
	int idx = GetIndexOfHad(PDG);
	return GetVecElem("mc_FSPartE",idx);
    }

    double GetEhadTrueGeV() const
    {
	return GetEhadTrue() / 1e+3;
    }

    double GetPhadTrue() const // MeV
    {
	const int PDG_pi = 211;
	const int PDG_p  = 2212;
	const int PDG = (GetNPiPlusTrue() > GetNProtonsTrue()) ? PDG_pi : PDG_p;
	int idx = GetIndexOfHad(PDG);
	
	double px = GetVecElem("mc_FSPartPx",idx);
	double py = GetVecElem("mc_FSPartPy",idx);
	double pz = GetVecElem("mc_FSPartPz",idx);
	return sqrt(px*px + py*py + pz*pz);
    }

    double GetPhadTrueGeV() const
    {
	return GetPhadTrue() / 1e+3;
    }

    double GetHadronThetaTrueDeg() const { //for variable declaration
	const int PDG_pi = 211;
	const int PDG_p  = 2212;
	const int PDG = (GetNPiPlusTrue() > GetNProtonsTrue()) ? PDG_pi : PDG_p;
	return GetThetaHadTrue(PDG) * rad2deg;
    }

    // truth theta calcs consistent with TruthFunctions.h
    double GetThetaHadTrue(const int PDG = 211) const
    {
	int idx = GetIndexOfHad(PDG);
	double px = GetVecElem("mc_FSPartPx", idx);
	double py = GetVecElem("mc_FSPartPy", idx);
	double pz = GetVecElem("mc_FSPartPz", idx);
	TVector3 p3had(px,py,pz);
	p3had.RotateX(MinervaUnits::numi_beam_angle_rad);
	return p3had.Theta();
    }

    double GetPhiHadTrue(const int PDG = 211) const
    {
	int idx = GetIndexOfHad(PDG);
	double px = GetVecElem("mc_FSPartPx", idx);
	double py = GetVecElem("mc_FSPartPy", idx);
	double pz = GetVecElem("mc_FSPartPz", idx);
	TVector3 p3had(px,py,pz);
	p3had.RotateX(MinervaUnits::numi_beam_angle_rad);
	return p3had.Phi();
    }

    double GetMuHadAngleTrue(const int PDG = 211) const
    {
	double mupx = GetVecElem("mc_primFSLepton", 0);
	double mupy = GetVecElem("mc_primFSLepton", 1);
	double mupz = GetVecElem("mc_primFSLepton", 2);
	double Pmu  = GetPlepTrue();

	int idx = GetIndexOfHad(PDG);
	double hadpx = GetVecElem("mc_FSPartPx", idx);
	double hadpy = GetVecElem("mc_FSPartPy", idx);
	double hadpz = GetVecElem("mc_FSPartPz", idx);
	double Phad  = GetPhadTrue();

	double num = mupx * hadpx + mupy * hadpy + mupz * hadpz;
	double den = Pmu * Phad;

	return std::acos(num / den);
    }

    double GetThetaMuPi() const {
	const double pxm = GetPxmuWrtNuBeam();
	const double pym = GetPymuWrtNuBeam();
	const double pzm = GetPzmuWrtNuBeam();
	const double pm  = std::sqrt( pxm*pxm + pym*pym + pzm*pzm );

	const double pxp = GetPxpiWrtBeam();
	const double pyp = GetPypiWrtBeam();
	const double pzp = GetPzpiWrtBeam();
	const double pp  = std::sqrt( pxp*pxp + pyp*pyp + pzp*pzp );

	const double num = pxm*pxp + pym*pyp + pzm*pzp;
	const double den = pm*pp;

	return TMath::ACos( num / den );
    }

    double GetThetaMuPiDeg() const {
	return GetThetaMuPi() * rad2deg;
    }

    double GetThetaMuPiTrue() const {
	const double pxm = GetPxmuTrueWrtNuBeam();
	const double pym = GetPymuTrueWrtNuBeam();
	const double pzm = GetPzmuTrueWrtNuBeam();
	const double pm  = std::sqrt( pxm*pxm + pym*pym + pzm*pzm );

	const double pxp = GetPxpiTrueWrtBeam();
	const double pyp = GetPypiTrueWrtBeam();
	const double pzp = GetPzpiTrueWrtBeam();
	const double pp  = std::sqrt( pxp*pxp + pyp*pyp + pzp*pzp );

	const double num = pxm*pxp + pym*pyp + pzm*pzp;
	const double den = pm*pp;

	return TMath::ACos( num / den );
    }

    double GetThetaMuPiTrueDeg() const {
	return GetThetaMuPiTrue() * rad2deg;
    }

    double GetMMuPi() const { // MeV
	const double pxm = GetPxmuWrtNuBeam();
	const double pym = GetPymuWrtNuBeam();
	const double pzm = GetPzmuWrtNuBeam();
	const double Em  = GetEmuMAD();

	const double pxp = GetPxpiWrtBeam();
	const double pyp = GetPypiWrtBeam();
	const double pzp = GetPzpiWrtBeam();
	const double Ep  = GetEpi();

	ROOT::Math::PxPyPzEVector mu4V( pxm, pym, pzm, Em );
	ROOT::Math::PxPyPzEVector pi4V( pxp, pyp, pzp, Ep );
	ROOT::Math::PxPyPzEVector mp4V = mu4V + pi4V;
	return mp4V.M();
    }

    double GetMMuPiTrue() const { // MeV
	const double pxm = GetPxmuTrueWrtNuBeam();
	const double pym = GetPymuTrueWrtNuBeam();
	const double pzm = GetPzmuTrueWrtNuBeam();
	const double Em  = GetEmuMADTrue();

	const double pxp = GetPxpiTrueWrtBeam();
	const double pyp = GetPypiTrueWrtBeam();
	const double pzp = GetPzpiTrueWrtBeam();
	const double Ep  = GetEpiTrue();

	ROOT::Math::PxPyPzEVector mu4V( pxm, pym, pzm, Em );
	ROOT::Math::PxPyPzEVector pi4V( pxp, pyp, pzp, Ep );
	ROOT::Math::PxPyPzEVector mp4V = mu4V + pi4V;
	return mp4V.M();
    }

    //calc |t|
    /*
    double GetAbsTTrue() const { // MeV^2
	ROOT::Math::PxPyPzEVector mu4V = GetLep4VTrue();
	const int PDG_pi = 211;
	const int PDG_p  = 2212;
	const int PDG = (GetNPiPlusTrue() > GetNProtonsTrue()) ? PDG_pi : PDG_p;
	ROOT::Math::PxPyPzEVector had4V = GetHadron4VTrue(PDG);
	const double Enu = GetEnuTrue();
	//const double py = Enu * TMath::Sin( MinervaUnits::numi_beam_angle_rad ); // beam points down by 58 mrad
	//const double pz = Enu * TMath::Cos( MinervaUnits::numi_beam_angle_rad );
	const double py = 0.0; const double pz = Enu;
	ROOT::Math::PxPyPzEVector nu4V(0.,py,pz,Enu);

	ROOT::Math::PxPyPzEVector sys4V = nu4V - mu4V - had4V;
	return std::abs(sys4V.M() * sys4V.M());
    }
    */

    double GetAbsTTrue() const { // MeV^2
	//ROOT::Math::PxPyPzEVector mu4V = GetLep4VTrue();
	const double Emu = GetEmuMADTrue();
	const double pmx = GetPxmuTrueWrtNuBeam();
	const double pmy = GetPymuTrueWrtNuBeam();
	const double pmz = GetPzmuTrueWrtNuBeam();
	ROOT::Math::PxPyPzEVector mu4V( pmx, pmy, pmz, Emu );
	
	const double Enu = GetEnuTrue();
	ROOT::Math::PxPyPzEVector nu4V(0.,0.,Enu,Enu);

	const double Epi = GetEpiTrue();
	const double px = GetPxpiTrueWrtBeam();
	const double py = GetPypiTrueWrtBeam();
	const double pz = GetPzpiTrueWrtBeam();
	ROOT::Math::PxPyPzEVector pi4V( px, py, pz, Epi );

	ROOT::Math::PxPyPzEVector sys4V = nu4V - mu4V - pi4V;
	return std::abs( sys4V.M() * sys4V.M() );
    }

    // some diagnostics for sys4V
    double GetSys4VETrue() const { // MeV
	const double muE = GetEmuMADTrue();
	const double piE = GetEpiTrue();
	const double nuE = GetEnuTrue();
	return nuE - muE - piE;
    }

    double GetSys4VPxTrue() const { // MeV
	const double muPx = GetPxmuTrueWrtNuBeam();
	const double piPx = GetPxpiTrueWrtBeam();
	const double nuPx = 0.0;
	return nuPx - muPx - piPx;
    }

    double GetSys4VPyTrue() const { // MeV
	const double muPy = GetPymuTrueWrtNuBeam();
	const double piPy = GetPypiTrueWrtBeam();
	const double nuPy = 0.0;
	return nuPy - muPy - piPy;
    }

    double GetSys4VPzTrue() const { // MeV
	const double muPz = GetPzmuTrueWrtNuBeam();
	const double piPz = GetPzpiTrueWrtBeam();
	const double nuPz = GetEnuTrue();
	return nuPz - muPz - piPz;
    }

    double GetAbsTTrueGeV() const {
	return 1e-6 * GetAbsTTrue(); // truth is 100x larger than reco!!!!
    }

    double GetAlexAbsTTrue() const { //MeV^2
	//const int PDG_pi = 211; //need to declare within function scope
	//const int PDG_p  = 2212;
	//const int PDG = (GetNPiPlusTrue() > GetNProtonsTrue()) ? PDG_pi : PDG_p;

	/*
	double Emu        = GetElepTrue();
	double Ehad       = GetEhadTrue();
	double Pmu        = GetPlepTrue();
	double Phad       = GetPhadTrue();
	double thetamu    = GetThetalepTrue();
	double thetahad   = GetThetaHadTrue(PDG);
	double thetaMuHad = GetMuHadAngleTrue(PDG); 
	*/

	const double Emu = GetEmuMADTrue();
	const double pmx = GetPxmuTrueWrtNuBeam();
	const double pmy = GetPymuTrueWrtNuBeam();
	const double pmz = GetPzmuTrueWrtNuBeam();
	const double pm  = std::sqrt( pmx*pmx + pmy*pmy + pmz*pmz );

	const double Epi = GetEpiTrue();
	const double px  = GetPxpiTrueWrtBeam();
	const double py  = GetPypiTrueWrtBeam();
	const double pz  = GetPzpiTrueWrtBeam();
	const double ppi = std::sqrt( px*px + py*py + pz*pz );

	const double inProdNum  = px*pmx + py*pmy + pz*pmz;
	const double inProdDen  = pm*ppi;
	const double thetaMuHad = TMath::ACos( inProdNum / inProdDen );

	const double Mmu = MinervaUnits::M_mu;
	const double Mpi = MinervaUnits::M_pion;

	double sumE = Emu + Epi;
	double term1 = -2. * sumE * ( Emu - pmz );
	double term2 = Mmu * Mmu;
	double term3 = -2. * ( Epi*Epi - sumE*pz + inProdNum );
	double term4 = Mpi * Mpi;

	return std::abs(term1 + term2 + term3 + term4);
    }

    double GetAlexAbsTTrueGeV() const {
	return 1.0e-6 * GetAlexAbsTTrue();
    }

    double GetWTrue() const {
	//return GetTrueExperimentersW();
	return GetDouble("mc_w");
    }

    double GetWTrueGeV() const {
	return GetWTrue() * 1e-3;
    }

    //and mu-had invariant mass
    double GetMuHadW() const { //MeV
	ROOT::Math::PxPyPzEVector mu4V = GetMuon4V();
	const double Mhad = ( GetHadronLLR() >= 0.0 ) ? MinervaUnits::M_pion : MinervaUnits::M_p ;
	ROOT::Math::PxPyPzEVector had4V = GetHadron4V(Mhad);

	ROOT::Math::PxPyPzEVector sys4V = mu4V + had4V;
	return sys4V.M();
    }

    double GetMuHadWGeV() const {
	return GetMuHadW() * 1e-3;
    }

    double GetMuHadWTrue() const {
	ROOT::Math::PxPyPzEVector mu4V = GetLep4VTrue();
	const int PDG_pi = 211;
	const int PDG_p  = 2212;
        int PDG = 0;
	if(GetNPiPlusTrue() == 1) PDG = PDG_pi;
	else if(GetNProtonsTrue() == 1) PDG = PDG_p;
	else return 0;
	
	ROOT::Math::PxPyPzEVector had4V = GetHadron4VTrue(PDG);
	ROOT::Math::PxPyPzEVector sys4V = mu4V + had4V;
	return sys4V.M();
    }

    double GetMuHadWTrueGeV() const {
	return GetMuHadWTrue() * 1e-3;
    }

    //============================================
    // Here be electron reco stuff
    // for non-muon events
    //============================================

    // RETHERE: Want to look at blob_nuefuzz_* branches?
    
    double GetConeEnergyVis() const {
	return GetDouble( "ConeEnergyVis" );
    }

    double GetExtraEnergyVis() const {
	return GetDouble( "ExtraEnergyVis" );
    }

    double GetPsi() const {
	return GetDouble( "Psi" );
    }

    bool GetIsEVertexMismatched() const {
	return GetInt( "HasNoVertexMismatch" ) < 1 ; // encapsulates non-nue candidates @ -1
    }
    
    //============================================
    // Vanilla Mnv101
    //============================================

    // Quantities only needed for cuts
    // Although unlikely, in principle these quanties could be shifted by a
    // systematic. And when they are, they'll only be shifted correctly if we
    // write these accessor functions.
  
    //Muon kinematics
    double GetMuonPT() const //GeV/c
    {
	return GetPmu()/1000. * sin(GetThetamu());
    }

    double GetMuonPz() const //GeV/c
    {
	return GetPmu()/1000. * cos(GetThetamu());
    }

    double GetMuonPTTrue() const //GeV/c
    {
	return GetPlepTrue()/1000. * sin(GetThetalepTrue());
    }

    double GetMuonPzTrue() const //GeV/c
    {
	return GetPlepTrue()/1000. * cos(GetThetalepTrue());
    }

    double GetEmuGeV() const //GeV
    {
	return GetEmuMAD()/1000.;
    }

    double GetElepTrueGeV() const //GeV
    {
	return GetElepTrue()/1000.;
    }

    int GetInteractionType() const {
	return GetInt("mc_intType");
    }

    int GetTargetNucleon() const {
	return GetInt("mc_targetNucleon");
    }
  
    double GetBjorkenXTrue() const {
	return GetDouble("mc_Bjorkenx");
    }

    double GetBjorkenYTrue() const {
	return GetDouble("mc_Bjorkeny");
    }

    virtual bool IsMinosMatchMuon() const {
	return GetInt("has_interaction_vertex") == 1;
    }

    // same as above but more transparent nomenclature
    virtual bool IsSomeInteraction() const {
	return GetInt("has_interaction_vertex") == 1;
    }
  
    ROOT::Math::XYZTVector GetVertex() const
    {
	ROOT::Math::XYZTVector result;
	result.SetCoordinates(GetVec<double>("vtx").data());
	return result;
    }

    ROOT::Math::XYZTVector GetTrueVertex() const
    {
	ROOT::Math::XYZTVector result;
	result.SetCoordinates(GetVec<double>("mc_vtx").data());
	return result;
    }

    virtual int GetTDead() const {
	return GetInt("phys_n_dead_discr_pair_upstream_prim_track_proj");
    }
  
    //TODO: If there was a spline correcting Eavail, it might not really be Eavail.
    //      Our energy correction spline, one of at least 2 I know of, corrects q0
    //      so that we get the right neutrino energy in an inclusive sample.  So,
    //      this function could be correcting for neutron energy which Eavail should
    //      not do.
    virtual double GetEavail() const {
	return GetDouble("recoilE_SplineCorrected");
    }
  
    virtual double GetQ2Reco() const{
	return GetDouble("qsquared_recoil");
    }

    //GetRecoilE is designed to match the NSF validation suite
    virtual double GetRecoilE() const {
	return GetVecElem("recoil_summed_energy", 0);
    }
  
    virtual double Getq3() const{
	double eavail = GetEavail()/pow(10,3);
	double q2 = GetQ2Reco() / pow(10,6);
	double q3mec = sqrt(eavail*eavail + q2);
	return q3mec;
    }
   
    virtual int GetCurrent() const { return GetInt("mc_current"); }

    virtual int GetTruthNuPDG() const { return GetInt("mc_incoming"); }

    virtual double GetMuonQP() const {
	return GetDouble((GetAnaToolName() + "_minos_trk_qp").c_str());
    }

    //Some functions to match CCQENuInclusive treatment of DIS weighting. Name matches same Dan area as before.
    virtual double GetTrueExperimentersQ2() const {
	double Enu = GetEnuTrue(); //MeV
	double Emu = GetElepTrue(); //MeV
	double thetaMu = GetThetalepTrue();
	return 4.0*Enu*Emu*pow(sin(thetaMu/2.0),2.0);//MeV^2
    }

    virtual double CalcTrueExperimentersQ2(double Enu, double Emu, double thetaMu) const{
	return 4.0*Enu*Emu*pow(sin(thetaMu/2.0),2.0);//MeV^2
    }

    virtual double GetTrueExperimentersW() const {
	double nuclMass = M_nucleon;
	int struckNucl = GetTargetNucleon();
	if (struckNucl == PDG_n){
	    nuclMass=M_n;
	}
	else if (struckNucl == PDG_p){
	    nuclMass=M_p;
	}
	double Enu = GetEnuTrue();
	double Emu = GetElepTrue();
	double thetaMu = GetThetalepTrue();
	double Q2 = CalcTrueExperimentersQ2(Enu, Emu, thetaMu);
	return TMath::Sqrt(pow(nuclMass,2) + 2.0*(Enu-Emu)*nuclMass - Q2);
    }

    //Still needed for some systematics to compile, but shouldn't be used for reweighting anymore.
protected:
#include "PlotUtils/WeightFunctions.h" // Get*Weight
};

#endif
