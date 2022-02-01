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

    // helpers

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

    // where is primary hadron going

    double HadronIsExitingID() const
    {
	return GetInt( (GetAnaToolName()+"_hadron_isExiting").c_str() );
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
	return GetDouble((GetAnaToolName()+"_hadron_piFit_scoreLLR").c_str());
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
	const int kPiPlusPDG = 211;
	int npi = 0;
	
	std::string pdgBranch = "mc_FSPartPDG";
	for(int i = 0; i < GetNFSPart(); i++){
	    if(GetVecElem(pdgBranch.c_str(), i) == kPiPlusPDG) npi++;
	}
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

    //============================================
    // I want to leave HadronFunctions for
    // MAD 1.49. Older branches for Epi etc
    // I'll put here. RETHERE REMOVE (?)
    //============================================

    // reco

    double GetEVtx() const // MeV
    {
	double EVTX = GetDouble("vtx_blobs_energy");
	double EISO = GetDouble("vtx_iso_blobs_energy_outside_radius"); // 150 mm by default
	return EVTX - EISO;
    }

    double GetERecoilVtx150mm() const // MeV
    {
	double EVNM = GetDouble("recoil_energy_nonmuon_vtx150mm");
	return EVNM;
    }

    double GetEhad() const //MeV
    {
	//double hadE = GetVecElem((GetAnaToolName()+"_HadronE").c_str(), 3);
	// VITAL: THIS IS NEUTRON ENERGY. NOT PROTON!!!!!
	double hadE = GetDouble((GetAnaToolName()+"_proton_E_fromdEdx").c_str());
	//double hadE = GetDouble((GetAnaToolName()+"_proton_calib_energy").c_str());
	//correction in v1.34 to get *pion* energy! Seems to all be proton stuff
	const double pionScore = GetPionScore();
	const double protonScore = GetProtonScore();
        return (pionScore >= protonScore) ? hadE - MinervaUnits::M_p + MinervaUnits::M_pion : hadE;
	//return hadE;
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

    ROOT::Math::PxPyPzEVector GetHadron4VTrue(const int PDG) const {
	int idx = GetIndexOfHad(PDG);
	double Px = GetVecElem("mc_FSPartPx",idx);
	double Py = GetVecElem("mc_FSPartPy",idx);
	double Pz = GetVecElem("mc_FSPartPz",idx);
	double E = GetVecElem("mc_FSPartE",idx);

	return ROOT::Math::PxPyPzEVector(Px,Py,Pz,E);
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
	const double piScore = GetPionScore();
	const double prScore = GetProtonScore();
	const double Mhad = (piScore >= prScore) ? MinervaUnits::M_pion : MinervaUnits::M_p;
	return GetMuHadAngle(Mhad) * 180. / TMath::Pi();
    }

    // simpler |t| calc using A. Mislivec's Eq. (8.6)
    // |t| = |(p_\nu - p_\mu - p_\pi)^2|

    double GetAbsT() const { //MeV^2
	ROOT::Math::PxPyPzEVector mu4V = GetMuon4V();
	const double pion_score = GetPionScore();
	const double proton_score = GetProtonScore();
	const double Mhad = (pion_score >= proton_score) ? MinervaUnits::M_pion : MinervaUnits::M_p ;
	ROOT::Math::PxPyPzEVector had4V = GetHadron4V(Mhad);
	const double Enu = GetEnu();
	ROOT::Math::PxPyPzEVector nu4V(0.,0.,Enu,Enu);

	ROOT::Math::PxPyPzEVector sys4V = nu4V - mu4V - had4V;
	return std::abs(sys4V.M() * sys4V.M());
    }

    double GetAbsTGeV() const {
	return GetAbsT() / 1e+6;
    }
    
    // now construct |t| using Alex's Eq. (6.1)
    double GetAlexAbsT() const { //MeV^2
	const double piScore = GetPionScore();
	const double prScore = GetProtonScore();
	const double Mhad = (piScore >= prScore) ?  MinervaUnits::M_pion : MinervaUnits::M_proton; //need to declare in func
	    
	double Emu      = GetEmu();
	double Ehad     = GetEhad();
	double Pmu      = GetPmu();
	double Phad     = GetPhad(Mhad);
	double thetamu  = GetThetamu();
	double thetahad = GetHadronTheta();
	double thetaMuHad = GetMuHadAngle(Mhad);
	
	const double Mmu = MinervaUnits::M_mu;

	double sumE = Emu + Ehad;
	double term1 = -2. * sumE * (Emu - Pmu * std::cos(thetamu));
	double term2 = Mmu * Mmu;
	double term3 = -2. * (Ehad*Ehad - sumE*Phad*std::cos(thetahad) + Pmu*Phad*std::cos(thetaMuHad));
	double term4 = Mhad * Mhad;

	return std::abs(term1 + term2 + term3 + term4);
    }

    double GetAlexAbsTGeV() const {
	return pow(1e-3, 2) * GetAlexAbsT();
    }

    // calculator function for Enu under COH
    // Equation (3.7) in Alex's thesis

    double GetCOHEnu() const {
	return GetEmu() + GetEhad();
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
	double Emu     = GetEmuGeV();
	double Pmu     = GetPmu() * 1e-3;
	double thetamu = GetThetamu();
	double Enu     = GetCOHEnuGeV();

	return 2. * Enu * (Emu - Pmu * std::cos(thetamu)) - mmu * mmu;
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
	return GetEmu() + GetEhad(); //good assuming measurement of Ehad is good
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

    // and a version for Variable declaration
    double GetThetaMuHadTrueVar() const {
	const int PDG_pi = 211; //need to declare within function scope
	const int PDG_p  = 2212;
	const int PDG = (GetNPiPlusTrue() > GetNProtonsTrue()) ? PDG_pi : PDG_p;
	return GetMuHadAngleTrue(PDG) * 180. / TMath::Pi();
    }

    //calc |t|
    double GetAbsTTrue() const { // MeV^2
	ROOT::Math::PxPyPzEVector mu4V = GetLep4VTrue();
	const int PDG_pi = 211;
	const int PDG_p  = 2212;
	const int PDG = (GetNPiPlusTrue() > GetNProtonsTrue()) ? PDG_pi : PDG_p;
	ROOT::Math::PxPyPzEVector had4V = GetHadron4VTrue(PDG);
	const double Enu = GetEnuTrue();
	ROOT::Math::PxPyPzEVector nu4V(0.,0.,Enu,Enu);

	ROOT::Math::PxPyPzEVector sys4V = nu4V - mu4V - had4V;
	return std::abs(sys4V.M() * sys4V.M());
    }

    double GetAbsTTrueGeV() const {
	return 1e-6 * GetAbsTTrue(); // truth is 100x larger than reco!!!!
    }

    double GetAlexAbsTTrue() const { //MeV^2
	const int PDG_pi = 211; //need to declare within function scope
	const int PDG_p  = 2212;
	const int PDG = (GetNPiPlusTrue() > GetNProtonsTrue()) ? PDG_pi : PDG_p;
	
	double Emu        = GetElepTrue();
	double Ehad       = GetEhadTrue();
	double Pmu        = GetPlepTrue();
	double Phad       = GetPhadTrue();
	double thetamu    = GetThetalepTrue();
	double thetahad   = GetThetaHadTrue(PDG);
	double thetaMuHad = GetMuHadAngleTrue(PDG); 

	const double Mmu = MinervaUnits::M_mu;
	double Mhad = 0;
	switch(PDG){
	    case 211:  Mhad = MinervaUnits::M_pion; break;
	    case 2212: Mhad = MinervaUnits::M_p;    break;
	    case 2112: Mhad = MinervaUnits::M_n;    break;
	    default:   Mhad = 0;
	}

	double sumE = Emu + Ehad;
	double term1 = -2. * sumE * (Emu - Pmu * std::cos(thetamu));
	double term2 = Mmu * Mmu;
	double term3 = -2. * (Ehad*Ehad - sumE*Phad*std::cos(thetahad) + Pmu*Phad*std::cos(thetaMuHad));
	double term4 = Mhad * Mhad;

	return std::abs(term1 + term2 + term3 + term4);
    }

    double GetAlexAbsTTrueGeV() const {
	return pow(1.0e-3, 2.0) * GetAlexAbsTTrue();
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
	const double piScore = GetPionScore();
	const double prScore = GetProtonScore();
	const double Mhad = (piScore > prScore) ? MinervaUnits::M_pion : MinervaUnits::M_p ;
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
    // Get vertex energy
    // Get this from all Prongs assoc with vtx
    //============================================
    // make and fill my branches with this
    
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
	return GetEmu()/1000.;
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
