//===============================================
// Signal definition for COH
// corresponding to `JohnsCuts.h`
//===============================================

#include "PlotUtils/Cut.h"
#include "event/CVUniverse.h"
#ifndef JOHNSSIGNALDEFINITION_H
#define JOHNSSIGNALDEFINITION_H

namespace Jtruth
{
    template <class UNIVERSE>
    class IsNumu: public PlotUtils::SignalConstraint<UNIVERSE>
    {
    public:
	IsNumu(): PlotUtils::SignalConstraint<UNIVERSE>("IsNumu") {}
	
    private:
	bool checkConstraint(const UNIVERSE& univ) const override
	{
	    return univ.GetTruthNuPDG() == 14;
	}
    };
    
    template <class UNIVERSE>
    class IsNumubar: public PlotUtils::SignalConstraint<UNIVERSE>
    {
    public:
	IsNumubar(): PlotUtils::SignalConstraint<UNIVERSE>("IsNumubar") {}

    private:
	bool checkConstraint(const UNIVERSE& univ) const override
	{
	    return univ.GetTruthNuPDG() == 14;
	}
    };
    
    template <class UNIVERSE>
    class IsNue: public PlotUtils::SignalConstraint<UNIVERSE>
    {
    public:
	IsNue(): PlotUtils::SignalConstraint<UNIVERSE>("IsNue") {}
	
    private:
	bool checkConstraint(const UNIVERSE& univ) const override
	{
	    return univ.GetTruthNuPDG() == 12;
	}
    };
    
    template <class UNIVERSE>
    class IsNuebar: public PlotUtils::SignalConstraint<UNIVERSE>
    {
    public:
	IsNuebar(): PlotUtils::SignalConstraint<UNIVERSE>("IsNuebar") {}

    private:
	bool checkConstraint(const UNIVERSE& univ) const override
	{
	    return univ.GetTruthNuPDG() == 12;
	}
    };
    
    template <class UNIVERSE>
    class IsCC: public PlotUtils::SignalConstraint<UNIVERSE>
    {
    public:
	IsCC(): PlotUtils::SignalConstraint<UNIVERSE>("IsCC") {}
	
    private:
	bool checkConstraint(const UNIVERSE& univ) const override
	{
	    return univ.GetCurrent() == 1;
	}
    };

    template <class UNIVERSE>
    class ThreeFSParticles: public PlotUtils::SignalConstraint<UNIVERSE>
    {
    public:
	ThreeFSParticles(): PlotUtils::SignalConstraint<UNIVERSE>("2 FS particles + nucleus") {}
	
    private:
	bool checkConstraint(const UNIVERSE& univ) const override
	{
	    // need to query from "truth" or "mc" branches.
	    //return univ.GetNFSPart() == 3; //NOT STRONG ENOUGH! Sometimes recoil nucleus doesn't get written and you get e.g. mu + p + pi+ RES.
	    bool passesNFS = (univ.GetNFSPart() == 3);
	    // want ONE charged lepton in position 0
	    bool passesLepton = (univ.GetPDGFromRecord(0) == 13) || (univ.GetPDGFromRecord(0) == 11);
	    // want ONE nucleus in positions 1 or 2
	    bool passesNucleus = (univ.GetPDGFromRecord(1) >= 1e+5 || univ.GetPDGFromRecord(2) >= 1e+5);
	    return passesNFS && passesLepton && passesNucleus;
	    // something + lepton + recoil nucleus
	}
    };

    template <class UNIVERSE>
    class AtLeastFourFSParticles: public PlotUtils::SignalConstraint<UNIVERSE>
    {
    public:
	AtLeastFourFSParticles(): PlotUtils::SignalConstraint<UNIVERSE>(">=4 FS particles OR no nucleus") {}
	
    private:
	bool checkConstraint(const UNIVERSE& univ) const override
	{
	    if(univ.GetNFSPart() >= 4) return true;
	    bool crossesNucleus = true;
	    for(int i = 0; i < univ.GetNFSPart(); i++){
		if(univ.GetPDGFromRecord(i) >= 1e+5) crossesNucleus = false;
	    }
	    return crossesNucleus; // NFSPart == 3 but none is a nucleus
	}
    };
    
    template <class UNIVERSE>
    class OnePiPlus: public PlotUtils::SignalConstraint<UNIVERSE>
    {
    public:
	OnePiPlus(): PlotUtils::SignalConstraint<UNIVERSE>("One pi+") {}
	
    private:
	bool checkConstraint(const UNIVERSE& univ) const override
	{
	    // need to query from "truth" or "mc" branches.
	    return univ.GetNPiPlusTrue() == 1;
	    // pi+ + lepton + recoil nucleus
	}
    };

    template <class UNIVERSE>
    class NTrueParticles: public PlotUtils::SignalConstraint<UNIVERSE>
    {
    public:
	NTrueParticles(const std::string& name, const int nPart, const int PDG): PlotUtils::SignalConstraint<UNIVERSE>(name), fnum(nPart), fPDG(PDG) {}

    private:
	bool checkConstraint(const UNIVERSE& univ) const override
	{
	    return univ.GetNPartTrue(fPDG) == fnum;
	}

	const int fnum;
	const int fPDG;
    };

    template <class UNIVERSE>
    class JunkParticlesAllowed: public PlotUtils::SignalConstraint<UNIVERSE>
    { // defined as kaons + neutrinos + gammas + electrons
    public:
	JunkParticlesAllowed(const std::string& name, const bool isAllowed): PlotUtils::SignalConstraint<UNIVERSE>(name), fIsAllowed(isAllowed) {}

    private:
	bool checkConstraint(const UNIVERSE& univ) const override
	{
	    const int KpPDG = 321, KmPDG = -321, K0PDG = 311, K0LPDG = 130, K0SPDG = 310, gamPDG = 22, nuePDG = 12, numuPDG = 14, nutauPDG = 16, ePDG = 11, tauPDG = 15, lambdaPDG = 3122, SigmapPDG = 3222, Sigma0PDG = 3212;
	    int PDGToCheck[] = {KpPDG, KmPDG, K0PDG, K0LPDG, K0SPDG, gamPDG, nuePDG, numuPDG, nutauPDG, ePDG, tauPDG, lambdaPDG, SigmapPDG, Sigma0PDG};

	    int njunk = 0;
	    for(int i = 0 ; i < sizeof(PDGToCheck)/sizeof(int) ; i++){
		int tmpPDG = PDGToCheck[i];
		njunk += univ.GetNPartTrue(tmpPDG);
	    }
	    
	    return fIsAllowed ? njunk != 0 : njunk == 0;
	}

	const bool fIsAllowed;
    };

    template <class UNIVERSE>
    class FormalNegationOfJunkNotAllowed: public PlotUtils::SignalConstraint<UNIVERSE>{
    public:
	FormalNegationOfJunkNotAllowed(): PlotUtils::SignalConstraint<UNIVERSE>("Events with at least 1 junk particle") {}

    private:
	bool checkConstraint(const UNIVERSE& univ) const override
	{
	    const int KpPDG = 321, KmPDG = -321, K0PDG = 311, K0LPDG = 130, K0SPDG = 310, gamPDG = 22, nuePDG = 12, numuPDG = 14, nutauPDG = 16, ePDG = 11, tauPDG = 15, lambdaPDG = 3122, SigmapPDG = 3222, Sigma0PDG = 3212;
	    int PDGToCheck[] = {KpPDG, KmPDG, K0PDG, K0LPDG, K0SPDG, gamPDG, nuePDG, numuPDG, nutauPDG, ePDG, tauPDG, lambdaPDG, SigmapPDG, Sigma0PDG};

	    int njunk = 0;
	    for(int i = 0 ; i < sizeof(PDGToCheck)/sizeof(int) ; i++){
		int tmpPDG = PDGToCheck[i];
		njunk += univ.GetNPartTrue(tmpPDG);
	    }
	    
	    return njunk > 0;
	}
    };

    template <class UNIVERSE>
    class NumNucleons: public PlotUtils::SignalConstraint<UNIVERSE>{
    public:
	NumNucleons(const std::string& name, const int nnuc): PlotUtils::SignalConstraint<UNIVERSE>(name), fNuc(nnuc) {}
	
    private:
	bool checkConstraint(const UNIVERSE& univ) const override
	{
	    const int pPDG = 2212, nPDG = 2112;
	    int ntnuc = 0;
	    ntnuc += univ.GetNPartTrue(pPDG);
	    ntnuc += univ.GetNPartTrue(nPDG);

	    return ntnuc == fNuc;
	}

	const int fNuc;
    };

    template <class UNIVERSE>
    class AtLeastNParticles: public PlotUtils::SignalConstraint<UNIVERSE>
    {
    public:
	AtLeastNParticles(const std::string& name, const int pdg, const int num): PlotUtils::SignalConstraint<UNIVERSE>(name), fPDG(pdg), fNum(num) {}

    private:
	bool checkConstraint(const UNIVERSE& univ) const override
	{   
	    return univ.GetNPartTrue(fPDG) >= fNum;
	}

	const int fPDG;
	const int fNum;
    };

    template <class UNIVERSE>
    class FormalNegationOfPionBand: public PlotUtils::SignalConstraint<UNIVERSE>
    {
    public:
	FormalNegationOfPionBand(): PlotUtils::SignalConstraint<UNIVERSE>("!{1,0} && !{1,1} && !{2,0} pi(+,0) && > 2pi(+,0)") {}

    private:
	bool checkConstraint(const UNIVERSE& univ) const override
	{
	    // note !(A&&B) == (!A) || (!B) !!!!
	    //bool not10 = (univ.GetNPartTrue(211) != 1 || univ.GetNPartTrue(111) != 0);
	    //bool not11 = (univ.GetNPartTrue(211) != 1 || univ.GetNPartTrue(111) != 1);
	    //bool not20 = (univ.GetNPartTrue(211) != 2 || univ.GetNPartTrue(111) != 0);
	    int npi = univ.GetNPartTrue(211) + univ.GetNPartTrue(111);
	    return (npi > 2);
	}
    };

    template <class UNIVERSE>
    class NNucleons: public PlotUtils::SignalConstraint<UNIVERSE>{
    public:
	NNucleons(): PlotUtils::SignalConstraint<UNIVERSE>(">1 n/p") {}
	
    private:
	bool checkConstraint(const UNIVERSE& univ) const override
	{
	    const int pPDG = 2212, nPDG = 2112;
	    int ntnuc = 0;
	    ntnuc += univ.GetNPartTrue(pPDG);
	    ntnuc += univ.GetNPartTrue(nPDG);

	    return ntnuc > 1;
	}
    };

    template <class UNIVERSE>
    class OnePiOrP: public PlotUtils::SignalConstraint<UNIVERSE>
    {
    public:
	OnePiOrP(): PlotUtils::SignalConstraint<UNIVERSE>("One pi+ or p") {}
	
    private:
	bool checkConstraint(const UNIVERSE& univ) const override
	{
	    //return (univ.GetNFSPart() <= 3 && (univ.GetNProtonsTrue() == 1 || univ.GetNPiPlusTrue() == 1) && univ.GetNProtonsTrue() != univ.GetNPiPlusTrue());
	    return (univ.GetNProtonsTrue() + univ.GetNPiPlusTrue() == 1);
	}
    };
    
    template <class UNIVERSE>
    class NoNeutrals: public PlotUtils::SignalConstraint<UNIVERSE>
    {
    public:
	NoNeutrals(): PlotUtils::SignalConstraint<UNIVERSE>("No neutral particles") {}
	
    private:
	bool checkConstraint(const UNIVERSE& univ) const override
	{
	    return univ.GetNNeutral() == 0;
	}
    };

    template <class UNIVERSE>
    class NoNegatives: public PlotUtils::SignalConstraint<UNIVERSE>
    {
    public:
	NoNegatives(): PlotUtils::SignalConstraint<UNIVERSE>("No negative hadrons") {}
	
    private:
	bool checkConstraint(const UNIVERSE& univ) const override
	{
	    return univ.GetNNegativeHadrons() == 0;
	}
    };

    template <class UNIVERSE>
    class MuonAngle: public PlotUtils::SignalConstraint<UNIVERSE>
    {
    public:
	MuonAngle(const double angleMax): PlotUtils::SignalConstraint<UNIVERSE>(std::string("Muon Angle ") + std::to_string(angleMax)), fMax(angleMax*M_PI/180.) {}
	
    private:
	bool checkConstraint(const UNIVERSE& univ) const override
	{
	    return univ.GetThetalepTrue() <= fMax;
	}
	
	const double fMax;
    };

    template <class UNIVERSE>
    class Apothem: public PlotUtils::SignalConstraint<UNIVERSE>
    {
    public:
	Apothem(const double apothem): PlotUtils::SignalConstraint<UNIVERSE>(std::string("Apothem ") + std::to_string(apothem)), fApothem(apothem) {}
	
    private:
      bool checkConstraint(const UNIVERSE& univ) const override
	{
	    const auto vertex = univ.GetTrueVertex();
	    return (fabs(vertex.y()) < fSlope*fabs(vertex.x()) + 2.*fApothem/sqrt(3.))
		&& (fabs(vertex.x()) < fApothem);
	}
	
	const double fApothem;
	const double fSlope = 1./sqrt(3.); //A regular hexagon has angles of 2*M_PI/3, so I can find this is 1/tan(M_PI/3.)
    };

    template <class UNIVERSE>
    class ZRange: public PlotUtils::SignalConstraint<UNIVERSE>
    {
    public:
	ZRange(const std::string& name, const double zMin, const double zMax): PlotUtils::SignalConstraint<UNIVERSE>(name), fMin(zMin), fMax(zMax) {}
	
    private:
	bool checkConstraint(const UNIVERSE& univ) const override
	{
	    return univ.GetTrueVertex().z() >= fMin && univ.GetTrueVertex().z() <= fMax;
	}

	const double fMin;
	const double fMax;
    };

    template <class UNIVERSE>
    class WRange: public PlotUtils::SignalConstraint<UNIVERSE>
    {
    public:
	WRange(const std::string& name, const double WMin, const double WMax): PlotUtils::SignalConstraint<UNIVERSE>(name), fMin(WMin), fMax(WMax) {}
	
    private:
	bool checkConstraint(const UNIVERSE& univ) const override
	{
	    return univ.GetWTrueGeV() >= fMin && univ.GetWTrueGeV() < fMax;
	}

	const double fMin;
	const double fMax;
    };

    template <class UNIVERSE>
    class EnuRange: public PlotUtils::SignalConstraint<UNIVERSE>
    {
    public:
	EnuRange(const std::string& name, const double EnuMin, const double EnuMax): PlotUtils::SignalConstraint<UNIVERSE>(name), fMin(EnuMin), fMax(EnuMax) {}
	
    private:
	bool checkConstraint(const UNIVERSE& univ) const override
	{
	    return univ.GetEnuTrueGeV() >= fMin && univ.GetEnuTrueGeV() <= fMax;
	}

	const double fMin;
	const double fMax;
    };

    template <class UNIVERSE>
    class EpiRange: public PlotUtils::SignalConstraint<UNIVERSE>
    {
    public:
	EpiRange(const std::string& name, const double EpiMin, const double EpiMax): PlotUtils::SignalConstraint<UNIVERSE>(name), fMin(EpiMin), fMax(EpiMax) {}
	
    private:
	bool checkConstraint(const UNIVERSE& univ) const override
	{
	    return univ.GetEpiTrueGeV() >= fMin && univ.GetEpiTrueGeV() < fMax;
	}

	const double fMin;
	const double fMax;
    };

    template <class UNIVERSE>
    class IsQE: public PlotUtils::SignalConstraint<UNIVERSE>
    {
    public:
	IsQE(): PlotUtils::SignalConstraint<UNIVERSE>("IsQE") {}

    private:
	bool checkConstraint(const UNIVERSE& univ) const override
	{
	    return univ.GetInteractionType() == 1;
	}
    };

    template <class UNIVERSE>
    class IsCOH: public PlotUtils::SignalConstraint<UNIVERSE>
    {
    public:
	IsCOH(): PlotUtils::SignalConstraint<UNIVERSE>("IsCOH") {}
	
    private:
	bool checkConstraint(const UNIVERSE& univ) const override
	{
	    return univ.GetInteractionType() == 4;
	}
    };

    template <class UNIVERSE>
    class IsNotQEOrCOH: public PlotUtils::SignalConstraint<UNIVERSE>
    {
    public:
	IsNotQEOrCOH(): PlotUtils::SignalConstraint<UNIVERSE>("IsNotQEOrCOH") {}

    private:
	bool checkConstraint(const UNIVERSE& univ) const override
	{
	    return univ.GetInteractionType() != 1 && univ.GetInteractionType() != 4;
	}
    };

    // notnumu OR notCC
    template <class UNIVERSE>
    class IsOther: public PlotUtils::SignalConstraint<UNIVERSE>
    {
    public:
	IsOther(): PlotUtils::SignalConstraint<UNIVERSE>("NC or not numu") {}

    private:
	bool checkConstraint(const UNIVERSE& univ) const override
	{
	    bool notCCNumu  = (univ.GetTruthNuPDG() != 14 || univ.GetCurrent() != 1);
	    return notCCNumu;
	}
    };
    
}

#endif // #ifndef JOHNSSIGNALDEFINITION_H
