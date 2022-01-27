#include "PlotUtils/Cut.h"
#include "event/CVUniverse.h"
#include <sstream>
#ifndef JOHNSCUTS_H
#define JOHNSCUTS_H

#ifdef __GCCXML__
#define override
namespace std {
    std::string to_string(double x) {
	std::stringstream ss;
	ss << x;
	return ss.str();
    }
}
#endif //#ifndef __GCCXML__

//===========================================================
// Simple file containing all my cuts.
// Author: John Plows, U. of Oxford
// komninosjohn.plows@physics.ox.ac.uk
//===========================================================

//===========================================================
// Neutrino, antineutrino
// MINOS QP significance gives us QP < 0 ? Nu : Nubar
// I presume this is just the Q/P ratio so checks out
//===========================================================

namespace Jreco
{
    template <class UNIVERSE, class EVENT = PlotUtils::detail::empty>
    class MINOSNumu: public PlotUtils::Cut<UNIVERSE, EVENT>
    {
    public:
	MINOSNumu(): PlotUtils::Cut<UNIVERSE, EVENT>("MINOS Q/P < 0") {}
	
    private:
	bool checkCut(const UNIVERSE& univ, EVENT& /*evt*/) const override
	{
	    return univ.GetMuonQP() < 0; 
	}
    };
    
    template <class UNIVERSE, class EVENT = PlotUtils::detail::empty>
    class MINOSNumubar: public PlotUtils::Cut<UNIVERSE, EVENT>
    {
    public:
	MINOSNumubar(): PlotUtils::Cut<UNIVERSE, EVENT>("MINOS Q/P > 0") {}
    private:
	bool checkCut(const UNIVERSE& univ, EVENT& /*evt*/) const override
	{
	    return univ.GetMuonQP() > 0; 
	}
    };
    
    //===========================================================
    // Need MINOS match!
    //===========================================================
    template <class UNIVERSE, class EVENT = PlotUtils::detail::empty>
    class HasMINOSMatch: public PlotUtils::Cut<UNIVERSE, EVENT>
    {
    public:
	HasMINOSMatch(): PlotUtils::Cut<UNIVERSE, EVENT>("Has MINOS match") {}
	  
    private:
	bool checkCut(const UNIVERSE& univ, EVENT& /*evt*/) const override
	{
	    return univ.IsMinosMatchMuon();
	}
    };
    
    //============================================================
    // Stuff related to MINOS acceptance etc.
    //============================================================
    
    //muon angle input DEG output RAD
    template <class UNIVERSE, class EVENT = PlotUtils::detail::empty>
    class MaxMuonAngle: public PlotUtils::Cut<UNIVERSE, EVENT>
    {
    public:
	MaxMuonAngle(const double max): PlotUtils::Cut<UNIVERSE, EVENT>(std::string("Muon Theta < ") + std::to_string(max)), fMax(max*M_PI/180.) {}
	  
    private:
	bool checkCut(const UNIVERSE& univ, EVENT& /*evt*/) const override
	{
	    return univ.GetThetamu() < fMax;
	}
	  
	const double fMax; //radians
    };

    //===========================================================
    // Reco needs 1 pion
    //===========================================================
   
    template <class UNIVERSE, class EVENT = PlotUtils::detail::empty>
    class Has1Pion: public PlotUtils::Cut<UNIVERSE, EVENT>
    {
    public:
	// C'tor
	Has1Pion(): PlotUtils::Cut<UNIVERSE, EVENT>("Has 1 pion") {}
      
    private:
	// 1 hadron track, no secondaries, pion score > proton score
	bool checkCut(const UNIVERSE& univ, EVENT& /*evt*/) const override
	{
	    double NONSENSE = -999.;
	    
	    // old selection based on dEdXTool, not good - use MAD's LLR when available!
	    // pion score and proton score
	    double piScore = univ.GetPionScore();
	    double prScore = univ.GetProtonScore();
	    // need to be not nonsense
	    if(piScore <= NONSENSE) return false;
	    // want no secondaries
	    int nSec = 0;
	    nSec = univ.GetNSecondaryHadrons();
	    if(nSec > 0) return false;
	    // now only 1 pion. Return TRUE if pion score > proton score
	    return piScore >= prScore;
	}
    };

    //isolate CCQE as well!
    template <class UNIVERSE, class EVENT = PlotUtils::detail::empty>
    class Has1Proton: public PlotUtils::Cut<UNIVERSE, EVENT>
    {
    public:
	// C'tor
	Has1Proton(): PlotUtils::Cut<UNIVERSE, EVENT>("Has 1 proton") {}
      
    private:
	// 1 hadron track, no secondaries, pion score > proton score
	bool checkCut(const UNIVERSE& univ, EVENT& /*evt*/) const override
	{
	    double NONSENSE = -999.;

	    // old selection based on dEdXTool, not good  use MAD's LLR when available!
	    // pion score and proton score
	    double piScore = univ.GetPionScore();
	    double prScore = univ.GetProtonScore();
	    // need to be not nonsense
	    if(piScore <= NONSENSE) return false;
	    // want no secondaries
	    int nSec = 0;
	    nSec = univ.GetNSecondaryHadrons();
	    if(nSec > 0) return false;
	    // now only 1 pion. Return TRUE if pion score > proton score
	    return piScore <= prScore;
	}
    };

    template <class UNIVERSE, class EVENT = PlotUtils::detail::empty>
    class Has1PiOrP: public PlotUtils::Cut<UNIVERSE, EVENT>
    {
    public:
	// C'tor
	Has1PiOrP(): PlotUtils::Cut<UNIVERSE, EVENT>("Has 1 pi+ or p") {}
      
    private:
	// 1 hadron (= pion or proton? what about K?) track, no secondaries
	bool checkCut(const UNIVERSE& univ, EVENT& /*evt*/) const override
	{
	    double NONSENSE = -999.;
	    //double NONSENSE = -101.;

	    // pion score and proton score
	    double piScore = univ.GetPionScore();
	    // need to be not nonsense
	    if(piScore <= NONSENSE) return false;
	    /*
	    double hadScore = univ.GetHadronLLR();
	    if(hadScore == NONSENSE) return false;
	    */
	    // want no secondaries
	    int nSec = univ.GetNSecondaryHadrons();
	    return nSec == 0;
	}
    };

    
    template <class UNIVERSE, class EVENT = PlotUtils::detail::empty>
    class ZRange: public PlotUtils::Cut<UNIVERSE, EVENT>
    {
    public:
	ZRange(const std::string& name, const double zMin, const double zMax): PlotUtils::Cut<UNIVERSE, EVENT>(name), fMin(zMin), fMax(zMax) {}
	  
    private:
	bool checkCut(const UNIVERSE& univ, EVENT& /*evt*/) const override
	{
	    return univ.GetVertexZ() >= fMin && univ.GetVertexZ() <= fMax;
	}
	  
	const double fMin;
	const double fMax;
    };
    
    template <class UNIVERSE, class EVENT = PlotUtils::detail::empty>
    class Apothem: public PlotUtils::Cut<UNIVERSE, EVENT>
    {
    public:
	Apothem(const double apothem): PlotUtils::Cut<UNIVERSE, EVENT>(std::string("Apothem ") + std::to_string(apothem)), fApothem(apothem), fSlope(1./sqrt(3.)) {} //A regular hexagon has angles of 2*M_PI/3, so I can find this is 1/tan(M_PI/3.) 
	
    private:
	bool checkCut(const UNIVERSE& univ, EVENT& /*evt*/) const override
	{
	    const ROOT::Math::XYZTVector vertex = univ.GetVertex();
	    return (fabs(vertex.y()) < fSlope*fabs(vertex.x()) + 2.*fApothem/sqrt(3.))
		&& (fabs(vertex.x()) < fApothem);
	}
	  
	const double fApothem;
	const double fSlope; 
    };
    
    //define W sidebands. Use min and max (GeV!)
    template <class UNIVERSE, class EVENT = PlotUtils::detail::empty>
    class WRange: public PlotUtils::Cut<UNIVERSE, EVENT>
    {
    public:
	WRange(const std::string& name, const double WMin, const double WMax): PlotUtils::Cut<UNIVERSE, EVENT>(name), fMin(WMin), fMax(WMax) {}
	
    private:
	bool checkCut(const UNIVERSE& univ, EVENT& /*evt*/) const override
	{
	    return univ.GetWGeV() >= fMin && univ.GetWGeV() < fMax;
	}
	  
	const double fMin;
	const double fMax;
    };
    
    template <class UNIVERSE, class EVENT = PlotUtils::detail::empty>
    class EnuRange: public PlotUtils::Cut<UNIVERSE, EVENT>
    {
    public:
	EnuRange(const std::string& name, const double EnuMin, const double EnuMax): PlotUtils::Cut<UNIVERSE, EVENT>(name), fMin(EnuMin), fMax(EnuMax) {}
	  
    private:
	bool checkCut(const UNIVERSE& univ, EVENT& /*evt*/) const override
	{
	    return univ.GetEnuGeV() >= fMin && univ.GetEnuGeV() <= fMax;
	}
	  
	const double fMin;
	const double fMax;
    };

#ifndef __GCCXML__
    template <class UNIVERSE, class EVENT = PlotUtils::detail::empty>
    using NoDeadtime = PlotUtils::Maximum<UNIVERSE, int, &UNIVERSE::GetTDead, EVENT>; //Andrew's fancy template notation
#endif
}

#endif //#ifndef JOHNSCUTS_H
