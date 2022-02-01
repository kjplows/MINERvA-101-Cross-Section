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
// Reorganised to follow Alex Ramirez's COH analysis
//===========================================================

namespace Jreco
{

    // -- Fiducial volume selection
    
    // -- ZRange gives allowed z in MINERvA coords

    template <class UNIVERSE, class EVENT = PlotUtils::detail::empty>
    class ZRange: public PlotUtils::Cut<UNIVERSE, EVENT>
    {
	
    public:
	ZRange(const std::string& name, const double zMin, const double zMax): PlotUtils::Cut<UNIVERSE, EVENT>(name), fMin(zMin), fMax(zMax) {}
	
    private:
	bool checkCut(const UNIVERSE& univ, EVENT& /*evt*/) const override
	{ return univ.GetVertexZ() >= fMin && univ.GetVertexZ() <= fMax; }
	
	const double fMin;
	const double fMax;
	
    }; // class ZRange

    // -- Apothem gives the hexagonal region
    
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
	
    }; // class Apothem
    
    // -- MINOS acceptance

    template <class UNIVERSE, class EVENT = PlotUtils::detail::empty>
    class HasMINOSMatch: public PlotUtils::Cut<UNIVERSE, EVENT>
    {
	
    public:
	HasMINOSMatch(): PlotUtils::Cut<UNIVERSE, EVENT>("Has MINOS match") {}
	  
    private:
	bool checkCut(const UNIVERSE& univ, EVENT& /*evt*/) const override
	{ return univ.IsMinosMatchMuon(); }
	
    }; // class HasMINOSMatch

    // -- No dead time
    
#ifndef __GCCXML__
    template <class UNIVERSE, class EVENT = PlotUtils::detail::empty>
    using NoDeadtime = PlotUtils::Maximum<UNIVERSE, int, &UNIVERSE::GetTDead, EVENT>; //Andrew's fancy template notation
#endif

    // -- Neutrino mode selection

    // -- mu ?
    
    template <class UNIVERSE, class EVENT = PlotUtils::detail::empty>
    class MINOSNumu: public PlotUtils::Cut<UNIVERSE, EVENT>
    {
    public:
	MINOSNumu(): PlotUtils::Cut<UNIVERSE, EVENT>("MINOS Numu") {}
	
    private:
	bool checkCut(const UNIVERSE& univ, EVENT& /*evt*/) const override
	{ return univ.GetMuonQP() < 0; }
	
    }; // class MINOSNumu

    // -- mubar ?
    
    template <class UNIVERSE, class EVENT = PlotUtils::detail::empty>
    class MINOSNumubar: public PlotUtils::Cut<UNIVERSE, EVENT>
    {
    public:
	MINOSNumubar(): PlotUtils::Cut<UNIVERSE, EVENT>("MINOS Numubar") {}
	
    private:
	bool checkCut(const UNIVERSE& univ, EVENT& /*evt*/) const override
	{ return univ.GetMuonQP() > 0; }
	
    }; // class MINOSNumubar

    // -- One hadron only please!

    template <class UNIVERSE, class EVENT = PlotUtils::detail::empty>
    class OneHadron: public PlotUtils::Cut<UNIVERSE, EVENT>
    {
    public:
	OneHadron(): PlotUtils::Cut<UNIVERSE, EVENT>("One hadron") {}

    private:
	bool checkCut(const UNIVERSE& univ, EVENT& /*evt*/) const override
	{ return univ.GetNHadrons() == 1; }

    }; // class OneHadron
    
    // -- Hadron containment in tracker
    
    template <class UNIVERSE, class EVENT = PlotUtils::detail::empty>
    class ContainedHadron: public PlotUtils::Cut<UNIVERSE, EVENT>
    {
    public:
	ContainedHadron(): PlotUtils::Cut<UNIVERSE, EVENT>("ID-Contained hadron") {}
	
    private:
	bool checkCut(const UNIVERSE& univ, EVENT& /*evt*/) const override
	{ return univ.HadronIsExitingID() == 0; }
	
    }; // class ContainedHadron

    // -- Neutrino energy

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
	
    }; // class EnuRange

    // -- PID selection

    // -- Is pion: log(L_pi) >= log(L_p) ==> L_pi >= L_p
    
    template <class UNIVERSE, class EVENT = PlotUtils::detail::empty>
    class IsPion: public PlotUtils::Cut<UNIVERSE, EVENT>
    {
    public:
	IsPion(): PlotUtils::Cut<UNIVERSE, EVENT>("Hadron is pion") {}

    private:
	bool checkCut(const UNIVERSE& univ, EVENT& /*evt*/) const override
	{
	    return univ.GetHadronLLR() >= 0.0;
	}
	
    }; // class IsPion

    // -- E-vtx cut pending distributions of COH evts & cuts up to here

    // -- Testing various Erecoil branch cuts @ different values here.

    template <class UNIVERSE, class EVENT = PlotUtils::detail::empty>
    class VtxECut: public PlotUtils::Cut<UNIVERSE, EVENT>
    {
    public:
	VtxECut( const std::string & name,
		 double Elow, double Ehigh ):
	    PlotUtils::Cut<UNIVERSE, EVENT>(name), fElow(Elow), fEhigh(Ehigh) {}

    private:
	bool checkCut(const UNIVERSE& univ, EVENT & /* evt */) const override
	{
	    const double EVTX = univ.GetERecoilVtx200mm();
	    return ( fElow <= EVTX && fEhigh >= EVTX );
	}

	double fElow, fEhigh;
	
    }; // class VtxECut

    // -- |t| cut pending distributions
    
} // namespace Jreco

#endif //#ifndef JOHNSCUTS_H
