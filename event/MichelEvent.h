#ifndef MichelEvent_h
#define MichelEvent_h

#include "event/CVUniverse.h"

struct MichelEvent {
    int m_idx; // Index for Best Michel in nmichels
    double m_bestdist; // in mm 
    std::vector<double> m_best2D; //0: XZ, 1: UZ, 2:VZ   
    double m_best_XZ;
    double m_best_UZ;
    double m_best_VZ;
    //std::vector<Michel*> m_nmichels; //nmatched michels
};

// Need to match Michels to pions in order to ensure good reco.
// I will mirror the way Ben does it : 
// At least one (endpoint) Michel first, then LLR cut, then leading pion
// See https://github.com/MinervaExpt/CC-CH-pip-ana/ --> Cuts.cxx, CCPiEvent.cxx, Michel.h

// Also see Aaron Bercellie's event selection
// https://minerva-docdb.fnal.gov/cgi-bin/sso/RetrieveFile?docid=29445&filename=ME-Pions-Plots-Results.pdf&version=1 page 8

class Michel;

typedef std::map< int, Michel > MichelMap;

bool IsQualityMatchedMichel_Fit(double fit_dist, double fit_cut);
bool IsQualityMatchedMichel_NoFit(double nofit_dist, double nofit_cut);
bool IsQualityMatchedMichel_OneView(double ov_dist, double ov_cut);

class Michel {
 public:
    enum EMatchCategory { kNoMatch, kFit, kNoFit, kOV, kNMatchCategories };
    
    // constructors
    Michel(const CVUniverse& univ, int i, int v);
  Michel()
      : idx(-105),
        vtx(-106),
        had_idx(-107),
        match_category(kNoMatch),
        fit_distance(-1.) {}
    
    // integer that uniquely identifies cluster of hits
    int idx;
    
    // integer corresponding to vertex to which cluster was matched.
    // vtx == 0 --> interaction vertex
    // vtx == 1 --> "first" track endpoint, corresponds to hadron index 0.
    //              i.e. hadron index = michel vtx - 1.
    int vtx;
    
    // hadron index to which this michel is matched (had_idx = vtx - 1)
    int had_idx;
    
    EMatchCategory match_category;
    double fit_distance;
};

Michel::Michel(const CVUniverse& univ, int i, int v)
: idx(i), vtx(v), had_idx(v - 1) {
    bool isIntVtx = (vtx == 0);
    // distances for fitted, 2/3-view nofit, and 1-view michels, respectively
    double mm_fit_dist = univ.GetVecElem("matched_michel_end_dist", idx) / 10.;
    //double mm_nofit_dist = univ.GetVecElem("matched_michel_avg_dist", idx) / 10.;
    double mm_ov_dist = univ.GetVecElem("matched_michel_ov_dist", idx) / 10.;
    // bad heresy here, RETHER
    double mm_nofit_dist = ( mm_fit_dist + mm_ov_dist ) / 2.0;
    // NEW
    const double FIT_CUT = isIntVtx ? 9.0 : 7.5;     // cm
    const double NOFIT_CUT = isIntVtx ? 10.0 : 50.;  // cm
    // OLD
    // const double FIT_CUT   = isIntVtx ? 9.0  :  5.; // cm
    // const double NOFIT_CUT = isIntVtx ? 10.0 : 10.; // cm
    // TODO distinguish between 2/3 view and OV
    // const double FIT_CUT   = isIntVtx ? 9. :  15.0; // cm
    // const double NOFIT_CUT = isIntVtx ? 10. : 15.0; // cm
    
    if (IsQualityMatchedMichel_Fit(mm_fit_dist, FIT_CUT)) {
	match_category = kFit;
	fit_distance = mm_fit_dist;
	return;
    } else if (IsQualityMatchedMichel_NoFit(mm_nofit_dist, NOFIT_CUT)) {
	match_category = kNoFit;
	fit_distance = mm_nofit_dist;
	return;
    } else if (IsQualityMatchedMichel_OneView(mm_ov_dist, NOFIT_CUT)) {
	match_category = kOV;
	fit_distance = mm_ov_dist;
	return;
    } else {
	match_category = kNoMatch;
	fit_distance = -2.;
	return;
    }
}

//==============================================================================
// MICHEL QUALITY CUTS
//==============================================================================
bool IsQualityMatchedMichel_Fit(double fit_dist, double fit_cut) {
    if (fit_dist >= 0 && fit_dist < fit_cut)
	return true;
    else
	return false;
}

bool IsQualityMatchedMichel_NoFit(double nofit_dist, double nofit_cut) {
    if (nofit_dist >= 0 && nofit_dist < nofit_cut)
	return true;
    else
	return false;
}

bool IsQualityMatchedMichel_OneView(double ov_dist, double ov_cut) {
    if (ov_dist >= 0 && ov_dist * (2.0 / 3) < ov_cut)
	return true;
    else
	return false;
}

#endif
