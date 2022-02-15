/*

  Helper functions for Ben's cuts to work. 
  Migrated here to keep BensCuts.h clean.

 */

#ifndef BENSHELPERS_H
#define BENSHELPERS_H

#include "BensCuts.h"

// -- Given a single michel cluster matched to two vertices
//    return vertex with the better-matched michel.
Michel CompareMichels(Michel r, Michel c);

// Add michel to MichelMap. Check if this cluster has already been matched.
// Then use only the best match.
bool AddOrReplaceMichel(MichelMap& mm, Michel m);

// Get the pion candidate indexes from a map of michels
std::vector<int> GetHadIdxsFromMichels(MichelMap michels);

//============================================================================
//  Collect the good michels in this event
//============================================================================
// Get quality michels map<(int michel cluster ID, Michel)>
// * A Michel is uniquely ID-ed by its cluster(of hits) integer index.
//   * The michel tool has already vetted the clusters themselves.
//     Whereas here we evaluate the quality of the fit.
// * At most one interaction vertex michel (which MUST be "fitted")
// * An OV michel will only be kept if no better michels are in the event.
// * In the case of 2+ OV michels, only keep the best one.
// * required: one-to-one michel<->vertex matching:
//    * when one michel cluster is matched to multiple vertices, the vtx
//    with the better match is chosen.
//    * We don't have the problem in the other direction -- michel tool: a
//    vertex cannot be matched to more than one cluster.
// * At the moment, we can return a single michel that is a quality
//   interaction vertex michel. We'd like to call these signal, but we don't
//   know the pion energy...yet. In the meantime, cut them.
MichelMap GetQualityMichels(const CVUniverse& univ);

// implementation here

Michel CompareMichels(Michel r, Michel c) {
    if (r.match_category > c.match_category)
	return r;
    else if (r.match_category < c.match_category)
	return c;
    else {
	if (r.fit_distance < c.fit_distance)
	    return r;
	else if (r.fit_distance > c.fit_distance)
	    return c;
	else {
	    // This must mean we're comparing the same michels with the same fits.
	    // std::cout << "WEIRD COMPAREMICHELS PROBLEM" << std::endl;
	    return r;
	}
    }
}

bool AddOrReplaceMichel(MichelMap& mm, Michel m) {
    std::pair<MichelMap::iterator, bool> SC;
    if (mm.count(m.idx) == 0)
	SC = mm.insert(std::pair<int, Michel>(m.idx, m));
    else {
	Michel reigning_michel = mm.at(m.idx);
	Michel best_michel = CompareMichels(reigning_michel, m);
	mm[m.idx] = best_michel;
    }
    return true;
}

std::vector<int> GetHadIdxsFromMichels(MichelMap michels) {
    std::vector<int> ret;
    for (auto m : michels) ret.push_back(m.second.had_idx);
    return ret;
}

MichelMap GetQualityMichels(const CVUniverse& univ) {
    std::map<int, Michel> ret_michels;
    std::vector<int> matched_michel_idxs = univ.GetVec<int>("matched_michel_idx");

    // Loop vertices in the event, i.e. the indices of the michel index vector
    for (uint vtx = 0; vtx < matched_michel_idxs.size(); ++vtx) {
	int mm_idx = matched_michel_idxs[vtx];

	// NO MATCH -- GO TO NEXT VTX. No michel cluster matched to this vertex.
	if (mm_idx < 0) continue;

	// MICHEL CONSTRUCTOR
	// Set match category (e.g. fitted, one-view, etc), match distance.
	Michel current_michel = Michel(univ, mm_idx, vtx);

	// MICHEL MATCH QUALITY CUT -- NEXT VTX
	// A michel cluster was matched to this vertex, but the match doesn't
	// pass match quality cuts.
	if (current_michel.match_category == Michel::kNoMatch) continue;

	// VERTEX MICHEL -- must meet the gold standard for its fit
	bool isIntVtx = (vtx == 0);
	if (isIntVtx && current_michel.match_category <= Michel::kNoFit) {
	    continue;
	}
	// ENDPOINT MICHELS
	else {
	    // ZERO MICHELS FOUND SO FAR
	    if (ret_michels.size() == 0)
		;

	    // ONE MICHEL FOUND SO FAR
	    // -- If either the michel we have already or this michel is OV, pick
	    //    the better of the two. Only one will remain.
	    else if (ret_michels.size() == 1) {
		Michel reigning_michel = (ret_michels.begin())->second;
		if (reigning_michel.match_category == Michel::kOV ||
		    current_michel.match_category == Michel::kOV)
		    ret_michels.clear();
		current_michel = CompareMichels(reigning_michel, current_michel);
	    }

	    // 2+ MICHELS FOUND SO FAR
	    else {
		if (current_michel.match_category == Michel::kOV) continue;
	    }
	}

	// ADD THIS MICHEL
	// When a cluster is matched to two vertices, pick the vtx with the
	// better match.
	bool SC = AddOrReplaceMichel(ret_michels, current_michel);
    }  // end vtx loop
    return ret_michels;
}


#endif // #ifndef BENSHELPERS_H
