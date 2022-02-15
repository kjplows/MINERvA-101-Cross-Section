/*
  
  Cuts from Ben Messerly's CC-CH-pip analysis supporting Michel quality + matching cuts.
  
  All other cuts on pion candidates should come *after* these are done.

  https://github.com/MinervaExpt/CC-CH-pip-ana/

 */

#ifndef BENSCUTS_H
#define BENSCUTS_H

#include <tuple>
#include <vector>

#include "event/CVUniverse.h"
#include "event/MichelEvent.h"
#include "BensHelpers.h"

namespace Jreco { // add to cuts from JohnsCuts.h

    // ==== At Least One Michel ==== // modified to work with MAT
    // For now, we need at least one ENDPOINT michel (any # of vtx michels).
    // This cut fills our michel containers, which we use to ID pion tracks
    // and subsequently make track cuts (LLR, node).
    template< class UNIVERSE, class EVENT = PlotUtils::detail::empty >
    class AtLeastOneMichel: public PlotUtils::Cut< UNIVERSE, EVENT >
    {
    public:
	AtLeastOneMichel() : PlotUtils::Cut< UNIVERSE, EVENT >( "At least 1 endpoint Michel" ) {}

    private:
	bool checkCut( const UNIVERSE& univ, EVENT& /* evt */ ) const override
	{
	    MichelMap all_michels = GetQualityMichels(univ);

	    MichelMap vertex_michels, endpoint_michels;
	    
	    for ( auto m : all_michels ) {
		if (m.second.had_idx == -1){
		    vertex_michels.insert(m);
		    // fVertex_michels.insert(m) crashes - don't initialise these maps as members!
		}
		else
		    endpoint_michels.insert(m);
		    //fEndpoint_michels.insert(m); // near match, but no dice
	    }
	    return endpoint_michels.size() > 0;
	}
    };
	
} // namespace Jreco
    
#endif // #ifndef BENSCUTS_H
