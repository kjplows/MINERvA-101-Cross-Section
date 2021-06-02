//File: GetMnvTunev1.h
//Brief: Sets up MnvTunev1 for PlotUtils::Model.  TODO: Link to a paper or docdb entry explaining what's supposed to be in it.
//Author: Andrew Olivier aolivier@ur.rochester.edu

//PlotUtils/reweighters includes
#include "PlotUtils/reweighters/Model.h"
#include "PlotUtils/reweighters/AMUDISReweighter.h"
#include "PlotUtils/reweighters/FluxAndCVReweighter.h"
#include "PlotUtils/reweighters/GENIEReweighter.h"
#include "PlotUtils/reweighters/LowRecoil2p2hReweighter.h"
#include "PlotUtils/reweighters/MINOSEfficiencyReweighter.h"
#include "PlotUtils/reweighters/RPAReweighter.h"

template <class UNIVERSE, class EVENT = PlotUtils::detail::empty>
std::vector<std::unique_ptr<PlotUtils::Reweighter<UNIVERSE, EVENT>>> GetMnvTunev1()
{
  std::vector<std::unique_ptr<PlotUtils::Reweighter<UNIVERSE, EVENT>>> weighters;

  weighters.emplace_back(new PlotUtils::FluxAndCVReweighter<UNIVERSE, EVENT>());
  weighters.emplace_back(new PlotUtils::GENIEReweighter<UNIVERSE, EVENT>(true, false));
  weighters.emplace_back(new PlotUtils::LowRecoil2p2hReweighter<UNIVERSE, EVENT>());
  weighters.emplace_back(new PlotUtils::MINOSEfficiencyReweighter<UNIVERSE, EVENT>());
  weighters.emplace_back(new PlotUtils::RPAReweighter<UNIVERSE, EVENT>());
  weighters.emplace_back(new PlotUtils::AMUDISReweighter<UNIVERSE, EVENT>());

  return weighters;
}
