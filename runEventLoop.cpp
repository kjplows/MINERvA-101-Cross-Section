#define MC_OUT_FILE_NAME "runEventLoopMC.root"
#define DATA_OUT_FILE_NAME "runEventLoopData.root"

#define USAGE \
"\n*** USAGE ***\n"\
"runEventLoop <dataPlaylist.txt> <mcPlaylist.txt>\n\n"\
"*** Explanation ***\n"\
"MAD files reduced to event-selection histos\n"\
"that supports my analysis for HNL in MINERvA. MC backgrounds specified here.\n\n"\
"*** The Input Files ***\n"\
"Playlist files are plaintext files with 1 file name per line.  Filenames may be\n"\
"xrootd URLs or refer to the local filesystem.  The first playlist file's\n"\
"entries will be treated like data, and the second playlist's entries must\n"\
"have the \"Truth\" tree to use for calculating the efficiency denominator.\n\n"\
"*** Output ***\n"\
"Produces a two files, " MC_OUT_FILE_NAME " and " DATA_OUT_FILE_NAME ", with\n"\
"all histograms needed for the ExtractCrossSection program also built by this\n"\
"package.  You'll need a .rootlogon.C that loads ROOT object definitions from\n"\
"PlotUtils to access systematics information from these files.\n\n"\
"*** Environment Variables ***\n"\
"Setting up this package appends to PATH and LD_LIBRARY_PATH.  PLOTUTILSROOT,\n"\
"MPARAMFILESROOT, and MPARAMFILES must be set according to the setup scripts in\n"\
"those packages for systematics and flux reweighters to function.\n"\
"If MNV101_SKIP_SYST is defined at all, output histograms will have no error bands.\n"\
"This is useful for debugging the CV and running warping studies.\n\n"\
"*** Return Codes ***\n"\
"0 indicates success.  All histograms are valid only in this case.  Any other\n"\
"return code indicates that histograms should not be used.  Error messages\n"\
"about what went wrong will be printed to stderr.  So, they'll end up in your\n"\
"terminal, but you can separate them from everything else with something like:\n"\
"\"runEventLoop data.txt mc.txt 2> errors.txt\"\n"

enum ErrorCodes
{
  success = 0,
  badCmdLine = 1,
  badInputFile = 2,
  badFileRead = 3,
  badOutputFile = 4
};

//PlotUtils includes
//No junk from PlotUtils please!  I already
//know that MnvH1D does horrible horrible things.
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Woverloaded-virtual"

//Includes from this package
#include "event/CVUniverse.h"
#include "event/MichelEvent.h"
#include "systematics/Systematics.h"
#include "cuts/MaxPzMu.h"
#include "util/Variable.h"
#include "util/Variable2D.h"
#include "util/GetFluxIntegral.h"
#include "util/GetPlaylist.h"
#include "cuts/SignalDefinition.h"
#include "cuts/q3RecoCut.h"
#include "studies/Study.h"
//#include "Binning.h" //TODO: Fix me

//Custom cuts from me
#include "cuts/BensCuts.h"
#include "cuts/JohnsCuts.h"
#include "cuts/JohnsSignalDefinition.h"

//PlotUtils includes
#include "PlotUtils/makeChainWrapper.h"
#include "PlotUtils/HistWrapper.h"
#include "PlotUtils/Hist2DWrapper.h"
#include "PlotUtils/MacroUtil.h"
#include "PlotUtils/MnvPlotter.h"
#include "PlotUtils/CCInclusiveCuts.h"
#include "PlotUtils/CCInclusiveSignal.h"
#include "PlotUtils/CrashOnROOTMessage.h" //Sets up ROOT's debug callbacks by itself
#include "PlotUtils/Cutter.h"
#include "PlotUtils/Model.h"
#include "PlotUtils/FluxAndCVReweighter.h"
#include "PlotUtils/GENIEReweighter.h"
#include "PlotUtils/LowRecoil2p2hReweighter.h"
#include "PlotUtils/RPAReweighter.h"
#include "PlotUtils/MINOSEfficiencyReweighter.h"
#include "PlotUtils/TargetUtils.h"
#pragma GCC diagnostic pop

//ROOT includes
#include "TParameter.h"

//c++ includes
#include <iostream>
#include <cstdlib> //getenv()

//==============================================================================
// Loop and Fill
//==============================================================================
void LoopAndFillEventSelection(
    PlotUtils::ChainWrapper* chain,
    std::map<std::string, std::vector<CVUniverse*> > error_bands,
    std::vector<Variable*> vars,
    std::vector<Variable2D*> vars2D,
    std::vector<Variable*> resolutions,
    std::vector<Study*> studies,
    PlotUtils::Cutter<CVUniverse, MichelEvent>& michelcuts,
    PlotUtils::Model<CVUniverse, MichelEvent>& model)
{
  assert(!error_bands["cv"].empty() && "\"cv\" error band is empty!  Can't set Model weight.");
  auto& cvUniv = error_bands["cv"].front();

  std::cout << "Starting MC reco loop...\n";
  const int nEntries = chain->GetEntries();
  for (int i=0; i<nEntries; ++i)
  {
    if(i%1000==0) std::cout << i << " / " << nEntries << "\r" << std::flush;

      //if(i >= 20000){std::cout << "HIT CEILING AND EXITING" << std::endl; exit(0);}

    MichelEvent cvEvent;
    cvUniv->SetEntry(i);
    model.SetEntry(*cvUniv, cvEvent);
    const double cvWeight = model.GetWeight(*cvUniv, cvEvent);

    //=========================================
    // Systematics loop(s)
    //=========================================
    for (auto band : error_bands)
    {
      std::vector<CVUniverse*> error_band_universes = band.second;
      for (auto universe : error_band_universes)
      {
        MichelEvent myevent; // make sure your event is inside the error band loop. 
    
        // Tell the Event which entry in the TChain it's looking at
        universe->SetEntry(i);
        
        // This is where you would Access/create a Michel

        //weight is ignored in isMCSelected() for all but the CV Universe.
        if (!michelcuts.isMCSelected(*universe, myevent, cvWeight).all()) continue; //all is another function that will later help me with sidebands
        const double weight = model.GetWeight(*universe, myevent); //Only calculate the per-universe weight for events that will actually use it.

        for(auto& var: vars) var->selectedMCReco->FillUniverse(universe, var->GetRecoValue(*universe), weight); //"Fake data" for closure

        const bool isSignal = michelcuts.isSignal(*universe, weight);

	const bool isPhaseSpace = michelcuts.isPhaseSpace(*universe, weight);

        if(isSignal)
        {
	    
          for(auto& study: studies) study->SelectedSignal(*universe, myevent, weight);

          for(auto& var: vars)
          {
            //Cross section components
	      if(isPhaseSpace){
		  var->efficiencyNumerator->FillUniverse(universe, var->GetTrueValue(*universe), weight);
		  var->migration->FillUniverse(universe, var->GetRecoValue(*universe), var->GetTrueValue(*universe), weight);
		  var->selectedSignalReco->FillUniverse(universe, var->GetRecoValue(*universe), weight); //Efficiency numerator in reco variables.  Useful for warping studies.
	      }
          }

          for(auto& var: vars2D)
          {
	      if(isPhaseSpace){
		  var->efficiencyNumerator->FillUniverse(universe, var->GetTrueValueX(*universe), var->GetTrueValueY(*universe), weight);
		  //var->efficiencyNumerator->FillUniverse(universe, var->GetRecoValueX(*universe), var->GetRecoValueY(*universe), weight); // lazy!
		  //var->recoCorrelation->FillUniverse(universe, var->GetRecoValueX(*universe), var->GetRecoValueY(*universe), weight);
	      }
          }

	  for(auto& var: resolutions)
	  {
	      //Grab reco and true value and fill the resolution histo
	      const double vreco = var->GetRecoValue(*universe);
	      const double vtrue = var->GetTrueValue(*universe);
	      const double vden = ( vtrue != 0.0 ) ? vtrue : 1.0;
	      
	      const double vfill = ( vreco - vtrue );
	      if(isPhaseSpace){
		  var->resolution->FillUniverse(universe, vfill, weight );

		  if( vden != 1.0 ){ var->fractionalResolution->FillUniverse( universe, vfill / vden, weight ); }
		  var->resolution2D->FillUniverse(universe, vfill, vtrue, weight);
		  
	      } // if(isPhaseSpace)
	  } // for(auto & var: resolutions)
        }
        else
        {
          int bkgd_ID = -1;
          if (universe->GetCurrent()==2)bkgd_ID=0;
          else bkgd_ID=1;

          for(auto& var: vars) (*var->m_backgroundHists)[bkgd_ID].FillUniverse(universe, var->GetRecoValue(*universe), weight);
          for(auto& var: vars2D) (*var->m_backgroundHists)[bkgd_ID].FillUniverse(universe, var->GetRecoValueX(*universe), var->GetRecoValueY(*universe), weight);
        }
      } // End band's universe loop
    } // End Band loop
  } //End entries loop
  std::cout << "Finished MC reco loop.\n";
}

void LoopAndFillData( PlotUtils::ChainWrapper* data,
			        std::vector<CVUniverse*> data_band,
				std::vector<Variable*> vars,
                                std::vector<Variable2D*> vars2D,
                                std::vector<Study*> studies,
				PlotUtils::Cutter<CVUniverse, MichelEvent>& michelcuts)

{
  std::cout << "Starting data loop...\n";
  const int nEntries = data->GetEntries();
  for (int i=0; i<data->GetEntries(); ++i) {
    for (auto universe : data_band) {
      universe->SetEntry(i);
      if(i%1000==0) std::cout << i << " / " << nEntries << "\r" << std::flush;
      MichelEvent myevent; 
      if (!michelcuts.isDataSelected(*universe, myevent).all()) continue;

      for(auto& study: studies) study->Selected(*universe, myevent, 1); 

      for(auto& var: vars)
      {
        var->dataHist->FillUniverse(universe, var->GetRecoValue(*universe, myevent.m_idx), 1);
      }

      for(auto& var: vars2D)
      {
        var->dataHist->FillUniverse(universe, var->GetRecoValueX(*universe), var->GetRecoValueY(*universe), 1);
      }
    }
  }
  std::cout << "Finished data loop.\n";
}

void LoopAndFillEffDenom( PlotUtils::ChainWrapper* truth,
    				std::map<std::string, std::vector<CVUniverse*> > truth_bands,
    				std::vector<Variable*> vars,
                                std::vector<Variable2D*> vars2D,
    				PlotUtils::Cutter<CVUniverse, MichelEvent>& michelcuts,
                                PlotUtils::Model<CVUniverse, MichelEvent>& model)
{
  assert(!truth_bands["cv"].empty() && "\"cv\" error band is empty!  Could not set Model entry.");
  auto& cvUniv = truth_bands["cv"].front();

  std::cout << "Starting efficiency denominator loop...\n";
  const int nEntries = truth->GetEntries();
  for (int i=0; i<nEntries; ++i)
  {
    if(i%1000==0) std::cout << i << " / " << nEntries << "\r" << std::flush;

    MichelEvent cvEvent;
    cvUniv->SetEntry(i);
    model.SetEntry(*cvUniv, cvEvent);
    const double cvWeight = model.GetWeight(*cvUniv, cvEvent);

    //=========================================
    // Systematics loop(s)
    //=========================================
    for (auto band : truth_bands)
    {
      std::vector<CVUniverse*> truth_band_universes = band.second;
      for (auto universe : truth_band_universes)
      {
        MichelEvent myevent; //Only used to keep the Model happy

        // Tell the Event which entry in the TChain it's looking at
        universe->SetEntry(i);

        if (!michelcuts.isEfficiencyDenom(*universe, cvWeight)) continue; //Weight is ignored for isEfficiencyDenom() in all but the CV universe 
        const double weight = model.GetWeight(*universe, myevent); //Only calculate the weight for events that will use it

        //Fill efficiency denominator now: 
        for(auto var: vars)
        {
          var->efficiencyDenominator->FillUniverse(universe, var->GetTrueValue(*universe), weight);
        }

        for(auto var: vars2D)
        {
          var->efficiencyDenominator->FillUniverse(universe, var->GetTrueValueX(*universe), var->GetTrueValueY(*universe), weight);
        }
      }
    }
  }
  std::cout << "Finished efficiency denominator loop.\n";
}

//Returns false if recoTreeName could not be inferred
bool inferRecoTreeNameAndCheckTreeNames(const std::string& mcPlaylistName, const std::string& dataPlaylistName, std::string& recoTreeName)
{
  const std::vector<std::string> knownTreeNames = {"Truth", "Meta"};
  bool areFilesOK = false;

  std::ifstream playlist(mcPlaylistName);
  std::string firstFile = "";
  playlist >> firstFile;
  auto testFile = TFile::Open(firstFile.c_str());
  if(!testFile)
  {
    std::cerr << "Failed to open the first MC file at " << firstFile << "\n";
    return false;
  }

  //Does the MC playlist have the Truth tree?  This is needed for the efficiency denominator.
  const auto truthTree = testFile->Get("Truth");
  if(truthTree == nullptr || !truthTree->IsA()->InheritsFrom(TClass::GetClass("TTree")))
  {
    std::cerr << "Could not find the \"Truth\" tree in MC file named " << firstFile << "\n";
    return false;
  }

  //Figure out what the reco tree name is
  for(auto key: *testFile->GetListOfKeys())
  {
    if(static_cast<TKey*>(key)->ReadObj()->IsA()->InheritsFrom(TClass::GetClass("TTree"))
       && std::find(knownTreeNames.begin(), knownTreeNames.end(), key->GetName()) == knownTreeNames.end())
    {
      recoTreeName = key->GetName();
      areFilesOK = true;
    }
  }
  delete testFile;
  testFile = nullptr;

  //Make sure the data playlist's first file has the same reco tree
  playlist.open(dataPlaylistName);
  playlist >> firstFile;
  testFile = TFile::Open(firstFile.c_str());
  if(!testFile)
  {
    std::cerr << "Failed to open the first data file at " << firstFile << "\n";
    return false;
  }

  const auto recoTree = testFile->Get(recoTreeName.c_str());
  if(recoTree == nullptr || !recoTree->IsA()->InheritsFrom(TClass::GetClass("TTree")))
  {
    std::cerr << "Could not find the \"" << recoTreeName << "\" tree in data file named " << firstFile << "\n";
    return false;
  }

  return areFilesOK;
}

//==============================================================================
// Main
//==============================================================================
int main(const int argc, const char** argv)
{
  TH1::AddDirectory(false);

  //Validate input.
  //I expect a data playlist file name and an MC playlist file name which is exactly 2 arguments.
  const int nArgsExpected = 2;
  if(argc != nArgsExpected + 1) //argc is the size of argv.  I check for number of arguments + 1 because
                                //argv[0] is always the path to the executable.
  {
    std::cerr << "Expected " << nArgsExpected << " arguments, but got " << argc - 1 << "\n" << USAGE << "\n";
    return badCmdLine;
  }

  //One playlist must contain only MC files, and the other must contain only data files.
  //Only checking the first file in each playlist because opening each file an extra time
  //remotely (e.g. through xrootd) can get expensive.
  //TODO: Look in INSTALL_DIR if files not found?
  const std::string mc_file_list = argv[2],
                    data_file_list = argv[1];

  //Check that necessary TTrees exist in the first file of mc_file_list and data_file_list
  std::string reco_tree_name;
  if(!inferRecoTreeNameAndCheckTreeNames(mc_file_list, data_file_list, reco_tree_name))
  {
    std::cerr << "Failed to find required trees in MC playlist " << mc_file_list << " and/or data playlist " << data_file_list << ".\n" << USAGE << "\n";
    return badInputFile;
  }

  const bool doCCQENuValidation = (reco_tree_name == "CCQENu"); //Enables extra histograms and might influence which systematics I use.

  const bool is_grid = false;
  //PlotUtils::MacroUtil options(reco_tree_name, mc_file_list, data_file_list, "minervame1A", true, is_grid);
  PlotUtils::MacroUtil options(reco_tree_name, mc_file_list, data_file_list, "minervame1A", true);
  options.m_plist_string = util::GetPlaylist(*options.m_mc, true); //TODO: Put GetPlaylist into PlotUtils::MacroUtil

  // You're required to make some decisions
  PlotUtils::MinervaUniverse::SetNuEConstraint(true);
  PlotUtils::MinervaUniverse::SetPlaylist(options.m_plist_string); //TODO: Infer this from the files somehow?
  PlotUtils::MinervaUniverse::SetAnalysisNuPDG(14);
  PlotUtils::MinervaUniverse::SetNFluxUniverses(100);
  PlotUtils::MinervaUniverse::SetZExpansionFaReweight(false);

  //Now that we've defined what a cross section is, decide which sample and model
  //we're extracting a cross section for.
  PlotUtils::Cutter<CVUniverse, MichelEvent>::reco_t sidebands, preCuts;
  PlotUtils::Cutter<CVUniverse, MichelEvent>::truth_t signalDefinition, phaseSpace;

  //==========================================================
  // Declare cut & signal parameters here
  //==========================================================
  const double minZ = 5980, maxZ = 8422, apothem = 850; //All in mm
  const double maxMuTheta = 20.; //deg
  //const double minEnu = 2.0, maxEnu = 20.0; //GeV, see Alex's thesis
  const double minEnu = 2.0, maxEnu = 50.0, medEnu = 20.0; //GeV
  //const double minW = 0.0, maxW = 1.4; // GeV [LOW]
  //const double minW = 1.4, maxW = 2.0; // GeV [MED]
  //const double minW = 2.0, maxW = 1e+5; // GeV [HIGH]

  const int distmm = 200; // up to how far from vtx do we look for energy depositions
  // note this isn't an argument in the actual Cut, you'll have to modify manually
  const double vtxELow = 10.0, vtxEHigh = 60.0; // highest by eff * pur
  const double Thigh = 0.15; // GeV^2
  
  const int pPDG = 2212, nPDG = 2112, pipPDG = 211, pimPDG = -211, pi0PDG = 111;

  //==========================================================
  // The actual cuts. Use Jreco and Jtruth namespaces.
  //==========================================================
  // Data cuts actually *remove* events that don't pass from loop
  //==========================================================
  preCuts.emplace_back(new Jreco::ZRange<CVUniverse, MichelEvent>("Tracker", minZ, maxZ));
  preCuts.emplace_back(new Jreco::Apothem<CVUniverse, MichelEvent>(apothem));
  preCuts.emplace_back(new Jreco::NoDeadtime<CVUniverse, MichelEvent>(1, "Deadtime"));
  preCuts.emplace_back(new Jreco::HasMINOSMatch<CVUniverse, MichelEvent>());
  preCuts.emplace_back(new Jreco::MINOSNumu<CVUniverse, MichelEvent>());
  preCuts.emplace_back(new Jreco::OneHadron<CVUniverse, MichelEvent>());
  //preCuts.emplace_back(new Jreco::AtLeastOneMichel<CVUniverse, MichelEvent>()); // gotta fix?
  preCuts.emplace_back(new Jreco::ContainedHadron<CVUniverse, MichelEvent>());
  preCuts.emplace_back(new Jreco::IsPion<CVUniverse, MichelEvent>());
  preCuts.emplace_back(new Jreco::EnuRange<CVUniverse, MichelEvent>(Form("%1.1f <= Enu [GeV] <= %1.1f", minEnu, medEnu), minEnu, medEnu)); // 2 - 20 GeV
  // vtx energy cut
  preCuts.emplace_back(new Jreco::VtxECut<CVUniverse, MichelEvent>(Form("%1.1f <= E_vtx (%d mm) <= %1.1f [MeV]", vtxELow, distmm, vtxEHigh), vtxELow, vtxEHigh));
  //WIP: testing |t|
  //preCuts.emplace_back(new Jreco::TCut<CVUniverse, MichelEvent>(Form("|t| [GeV^{2}] <= %1.2f", Thigh), Thigh));
  //==========================================================
  // Sidebands are specific non-signal events that I wanna analyse
  //==========================================================
  
  //==========================================================
  // Signal definition sifts through data-OK evts
  //==========================================================
  signalDefinition.emplace_back(new Jtruth::IsNumu<CVUniverse>());
  signalDefinition.emplace_back(new Jtruth::IsCC<CVUniverse>());
  //signalDefinition.emplace_back(new Jtruth::IsOther<CVUniverse>());
  //------------------------------------------
  //signalDefinition.emplace_back(new Jtruth::ThreeFSParticles<CVUniverse>());
  //signalDefinition.emplace_back(new Jtruth::AtLeastFourFSParticles<CVUniverse>());
  signalDefinition.emplace_back(new Jtruth::IsCOH<CVUniverse>());
  //signalDefinition.emplace_back(new Jtruth::IsQE<CVUniverse>());
  //signalDefinition.emplace_back(new Jtruth::IsNotQEOrCOH<CVUniverse>());
  //------------------------------------------
  // -- pion bands
  //signalDefinition.emplace_back(new Jtruth::NTrueParticles<CVUniverse>("1 pi+", 1, pipPDG));
  //signalDefinition.emplace_back(new Jtruth::NTrueParticles<CVUniverse>("0 pi0", 0, pi0PDG));
  //signalDefinition.emplace_back(new Jtruth::AtLeastNParticles<CVUniverse>(">= 2pi+", 2, pipPDG));
  //signalDefinition.emplace_back(new Jtruth::FormalNegationOfPionBand<CVUniverse>());
  // -- nucleon bands
  //signalDefinition.emplace_back(new Jtruth::NTrueParticles<CVUniverse>("0 p", 0, pPDG));
  //signalDefinition.emplace_back(new Jtruth::NTrueParticles<CVUniverse>("0 n", 0, nPDG));
  //signalDefinition.emplace_back(new Jtruth::NumNucleons<CVUniverse>("1 nucleon", 1));
  //signalDefinition.emplace_back(new Jtruth::NNucleons<CVUniverse>());
  // -- junk?
  //signalDefinition.emplace_back(new Jtruth::JunkParticlesAllowed<CVUniverse>("Allow junk particles: FALSE", false));
  //signalDefinition.emplace_back(new Jtruth::FormalNegationOfJunkNotAllowed<CVUniverse>());
  
  //==========================================================
  // Phase space rejects events that are likely poorly reco'd
  //==========================================================
  phaseSpace.emplace_back(new Jtruth::ZRange<CVUniverse>("Tracker", minZ, maxZ));
  phaseSpace.emplace_back(new Jtruth::Apothem<CVUniverse>(apothem));
  //phaseSpace.emplace_back(new Jtruth::MuonAngle<CVUniverse>(maxMuTheta));
  phaseSpace.emplace_back(new Jtruth::EnuRange<CVUniverse>(Form("%1.1f <= Enu [GeV] <= %1.1f", minEnu, medEnu), minEnu, medEnu)); // 2 - 20 GeV
  //phaseSpace.emplace_back(new Jtruth::WRange<CVUniverse>(Form("%1.1f <= W [GeV] < %1.1f", minW, maxW), minW, maxW));
  //==========================================================
  
  PlotUtils::Cutter<CVUniverse, MichelEvent> mycuts(std::move(preCuts), std::move(sidebands) , std::move(signalDefinition),std::move(phaseSpace));

  std::vector<std::unique_ptr<PlotUtils::Reweighter<CVUniverse, MichelEvent>>> MnvTunev1;
  MnvTunev1.emplace_back(new PlotUtils::FluxAndCVReweighter<CVUniverse, MichelEvent>());
  MnvTunev1.emplace_back(new PlotUtils::GENIEReweighter<CVUniverse, MichelEvent>(true, false));
  MnvTunev1.emplace_back(new PlotUtils::LowRecoil2p2hReweighter<CVUniverse, MichelEvent>());
  MnvTunev1.emplace_back(new PlotUtils::MINOSEfficiencyReweighter<CVUniverse, MichelEvent>());
  MnvTunev1.emplace_back(new PlotUtils::RPAReweighter<CVUniverse, MichelEvent>());

  PlotUtils::Model<CVUniverse, MichelEvent> model(std::move(MnvTunev1));

  // Make a map of systematic universes
  // Leave out systematics when making validation histograms
  const bool doSystematics = (getenv("MNV101_SKIP_SYST") == nullptr);
  if(!doSystematics){
    std::cout << "Skipping systematics (except 1 flux universe) because environment variable MNV101_SKIP_SYST is set.\n";
    PlotUtils::MinervaUniverse::SetNFluxUniverses(2); //Necessary to get Flux integral later...  Doesn't work with just 1 flux universe though because _that_ triggers "spread errors".
  }

  std::map< std::string, std::vector<CVUniverse*> > error_bands;
  if(doSystematics) error_bands = GetStandardSystematics(options.m_mc);
  else{
    std::map<std::string, std::vector<CVUniverse*> > band_flux = PlotUtils::GetFluxSystematicsMap<CVUniverse>(options.m_mc, CVUniverse::GetNFluxUniverses());
    error_bands.insert(band_flux.begin(), band_flux.end()); //Necessary to get flux integral later...
  }
  error_bands["cv"] = {new CVUniverse(options.m_mc)};
  std::map< std::string, std::vector<CVUniverse*> > truth_bands;
  if(doSystematics) truth_bands = GetStandardSystematics(options.m_truth);
  truth_bands["cv"] = {new CVUniverse(options.m_truth)};

  std::vector<double> dansPTBins = {0, 0.075, 0.15, 0.25, 0.325, 0.4, 0.475, 0.55, 0.7, 0.85, 1, 1.25, 1.5, 2.5, 4.5},
      dansPzBins = {1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 6, 7, 8, 9, 10, 15, 20, 40, 60},
      robsEmuBins = {0,1,2,3,4,5,7,9,12,15,18,22,36,50,75,100,120},
      robsRecoilBins,
      johnsQ2Bins = {0., 0.02, 0.04, 0.06, 0.08, 0.1, 0.12, 0.14, 0.16, 0.18, 0.2, 0.22, 0.24, 0.26, 0.28, 0.3, 0.35, 0.4, 0.45, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0},
      //johnsTBins = {0., 1.e+2, 2.e+2, 3.e+2, 4.e+2, 5.e+2, 6.e+2, 7.e+2, 8.e+2, 9.e+2, 1.e+3, 1.5e+3, 2.e+3, 2.5e+3, 3.e+3}, // NO THAT WAS FOR SQRT(|t|)!!
      johnsTBins = {0., 0.025, 0.05, 0.075, 0.1, 0.125, 0.15, 0.175, 0.2, 0.25, 0.3, 0.35, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.},
      johnsTMeVBins = {0., 2.5e+4, 5.0e+4, 7.5e+4, 1.0e+5, 1.25e+5, 1.5e+5, 1.75e+5, 2.0e+5, 2.5e+5, 3.0e+5, 3.5e+5, 4.0e+5, 5.0e+5, 6.0e+5, 7.0e+5, 8.0e+5, 9.0e+5, 1.0e+6},
      johnsFineTBins = {0., 0.025, 0.05, 0.075, 0.1, 0.125, 0.15, 0.175, 0.2, 0.225, 0.25, 0.275, 0.3, 0.325, 0.35, 0.375, 0.4, 0.425, 0.45, 0.475, 0.5, 0.525, 0.55, 0.575, 0.6, 0.625, 0.65, 0.675, 0.7, 0.725, 0.75, 0.775, 0.8, 0.825, 0.85, 0.875, 0.9, 0.925, 0.95, 0.975, 1.},
      johnsEmuBins = {0., 0.5, 1., 1.5, 2., 2.5, 3., 3.5, 4., 4.5, 5., 5.5, 6., 6.5, 7., 7.5, 8., 8.5, 9., 9.5, 10., 11., 12., 13., 14., 15., 16., 17., 18., 19., 20.},
      johnsEhadBins = {0., 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0},
      johnsEnuBins = {0., 0.5, 1., 1.5, 2., 2.5, 3., 3.5, 4., 4.5, 5., 5.5, 6., 6.5, 7., 7.5, 8., 8.5, 9., 9.5, 10., 10.5, 11., 11.5, 12., 12.5, 13., 13.5, 14., 14.5, 15., 15.5, 16., 16.5, 17., 17.5, 18., 18.5, 19., 19.5, 20.},
      johnsThetaMuPiBins = {0., 5., 10., 12., 14., 16., 18., 20., 22., 24., 26., 28., 30., 32., 34., 36., 38., 40., 42., 44., 46., 48., 50., 55., 60.},
      johnsThetaMuBins = {0.,2.,4.,6.,8.,10.,12.,14.,16.,18.,20.},
      johnsThetaHadBins = {0.,2.,4.,6.,8.,10.,12.,14.,16.,18.,20.,30.,40.,55.,70.},
      johnsMMuPiBins = {0., 100., 200., 300., 400., 500., 600., 700., 800., 900., 1000., 1100., 1200., 1300., 1400., 1500., 1600., 1700., 1800., 1900., 2000.},
      johnsWBins = {0., 500., 600., 700., 800., 900., 1000., 1100., 1200., 1300., 1400., 1500., 1600., 1700., 1800., 1900., 2000., 2200., 2400., 2600., 2800., 3000.},
      johnsSysWBins = {100., 150., 200., 250., 300., 350., 400., 450., 500., 600., 700., 800., 900., 1000., 1250., 1500., 2000., 2500., 3000.},
      // now direct branches
      johnsScoreBins = {0., 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0};

  std::vector<double> johnsEvtxBins, johnsERecoilBins,
      johnsPTPiBins, johnsPzPiBins, johnsEPiBins,
      johnsPTMuBins, johnsPzMuBins, johnsEMuBins,
      johnsPTSysBins, johnsPzSysBins, johnsESysBins;
  
  for(Int_t j = 0; j < 101; j++){
      double elo = j * 5.0; johnsEvtxBins.push_back(elo); johnsERecoilBins.push_back(elo);
  }

  for(Int_t j = 0; j < 21; j++){
      double elo = j * 100.0; johnsEPiBins.push_back(elo);
  }

  for(Int_t j = 0; j < 31; j++){
      double elo = -1000.0 + j * 100.0; johnsPzPiBins.push_back(elo);
  }

  for(Int_t j = 0; j < 9; j++){
      double elo = -200.0 + j * 50.0; johnsPTPiBins.push_back(elo);
  }

  for(Int_t j = 0; j < 26; j++){
      double elo = j * 400.0; johnsEMuBins.push_back(elo);
  }

  for(Int_t j = 0; j < 22; j++){
      double elo = -400.0 + j * 400.0; johnsPzMuBins.push_back(elo);
  }

  for(Int_t j = 0; j < 17; j++){
      double elo = -400.0 + j * 50.0; johnsPTMuBins.push_back(elo);
  }

  for(Int_t j = 0; j < 21; j++){
      double elo = -500.0 + j * 50.0; johnsPTSysBins.push_back(elo);
  }

  for(Int_t j = 0; j < 11; j++){
      double elo = -100.0 + j * 50.0; johnsPzSysBins.push_back(elo);
  }

  for(Int_t j = 0; j < 9; j++){
      double elo = -100.0 + j * 25.0; johnsESysBins.push_back(elo);
  }

  const double robsRecoilBinWidth = 50; //MeV
  for(int whichBin = 0; whichBin < 100 + 1; ++whichBin) robsRecoilBins.push_back(robsRecoilBinWidth * whichBin);

  std::vector<Variable*> vars = {
      new Variable("pTmu", "p_{T, #mu} [GeV/c]", dansPTBins, &CVUniverse::GetMuonPT, &CVUniverse::GetMuonPTTrue), // 0
      new Variable("pionE", "E_{#pi} [MeV]", johnsEPiBins, &CVUniverse::GetEpi, &CVUniverse::GetEpiTrue),
      new Variable("pionP", "p_{#pi} [MeV]", johnsEPiBins, &CVUniverse::GetPpi, &CVUniverse::GetPpiTrue),
      new Variable("pionPx", "p_{x;#pi} [MeV]", johnsPTPiBins, &CVUniverse::GetPxpi, &CVUniverse::GetPxpiTrue),
      new Variable("pionPy", "p_{y;#pi} [MeV]", johnsPTPiBins, &CVUniverse::GetPypi, &CVUniverse::GetPypiTrue), // 4
      new Variable("pionPz", "p_{z;#pi} [MeV]", johnsPzPiBins, &CVUniverse::GetPzpi, &CVUniverse::GetPzpiTrue),
      new Variable("muonE", "E_{#mu} [MeV]", johnsEMuBins, &CVUniverse::GetEmuMAD, &CVUniverse::GetEmuMADTrue),
      new Variable("muonPx", "p_{x;#mu} [MeV]", johnsPTMuBins, &CVUniverse::GetPxmuWrtNuBeam, &CVUniverse::GetPxmuTrueWrtNuBeam),
      new Variable("muonPy", "p_{y;#mu} [MeV]", johnsPTMuBins, &CVUniverse::GetPymuWrtNuBeam, &CVUniverse::GetPymuTrueWrtNuBeam),
      new Variable("muonPz", "p_{z;#mu} [MeV]", johnsPzMuBins, &CVUniverse::GetPzmuWrtNuBeam, &CVUniverse::GetPzmuTrueWrtNuBeam), // 9
      new Variable("Enu", "E_{#nu} [MeV]", johnsEMuBins, &CVUniverse::GetEnu, &CVUniverse::GetEnuTrue),
      new Variable("sysE", "E_{sys} [MeV]", johnsESysBins, &CVUniverse::GetSys4VE, &CVUniverse::GetSys4VETrue),
      new Variable("sysPx", "p_{x;sys} [MeV]", johnsPTSysBins, &CVUniverse::GetSys4VPx, &CVUniverse::GetSys4VPxTrue),
      new Variable("sysPy", "p_{y;sys} [MeV]", johnsPTSysBins, &CVUniverse::GetSys4VPy, &CVUniverse::GetSys4VPyTrue),
      new Variable("sysPz", "p_{z;sys} [MeV]", johnsPzSysBins, &CVUniverse::GetSys4VPz, &CVUniverse::GetSys4VPzTrue), // 14
      new Variable("AbsT", "|t| [(GeV)^{2}]", johnsTBins, &CVUniverse::GetAbsTGeV, &CVUniverse::GetAbsTTrueGeV),
      new Variable("AlexAbsT", "|t|_{COH} [(GeV)^{2}]", johnsTBins, &CVUniverse::GetAlexAbsTGeV, &CVUniverse::GetAlexAbsTTrueGeV),
      new Variable("Q2", "Q^{2} [(GeV)^{2}]", johnsQ2Bins, &CVUniverse::GetCOHQ2GeV, &CVUniverse::GetCOHQ2TrueGeV),
      new Variable("thetaMuPi", "#theta_{#mu#pi} [#circ]", johnsThetaMuPiBins, &CVUniverse::GetThetaMuPiDeg, &CVUniverse::GetThetaMuPiTrueDeg),
      new Variable("MMuPi", "M_{#mu#pi} [MeV]", johnsMMuPiBins, &CVUniverse::GetMMuPi, &CVUniverse::GetMMuPiTrue), // 19
      new Variable("Evtx", "E_{vtx} [MeV]", johnsEvtxBins, &CVUniverse::GetERecoilVtx200mm, &CVUniverse::GetEpiTrue),
  };

  // bins for resolution studies

  const double PTMuHigh = std::max( std::abs( dansPTBins.at(0) ), std::abs( dansPTBins.at( dansPTBins.size() - 1 ) ) ) / 4.0;
  const double PzMuHigh = std::max( std::abs( johnsPzMuBins.at(0) ), std::abs( johnsPzMuBins.at( johnsPzMuBins.size() - 1 ) ) ) / 2.0;
  const double EMuHigh = std::max( std::abs( johnsEMuBins.at(0) ), std::abs( johnsEMuBins.at( johnsEMuBins.size() - 1 ) ) ) / 2.0;
  const double PTPiHigh = std::max( std::abs( johnsPTPiBins.at(0) ), std::abs( johnsPTPiBins.at( johnsPTPiBins.size() - 1 ) ) );
  const double PzPiHigh = std::max( std::abs( johnsPzPiBins.at(0) ), std::abs( johnsPzPiBins.at( johnsPzPiBins.size() - 1 ) ) );
  const double EPiHigh = std::max( std::abs( johnsEPiBins.at(0) ), std::abs( johnsEPiBins.at( johnsEPiBins.size() - 1 ) ) );
  const double PTSysHigh = std::max( std::abs( johnsPTSysBins.at(0) ), std::abs( johnsPTSysBins.at( johnsPTSysBins.size() - 1 ) ) );
  const double PzSysHigh = std::max( std::abs( johnsPzSysBins.at(0) ), std::abs( johnsPzSysBins.at( johnsPzSysBins.size() - 1 ) ) ) / 2.0;
  const double ESysHigh = std::max( std::abs( johnsESysBins.at(0) ), std::abs( johnsESysBins.at( johnsESysBins.size() - 1 ) ) ) / 4.0;
  const double THigh = std::max( std::abs( johnsTBins.at(0) ), std::abs( johnsTBins.at( johnsTBins.size() - 1 ) ) ) / 4.0;
  const double Q2High = std::max( std::abs( johnsQ2Bins.at(0) ), std::abs( johnsQ2Bins.at( johnsQ2Bins.size() - 1 ) ) ) / 4.0;
  const double thetaMuPiHigh = std::max( std::abs( johnsThetaMuPiBins.at(0) ), std::abs( johnsThetaMuPiBins.at( johnsThetaMuPiBins.size() - 1 ) ) ) / 2.0;
  const double MMuPiHigh = std::max( std::abs( johnsMMuPiBins.at(0) ), std::abs( johnsMMuPiBins.at( johnsMMuPiBins.size() - 1 ) ) ) / 2.0;

  std::vector<double> resPTMuBins, resPzMuBins, resEMuBins, resPTPiBins, resPzPiBins, resEPiBins, resPTSysBins, resPzSysBins, resESysBins, resTBins, resQ2Bins, resThetaMuPiBins, resMMuPiBins;

  for( auto j = 0; j < 101; j++ ) {
      resPTMuBins.emplace_back( PTMuHigh * ( -1.0 + 0.02 * j ) );
      resPzMuBins.emplace_back( PzMuHigh * ( -1.0 + 0.02 * j ) );
      resEMuBins.emplace_back( EMuHigh * ( -1.0 + 0.02 * j ) );
      resPTPiBins.emplace_back( PTPiHigh * ( -1.0 + 0.02 * j ) );
      resPzPiBins.emplace_back( PzPiHigh * ( -1.0 + 0.02 * j ) );
      resEPiBins.emplace_back( EPiHigh * ( -1.0 + 0.02 * j ) );
      resPTSysBins.emplace_back( PTSysHigh * ( -1.0 + 0.02 * j ) );
      resPzSysBins.emplace_back( PzSysHigh * ( -1.0 + 0.02 * j ) );
      resESysBins.emplace_back( ESysHigh * ( -1.0 + 0.02 * j ) );
      resTBins.emplace_back( THigh * ( -1.0 + 0.02 * j ) );
      resQ2Bins.emplace_back( Q2High * ( -1.0 + 0.02 * j ) );
      resThetaMuPiBins.emplace_back( thetaMuPiHigh * ( -1.0 + 0.02 * j ) );
      resMMuPiBins.emplace_back( MMuPiHigh * ( -1.0 + 0.02 * j ) );
  }
  
  std::vector<Variable*> res = {
      new Variable("pTmu", "p_{T, #mu} [GeV/c]", resPTMuBins, &CVUniverse::GetMuonPT, &CVUniverse::GetMuonPTTrue), // 0
      new Variable("pionE", "E_{#pi} [MeV]", resEPiBins, &CVUniverse::GetEpi, &CVUniverse::GetEpiTrue),
      new Variable("pionP", "p_{#pi} [MeV]", resEPiBins, &CVUniverse::GetPpi, &CVUniverse::GetPpiTrue),
      new Variable("pionPx", "p_{x;#pi} [MeV]", resPTPiBins, &CVUniverse::GetPxpi, &CVUniverse::GetPxpiTrue),
      new Variable("pionPy", "p_{y;#pi} [MeV]", resPTPiBins, &CVUniverse::GetPypi, &CVUniverse::GetPypiTrue), // 4
      new Variable("pionPz", "p_{z;#pi} [MeV]", resPzPiBins, &CVUniverse::GetPzpi, &CVUniverse::GetPzpiTrue),
      new Variable("muonE", "E_{#mu} [MeV]", resEMuBins, &CVUniverse::GetEmuMAD, &CVUniverse::GetEmuMADTrue),
      new Variable("muonPx", "p_{x;#mu} [MeV]", resPTMuBins, &CVUniverse::GetPxmuWrtNuBeam, &CVUniverse::GetPxmuTrueWrtNuBeam),
      new Variable("muonPy", "p_{y;#mu} [MeV]", resPTMuBins, &CVUniverse::GetPymuWrtNuBeam, &CVUniverse::GetPymuTrueWrtNuBeam),
      new Variable("muonPz", "p_{z;#mu} [MeV]", resPzMuBins, &CVUniverse::GetPzmuWrtNuBeam, &CVUniverse::GetPzmuTrueWrtNuBeam), // 9
      new Variable("Enu", "E_{#nu} [MeV]", resEMuBins, &CVUniverse::GetEnu, &CVUniverse::GetEnuTrue),
      new Variable("sysE", "E_{sys} [MeV]", resESysBins, &CVUniverse::GetSys4VE, &CVUniverse::GetSys4VETrue),
      new Variable("sysPx", "p_{x;sys} [MeV]", resPTSysBins, &CVUniverse::GetSys4VPx, &CVUniverse::GetSys4VPxTrue),
      new Variable("sysPy", "p_{y;sys} [MeV]", resPTSysBins, &CVUniverse::GetSys4VPy, &CVUniverse::GetSys4VPyTrue),
      new Variable("sysPz", "p_{z;sys} [MeV]", resPzSysBins, &CVUniverse::GetSys4VPz, &CVUniverse::GetSys4VPzTrue), // 14
      new Variable("AbsT", "|t| [(GeV)^{2}]", resTBins, &CVUniverse::GetAbsTGeV, &CVUniverse::GetAbsTTrueGeV),
      new Variable("AlexAbsT", "|t|_{COH} [(GeV)^{2}]", resTBins, &CVUniverse::GetAlexAbsTGeV, &CVUniverse::GetAlexAbsTTrueGeV),
      new Variable("Q2", "Q^{2} [(GeV)^{2}]", resQ2Bins, &CVUniverse::GetCOHQ2GeV, &CVUniverse::GetCOHQ2TrueGeV),
      new Variable("thetaMuPi", "#theta_{#mu#pi} [#circ]", resThetaMuPiBins, &CVUniverse::GetThetaMuPiDeg, &CVUniverse::GetThetaMuPiTrueDeg),
      new Variable("MMuPi", "M_{#mu#pi} [MeV]", resMMuPiBins, &CVUniverse::GetMMuPi, &CVUniverse::GetMMuPiTrue), // 19
      };

  std::vector<Variable2D*> vars2D = {
      new Variable2D(*vars[1], *vars[17]),
      new Variable2D(*vars[1], *vars[6]),
      new Variable2D(*vars[1], *vars[20]),
      /*
      new Variable2D(*vars[15], *vars[16]),
      new Variable2D(*vars[15], *vars[17]),
      new Variable2D(*vars[15], *vars[18]),
      new Variable2D(*vars[15], *vars[19]),
      new Variable2D(*vars[17], *vars[18]),
      new Variable2D(*vars[17], *vars[19]),
      new Variable2D(*vars[18], *vars[19]),
      */
  };
  /* = {
      new Variable2D(*vars[11], *vars[10])
      }; */
  
  if(doCCQENuValidation)
  {
    std::cerr << "Detected that tree name is CCQENu.  Making validation histograms.\n";
    vars.push_back(new Variable("pzmu", "p_{||, #mu} [GeV/c]", dansPzBins, &CVUniverse::GetMuonPz, &CVUniverse::GetMuonPzTrue));
    vars.push_back(new Variable("Emu", "E_{#mu} [GeV]", robsEmuBins, &CVUniverse::GetEmuGeV, &CVUniverse::GetElepTrueGeV));
    vars.push_back(new Variable("Erecoil", "E_{recoil}", robsRecoilBins, &CVUniverse::GetRecoilE, &CVUniverse::Getq0True)); //TODO: q0 is not the same as recoil energy without a spline correction
    vars2D.push_back(new Variable2D(*vars[1], *vars[0]));
  }

  std::vector<Study*> studies;

  CVUniverse* data_universe = new CVUniverse(options.m_data);
  std::vector<CVUniverse*> data_band = {data_universe};
  std::map<std::string, std::vector<CVUniverse*> > data_error_bands;
  data_error_bands["cv"] = data_band;
  
  std::vector<Study*> data_studies;

  //remind myself which variables I am using
  std::cout << "====================================================\n\tUsing the following variables:\n";
  for(auto& var: vars) std::cout << (var->GetName()).c_str() << " ";
  std::cout << "\n====================================================\n\n" << std::endl;

  for(auto& var: vars) var->InitializeMCHists(error_bands, truth_bands);
  for(auto& var: vars) var->InitializeDATAHists(data_band);
  for(auto& var: res) var->InitializeMCHists(error_bands, truth_bands);

  for(auto& var: vars2D) var->InitializeMCHists(error_bands, truth_bands);
  for(auto& var: vars2D) var->InitializeDATAHists(data_band);

  // Loop entries and fill
  try
  {
    CVUniverse::SetTruth(false);
    LoopAndFillEventSelection(options.m_mc, error_bands, vars, vars2D, res, studies, mycuts, model);
    CVUniverse::SetTruth(true);
    LoopAndFillEffDenom(options.m_truth, truth_bands, vars, vars2D, mycuts, model);
    options.PrintMacroConfiguration(argv[0]);
    std::cout << "MC cut summary:\n" << mycuts << "\n";
    mycuts.resetStats();

    CVUniverse::SetTruth(false);
    LoopAndFillData(options.m_data, data_band, vars, vars2D, data_studies, mycuts);
    std::cout << "Data cut summary:\n" << mycuts << "\n";

    //Write MC results
    TFile* mcOutDir = TFile::Open(MC_OUT_FILE_NAME, "RECREATE");
    if(!mcOutDir)
    {
      std::cerr << "Failed to open a file named " << MC_OUT_FILE_NAME << " in the current directory for writing histograms.\n";
      return badOutputFile;
    }

    for(auto& study: studies) study->SaveOrDraw(*mcOutDir);
    for(auto& var: vars) var->WriteMC(*mcOutDir);
    for(auto& var: res){ // write resolution histos
	if(var->resolution && var->resolution->hist->GetEntries() != 0)
	{
	    mcOutDir->cd();
	    var->resolution->hist->SetDirectory(&(*mcOutDir));
	    var->resolution->hist->Write((var->GetName() + "_resolution").c_str());
	}

	if(var->fractionalResolution && var->fractionalResolution->hist->GetEntries() != 0)
	{
	    mcOutDir->cd();
	    var->fractionalResolution->hist->SetDirectory(&(*mcOutDir));
	    var->fractionalResolution->hist->Write((var->GetName() + "_fractional_resolution").c_str());
	}

	if(var->resolution2D && var->resolution2D->hist->GetEntries() != 0)
	{
	    mcOutDir->cd();
	    var->resolution2D->hist->SetDirectory(&(*mcOutDir));
	    var->resolution2D->hist->Write((var->GetName() + "_resolution2D").c_str());
	}
    }
    for(auto& var: vars2D) var->Write(*mcOutDir);

    //Protons On Target
    auto mcPOT = new TParameter<double>("POTUsed", options.m_mc_pot);
    mcPOT->Write();

    PlotUtils::TargetUtils targetInfo;
    assert(error_bands["cv"].size() == 1 && "List of error bands must contain a universe named \"cv\" for the flux integral.");

    for(const auto& var: vars)
    {
      //Flux integral only if systematics are being done (temporary solution)
      util::GetFluxIntegral(*error_bands["cv"].front(), var->efficiencyNumerator->hist)->Write((var->GetName() + "_reweightedflux_integrated").c_str());
      //Always use MC number of nucleons for cross section
      auto nNucleons = new TParameter<double>((var->GetName() + "_fiducial_nucleons").c_str(), targetInfo.GetTrackerNNucleons(minZ, maxZ, true, apothem));
      nNucleons->Write();
    }

    //Write data results
    TFile* dataOutDir = TFile::Open(DATA_OUT_FILE_NAME, "RECREATE");
    if(!dataOutDir)
    {
      std::cerr << "Failed to open a file named " << DATA_OUT_FILE_NAME << " in the current directory for writing histograms.\n";
      return badOutputFile;
    }

    for(auto& var: vars) var->WriteData(*dataOutDir);

    //Protons On Target
    auto dataPOT = new TParameter<double>("POTUsed", options.m_data_pot);
    dataPOT->Write();

    std::cout << "Success" << std::endl;
  }
  catch(const ROOT::exception& e)
  {
    std::cerr << "Ending on a ROOT error message.  No histograms will be produced.\n"
              << "If the message talks about \"TNetXNGFile\", this could be a problem with dCache.  The message is:\n"
              << e.what() << "\n" << USAGE << "\n";
    return badFileRead;
  }

  return success;
}
