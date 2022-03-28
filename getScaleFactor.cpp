#define ALPHA_OUT_FILE_NAME "scaleFactor.root"

#define USAGE \
    "\n*** USAGE ***\n" \
    "getScaleFactor <dataPlaylist.txt> <mcPlaylist.txt>\n\n" \
    "*** Explanation ***\n" \
    "Study of pion energy scale factor in tracker\n" \
    "Eats MC files, compares truth to reco Epi, and minimises a function\n" \
    "with a Legendre multiplier alpha. See docDB 18494.\n"

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

//ROOT includes
#include "TParameter.h"
#include "TMinuit.h"
#include "TF1.h"
#include "TGraph.h"

//c++ includes
#include <iostream>
#include <cstdlib> //getenv()

// global vector to hold all the ratios
std::vector< double > ratioVec;

//==============================================================================
// Plot the error function wrt alpha
//==============================================================================

TGraph *alphafun( ){

    Int_t np = ratioVec.size();
    Int_t ns = 10000;

    // we gotta create a loss function that accepts alpha as param, set a value, eval,
    // and get alpha ==> x, eval ==> y.

    TF1 *ftemp = new TF1("ftemp", "( TMath::ATan( [0] * x ) - TMath::Pi() / 4 ) * ( TMath::ATan( [0] * x ) - TMath::Pi() / 4 )", -1.0, 4.0);

    Double_t x[ns], y[ns];
    Double_t amin = -1.0, amax = 4.0;

    for( Int_t i = 0; i < ns; i++ ){
	Double_t atemp = amin + i * ( amax - amin ) / ns;
	ftemp->SetParameter(0, atemp);

	x[i] = 0.0;
	y[i] = 0.0;
	
	x[i] = atemp;
	for( Int_t j = 0; j < np; j++ ){
	    Double_t val = ratioVec.at(j);
	    y[i] += ftemp->Eval(val);
	}

	if( i % 500 == 0 || i == ns - 1 ){
	    std::cout << i << " / " << ns << "\r" << std::flush;
	}
    }

    TGraph *graph = new TGraph( ns, x, y );
    graph->SetTitle( "Calorimetric scale minimisation; #alpha; f(#alpha)" );

    delete ftemp;
    return graph;
    
}

//==============================================================================
// Minimisation function
//==============================================================================

// f is the loss function to be minimised
// ffit takes the ratio Ereco / Etrue as argument
void funcform( Int_t &npar, double *gin, double &f, double *par, int iflag ){ // why gin?
    Double_t alpha = par[0];
    f = 0.0;
    TF1 *ffit = new TF1("minfun", "( TMath::ATan( [0] * x ) - TMath::Pi() / 4 ) * ( TMath::ATan( [0] * x ) - TMath::Pi() / 4 )", -1.0, 4.0);
    ffit->SetParameter( 0, alpha );

    // I will pass the ratios as a global vector
    Int_t Npoints = ratioVec.size();
    for( Int_t k = 0; k < Npoints; k++ ){
	Double_t rat = ratioVec.at(k);
	Double_t z = ffit->Eval(rat);
	f += z; // simple loss function
    }
}

// minuit wants this wrapper
TF1 *minfun( Double_t *sv, Double_t *par, Double_t *err, Double_t *chisqr, Int_t *NDF ){
    TMinuit minuit(1); // 1 parameter
    TF1 *ffitold = (TF1*) gROOT->GetListOfFunctions()->FindObject("minfun");
    if( ffitold ) delete ffitold;

    minuit.SetFCN( funcform );
    Int_t ierflag = 0; Double_t arglist[10];
    minuit.mnparm( 0, "Alpha", sv[0], 0.001, 0.5, 20.0, ierflag );
    arglist[0] = 2; // minimisation strategy = "ACCURATE"
    minuit.mnexcm("SET STR", arglist, 1, ierflag);
    minuit.SetErrorDef(1.); // 1 sigma tolerance
    minuit.SetPrintLevel(0);
    Double_t maxcalls = 10000; Double_t tol = 1.0e-5;
    arglist[0] = maxcalls; arglist[1] = tol;
    minuit.mnexcm("MINIMIZE",arglist,2,ierflag);

    // fit done, now obtain parameters and the outgoing fitted function
    TF1 *ffit = ( TF1* ) gROOT->GetListOfFunctions()->FindObject( "minfun" );
    Double_t finalpar[1], finalerr[1];
    Double_t fmin, fedm, errdef; Int_t npari, nparx, istat;
    
    minuit.GetParameter( 0, finalpar[0], finalerr[0] );
    minuit.mnstat( fmin, fedm, errdef, npari, nparx, istat );
    ffit->SetParameters(finalpar); ffit->SetParErrors(finalerr);
    ffit->GetParameters(finalpar);
    err[0] = ffit->GetParError(0);

    chisqr[0] = fmin; ffit->SetChisquare(fmin);
    NDF[0] = npari; ffit->SetNDF(npari);

    return ffit;
}

//==============================================================================
// Loop and Fill
//==============================================================================
void LoopAndFillScaleFactor(
    PlotUtils::ChainWrapper* chain,
    std::map<std::string, std::vector<CVUniverse*> > error_bands,
    std::vector<Variable*> vars,
    std::vector<Study*> studies,
    PlotUtils::Cutter<CVUniverse, MichelEvent>& michelcuts,
    PlotUtils::Model<CVUniverse, MichelEvent>& model,
    double scaleFactor = 1.0)
{
  assert(!error_bands["cv"].empty() && "\"cv\" error band is empty!  Can't set Model weight.");
  auto& cvUniv = error_bands["cv"].front();

  std::cout << "Starting scale factor loop...\n";
  const int nEntries = chain->GetEntries();
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

		  // obtain reco over truth pion energy
		  const double trueVar = var->GetTrueValue(*universe);
		  const double recoVar = var->GetRecoValue(*universe);

		  if( trueVar == 0.0 ) continue;

		  const double ratio = recoVar / trueVar;

		  /*
		  std::cerr << "\n\n True | Reco | alpha | a*ratio = " << trueVar << " | " << recoVar << " | " << scaleFactor << " | " << scaleFactor*ratio << std::endl;
		  std::cerr << "KILLING THE PROCESS NOW." << std::endl; exit(0);
		  */
		  
		  ratioVec.emplace_back( ratio );
		  var->selectedSignalReco->FillUniverse(universe, ratio, weight);
	      }
          }
        }
        else
        {
          int bkgd_ID = -1;
          if (universe->GetCurrent()==2)bkgd_ID=0;
          else bkgd_ID=1;

          for(auto& var: vars) (*var->m_backgroundHists)[bkgd_ID].FillUniverse(universe, var->GetRecoValue(*universe), weight);
        }
      } // End band's universe loop
    } // End Band loop
    
  } //End entries loop
  std::cout << "Finished MC reco loop.\n";
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
  const double minEnu = 2.0, maxEnu = 50.0, medEnu = 20.0; //GeV

  const int distmm = 200;
  const double vtxELow = 10.0, vtxEHigh = 60.0;

  //==========================================================
  // The actual cuts. Use Jreco and Jtruth namespaces.
  //==========================================================

  preCuts.emplace_back(new Jreco::ZRange<CVUniverse, MichelEvent>("Tracker", minZ, maxZ));
  preCuts.emplace_back(new Jreco::Apothem<CVUniverse, MichelEvent>(apothem));
  preCuts.emplace_back(new Jreco::NoDeadtime<CVUniverse, MichelEvent>(1, "Deadtime"));
  preCuts.emplace_back(new Jreco::HasMINOSMatch<CVUniverse, MichelEvent>());
  preCuts.emplace_back(new Jreco::MINOSNumu<CVUniverse, MichelEvent>());
  preCuts.emplace_back(new Jreco::OneHadron<CVUniverse, MichelEvent>());
  //preCuts.emplace_back(new Jreco::ContainedHadron<CVUniverse, MichelEvent>());
  preCuts.emplace_back(new Jreco::IsPion<CVUniverse, MichelEvent>());
  preCuts.emplace_back(new Jreco::EnuRange<CVUniverse, MichelEvent>(Form("%1.1f <= Enu [GeV] <= %1.1f", minEnu, medEnu), minEnu, medEnu)); // 2 - 20 GeV
  preCuts.emplace_back(new Jreco::VtxECut<CVUniverse, MichelEvent>(Form("%1.1f <= E_vtx (%d mm) <= %1.1f [MeV]", vtxELow, distmm, vtxEHigh), vtxELow, vtxEHigh));

  //==========================================================
  
  signalDefinition.emplace_back(new Jtruth::IsNumu<CVUniverse>());
  signalDefinition.emplace_back(new Jtruth::IsCC<CVUniverse>());
  //------------------------------------------
  signalDefinition.emplace_back(new Jtruth::IsCOH<CVUniverse>());
  //------------------------------------------
  
  //==========================================================
  phaseSpace.emplace_back(new Jtruth::ZRange<CVUniverse>("Tracker", minZ, maxZ));
  phaseSpace.emplace_back(new Jtruth::Apothem<CVUniverse>(apothem));
  //phaseSpace.emplace_back(new Jtruth::MuonAngle<CVUniverse>(maxMuTheta));
  phaseSpace.emplace_back(new Jtruth::EnuRange<CVUniverse>(Form("%1.1f <= Enu [GeV] <= %1.1f", minEnu, medEnu), minEnu, medEnu)); // 2 - 20 GeV
  //phaseSpace.emplace_back(new Jtruth::WRange<CVUniverse>(Form("%1.1f <= W [GeV] < %1.1f", minW, maxW), minW, maxW));
  //==========================================================
  
  PlotUtils::Cutter<CVUniverse, MichelEvent> mycuts(std::move(preCuts), std::move(sidebands) , std::move(signalDefinition),std::move(phaseSpace));

  // need these for event loop
  std::vector<std::unique_ptr<PlotUtils::Reweighter<CVUniverse, MichelEvent>>> MnvTunev1;
  MnvTunev1.emplace_back(new PlotUtils::FluxAndCVReweighter<CVUniverse, MichelEvent>());
  MnvTunev1.emplace_back(new PlotUtils::GENIEReweighter<CVUniverse, MichelEvent>(true, false));
  MnvTunev1.emplace_back(new PlotUtils::LowRecoil2p2hReweighter<CVUniverse, MichelEvent>());
  MnvTunev1.emplace_back(new PlotUtils::MINOSEfficiencyReweighter<CVUniverse, MichelEvent>());
  MnvTunev1.emplace_back(new PlotUtils::RPAReweighter<CVUniverse, MichelEvent>());

  PlotUtils::Model<CVUniverse, MichelEvent> model(std::move(MnvTunev1));
  
  std::map< std::string, std::vector<CVUniverse*> > error_bands;
  if(false) error_bands = GetStandardSystematics(options.m_mc);
  else{
    std::map<std::string, std::vector<CVUniverse*> > band_flux = PlotUtils::GetFluxSystematicsMap<CVUniverse>(options.m_mc, CVUniverse::GetNFluxUniverses());
    error_bands.insert(band_flux.begin(), band_flux.end()); //Necessary to get flux integral later...
  }
  error_bands["cv"] = {new CVUniverse(options.m_mc)};
  std::map< std::string, std::vector<CVUniverse*> > truth_bands;
  if(false) truth_bands = GetStandardSystematics(options.m_truth);
  truth_bands["cv"] = {new CVUniverse(options.m_truth)};

  // I am just getting pion energy and comparing to truth.

  std::vector<double> johnsEPiBins, johnsRatioBins;
  for(Int_t j = 0; j < 21; j++){
      double elo = j * 100.0; johnsEPiBins.push_back(elo);
      double ero = j * 0.10; johnsRatioBins.push_back(ero);
  }

  std::vector<Variable*> vars = {
      new Variable("Epi", "Epi", johnsRatioBins, &CVUniverse::GetEpi, &CVUniverse::GetEpiTrue)
  };

  for(auto& var: vars) var->InitializeMCHists(error_bands, truth_bands);

  std::vector<Study*> studies;

  // Loop entries and fill
  try
  {
    CVUniverse::SetTruth(false);
    // pass over all signal events and fill vector of ratios
    LoopAndFillScaleFactor(options.m_mc, error_bands, vars, studies, mycuts, model);
    // now minimise the function
    Double_t par[1], err[1], chisqr[1]; Int_t ndf[1];
    Double_t sv[] = {1.0};
    TF1 *fun = minfun( sv, par, err, chisqr, ndf );
    std::cout << "\n\n ALPHA = " << fun->GetParameter(0) << " \u00b1 " << fun->GetParError(0) << "\n\n";

    // let's plot this too!
    TGraph *afun = alphafun();
    
    options.PrintMacroConfiguration(argv[0]);
    std::cout << "MC cut summary:\n" << mycuts << "\n";
    mycuts.resetStats();

    //Write MC results
    TFile* mcOutDir = TFile::Open(ALPHA_OUT_FILE_NAME, "RECREATE");
    if(!mcOutDir)
    {
      std::cerr << "Failed to open a file named " << ALPHA_OUT_FILE_NAME << " in the current directory for writing histograms.\n";
      return badOutputFile;
    }
    
    for( auto &var: vars ){
	// write directly
	if( var->selectedSignalReco ){
	    var->selectedSignalReco->SyncCVHistos();
	    mcOutDir->cd();
	    
	    var->selectedSignalReco->hist->SetDirectory( mcOutDir );
	    var->selectedSignalReco->hist->Write( (var->GetName() + "_scale").c_str() );
	}
    }

    if(afun){
	afun->Write( "Epi_minimisation" );
    }

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
