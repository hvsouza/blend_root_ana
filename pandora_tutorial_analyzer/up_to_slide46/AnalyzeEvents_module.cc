////////////////////////////////////////////////////////////////////////
// Class:       AnalyzeEvents
// Plugin Type: analyzer (Unknown Unknown)
// File:        AnalyzeEvents_module.cc
//
// Generated at Thu Apr  6 12:38:09 2023 by Henrique Vieira de Souza using cetskelgen
// from  version .
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// Additional framwork includes
#include "art_root_io/TFileService.h"
#include "canvas/Persistency/Common/FindManyP.h"

// Additional LArSoft includes
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardataobj/AnalysisBase/T0.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Slice.h"

// ROOT includes
#include <TH1F.h>
#include "TTree.h"

// STL includes
#include <string>
#include <vector>

namespace test {
  class AnalyzeEvents;
}


class test::AnalyzeEvents : public art::EDAnalyzer {
public:
  explicit AnalyzeEvents(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  AnalyzeEvents(AnalyzeEvents const&) = delete;
  AnalyzeEvents(AnalyzeEvents&&) = delete;
  AnalyzeEvents& operator=(AnalyzeEvents const&) = delete;
  AnalyzeEvents& operator=(AnalyzeEvents&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

private:

  // Declare member data here.

  //Create output TTree
  TTree *fTree = nullptr;
  TH1F *fTrackLengthHist = nullptr;

  // Tree Variables
  unsigned int fEventID;
  unsigned int fNPFParticles;
  unsigned int fNPrimaries;
  int fNPrimaryDaughters;
  float fT0;


  std::vector<float> fDaughterTrackLengths;
  /* std::vector<bool> fDaughterLongestTrack; */

  /* std::vector<std::vector<float>> fDaughterTrackdEdx; */
  /* std::vector<std::vector<float>> fDaughterTrackResidualRange; */

  // Declare input labels
  const std::string fPFParticleLabel;
  const std::string fTrackLabel;
  /* const std::string fCaloLabel; */
  /* const std::string fSliceLabel; */
  /* const std::string fOptLabel; */


};


test::AnalyzeEvents::AnalyzeEvents(fhicl::ParameterSet const& p)
  : EDAnalyzer{p},
  fPFParticleLabel(p.get<std::string>("PFParticleLabel")),
  fTrackLabel(p.get<std::string>("TrackLabel"))
  // More initializers here.
{
  // Call appropriate consumes<>() for any products to be retrieved by this module.
}

void test::AnalyzeEvents::analyze(art::Event const& e)
{
  // Implementation of required member function here.
  // Increment the event ID
  fEventID = e.id().event();


  // Reset all of our variables to 0 or empty vectors
  // This ensures things are not kept from the previous event
  fNPFParticles      = 0;
  fNPrimaries        = 0;
  fNPrimaryDaughters = 0;
  fDaughterTrackLengths.clear();

  // Load the PFParticles from pandora
  art::Handle<std::vector<recob::PFParticle>> pfpHandle;
  std::vector<art::Ptr<recob::PFParticle>> pfpVec;
  if (e.getByLabel(fPFParticleLabel, pfpHandle))
    art::fill_ptr_vector(pfpVec, pfpHandle);

  // If there are no PFParticles then give up and skip the event
  if (pfpVec.empty())
    return;
  

  // Initialise neutrino ID
  size_t neutrinoID(std::numeric_limits<size_t>::max());
  
  // Loop over the PFParticles and find the neutrino
  for (const art::Ptr<recob::PFParticle>& pfp : pfpVec) {
    fNPFParticles++;

    // Check that we are looking at a primary that has a neutrino PDG code, if not move on to the next PFP
    if (!(pfp->IsPrimary() && (std::abs(pfp->PdgCode()) == 14 || std::abs(pfp->PdgCode()) == 12)))
      continue;

    neutrinoID = pfp->Self();
    fNPrimaryDaughters = pfp->NumDaughters();
    fNPrimaries++;
  }

  // Check that we found a reconstructed neutrino, if not skip the event
  if (neutrinoID == std::numeric_limits<size_t>::max())
    return;

  // Load the associations between PFPs, Tracks and Calorimetries
  art::FindManyP<recob::Track> pfpTrackAssns(pfpVec, e, fTrackLabel);
  /* art::FindManyP<anab::Calorimetry> trackCaloAssns(trackVec, e, fCaloLabel); */


  // Now access the slices and corresponding timing information
  for (const art::Ptr<recob::PFParticle>& pfp : pfpVec) {
    // Start by assessing the neutrino PFParticle itself
    /* if(pfp->Self() != neutrinoID) continue; */
    // When interested in neutrino's daughters 
    if(pfp->Parent() != neutrinoID) continue;

    // Get the tracks associated with this PFParticle
    const std::vector<art::Ptr<recob::Track>> pfpTracks(pfpTrackAssns.at(pfp.key()));

    // There should only ever be 0 or 1 tracks associated with a sigle PFParticle
    if (pfpTracks.size() == 1) {
      // Get the first (only) element of the vector
      const art::Ptr<recob::Track>& pfpTrack(pfpTracks.front());

      // Add parameters from the track to the branch vector
      fDaughterTrackLengths.push_back(pfpTrack->Length());

      // Fill the histogram with the track length
      fTrackLengthHist->Fill(pfpTrack->Length());



    } // PFParticle Track
  } // PFParticles


  // Store the outputs in the TTree
  fTree->Fill();
}

void test::AnalyzeEvents::beginJob()
{
  // Implementation of optional member function here.
  // Get the TFileService to create the output TTree for us
  art::ServiceHandle<art::TFileService> tfs;
  fTree = tfs->make<TTree>("tree","Output TTree");
  fTrackLengthHist = tfs->make<TH1F>("trackLegthHist", "Reconstruced track lenghts; Trach Length [cm]", 20, 0, 350);

  // Add branches to TTree
  fTree->Branch("eventID", &fEventID);
  fTree->Branch("nPFParticles", &fNPFParticles);
  fTree->Branch("nPrimaries", &fNPrimaries);
  fTree->Branch("nPrimaryDaughters", &fNPrimaryDaughters);
  fTree->Branch("daughterTrackLengths", &fDaughterTrackLengths);
}

void test::AnalyzeEvents::endJob()
{
  // Implementation of optional member function here.
}

DEFINE_ART_MODULE(test::AnalyzeEvents)
