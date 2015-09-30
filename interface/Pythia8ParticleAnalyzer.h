#ifndef Pythia8ParticleAnalyzer_h
#define Pythia8ParticleAnalyzer_h

// system include files
#include <memory>
#include <algorithm>
#include <string>
#include <map>
#include <vector>
#include <iostream>
#include <sstream>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "SimGeneral/HepPDTRecord/interface/ParticleDataTable.h"

#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"

// ROOT unclude files
#include "TFile.h"
#include "TH1F.h"

//
// class declaration
//

class Pythia8ParticleAnalyzer : public edm::EDAnalyzer {
  public:
    explicit Pythia8ParticleAnalyzer(const edm::ParameterSet&);
    ~Pythia8ParticleAnalyzer();

    void analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup);
    //void analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) override;
    void endJob();


  private:
    // Return the name of the particle given a pdgId
    std::string getParticleName(int id) const;

    // Return true if the pdgId corresponds to any lepton
    bool isLepton(int pdgId) const;

    // Return true if the pdgId corresponds to any lepton or photon
    bool isLeptonOrPhoton(int pdgId) const;

    // Return true if the pdgId corresponds to any quark or gluon
    bool isQuarkOrGluon(int pdgId) const;

    // Return true if the genParticle has no lepton or photon parent
    bool py8NoLeptonOrPhotonParent(const reco::Candidate* cand) const;

    // Return true if the genParticle is a Pythia8 parton in the hard process
    bool py8Parton(const reco::Candidate* cand) const;

    // Return the particle index as used by Pythia8
    int py8Index(int idx) const;

    // Return the number of particles with pT > min_pT, |eta| < maxAbsEta
    int countParticles(std::vector<const reco::Candidate *>::const_iterator cands_begin,
            std::vector<const reco::Candidate *>::const_iterator cands_end,
            double minPt, double maxAbsEta);

    // Print the particle info to the 'out' stream
    void listParticles(const edm::Handle<reco::CandidateView>& particles, bool doPartonLevel, std::ostringstream& out) const;

    // Make plots
    void plotParticles(const edm::Handle<reco::CandidateView>& particles, bool doPartonLevel);
    void plotParticlesImpl(const std::vector<const reco::Candidate *>& cands, int ih);

    // Configurable parameters
    edm::InputTag src_;
    int numberShowInfo_;
    int numberShowProcess_;
    int numberShowEvent_;
    bool useMessageLogger_;
    std::string histFile_;
    bool makePlots_;

    // For event lisitng
    edm::ESHandle<ParticleDataTable> pdt_;
    int nProcessesAnalyzed_;
    int nEventsAnalyzed_;

    // For histograms
    TFile* tfile_;
    std::map<TString, TH1F*>  histos_;
};

#endif  // Pythia8ParticleAnalyzer_h
