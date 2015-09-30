#ifndef Pythia8ParticleAnalyzer_h
#define Pythia8ParticleAnalyzer_h

// system include files
#include <memory>
#include <algorithm>
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

//
// class declaration
//

class Pythia8ParticleAnalyzer : public edm::EDAnalyzer {
  public:
    explicit Pythia8ParticleAnalyzer(const edm::ParameterSet&);
    ~Pythia8ParticleAnalyzer();
    void analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup);
    //void analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) override;

  private:
    std::string getParticleName(int id) const;
    bool py8Parton(const reco::Candidate* cand) const;
    int py8Index(int idx) const;
    void listParticles(const edm::Handle<reco::CandidateView>& particles, bool doPartonLevel, std::ostringstream& out) const;

    edm::ESHandle<ParticleDataTable> pdt_;
    int nProcessesAnalyzed_;
    int nEventsAnalyzed_;

    edm::InputTag src_;
    int numberShowInfo_;
    int numberShowProcess_;
    int numberShowEvent_;
    bool useMessageLogger_;
};

#endif  // Pythia8ParticleAnalyzer_h
