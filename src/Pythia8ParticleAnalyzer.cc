#include "NWU/GenAnalyzers/interface/Pythia8ParticleAnalyzer.h"


Pythia8ParticleAnalyzer::Pythia8ParticleAnalyzer(const edm::ParameterSet& iConfig) :
  src_(iConfig.getParameter<edm::InputTag>("src")),
  numberShowInfo_(iConfig.getUntrackedParameter<int>("numberShowInfo",1)),
  numberShowProcess_(iConfig.getUntrackedParameter<int>("numberShowProcess",1)),
  numberShowEvent_(iConfig.getUntrackedParameter<int>("numberShowEvent",1)),
  useMessageLogger_(iConfig.getUntrackedParameter<bool>("useMessageLogger", false))
{
    nProcessesAnalyzed_ = 0;
    nEventsAnalyzed_ = 0;
}

Pythia8ParticleAnalyzer::~Pythia8ParticleAnalyzer()
{

}

// Stolen from PhysicsTools/HepMCCandAlgos/plugins/ParticleListDrawer.cc
std::string Pythia8ParticleAnalyzer::getParticleName(int id) const
{
    const ParticleData * pd = pdt_->particle( id );
    if (!pd) {
        std::ostringstream ss;
        ss << "P" << id;
        return ss.str();
    } else {
        return pd->name();
    }
}

bool Pythia8ParticleAnalyzer::py8Parton(const reco::Candidate* cand) const {
    int pdgId = cand->pdgId();
    int status = cand->status();
    if (status > 2 && status < 30)
        return true;
    if ((11 <= std::abs(pdgId) && std::abs(pdgId) <= 16) || std::abs(pdgId) == 21)
        if (status == 1 || status == 2)
            if (cand->numberOfMothers() > 0 && cand->mother(0)->pdgId() != pdgId &&
                cand->mother(0)->status() > 4 && cand->mother(0)->status() < 200)
                return true;
    return false;
}

int Pythia8ParticleAnalyzer::py8Index(int index) const {
    return index+1;
}

// Stolen from Pythia8 source in src/Event.cc
void Pythia8ParticleAnalyzer::listParticles(const edm::Handle<reco::CandidateView>& particles,
                                            bool doPartonLevel, std::ostringstream& out) const
{
    using namespace std;

    std::string headerIn = doPartonLevel ? "(hard process)" : "(complete event)";
    std::string headerList = "----------------------------------------";
    headerList.replace(0, headerIn.length() + 2, headerIn + "  ");

    int pdgId = 0;
    std::string particleName = "";

    int idx  = -1;
    int iMo1 = -1;
    int iMo2 = -1;
    int iDa1 = -1;
    int iDa2 = -1;
    int nMo = -1;
    int nDa = -1;

    std::map<const reco::Candidate *, int> candmap;
    std::map<const reco::Candidate *, int>::const_iterator found = candmap.begin();

    reco::Candidate::LorentzVector pSum;
    float chargeSum = 0.;


    // Header
    out << "\n --------  PYTHIA Event Listing  " << headerList << "----------"
        << "-------------------------------------------------\n \n    no    "
        << "    id   name            status     mothers   daughters     colou"
        << "rs      p_x        p_y        p_z         e          m \n";
    out << "     0        90   (system)           -11     0     0     0     0     0     0      0.000      0.000      0.000   8000.000   8000.000\n";  // FIXME: temporary fix

    // At high energy switch to scientific format for momenta.
    bool useFixed = (particles->empty() || particles->at(0).energy() < 1e5);

    // Build a map for quick find
    for (reco::CandidateView::const_iterator p = particles->begin();
        p != particles->end(); ++p) {

        if (doPartonLevel && !py8Parton(&*p))
            continue;

        candmap[&*p] = p - particles->begin();
    }

    for (reco::CandidateView::const_iterator p = particles->begin();
        p != particles->end(); ++p) {

        if (doPartonLevel && !py8Parton(&*p))
            continue;

        // Particle name
        pdgId = p->pdgId();
        particleName = getParticleName(pdgId);

        // Particle index
        idx = p - particles->begin();

        // Particle mothers and daughters
        nMo  = p->numberOfMothers();
        nDa  = p->numberOfDaughters();

        found = candmap.find(p->mother(0));
        iMo1 = (found != candmap.end()) ? found->second : -1;

        found = candmap.find(p->mother(nMo-1));
        iMo2 = (found != candmap.end()) ? found->second : -1;

        found = candmap.find(p->daughter(0));
        iDa1 = (found != candmap.end()) ? found->second : -1;

        found = candmap.find(p->daughter(nDa-1));
        iDa2 = (found != candmap.end()) ? found->second : -1;

        const int col=0, acol=0;  // colors information not stored

        // Basic line for a particle, always printed.
        out << setw(6) << py8Index(idx) << setw(10) << pdgId << "   " << left
            << setw(18) << particleName << right << setw(4)
            << p->status() << setw(6) << py8Index(iMo1) << setw(6)
            << py8Index(iMo2) << setw(6) << py8Index(iDa1) << setw(6)
            << py8Index(iDa2) << setw(6) << col << setw(6) << acol
            << ( (useFixed) ? fixed : scientific ) << setprecision(3)
            << setw(11) << p->px() << setw(11) << p->py() << setw(11)
            << p->pz() << setw(11) << p->energy() << setw(11) << p->mass() << "\n";

        // Statistics on momentum and charge.
        if (p->status() == 1) {
            pSum += p->p4();
            chargeSum += float(p->threeCharge())/3;
        }
    }

    // Line with sum charge, momentum, energy and invariant mass.
    out << fixed << setprecision(3) << "                                   "
        << "Charge sum:" << setw(7) << chargeSum << "           Momentum sum:"
        << ( (useFixed) ? fixed : scientific ) << setprecision(3)
        << setw(11) << pSum.px() << setw(11) << pSum.py() << setw(11)
        << pSum.pz() << setw(11) << pSum.energy() << setw(11) << pSum.mass()
        << "\n";

    // Listing finished.
    out << "\n --------  End PYTHIA Event Listing  ----------------------------"
        << "-------------------------------------------------------------------"
        << "\n";
}

void Pythia8ParticleAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    // Retrieve genParticles
    edm::Handle<reco::CandidateView> particles;
    iEvent.getByLabel(src_, particles);
    iSetup.getData(pdt_);

    std::ostringstream out;

    if (numberShowProcess_ < 0 || nProcessesAnalyzed_ < numberShowProcess_) {
        listParticles(particles, 1, out);

        nProcessesAnalyzed_++;
    }

    if (numberShowEvent_ < 0 || nEventsAnalyzed_ < numberShowEvent_) {
        listParticles(particles, 0, out);

        nEventsAnalyzed_++;
    }

    if (useMessageLogger_)
        edm::LogVerbatim("Pythia8ParticleAnalyzer") << out.str();
    else
        std::cout << out.str();
}

//define this as a plug-in
DEFINE_FWK_MODULE(Pythia8ParticleAnalyzer);
