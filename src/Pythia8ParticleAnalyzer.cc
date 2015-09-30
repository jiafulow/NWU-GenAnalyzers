#include "NWU/GenAnalyzers/interface/Pythia8ParticleAnalyzer.h"

#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/deltaPhi.h"

//
// See header file for documentation.
//

Pythia8ParticleAnalyzer::Pythia8ParticleAnalyzer(const edm::ParameterSet& iConfig) :
  src_(iConfig.getParameter<edm::InputTag>("src")),
  numberShowInfo_(iConfig.getUntrackedParameter<int>("numberShowInfo",1)),
  numberShowProcess_(iConfig.getUntrackedParameter<int>("numberShowProcess",1)),
  numberShowEvent_(iConfig.getUntrackedParameter<int>("numberShowEvent",1)),
  useMessageLogger_(iConfig.getUntrackedParameter<bool>("useMessageLogger", false)),
  histFile_(iConfig.getUntrackedParameter<std::string>("histFile", "histograms.root")),
  makePlots_(iConfig.getUntrackedParameter<bool>("makePlots", true))
{
    nProcessesAnalyzed_ = 0;
    nEventsAnalyzed_ = 0;

    tfile_ = 0;
    if (makePlots_) {
        tfile_ = TFile::Open(histFile_.c_str(), "RECREATE");
        assert(tfile_);

        TH1::AddDirectory(kFALSE);

        TString hname = "";
        for (int ih=0; ih<20; ih++) {
            hname = Form("pt_top_%i", ih);
            histos_[hname] = new TH1F(hname, "; p_{T}(top) [GeV]", 20, 0, 200);
            hname = Form("pt_higgs_%i", ih);
            histos_[hname] = new TH1F(hname, "; p_{T}(higgs) [GeV]", 20, 0, 200);
            hname = Form("pt_mu1_%i", ih);
            histos_[hname] = new TH1F(hname, "; p_{T}(#mu_{1}) [GeV]", 20, 0, 200);
            hname = Form("pt_mu2_%i", ih);
            histos_[hname] = new TH1F(hname, "; p_{T}(#mu_{2}) [GeV]", 20, 0, 200);
            hname = Form("pt_mumu_%i", ih);
            histos_[hname] = new TH1F(hname, "; p_{T}(#mu#mu) [GeV]", 20, 0, 200);
            hname = Form("pt_j1_%i", ih);
            histos_[hname] = new TH1F(hname, "; p_{T}(j_{1}) [GeV]", 20, 0, 200);
            hname = Form("pt_j2_%i", ih);
            histos_[hname] = new TH1F(hname, "; p_{T}(j_{2}) [GeV]", 20, 0, 200);
            hname = Form("pt_jj_%i", ih);
            histos_[hname] = new TH1F(hname, "; p_{T}(jj) [GeV]", 20, 0, 200);

            hname = Form("eta_top_%i", ih);
            histos_[hname] = new TH1F(hname, "; #eta(top)", 40, -2.4, 2.4);
            hname = Form("eta_higgs_%i", ih);
            histos_[hname] = new TH1F(hname, "; #eta(higgs)", 40, -2.4, 2.4);
            hname = Form("eta_mu1_%i", ih);
            histos_[hname] = new TH1F(hname, "; #eta(#mu_{1})", 40, -2.4, 2.4);
            hname = Form("eta_mu2_%i", ih);
            histos_[hname] = new TH1F(hname, "; #eta(#mu_{2})", 40, -2.4, 2.4);
            hname = Form("eta_mumu_%i", ih);
            histos_[hname] = new TH1F(hname, "; #eta(#mu#mu)", 40, -2.4, 2.4);
            hname = Form("eta_j1_%i", ih);
            histos_[hname] = new TH1F(hname, "; #eta(j_{1})", 40, -4.8, 4.8);
            hname = Form("eta_j2_%i", ih);
            histos_[hname] = new TH1F(hname, "; #eta(j_{2})", 40, -4.8, 4.8);
            hname = Form("eta_jj_%i", ih);
            histos_[hname] = new TH1F(hname, "; #eta(jj)", 40, -4.8, 4.8);

            hname = Form("pdgId_top_%i", ih);
            histos_[hname] = new TH1F(hname, "; pdgId(top)", 61, -30.5, 30.5);
            hname = Form("pdgId_j1_%i", ih);
            histos_[hname] = new TH1F(hname, "; pdgId(j_{1})", 61, -30.5, 30.5);
            hname = Form("pdgId_j2_%i", ih);
            histos_[hname] = new TH1F(hname, "; pdgId(j_{2})", 61, -30.5, 30.5);

            hname = Form("mass_mumu_%i", ih);
            histos_[hname] = new TH1F(hname, "; M(#mu#mu) [GeV]", 29, 12, 70);
            hname = Form("mass_jj_%i", ih);
            histos_[hname] = new TH1F(hname, "; M(jj) [GeV]", 15, 0, 600);
            hname = Form("mass_mumujj_%i", ih);
            histos_[hname] = new TH1F(hname, "; M(#mu#mujj) [GeV]", 20, 0, 800);
            hname = Form("mass_mumuj1_%i", ih);
            histos_[hname] = new TH1F(hname, "; M(#mu#muj_{1}) [GeV]", 20, 0, 500);
            hname = Form("mass_mumuj2_%i", ih);
            histos_[hname] = new TH1F(hname, "; M(#mu#muj_{2}) [GeV]", 20, 0, 500);

            hname = Form("dR_mu1_mu2_%i", ih);
            histos_[hname] = new TH1F(hname, "; #DeltaR(#mu_{1},#mu_{2})", 12, 0., 6./3);
            hname = Form("dR_j1_j2_%i", ih);
            histos_[hname] = new TH1F(hname, "; #DeltaR(j_{1},j_{2})", 12, 0., 6.);
            hname = Form("dR_mumu_jj_%i", ih);
            histos_[hname] = new TH1F(hname, "; #DeltaR(#mu#mu,jj)", 12, 0., 6.);
            hname = Form("dR_mumu_j1_%i", ih);
            histos_[hname] = new TH1F(hname, "; #DeltaR(#mu#mu,j_{1})", 12, 0., 6.);
            hname = Form("dR_mumu_j2_%i", ih);
            histos_[hname] = new TH1F(hname, "; #DeltaR(#mu#mu,j_{2})", 12, 0., 6.);

            hname = Form("njets_pt20_%i", ih);
            histos_[hname] = new TH1F(hname, "; # of jets", 12, 0, 12);
            hname = Form("njets_pt25_%i", ih);
            histos_[hname] = new TH1F(hname, "; # of jets", 12, 0, 12);
            hname = Form("njets_pt30_%i", ih);
            histos_[hname] = new TH1F(hname, "; # of jets", 12, 0, 12);
            hname = Form("njets_pt40_%i", ih);
            histos_[hname] = new TH1F(hname, "; # of jets", 12, 0, 12);
            hname = Form("njets_pt50_%i", ih);
            histos_[hname] = new TH1F(hname, "; # of jets", 12, 0, 12);

            hname = Form("ncjets_pt20_%i", ih);
            histos_[hname] = new TH1F(hname, "; # of jets in |#eta|<2.4", 12, 0, 12);
            hname = Form("ncjets_pt25_%i", ih);
            histos_[hname] = new TH1F(hname, "; # of jets in |#eta|<2.4", 12, 0, 12);
            hname = Form("ncjets_pt30_%i", ih);
            histos_[hname] = new TH1F(hname, "; # of jets in |#eta|<2.4", 12, 0, 12);
            hname = Form("ncjets_pt40_%i", ih);
            histos_[hname] = new TH1F(hname, "; # of jets in |#eta|<2.4", 12, 0, 12);
            hname = Form("ncjets_pt50_%i", ih);
            histos_[hname] = new TH1F(hname, "; # of jets in |#eta|<2.4", 12, 0, 12);
        }
    }  // end if
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

bool Pythia8ParticleAnalyzer::isLepton(int pdgId) const
{
    return ((11 <= std::abs(pdgId) && std::abs(pdgId) <= 16));
}

bool Pythia8ParticleAnalyzer::isLeptonOrPhoton(int pdgId) const
{
    return ((11 <= std::abs(pdgId) && std::abs(pdgId) <= 16) || pdgId == 22);
}

bool Pythia8ParticleAnalyzer::isQuarkOrGluon(int pdgId) const
{
    return ((1 <= std::abs(pdgId) && std::abs(pdgId) <= 5) || pdgId == 21);  // quarks 1-5 only
}

bool Pythia8ParticleAnalyzer::py8NoLeptonOrPhotonParent(const reco::Candidate* cand) const
{
    while (cand->numberOfMothers() > 0) {
        assert(cand->mother(0));
        int momPdgId = cand->mother(0)->pdgId();
        if (isLeptonOrPhoton(momPdgId))
            return false;
        cand = cand->mother(0);
    }
    return true;
}

bool Pythia8ParticleAnalyzer::py8Parton(const reco::Candidate* cand) const
{
    int pdgId = cand->pdgId();
    int status = cand->status();
    if (status > 2 && status < 30)
        return true;
    if (isLepton(pdgId))
        if (status == 1 || status == 2)
            if (cand->numberOfMothers() > 0 && py8NoLeptonOrPhotonParent(cand) &&
                cand->mother(0)->status() > 4 && cand->mother(0)->status() < 200)
                return true;
    return false;
}

int Pythia8ParticleAnalyzer::py8Index(int index) const
{
    return index+1;
}

int Pythia8ParticleAnalyzer::countParticles(std::vector<const reco::Candidate *>::const_iterator cands_begin,
            std::vector<const reco::Candidate *>::const_iterator cands_end,
            double minPt, double maxAbsEta)
{
    int n = 0;

    for (std::vector<const reco::Candidate *>::const_iterator cands_it = cands_begin;
        cands_it != cands_end; ++cands_it) {
        const reco::Candidate * cand = *cands_it;
        if (cand->pt() > minPt && std::abs(cand->eta()) < maxAbsEta)
            n++;
    }
    return n;
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

    // Loop over the particles and list their info
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

void Pythia8ParticleAnalyzer::plotParticles(const edm::Handle<reco::CandidateView>& particles, bool doPartonLevel)
{
    int pdgId = 0;
    int status = 0;

    // A container for the following 6 particles:
    //   0: top
    //   1: higgs
    //   2: mu_1
    //   3: mu_2
    //   4: j_1
    //   5: j_2
    std::vector<const reco::Candidate *> cands(6, NULL);

    int nStatus22 = 0;
    int nStatus23 = 0;
    for (reco::CandidateView::const_iterator p = particles->begin();
        p != particles->end(); ++p) {

        if (doPartonLevel && !py8Parton(&*p))
            continue;

        pdgId = p->pdgId();
        status = p->status();

        if (status == 22)
            nStatus22++;
        if (status == 23 || status == 1)
            nStatus23++;

        // Grab the 6 particles
        if (status == 22) {
            if (std::abs(pdgId) == 6) {  // top
                if (!cands.at(0))
                    cands.at(0) = &*p;
                else
                    assert(false);

            } else if (pdgId == 25) {  // higgs
                if (!cands.at(1))
                    cands.at(1) = &*p;
                else
                    assert(false);
            }

        } else if (status == 23 || status == 1) {
            if (std::abs(pdgId) == 13) {  // muon
                if (!cands.at(2))
                    cands.at(2) = &*p;
                else if (!cands.at(3))
                    cands.at(3) = &*p;
                else
                    assert(false);

            } else if (isQuarkOrGluon(pdgId)) {
                if (!cands.at(4))
                    cands.at(4) = &*p;
                else if (!cands.at(5))
                    cands.at(5) = &*p;
                else
                    assert(false);
            }
        }
    }  // end loop

    // Check if the number of particles make sense
    if (doPartonLevel) {
        if (nStatus22 != 2)
            std::cout << "[WARNING] Unexpected number of particles with status = 22: " << nStatus22 << std::endl;
        if (nStatus23 != 4)
            std::cout << "[WARNING] Unexpected number of particles with status = 23: " << nStatus23 << std::endl;

        for (const auto cand: cands)
            assert(cand);

        // Order by pT
        if (cands.at(2)->pt() < cands.at(3)->pt())
            std::swap(cands.at(2), cands.at(3));

        if (cands.at(4)->pt() < cands.at(5)->pt())
            std::swap(cands.at(4), cands.at(5));

        assert(cands.at(2)->pt() >= cands.at(3)->pt());
        assert(cands.at(4)->pt() >= cands.at(5)->pt());
    }

    const int nmuons_pt25 = countParticles(cands.begin()+2, cands.begin()+3+1, 25, 2.4);
    const int nmuons_pt20 = countParticles(cands.begin()+2, cands.begin()+3+1, 20, 2.4);
    const int nmuons_pt10 = countParticles(cands.begin()+2, cands.begin()+3+1, 10, 2.4);
    const int njets_pt30  = countParticles(cands.begin()+4, cands.begin()+5+1, 30, 5.0);
    const int ncjets_pt30 = countParticles(cands.begin()+4, cands.begin()+5+1, 30, 2.4);

    if (true)
        plotParticlesImpl(cands, 0);

    if (nmuons_pt25 >= 1)
        plotParticlesImpl(cands, 1);

    if (nmuons_pt25 >= 1 && njets_pt30 >= 1)
        plotParticlesImpl(cands, 2);

    if (nmuons_pt25 >= 1 && ncjets_pt30 >= 1)
        plotParticlesImpl(cands, 3);

    if (nmuons_pt20 >= 1 && nmuons_pt10 >= 1)
        plotParticlesImpl(cands, 11);

    if (nmuons_pt20 >= 1 && nmuons_pt10 >= 1 && njets_pt30 >= 1)
        plotParticlesImpl(cands, 12);

    if (nmuons_pt20 >= 1 && nmuons_pt10 >= 1 && ncjets_pt30 >= 1)
        plotParticlesImpl(cands, 13);
}

void Pythia8ParticleAnalyzer::plotParticlesImpl(const std::vector<const reco::Candidate *>& cands, int ih)
{
    const reco::Candidate* top    = cands.at(0);
    const reco::Candidate* higgs  = cands.at(1);
    const reco::Candidate* mu1    = cands.at(2);
    const reco::Candidate* mu2    = cands.at(3);
    const reco::Candidate* j1     = cands.at(4);
    const reco::Candidate* j2     = cands.at(5);
    assert(top && higgs && mu1 && mu2 && j1 && j2);

    const reco::Candidate::LorentzVector mumuP4 = mu1->p4() + mu2->p4();
    const reco::Candidate::LorentzVector jjP4 = j1->p4() + j2->p4();
    const reco::Candidate::LorentzVector mumujjP4 = mumuP4 + jjP4;
    const reco::Candidate::LorentzVector mumuj1P4 = mumuP4 + j1->p4();
    const reco::Candidate::LorentzVector mumuj2P4 = mumuP4 + j2->p4();

    TString hname = "";

    // Fill
    hname = Form("pt_top_%i", ih);
    histos_[hname]->Fill(top->pt());
    hname = Form("pt_higgs_%i", ih);
    histos_[hname]->Fill(higgs->pt());
    hname = Form("pt_mu1_%i", ih);
    histos_[hname]->Fill(mu1->pt());
    hname = Form("pt_mu2_%i", ih);
    histos_[hname]->Fill(mu2->pt());
    hname = Form("pt_mumu_%i", ih);
    histos_[hname]->Fill(mumuP4.pt());
    hname = Form("pt_j1_%i", ih);
    histos_[hname]->Fill(j1->pt());
    hname = Form("pt_j2_%i", ih);
    histos_[hname]->Fill(j2->pt());
    hname = Form("pt_jj_%i", ih);
    histos_[hname]->Fill(jjP4.pt());

    hname = Form("eta_top_%i", ih);
    histos_[hname]->Fill(top->eta());
    hname = Form("eta_higgs_%i", ih);
    histos_[hname]->Fill(higgs->eta());
    hname = Form("eta_mu1_%i", ih);
    histos_[hname]->Fill(mu1->eta());
    hname = Form("eta_mu2_%i", ih);
    histos_[hname]->Fill(mu2->eta());
    hname = Form("eta_mumu_%i", ih);
    histos_[hname]->Fill(mumuP4.eta());
    hname = Form("eta_j1_%i", ih);
    histos_[hname]->Fill(j1->eta());
    hname = Form("eta_j2_%i", ih);
    histos_[hname]->Fill(j2->eta());
    hname = Form("eta_jj_%i", ih);
    histos_[hname]->Fill(jjP4.eta());

    hname = Form("pdgId_top_%i", ih);
    histos_[hname]->Fill(top->pdgId());
    hname = Form("pdgId_j1_%i", ih);
    histos_[hname]->Fill(j1->pdgId());
    hname = Form("pdgId_j2_%i", ih);
    histos_[hname]->Fill(j2->pdgId());

    hname = Form("mass_mumu_%i", ih);
    histos_[hname]->Fill(mumuP4.mass());
    hname = Form("mass_jj_%i", ih);
    histos_[hname]->Fill(jjP4.mass());
    hname = Form("mass_mumujj_%i", ih);
    histos_[hname]->Fill(mumujjP4.mass());
    hname = Form("mass_mumuj1_%i", ih);
    histos_[hname]->Fill(mumuj1P4.mass());
    hname = Form("mass_mumuj2_%i", ih);
    histos_[hname]->Fill(mumuj2P4.mass());

    hname = Form("dR_mu1_mu2_%i", ih);
    histos_[hname]->Fill(reco::deltaR(mu1->p4(), mu2->p4()));
    hname = Form("dR_j1_j2_%i", ih);
    histos_[hname]->Fill(reco::deltaR(j1->p4(), j2->p4()));
    hname = Form("dR_mumu_jj_%i", ih);
    histos_[hname]->Fill(reco::deltaR(mumuP4, jjP4));
    hname = Form("dR_mumu_j1_%i", ih);
    histos_[hname]->Fill(reco::deltaR(mumuP4, j1->p4()));
    hname = Form("dR_mumu_j2_%i", ih);
    histos_[hname]->Fill(reco::deltaR(mumuP4, j2->p4()));

    const int njets_pt20 = countParticles(cands.begin()+4, cands.begin()+5+1, 20, 5.0);
    const int njets_pt25 = countParticles(cands.begin()+4, cands.begin()+5+1, 25, 5.0);
    const int njets_pt30 = countParticles(cands.begin()+4, cands.begin()+5+1, 30, 5.0);
    const int njets_pt40 = countParticles(cands.begin()+4, cands.begin()+5+1, 40, 5.0);
    const int njets_pt50 = countParticles(cands.begin()+4, cands.begin()+5+1, 50, 5.0);
    const int ncjets_pt20 = countParticles(cands.begin()+4, cands.begin()+5+1, 20, 2.4);
    const int ncjets_pt25 = countParticles(cands.begin()+4, cands.begin()+5+1, 25, 2.4);
    const int ncjets_pt30 = countParticles(cands.begin()+4, cands.begin()+5+1, 30, 2.4);
    const int ncjets_pt40 = countParticles(cands.begin()+4, cands.begin()+5+1, 40, 2.4);
    const int ncjets_pt50 = countParticles(cands.begin()+4, cands.begin()+5+1, 50, 2.4);

    hname = Form("njets_pt20_%i", ih);
    histos_[hname]->Fill(njets_pt20);
    hname = Form("njets_pt25_%i", ih);
    histos_[hname]->Fill(njets_pt25);
    hname = Form("njets_pt30_%i", ih);
    histos_[hname]->Fill(njets_pt30);
    hname = Form("njets_pt40_%i", ih);
    histos_[hname]->Fill(njets_pt40);
    hname = Form("njets_pt50_%i", ih);
    histos_[hname]->Fill(njets_pt50);
    hname = Form("ncjets_pt20_%i", ih);
    histos_[hname]->Fill(ncjets_pt20);
    hname = Form("ncjets_pt25_%i", ih);
    histos_[hname]->Fill(ncjets_pt25);
    hname = Form("ncjets_pt30_%i", ih);
    histos_[hname]->Fill(ncjets_pt30);
    hname = Form("ncjets_pt40_%i", ih);
    histos_[hname]->Fill(ncjets_pt40);
    hname = Form("ncjets_pt50_%i", ih);
    histos_[hname]->Fill(ncjets_pt50);
}

void Pythia8ParticleAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    // Retrieve genParticles
    edm::Handle<reco::CandidateView> particles;
    iEvent.getByLabel(src_, particles);
    iSetup.getData(pdt_);

    // Event listing
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

    // Histogram plotting
    if (makePlots_) {
        plotParticles(particles, 1);

    }

}

void Pythia8ParticleAnalyzer::endJob()
{
    if (makePlots_) {
        tfile_->cd();

        for (const auto kv: histos_) {
            kv.second->SetDirectory(gDirectory);
        }
        tfile_->Write();
        tfile_->Close();
    }
}

//define this as a plug-in
DEFINE_FWK_MODULE(Pythia8ParticleAnalyzer);
