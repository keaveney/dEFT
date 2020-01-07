// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/PartonicTops.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Tools/RivetMT2.hh"
#include <bitset>
#include <array>
#include <unordered_set>

namespace Rivet {

  /// @brief Add a short analysis description here
  class LHCb_TOP : public Analysis {
  private:

  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(LHCb_TOP);

    /// Book histograms and initialise projections before the run
    void init() {
     
      //      std::cout<<" in init 1  "<<std::endl;

      // Initialise and register projections
       declare(FinalState(Cuts::abseta < 6 && Cuts::pT > 0.*MeV), "FS");
      //      declare(FinalState(Cuts::abseta < 5 && Cuts::pT > 100*MeV), "FS");                                                                                                                              
      declare(PartonicTops(PartonicTops::E_MU, false), "LeptonicPartonTops");  

      Cut eta_full = Cuts::abseta < 6.0 and Cuts::pT > 0.; // james
      FinalState fs(eta_full);
     
      _h_t_pT  = bookHisto1D(1, 1, 1);
      _h_tbar_pT = bookHisto1D(2,1,1);
      _h_tlead_pT = bookHisto1D(3,1,1);
      _h_ttrail_pT = bookHisto1D(4,1,1);
      _h_tttrf_pT = bookHisto1D(5,1,1);

      _h_t_y = bookHisto1D(6,1,1);
      _h_tbar_y = bookHisto1D(7,1,1);
      _h_tlead_y = bookHisto1D(8,1,1);
      _h_ttrail_y = bookHisto1D(9,1,1);

      _h_tt_pT = bookHisto1D(10,1,1);
      _h_tt_y = bookHisto1D(11,1,1);
      _h_tt_mass = bookHisto1D(12,1,1);
      _h_tt_delphi = bookHisto1D(13,1,1);
      _h_tt_delrap = bookHisto1D(14,1,1);

      _h_l_pT = bookHisto1D(15,1,1);
      _h_lbar_pT = bookHisto1D(16,1,1);
      _h_l_eta = bookHisto1D(17,1,1);
      _h_lbar_eta = bookHisto1D(18,1,1);
      _h_llead_pT = bookHisto1D(19,1,1);
      _h_ltrail_pT = bookHisto1D(20,1,1);
      _h_llead_eta = bookHisto1D(21,1,1);
      _h_ltrail_eta = bookHisto1D(22,1,1);
      _h_llbar_pT = bookHisto1D(23,1,1);
      _h_llbar_mass = bookHisto1D(24,1,1);
      _h_llbar_delphi = bookHisto1D(25,1,1);
      _h_llbar_deleta = bookHisto1D(26,1,1);
      _h_blead_pT = bookHisto1D(27,1,1);
      _h_btrail_pT = bookHisto1D(28,1,1);
      _h_blead_eta = bookHisto1D(29,1,1);
      _h_btrail_eta = bookHisto1D(30,1,1);
      _h_bbbar_mass = bookHisto1D(31,1,1);
      _h_bbbar_pT = bookHisto1D(32,1,1);
      _h_njet = bookHisto1D(33,1,1);
      _h_llbar_absdeleta = bookHisto1D(34,1,1);


      //      std::cout<<" in init 3  "<<std::endl;
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {

          bool debug = false;

          if (debug) std::cout<<" in analyse  "<<std::endl;
          /// @todo Do the event by event analysis here
          const double weight = event.weight();

	  //  std::cout<<" weight =   "<< weight  <<std::endl;

	  _isValid = false;
	  //	  _theParticles.clear();
	  _wDecays1.clear();
	  _wDecays2.clear();
	  _jets.clear();
	  _leptons.clear();
	  _bjets.clear();
	  _ljets.clear();
	  _mode1 = _mode2 = CH_HADRON;

      //////******PRELIMINARY EVENT SELECTION *******/////////////
      // Get the parton level ttbar candidate
      // Do the analysis only for the ttbar full leptonic, without tau decay
      const Particles leptonicpartontops = apply<ParticleFinder>(event, "LeptonicPartonTops").particlesByPt();
      if (debug) std::cout<<" passed 1  "<<std::endl;
   
      //restrict ourselves to the dilepton channels
      if (leptonicpartontops.size() != 2) vetoEvent;
 
	  FourMomentum t1P4_parton;
	  FourMomentum t2P4_parton;

	  t1P4_parton = leptonicpartontops[0].momentum();
	  t2P4_parton = leptonicpartontops[1].momentum();

	  FourMomentum ttbarP4_parton = t1P4_parton + t2P4_parton;

          if (debug) std::cout<<" passed 2  "<<std::endl;

	  //////////////////////pseudo top 

      // Collect final state particles
	  Particles pForLep, pForJet;
	  Particles neutrinos; // Prompt neutrinos
	  foreach (const GenParticle* p, Rivet::particles(event.genEvent())) {
	    const int status = p->status();
	    const int pdgId = p->pdg_id();
	    if (status == 1) {
	      Particle rp(*p);
	      if (!PID::isHadron(pdgId) && !rp.fromHadron()) {
		// Collect particles not from hadron decay
		if (rp.isNeutrino()) {
		  // Prompt neutrinos are kept in separate collection
		  neutrinos.push_back(rp);
		} else if (pdgId == 22 || rp.isLepton()) {
		  // Leptons and photons for the dressing
		  pForLep.push_back(rp);
		}
	      } else if (!rp.isNeutrino()) {
		// Use all particles from hadron decay
		pForJet.push_back(rp);
	      }
	    } else if (PID::isHadron(pdgId) && PID::hasBottom(pdgId)) {
	      // NOTE: Consider B hadrons with pT > 5GeV - not in CMS proposal
	      //if ( p->momentum().perp() < 5 ) continue; 

	      // Do unstable particles, to be used in the ghost B clustering
	      // Use last B hadrons only
	      bool isLast = true;
	      foreach (GenParticle* pp, Rivet::particles(p->end_vertex(), HepMC::children)) {
		if (PID::hasBottom(pp->pdg_id())) {
		  isLast = false;
		  break;
		}
	      }
	      if (!isLast) continue;

	      // Rescale momentum by 10^-20
	      Particle ghost(pdgId, FourMomentum(p->momentum())*1e-20/p->momentum().rho());
	      pForJet.push_back(ghost);
	    }
	  }

	  // Start object building from trivial thing - prompt neutrinos
	  std::sort(neutrinos.begin(), neutrinos.end(), GreaterByPt());

	  // Proceed to lepton dressing
	  FastJets fjLep(FinalState(), FastJets::ANTIKT, _lepR);
	  fjLep.calc(pForLep);
	  //Jets leptons;
	  std::vector<int> leptonsId;
	  std::set<int> dressedIdxs;
	  foreach (const Jet& lep, fjLep.jetsByPt(_lepMinPt)) {
	    if (( lep.eta() < _lepMinEta)  ||  (lep.eta() > _lepMaxEta)) continue;

	    double leadingPt = -1;
	    int leptonId = 0;
	    foreach (const Particle& p, lep.particles()) {
	      dressedIdxs.insert(p.genParticle()->barcode());
	      if (p.isLepton() && p.momentum().pt() > leadingPt) {
		leadingPt = p.momentum().pt();
		leptonId = p.pdgId();
	      }
	    }
	    if (leptonId == 0) continue;
	    _leptons.push_back(lep);
	    leptonsId.push_back(leptonId);
	  }

	  // Re-use particles not used in lepton dressing
	  foreach (const Particle& rp, pForLep) {
	    const int barcode = rp.genParticle()->barcode();
	    // Skip if the particle is used in dressing
	    if (dressedIdxs.find(barcode) != dressedIdxs.end()) continue;

	    // Put back to be used in jet clustering
	    pForJet.push_back(rp);
	  }

	  // Then do the jet clustering
	  FastJets fjJet(FinalState(), FastJets::ANTIKT, _jetR);
	  //fjJet.useInvisibles(); // NOTE: CMS proposal to remove neutrinos
	  fjJet.calc(pForJet);
	  // Jets _bjets, _ljets; FIME
	  foreach (const Jet& jet, fjJet.jetsByPt(_jetMinPt)) {
	    if ((jet.eta() < _jetMinEta) || (jet.eta() > _jetMaxEta)) continue;
	    _jets.push_back(jet);

	    bool isBJet = false;
	    foreach (const Particle& rp, jet.particles()) {
	      if (PID::hasBottom(rp.pdgId())) {
		isBJet = true;
		break;
	      }
	    }

	    if ( isBJet ) _bjets.push_back(jet);
	    else _ljets.push_back(jet);
	  }
	  
	  if (debug) std::cout <<" n lepton  = "  << _leptons.size() <<"  n bjets " <<    _bjets.size() << "  n neutrinos " <<    neutrinos.size() << std::endl;

        
     //event selection
	  //if (_leptons.size() != 2) return;
	  //if (_bjets.size() < 1) return;
        
      // Every building blocks are ready. Continue to pseudo-W and pseudo-top combination

     if (debug) std::cout<<" passed 2+  "<<std::endl;

	  std::map<double, std::pair<size_t, size_t> > wLepCandIdxs;
	  std::map<double, std::pair<size_t, size_t> > wHadCandIdxs;

	  // Collect leptonic-decaying W's
	  for (size_t iLep = 0, nLep = _leptons.size(); iLep < nLep; ++iLep) {
	    const Jet& lep = _leptons.at(iLep);
	    for (size_t iNu = 0, nNu = neutrinos.size(); iNu < nNu; ++iNu) {
	      const Particle& nu = neutrinos.at(iNu);
	      const double m = (lep.momentum()+nu.momentum()).mass();
	      const double dm = std::abs(m-_wMass);
	      wLepCandIdxs[dm] = make_pair(iLep, iNu);
	    }
	  }

	  // Continue to hadronic decaying W's
	  for (size_t i = 0, nLjet = _ljets.size(); i < nLjet; ++i) {
	    const Jet& ljet1 = _ljets[i];
	    for (size_t j = i+1; j < nLjet; ++j) {
	      const Jet& ljet2 = _ljets[j];
	      const double m = (ljet1.momentum()+ljet2.momentum()).mass();
	      const double dm = std::abs(m-_wMass);
	      wHadCandIdxs[dm] = make_pair(i, j);
	    }
	  }

	  // Cleanup W candidate, choose pairs with minimum dm if they share decay products
	  // cleanup(wLepCandIdxs);
	  cleanup(wHadCandIdxs, true);
	  const size_t nWLepCand = wLepCandIdxs.size();
	  const size_t nWHadCand = wHadCandIdxs.size();

	  if (nWLepCand + nWHadCand < 2) return; // this cut should be superfluous, TBC

          if (debug) std::cout<<" passed 3  "<<std::endl;

	  int w1Q = 1, w2Q = -1;
	  int w1dau1Id = 1, w2dau1Id = -1;
	  FourMomentum w1dau1LVec, w1dau2LVec;
	  FourMomentum w2dau1LVec, w2dau2LVec;

          if (debug) std::cout<<" n leptons   "<< _leptons.size()  <<std::endl;

	  // Full leptonic case
	  const auto& idPair1 = wLepCandIdxs.begin()->second;
	  const auto& idPair2 = std::next(wLepCandIdxs.begin())->second;
	  const auto& w1dau1 = _leptons[idPair1.first];
	  const auto& w1dau2 = neutrinos[idPair1.second];
	  const auto& w2dau1 = _leptons[idPair2.first];
	  const auto& w2dau2 = neutrinos[idPair2.second];

	  w1dau1LVec = w1dau1.momentum();
	  w1dau2LVec = w1dau2.momentum();
	  w2dau1LVec = w2dau1.momentum();
	  w2dau2LVec = w2dau2.momentum();
	  w1dau1Id = leptonsId[idPair1.first];
	  w2dau1Id = leptonsId[idPair2.first];
	  w1Q = w1dau1Id > 0 ? -1 : 1;
	  w2Q = w2dau1Id > 0 ? -1 : 1;

          switch (w1dau1Id) {
	  case 13: case -13: _mode1 = CH_MUON; break;
	  case 11: case -11: _mode1 = CH_ELECTRON; break;
	   }
	  switch (w2dau1Id) {
	  case 13: case -13: _mode2 = CH_MUON; break;
	  case 11: case -11: _mode2 = CH_ELECTRON; break;
	   }

          if (debug) std::cout<<" passed 4  "<<std::endl;

	  const auto w1LVec = w1dau1LVec+w1dau2LVec;
	  const auto w2LVec = w2dau1LVec+w2dau2LVec;

	  // Combine b jets
	  double sumDm = 1e9;
	  FourMomentum b1LVec, b2LVec;
	  for (size_t i = 0, n = _bjets.size(); i < n; ++i) {
	    const Jet& bjet1 = _bjets[i];
	    const double mtop1 = (w1LVec+bjet1.momentum()).mass();
	    const double dmtop1 = std::abs(mtop1-_tMass);
	    for (size_t j=0; j<n; ++j) {
	      if (i == j) continue;
	      const Jet& bjet2 = _bjets[j];
	      const double mtop2 = (w2LVec+bjet2.momentum()).mass();
	      const double dmtop2 = std::abs(mtop2-_tMass);

	      if (sumDm <= dmtop1+dmtop2) continue;

	      sumDm = dmtop1+dmtop2;
	      b1LVec = bjet1.momentum();
	      b2LVec = bjet2.momentum();
	    }
	  }

	  if (sumDm >= 1e9) return; // Failed to make top, but this should not happen.

	  const auto t1LVec = w1LVec + b1LVec;
	  const auto t2LVec = w2LVec + b2LVec;
	  const auto ttLVec = t1LVec + t2LVec;
	  const auto llLVec = w1dau1LVec + w2dau1LVec;
          const auto bbLVec = b1LVec + b2LVec;

	  //apply lepton-bjet cleaning
	  if (( deltaR(w1dau1LVec,  b1LVec) < 0.4 )   ||  ( deltaR(w1dau1LVec,  b2LVec) < 0.4 ) || ( deltaR(w2dau1LVec,  b2LVec) < 0.4 ))   vetoEvent;

          if (debug) std::cout<<" passed 4.1  "<<std::endl;

	  // Put all of them into candidate collection
	  _t1 = Particle(w1Q*6, t1LVec);
	  _b1 = Particle(w1Q*5, b1LVec);
	  _w1 = Particle(w1Q*24, w1LVec);
	  _wDecays1.push_back(Particle(w1dau1Id, w1dau1LVec));
	  _wDecays1.push_back(Particle(-w1dau1Id+w1Q, w1dau2LVec));

          if (debug) std::cout<<" passed 4.11  "<<std::endl;

	  _t2 = Particle(w2Q*6, t2LVec);
	  _b2 = Particle(w2Q*5, b2LVec);
	  _w2 = Particle(w2Q*24, w2LVec);
	  _wDecays2.push_back(Particle(w2dau1Id, w2dau1LVec));
	  _wDecays2.push_back(Particle(-w2dau1Id+w2Q, w2dau2LVec));

	  _isValid = true;

      if (debug) std::cout<<" passed 4.2  "<<std::endl;

      Particle _t1P(6, FourMomentum(_t1.momentum()));
      Particle _t2P(-6, FourMomentum(_t2.momentum()));

      Particle _b1P(5, FourMomentum(_b1.momentum()));
      Particle _b2P(-5, FourMomentum(_b2.momentum()));

      Particle _lepton1P(w1dau1Id, w1dau1LVec);
      Particle _lepton2P(w2dau1Id, w2dau1LVec);

	  _tops.clear();
	  _bottoms.clear();
      _pleptons.clear();

	  _tops.push_back(_t1P);
	  _tops.push_back(_t2P);

      _bottoms.push_back(_b1P);
      _bottoms.push_back(_b2P);

      _pleptons.push_back(_lepton1P);
      _pleptons.push_back(_lepton2P);

      const FourMomentum t1_ttrf = LorentzTransform::mkObjTransformFromBeta( -1. * ttLVec.boostVector() ).transform( _t1 );
      const FourMomentum t2_ttrf = LorentzTransform::mkObjTransformFromBeta( -1. * ttLVec.boostVector() ).transform( _t2 );

	  std::sort(_tops.begin(), _tops.end(), GreaterByPt());
	  std::sort(_bottoms.begin(), _bottoms.end(), GreaterByPt());
	  std::sort(_pleptons.begin(), _pleptons.end(), GreaterByPt());

	  double tt_deltaphi = deltaPhi(t1LVec, t2LVec);
	  double tt_deltarap = fabs(t1LVec.rapidity()) - fabs(t2LVec.rapidity());

	  double ll_deltaphi = fabs(deltaPhi(w1dau1LVec, w2dau1LVec));
	  double ll_deltaeta = fabs(w1dau1LVec.pseudorapidity()) - fabs(w2dau1LVec.pseudorapidity());
      double ll_absdeltaeta = fabs(w1dau1LVec.pseudorapidity() - w2dau1LVec.pseudorapidity());

    
	  //fill histos 
	  _h_t_pT->fill(t1LVec.pT(), weight);
      _h_tbar_pT->fill(t2LVec.pT(), weight);
	  _h_tlead_pT->fill(_tops[0].pT(), weight);
	  _h_ttrail_pT->fill(_tops[1].pT(), weight);
	  _h_tttrf_pT->fill(t1_ttrf.pT(), weight);
	  _h_t_y->fill(t1LVec.rapidity(), weight);
	  _h_tbar_y->fill(t2LVec.rapidity(), weight);
	  _h_tlead_y->fill(_tops[0].rapidity(), weight);
	  _h_ttrail_y->fill(_tops[1].rapidity(), weight);
	  _h_tt_pT->fill(ttLVec.pT(), weight);
	  _h_tt_y->fill(ttLVec.rapidity(), weight);
	  _h_tt_mass->fill(ttLVec.mass(), weight);
	  _h_tt_delphi->fill(tt_deltaphi, weight);
      _h_tt_delrap->fill(tt_deltarap, weight);
	  _h_l_pT->fill(w1dau1LVec.pT(), weight);
	  _h_lbar_pT->fill(w2dau1LVec.pT(), weight);
	  _h_lbar_eta->fill(w1dau2LVec.pseudorapidity(), weight);
      _h_l_eta->fill(w2dau1LVec.pseudorapidity(), weight);
      _h_llead_pT->fill(_pleptons[0].pT(), weight);
      _h_ltrail_pT->fill(_pleptons[1].pT(), weight);
      _h_llead_eta->fill(_pleptons[0].pseudorapidity(), weight);
      _h_ltrail_eta->fill(_pleptons[1].pseudorapidity(), weight);
	  _h_llbar_pT->fill(llLVec.pT(), weight);
	  _h_llbar_mass->fill(llLVec.mass(), weight);
      _h_llbar_delphi->fill(ll_deltaphi, weight);
      _h_llbar_deleta->fill(ll_deltaeta, weight);
	  _h_blead_pT->fill(_bottoms[0].pT(), weight);
	  _h_btrail_pT->fill(_bottoms[1].pT(), weight);
	  _h_blead_eta->fill(_bottoms[0].rapidity(), weight);
	  _h_btrail_eta->fill(_bottoms[1].rapidity(), weight);
	  _h_bbbar_mass->fill(bbLVec.mass(), weight);
	  _h_bbbar_pT->fill(bbLVec.pT(), weight);
	  _h_njet->fill(_jets.size(), weight);
      _h_llbar_absdeleta->fill(ll_absdeltaeta, weight);

    }
    /// Normalise histograms etc., after the run
    void finalize() {
      // normalize(_h_YYYY); // normalize to unity
        
      scale(_h_t_pT, crossSection()/picobarn/sumOfWeights()); // norm to cross section
      scale(_h_tbar_pT, crossSection()/picobarn/sumOfWeights());
      scale(_h_tlead_pT, crossSection()/picobarn/sumOfWeights());
      scale(_h_ttrail_pT, crossSection()/picobarn/sumOfWeights());
      scale(_h_tttrf_pT, crossSection()/picobarn/sumOfWeights()); // norm to cross section
      scale(_h_t_y, crossSection()/picobarn/sumOfWeights());
      scale(_h_tbar_y, crossSection()/picobarn/sumOfWeights());
      scale(_h_tlead_y, crossSection()/picobarn/sumOfWeights());
      scale(_h_ttrail_y, crossSection()/picobarn/sumOfWeights());
      scale(_h_tt_pT, crossSection()/picobarn/sumOfWeights());
      scale(_h_tt_y, crossSection()/picobarn/sumOfWeights());
      scale(_h_tt_mass, crossSection()/picobarn/sumOfWeights());
      scale(_h_tt_delphi, crossSection()/picobarn/sumOfWeights());
      scale(_h_tt_delrap, crossSection()/picobarn/sumOfWeights());
      scale(_h_l_pT, crossSection()/picobarn/sumOfWeights());
      scale(_h_l_eta, crossSection()/picobarn/sumOfWeights());
      scale(_h_llead_eta, crossSection()/picobarn/sumOfWeights());
      scale(_h_ltrail_eta, crossSection()/picobarn/sumOfWeights());
      scale(_h_llbar_pT, crossSection()/picobarn/sumOfWeights());
      scale(_h_llbar_mass, crossSection()/picobarn/sumOfWeights());
      scale(_h_llbar_delphi, crossSection()/picobarn/sumOfWeights());
      scale(_h_llbar_deleta, crossSection()/picobarn/sumOfWeights());
      scale(_h_blead_pT, crossSection()/picobarn/sumOfWeights());
      scale(_h_btrail_pT, crossSection()/picobarn/sumOfWeights());
      scale(_h_blead_eta, crossSection()/picobarn/sumOfWeights());
      scale(_h_btrail_eta, crossSection()/picobarn/sumOfWeights());
      scale(_h_bbbar_mass, crossSection()/picobarn/sumOfWeights());
      scale(_h_bbbar_pT, crossSection()/picobarn/sumOfWeights());
      scale(_h_njet, crossSection()/picobarn/sumOfWeights());
      scale(_h_llbar_absdeleta, crossSection()/picobarn/sumOfWeights());
    
        
     /*   scale(_h_t_pT, 1.0/sumOfWeights());
         scale(_h_tbar_pT, 1.0/sumOfWeights());
         scale(_h_tlead_pT, 1.0/sumOfWeights());
         scale(_h_ttrail_pT, 1.0/sumOfWeights());
         scale(_h_tttrf_pT, 1.0/sumOfWeights());
         scale(_h_t_y, 1.0/sumOfWeights());
         scale(_h_tbar_y, 1.0/sumOfWeights());
         scale(_h_tlead_y, 1.0/sumOfWeights());
         scale(_h_ttrail_y, 1.0/sumOfWeights());
         scale(_h_tt_pT, 1.0/sumOfWeights());
         scale(_h_tt_y, 1.0/sumOfWeights());
         scale(_h_tt_mass, 1.0/sumOfWeights());
         scale(_h_tt_delphi, 1.0/sumOfWeights());
         scale(_h_tt_delrap, 1.0/sumOfWeights());
         scale(_h_l_pT, 1.0/sumOfWeights());
         scale(_h_l_eta, 1.0/sumOfWeights());
         scale(_h_llead_eta, 1.0/sumOfWeights());
         scale(_h_ltrail_eta, 1.0/sumOfWeights());
         scale(_h_llbar_pT, 1.0/sumOfWeights());
         scale(_h_llbar_mass, 1.0/sumOfWeights());
         scale(_h_llbar_delphi, 1.0/sumOfWeights());
         scale(_h_llbar_deleta, 1.0/sumOfWeights());
         scale(_h_blead_pT, 1.0/sumOfWeights());
         scale(_h_btrail_pT, 1.0/sumOfWeights());
         scale(_h_blead_eta, 1.0/sumOfWeights());
         scale(_h_btrail_eta, 1.0/sumOfWeights());
         scale(_h_bbbar_mass, 1.0/sumOfWeights());
         scale(_h_bbbar_pT, 1.0/sumOfWeights());
         scale(_h_njet, 1.0/sumOfWeights());
        */

          std::cout<<" cross section =   "<<  crossSection()   <<std::endl;
          std::cout<<" sumOfWeights =   "<<  sumOfWeights()   <<std::endl;


    }

    void cleanup(std::map<double, std::pair<size_t, size_t> >& v, const bool doCrossCleanup) const
    {
      std::vector<std::map<double, std::pair<size_t, size_t> >::const_iterator> toErase;
      std::set<size_t> usedLeg1, usedLeg2;
      if ( !doCrossCleanup ) {
	for (auto key = v.begin(); key != v.end(); ++key) {
	  const size_t leg1 = key->second.first;
	  const size_t leg2 = key->second.second;
	  if (usedLeg1.find(leg1) == usedLeg1.end() and
	      usedLeg2.find(leg2) == usedLeg2.end()) {
	    usedLeg1.insert(leg1);
	    usedLeg2.insert(leg2);
	  } else {
	    toErase.push_back(key);
	  }
	}
      }
      else {
	for (auto key = v.begin(); key != v.end(); ++key) {
	  const size_t leg1 = key->second.first;
	  const size_t leg2 = key->second.second;
	  if (usedLeg1.find(leg1) == usedLeg1.end() and
	      usedLeg1.find(leg2) == usedLeg1.end()) {
	    usedLeg1.insert(leg1);
	    usedLeg1.insert(leg2);
	  } else {
	    toErase.push_back(key);
	  }
	}
      }
      for (auto& key : toErase) v.erase(key);
    }


   struct GreaterByPt
   {
     bool operator()(const Particle& a, const Particle& b) {
       return a.pt() > b.pt();
     }
   };


     enum TTbarMode {CH_NONE=-1, CH_FULLHADRON = 0, CH_SEMILEPTON, CH_FULLLEPTON};
     enum DecayMode {CH_HADRON = 0, CH_MUON, CH_ELECTRON};

     TTbarMode mode() const {
        if (!_isValid) return CH_NONE;
        if (_mode1 == CH_HADRON && _mode2 == CH_HADRON) return CH_FULLHADRON;
        else if ( _mode1 != CH_HADRON && _mode2 != CH_HADRON) return CH_FULLLEPTON;
        else return CH_SEMILEPTON;
      }

     DecayMode mode1() const {return _mode1;}
     DecayMode mode2() const {return _mode2;}

    double _lepR = 0.1, _lepMinPt = 20, _lepMinEta = 2.0, _lepMaxEta = 5.0;
    double _jetR = 0.4 , _jetMinPt = 20 , _jetMinEta = 2.2 ,_jetMaxEta = 4.2;

    constexpr static double _tMass = 172.5;
    constexpr static double _wMass = 80.4;


    Profile1DPtr _h_XXXX;
    Histo1DPtr _h_YYYY;
    CounterPtr _h_ZZZZ;

    Histo1DPtr _h_t_pT,  _h_tbar_pT,  _h_tlead_pT, _h_ttrail_pT,  _h_tttrf_pT,  _h_t_y,  _h_tbar_y,  _h_tlead_y,  _h_ttrail_y,   _h_tt_pT,   _h_tt_y,   _h_tt_mass, _h_tt_delphi, _h_tt_delrap;
  
    Histo1DPtr _h_l_pT, _h_lbar_pT, _h_l_eta,  _h_lbar_eta,  _h_llead_pT,  _h_ltrail_pT,  _h_llead_eta,  _h_ltrail_eta,  _h_llbar_pT,  _h_llbar_mass,   _h_llbar_delphi,  _h_llbar_deleta,  _h_blead_pT,   _h_btrail_pT,  _h_blead_eta,   _h_btrail_eta,  _h_bbbar_mass,   _h_bbbar_pT,  _h_njet, _h_llbar_absdeleta ;
  
    bool _isValid;

    Particle _t1, _t2;
    Particle _b1, _b2;
    Particle _w1, _w2;
    ParticleVector _wDecays1, _wDecays2, _tops, _bottoms, _pleptons;
    Jets _jets, _bjets, _ljets, _leptons;

    DecayMode _mode1, _mode2;

  };

  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(LHCb_TOP);

}

