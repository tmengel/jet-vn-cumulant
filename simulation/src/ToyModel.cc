// ROOT includes
#include <TFile.h>
#include <TTree.h>
#include <TMath.h>
#include <TF1.h>
#include <TH1D.h>
#include <TRandom3.h>
#include <TTimeStamp.h>
 
// Pythia includes
#include <Pythia8/Pythia.h>

// FastJet includes
#include <fastjet/ClusterSequence.hh>
#include <fastjet/PseudoJet.hh>
#include <fastjet/ClusterSequenceArea.hh>
#include <fastjet/config.h>

// Custom includes
#include "Settings.h"
#include "BkgdFunctions.h"
#include "ToyModel.h"

// C++ includes
#include <iostream>
#include <vector>

namespace toymodel  
{

bool ToyModel::accept_pythia_event(std::vector<fastjet::PseudoJet> &jets, double pTmin, double pTmax)
{

        if(pTmin  < 0 ) return true;
        if(pTmax < 0) pTmax = 1000000; // if pTmax is negative, set to large number

        fastjet::JetDefinition jetdef(fastjet::antikt_algorithm, 0.4);
        fastjet::ClusterSequence clust_seq(jets, jetdef);
        
        std::vector<fastjet::PseudoJet> inclusive_jets = fastjet::sorted_by_pt(clust_seq.inclusive_jets(0));
        
        double leading_jet_pt = 0;
        for (auto jet : jets)
        {
            if (jet.pt() > leading_jet_pt) leading_jet_pt = jet.pt();
        }
        
        if(leading_jet_pt < pTmin || leading_jet_pt > pTmax)
        {
            return false;
        }


        return true;
}

void ToyModel::run()
{

    // ================================================================
    // ===================== CONFIGURATION ============================
    // ================================================================
    ToyModelSettings toymodel_settings;
    toymodel_settings.read(m_config_filename);

    if(toymodel_settings.Main.RandomSeed == 0 )
    {
        TTimeStamp *time = new TTimeStamp();
        toymodel_settings.Main.RandomSeed = static_cast<unsigned int>(time->GetNanoSec());
        delete time;
    }
    if (toymodel_settings.Main.RandomSeed > 900000000) toymodel_settings.Main.RandomSeed = toymodel_settings.Main.RandomSeed % 900000000;

    
    if(toymodel_settings.Main.Jet.MaxEta  < 0) toymodel_settings.Main.Jet.MaxEta  = toymodel_settings.Main.Particle.MaxEta - toymodel_settings.Main.Jet.R;
    if (toymodel_settings.Signal.Jet.MinPt < 0) toymodel_settings.Signal.Jet.MinPt = toymodel_settings.Main.Jet.MinPt;
    float dRMax = 0.75*toymodel_settings.Main.Jet.R;

    std::string pythia_config_file = toymodel_settings.Signal.Pythia.Config;
    std::vector<std::string> pythia_commands = toymodel_settings.Signal.Pythia.Commands;

    if(toymodel_settings.Main.Verbosity > 0)
    {
        toymodel_settings.print(std::cout);
    }


    double xsec_over_eventweight = 0.0;
    int n_accepted_events = 0;
    int n_generated_events = 0;
    int n_processed_events = 0;
    double sum_of_weights = 0.0;
    double integrated_luminosity = 0.0;

    double v2_ref_truth = toymodel::TruthRefVn(0, toymodel_settings.Main.Particle.MinPt);
    double v3_ref_truth = toymodel::TruthRefVn(1, toymodel_settings.Main.Particle.MinPt);
    double v4_ref_truth = toymodel::TruthRefVn(2, toymodel_settings.Main.Particle.MinPt);
    double v2_sum = 0.0;
    double v3_sum = 0.0;
    double v4_sum = 0.0;
    int n_tenngen_particles = 0;
    int max_nconst_truth = 0;
    double max_pt_reco = 0;
    double percentage_matched = 0;
    //============================================================
    //==================== OUTPUT CONFIGURATION ===================
    //============================================================
    // create output file
    TFile *outFile = new TFile(m_output_filename.c_str(),"RECREATE");
    
    // create output trees
    TTree * config_tree = new TTree("config", "config");
    TTree * event_tree = new TTree("tree", "tree");
   
    // config branches
    config_tree->Branch("nEvents", &toymodel_settings.Main.NumEvents, "nEvents/I");
    config_tree->Branch("seed", &toymodel_settings.Main.RandomSeed , "seed/i");
    config_tree->Branch("max_nconst_truth", &max_nconst_truth, "max_nconst_truth/I");
    config_tree->Branch("max_pt_reco", &max_pt_reco, "max_pt_reco/D");
    config_tree->Branch("percentage_matched", &percentage_matched, "percentage_matched/D");

    config_tree->Branch("max_particle_eta", & toymodel_settings.Main.Particle.MaxEta, "max_particle_eta/F");
    config_tree->Branch("min_particle_pt", & toymodel_settings.Main.Particle.MinPt, "min_particle_pt/F");
    
    config_tree->Branch("R", &toymodel_settings.Main.Jet.R, "R/F");
    config_tree->Branch("max_jet_eta", &toymodel_settings.Main.Jet.MaxEta, "max_jet_eta/F");
    config_tree->Branch("min_jet_pt", &toymodel_settings.Main.Jet.MinPt, "min_jet_pt/F");
    config_tree->Branch("min_pythia_jet_pt", &toymodel_settings.Signal.Jet.MinPt, "min_pythia_jet_pt/F");
    config_tree->Branch("min_jet_area", &toymodel_settings.Main.Jet.MinArea, "min_jet_area/F");
    config_tree->Branch("dRMax", &dRMax, "dRMax/F");

    config_tree->Branch("pT_hard_min", &toymodel_settings.Signal.Pythia.PtHardMin, "pT_hard_min/F");
    config_tree->Branch("pT_hard_max", &toymodel_settings.Signal.Pythia.PtHardMax, "pT_hard_max/F");

    config_tree->Branch("xsec_over_eventweight", &xsec_over_eventweight, "xsec_over_eventweight/D");
    config_tree->Branch("n_accepted_events", &n_accepted_events, "n_accepted_events/I");
    config_tree->Branch("n_generated_events", &n_generated_events, "n_generated_events/I");
    config_tree->Branch("n_processed_events", &n_processed_events, "n_processed_events/I");
    config_tree->Branch("sum_of_weights", &sum_of_weights, "sum_of_weights/D");
    config_tree->Branch("integrated_luminosity", &integrated_luminosity, "integrated_luminosity/D");

    config_tree->Branch("v2_ref_truth", &v2_ref_truth, "v2_ref_truth/D");
    config_tree->Branch("v3_ref_truth", &v3_ref_truth, "v3_ref_truth/D");
    config_tree->Branch("v4_ref_truth", &v4_ref_truth, "v4_ref_truth/D");
    config_tree->Branch("v2_sum", &v2_sum, "v2_sum/D");
    config_tree->Branch("v3_sum", &v3_sum, "v3_sum/D");
    config_tree->Branch("v4_sum", &v4_sum, "v4_sum/D");
    config_tree->Branch("n_tenngen_particles", &n_tenngen_particles, "n_tenngen_particles/I");

    // output variables (bkdg)
    int event_id;
    double weight;
    double psi1, psi2, psi3, psi4; // event plane angles    
    double average_pt;

    int n_forward_particles; // forward rapidity multiplicity
    double Q2_i, Q2_r; // q-vector real and imaginary parts
    double Q3_i, Q3_r; // q-vector real and imaginary parts
    double Q4_i, Q4_r; // q-vector real and imaginary parts
    double Q6_i, Q6_r; // q-vector real and imaginary parts
    double Q8_i, Q8_r; // q-vector real and imaginary parts
    
    // single truth jet variables
    double jet_pt_truth;
    int jet_nconst_truth;
    double jet_phi_truth;

    double jet_v2_truth;
    double jet_v3_truth;
    double jet_v4_truth;

    // single reco jet variables
    double jet_pt_reco;
    int jet_nconst_reco;
    double jet_phi_reco;

    // unmatched reco jet variables
    int n_unmatched_reco_jets;
    std::vector<double> unmatched_reco_pt;
    std::vector<int> unmatched_reco_nconst;
    std::vector<double> unmatched_reco_phi;

    // set up output branches
    event_tree->Branch("event_id", &event_id, "event_id/I");
    event_tree->Branch("weight", &weight, "weight/D");

    event_tree->Branch("psi1", &psi1, "psi1/D");
    event_tree->Branch("psi2", &psi2, "psi2/D");
    event_tree->Branch("psi3", &psi3, "psi3/D");
    event_tree->Branch("psi4", &psi4, "psi4/D");
    event_tree->Branch("average_pt", &average_pt, "average_pt/D");

    event_tree->Branch("n_forward_particles", &n_forward_particles, "n_forward_particles/I");
    event_tree->Branch("Q2_i", &Q2_i, "Q2_i/D");
    event_tree->Branch("Q2_r", &Q2_r, "Q2_r/D");
    event_tree->Branch("Q3_i", &Q3_i, "Q3_i/D");
    event_tree->Branch("Q3_r", &Q3_r, "Q3_r/D");
    event_tree->Branch("Q4_i", &Q4_i, "Q4_i/D");
    event_tree->Branch("Q4_r", &Q4_r, "Q4_r/D");
    event_tree->Branch("Q6_i", &Q6_i, "Q6_i/D");
    
    event_tree->Branch("Q6_r", &Q6_r, "Q6_r/D");
    event_tree->Branch("Q8_i", &Q8_i, "Q8_i/D");
    event_tree->Branch("Q8_r", &Q8_r, "Q8_r/D");
    
    event_tree->Branch("jet_pt_truth", &jet_pt_truth, "jet_pt_truth/D");
    event_tree->Branch("jet_phi_truth", &jet_phi_truth, "jet_phi_truth/D");
    event_tree->Branch("jet_nconst_truth", &jet_nconst_truth, "jet_nconst_truth/I");
    event_tree->Branch("jet_v2_truth", &jet_v2_truth, "jet_v2_truth/D");
    event_tree->Branch("jet_v3_truth", &jet_v3_truth, "jet_v3_truth/D");
    event_tree->Branch("jet_v4_truth", &jet_v4_truth, "jet_v4_truth/D");

    event_tree->Branch("jet_pt_reco", &jet_pt_reco, "jet_pt_reco/D");
    event_tree->Branch("jet_nconst_reco", &jet_nconst_reco, "jet_nconst_reco/I");
    event_tree->Branch("jet_phi_reco", &jet_phi_reco, "jet_phi_reco/D");

    event_tree->Branch("n_unmatched_reco_jets", &n_unmatched_reco_jets, "n_unmatched_reco_jets/I");
    event_tree->Branch("unmatched_reco_pt", &unmatched_reco_pt);
    event_tree->Branch("unmatched_reco_nconst", &unmatched_reco_nconst);
    event_tree->Branch("unmatched_reco_phi", &unmatched_reco_phi);


    
    //================================================================
    //==================== BKGD CONFIGURATION =======================
    //================================================================
    // initialize pT distributions
    if(toymodel_settings.Main.Verbosity > 0) std::cout << "Configuring Bkgd distributions" << std::endl;
    
    TF1 * f1_piPluspT = new TF1("f1_piPluspT", toymodel::BGBW,  toymodel_settings.Main.Particle.MinPt, 100, 5);
    TF1 * f1_piMinuspT = new TF1("f1_piMinuspT", toymodel::BGBW,  toymodel_settings.Main.Particle.MinPt, 100, 5);
    TF1 * f1_KPluspT = new TF1("f1_KPluspT", toymodel::BGBW,  toymodel_settings.Main.Particle.MinPt, 100, 5);
    TF1 * f1_KMinuspT = new TF1("f1_KMinuspT", toymodel::BGBW,  toymodel_settings.Main.Particle.MinPt, 100, 5);
    TF1 * f1_protonpT = new TF1("f1_protonpT", toymodel::BGBW,  toymodel_settings.Main.Particle.MinPt, 100, 5);
    TF1 * f1_pbarpT = new TF1("f1_pbarpT", toymodel::BGBW,  toymodel_settings.Main.Particle.MinPt, 100, 5);

    // initialize vn functions
    std::cout << "Configuring vn functions" << std::endl;
    TF1 * v2_pi =  new TF1("v2_pi", toymodel::VnFunction,  toymodel_settings.Main.Particle.MinPt, 100, 7);
    TF1 * v3_pi =  new TF1("v3_pi", toymodel::VnFunction,  toymodel_settings.Main.Particle.MinPt, 100 , 7);
    TF1 * v4_pi =  new TF1("v4_pi", toymodel::VnFunction,  toymodel_settings.Main.Particle.MinPt, 100 , 7);
    TF1 * v2_k =  new TF1("v2_k", toymodel::VnFunction,  toymodel_settings.Main.Particle.MinPt, 100 , 7);
    TF1 * v3_k =  new TF1("v3_k", toymodel::VnFunction,  toymodel_settings.Main.Particle.MinPt, 100 , 7);
    TF1 * v4_k =  new TF1("v4_k", toymodel::VnFunction,  toymodel_settings.Main.Particle.MinPt, 100 , 7);
    TF1 * v2_pro =  new TF1("v2_pro", toymodel::VnFunction,  toymodel_settings.Main.Particle.MinPt, 100 , 7);
    TF1 * v3_pro =  new TF1("v3_pro", toymodel::VnFunction,  toymodel_settings.Main.Particle.MinPt, 100 , 7);
    TF1 * v4_pro =  new TF1("v4_pro", toymodel::VnFunction,  toymodel_settings.Main.Particle.MinPt, 100 , 7);
    
    std::vector<TF1*> pT_distos = { f1_piPluspT,
                                    f1_piMinuspT, 
                                    f1_KPluspT, 
                                    f1_KMinuspT, 
                                    f1_protonpT, 
                                    f1_pbarpT};
    for (auto f : pT_distos) f->SetNpx(1000);

    for (int i=0;i<6;i++)
    {
        for (int j=0;j<5;j++)
        {
            pT_distos.at(i)->SetParameter(j,toymodel::blastwave_params[j][i]);
        }
    }
  
    std::vector<TF1*> v2_distos = {v2_pi, 
                                   v2_k, 
                                   v2_pro};
    std::vector<TF1*> v3_distos = {v3_pi, 
                                   v3_k, 
                                   v3_pro};
    std::vector<TF1*> v4_distos = {v4_pi, 
                                   v4_k, 
                                   v4_pro};
    for (auto f : v2_distos) f->SetNpx(1000);
    for (auto f : v3_distos) f->SetNpx(1000);
    for (auto f : v4_distos) f->SetNpx(1000);
    for (int i=0;i<3;i++)
    {
        if(i == 0)
        {
            for (int j=0;j<7;j++)
            {
                v2_distos.at(i)->SetParameter(j,toymodel::pi_vn_params[0][j]);
                v3_distos.at(i)->SetParameter(j,toymodel::pi_vn_params[1][j]);
                v4_distos.at(i)->SetParameter(j,toymodel::pi_vn_params[2][j]);
            }
        }
        else if(i == 1)
        {
            for (int j=0;j<7;j++)
            {
                v2_distos.at(i)->SetParameter(j,toymodel::k_vn_params[0][j]);
                v3_distos.at(i)->SetParameter(j,toymodel::k_vn_params[1][j]);
                v4_distos.at(i)->SetParameter(j,toymodel::k_vn_params[2][j]);
            }
        }
        else if(i == 2)
        {
            for (int j=0;j<7;j++)
            {
                v2_distos.at(i)->SetParameter(j,toymodel::pro_vn_params[0][j]);
                v3_distos.at(i)->SetParameter(j,toymodel::pro_vn_params[1][j]);
                v4_distos.at(i)->SetParameter(j,toymodel::pro_vn_params[2][j]);
            }
        }
       
    }

    // initialize random number generators
    TRandom3 *eta_rand = new TRandom3(static_cast<unsigned int>(toymodel_settings.Main.RandomSeed));
    TRandom3 *phi_rand = new TRandom3(static_cast<unsigned int>(toymodel_settings.Main.RandomSeed*7));
    TRandom3 *pt_rand = new TRandom3(static_cast<unsigned int>(toymodel_settings.Main.RandomSeed*3));
    TRandom3 *psi_rand = new TRandom3(static_cast<unsigned int>(toymodel_settings.Main.RandomSeed*5));


    //============================================================
    //==================== PYTHIA CONFIGURATION ==================
    //============================================================
    if(toymodel_settings.Main.Verbosity > 0) std::cout << "Configuring Pythia" << std::endl;

    Pythia8::Pythia pythia;
    Pythia8::Settings& settings = pythia.settings;
    const Pythia8::Info& info = pythia.info;
    Pythia8::Event& event = pythia.event;
    if(pythia_config_file != "")
    {
        pythia.readFile(pythia_config_file);
    }
    for (auto command : pythia_commands)
    {
        pythia.readString(command);
    }
    if(toymodel_settings.Signal.Pythia.PtHardMin > 0)
    {
        settings.parm("PhaseSpace:pTHatMin", toymodel_settings.Signal.Pythia.PtHardMin);
        if(toymodel_settings.Signal.Pythia.PtHardMax > 0)
        {
            settings.parm("PhaseSpace:pTHatMax", toymodel_settings.Signal.Pythia.PtHardMax);
        }
    }
    settings.readString("Random:setSeed = on");
    settings.parm("Random:seed", toymodel_settings.Main.RandomSeed);

    pythia.init();
 
    // ================================================================
    // ===================== FASTJET CONFIG  =============================
    // ================================================================
    
    // area definition
    fastjet::GhostedAreaSpec ghost_area_spec( toymodel_settings.Main.Particle.MaxEta, 1, 0.01);
    fastjet::AreaDefinition area_def(fastjet::active_area, ghost_area_spec);
    
    // jet definition
    fastjet::JetDefinition jet_def_antikt(fastjet::antikt_algorithm, toymodel_settings.Main.Jet.R, fastjet::E_scheme, fastjet::Best);

    // pythia jet selector
    fastjet::Selector pythia_jet_select = fastjet::SelectorNHardest(1) * (
            fastjet::SelectorAbsEtaMax(toymodel_settings.Main.Jet.MaxEta) &&
            fastjet::SelectorPtMin(toymodel_settings.Signal.Jet.MinPt));

    // reco jet selector
    fastjet::Selector reco_jet_select = (fastjet::SelectorAbsEtaMax(toymodel_settings.Main.Jet.MaxEta) && fastjet::SelectorPtMin(toymodel_settings.Main.Jet.MinPt));

    // selectors for jet constituents
    fastjet::Selector comp_select = !fastjet::SelectorIsPureGhost();


    // ================================================================
    // ===================== Jet Rotation =============================
    // ================================================================

    TF1 * jet_v2_func = nullptr;
    TF1 * jet_v3_func = nullptr;
    TF1 * jet_v4_func = nullptr;
    if(toymodel_settings.Signal.Jet.Rotate)
    {
        jet_v2_func = new TF1("jet_v2_func", toymodel_settings.Signal.Jet.v2Function.c_str(), 0, 250);
        jet_v3_func = new TF1("jet_v3_func", toymodel_settings.Signal.Jet.v3Function.c_str(), 0, 250);
        jet_v4_func = new TF1("jet_v4_func", toymodel_settings.Signal.Jet.v4Function.c_str(), 0, 250);
    }

    // ================================================================
    // ===================== Debug Histograms =========================
    // ================================================================
    TH1D * h1_pre_rotated_jet_const_deltaR  = new TH1D("h1_pre_rotated_jet_const_deltaR", "h1_pre_rotated_jet_const_deltaR", 100, 0, 1);
    TH1D * h1_post_rotated_jet_const_deltaR = new TH1D("h1_post_rotated_jet_const_deltaR", "h1_post_rotated_jet_const_deltaR", 100, 0, 1);
    
    TH1D * h1_tenngen_dNdpT = new TH1D("h1_tenngen_dNdpT", "h1_tenngen_dNdpT", 100, 0, 100);
    TH1D * h1_pythia_dNdpT = new TH1D("h1_pythia_dNdpT", "h1_pythia_dNdpT", 100, 0, 100);
    TH1D * h1_tenngen_dNdeta = new TH1D("h1_tenngen_dNdeta", "h1_tenngen_dNdeta", 100, -5, 5);
    TH1D * h1_pythia_dNdeta = new TH1D("h1_pythia_dNdeta", "h1_pythia_dNdeta", 100, -5, 5);
    TH1D * h1_tenngen_dNdphi = new TH1D("h1_tenngen_dNdphi", "h1_tenngen_dNdphi", 100, 0, 2*TMath::Pi());
    TH1D * h1_pythia_dNdphi = new TH1D("h1_pythia_dNdphi", "h1_pythia_dNdphi", 100, 0, 2*TMath::Pi());

    TH1D * h1_mixed_event_dNjetdpT = new TH1D("h1_mixed_event_dNjetdpT", "h1_mixed_event_dNjetdpT", 100, 0, 100);

    TH1D * h1_truth_jet_const_dNdpT = new TH1D("h1_truth_jet_const_dNdpT", "h1_truth_jet_const_dNdpT", 100, 0, 100);
    TH1D * h1_mixed_jet_const_dNdpT = new TH1D("h1_mixed_jet_const_dNdpT", "h1_mixed_jet_const_dNdpT", 100, 0, 100);


    // ================================================================
    // ===================== EVENT LOOP ===============================
    // ================================================================
    if(toymodel_settings.Main.Verbosity > 0) std::cout << "Starting event loop" << std::endl;
    int progess_update = int(toymodel_settings.Main.NumEvents / 10);
    int progress = 0;

    // loop over events
    int iEvent = 0;
    int n_merged_jets = 0;
    while(iEvent < toymodel_settings.Main.NumEvents)
    {
   
        // reset event variables
        event_id = iEvent;
        weight = 0.0;

        // get psi
        if(toymodel_settings.Bkgd.ConstEventPlane1 < 0){ psi1 = psi_rand->Uniform(0,2*TMath::Pi()); }
        else { psi1 = toymodel_settings.Bkgd.ConstEventPlane1; }
        if(toymodel_settings.Bkgd.ConstEventPlane2 < 0){ psi2 = psi_rand->Uniform(0,2*TMath::Pi()); }
        else { psi2 = toymodel_settings.Bkgd.ConstEventPlane2; }
        if(toymodel_settings.Bkgd.ConstEventPlane3 < 0){ psi3 = psi_rand->Uniform(0,2*TMath::Pi()); }
        else { psi3 = toymodel_settings.Bkgd.ConstEventPlane3; }
        if(toymodel_settings.Bkgd.ConstEventPlane4 < 0){ psi4 = psi_rand->Uniform(0,2*TMath::Pi()); }
        else { psi4 = toymodel_settings.Bkgd.ConstEventPlane4; }


        // clear variables
        average_pt = 0.0;
        n_forward_particles = 0;
        Q2_i = 0.0;
        Q2_r = 0.0;
        Q3_i = 0.0;
        Q3_r = 0.0;
        Q4_i = 0.0;
        Q4_r = 0.0;
        Q6_i = 0.0;
        Q6_r = 0.0;
        Q8_i = 0.0;
        Q8_r = 0.0;

        jet_pt_truth = 0.0;
        jet_nconst_truth = 0;
        jet_phi_truth = 0.0;
        jet_v2_truth = 0.0;
        jet_v3_truth = 0.0;
        jet_v4_truth = 0.0;

        jet_pt_reco = -1.0;  // -1 means no matched jet
        jet_nconst_reco = -1; // -1 means no matched jet
        jet_phi_reco = -1.0; // -1 means no matched jet

        n_unmatched_reco_jets = 0;
        unmatched_reco_pt.clear();
        unmatched_reco_nconst.clear();
        unmatched_reco_phi.clear();

        // start list of all particles
        std::vector<fastjet::PseudoJet> particles_all;
        std::vector<fastjet::PseudoJet> pythia_particles;
        // ================================================================
        // ===================== PYTHIA EVENT =============================
        // ================================================================
        bool acceptEvent = false;
        while(!acceptEvent)
        {   
            
            pythia_particles.clear();

            while (!pythia.next()) 
            {
                n_generated_events++;
            }
            
            n_generated_events++;
            
            // get event info
            weight = info.weight();

            for (unsigned int ipart = 0; ipart < event.size(); ipart++){
                
                // only consider final state charged particles
                if(!event[ipart].isFinal()) continue;
                
                // only consider charged particles
                if(!event[ipart].isCharged()) continue;

                // min pT cut
                if(event[ipart].pT() <  toymodel_settings.Main.Particle.MinPt) continue;

                // abs(eta) cut
                if(TMath::Abs(event[ipart].eta()) >  toymodel_settings.Main.Particle.MaxEta) continue;


                // push back particle
                pythia_particles.push_back(fastjet::PseudoJet(event[ipart].px(), event[ipart].py(), event[ipart].pz(), event[ipart].e()));
            }

            
            acceptEvent = accept_pythia_event(pythia_particles, toymodel_settings.Signal.Pythia.PtHardMin, toymodel_settings.Signal.Pythia.PtHardMax);
        }


        // ================================================================
        // ===================== PYTHIA JET ===============================
        // ================================================================

        fastjet::ClusterSequence clust_seq_pythia(pythia_particles, jet_def_antikt); 
        std::vector<fastjet::PseudoJet> pythia_jets = sorted_by_pt(pythia_jet_select(clust_seq_pythia.inclusive_jets()));
        if(pythia_jets.size() == 0){ continue; }

        // get leading jet 
        fastjet::PseudoJet pythia_jet = pythia_jets.at(0);
        std::vector<fastjet::PseudoJet> pythia_jet_constituents = pythia_jet.constituents();
        jet_pt_truth = pythia_jet.pt();
        jet_nconst_truth = pythia_jet_constituents.size();
        for (unsigned int iconst = 0; iconst < pythia_jet_constituents.size(); iconst++)
        {
            h1_pythia_dNdpT->Fill(pythia_jet_constituents.at(iconst).pt());
            h1_pythia_dNdeta->Fill(pythia_jet_constituents.at(iconst).eta());
            h1_pythia_dNdphi->Fill(pythia_jet_constituents.at(iconst).phi());
            h1_truth_jet_const_dNdpT->Fill(pythia_jet_constituents.at(iconst).pt());
        }
        if(jet_nconst_truth > max_nconst_truth) max_nconst_truth = jet_nconst_truth;

        // needed for matching later
        double jet_eta_truth = pythia_jet.eta();
        double jet_truth_phi_orig = pythia_jet.phi(); 
      
        particles_all.clear();
        if(toymodel_settings.Signal.Jet.Rotate)
        {
            jet_v2_truth = jet_v2_func->Eval(jet_pt_truth);
            jet_v3_truth = jet_v3_func->Eval(jet_pt_truth);
            jet_v4_truth = jet_v4_func->Eval(jet_pt_truth);
        
            // get truth jet phi
            double tmp_dndpi = 1.0 + 2.0*( TMath::Abs(jet_v2_truth - 0.02) + TMath::Abs(jet_v2_truth) + TMath::Abs(jet_v3_truth) + TMath::Abs(jet_v4_truth));
            double phi_tmp = 0.0;
            while(true)
            {
                double y1 = phi_rand->Uniform(0.0,tmp_dndpi);
                double x = phi_rand->Uniform(0.0,2*TMath::Pi());
                double y2 = 1.0
                    +2.0*(jet_v2_truth-0.02)*TMath::Cos(x-psi1)
                    +2.0*jet_v2_truth*TMath::Cos(2.0*(x-psi2))
                    +2.0*jet_v3_truth*TMath::Cos(3.0*(x-psi3))
                    +2.0*jet_v4_truth*TMath::Cos(4.0*(x-psi4));
                
                if(y1 < y2){ phi_tmp = x; break; }
            }

            jet_phi_truth = phi_tmp;

            // rotate jet constituents
            double dphi_rot = jet_phi_truth - jet_truth_phi_orig;
            if(dphi_rot > TMath::Pi()) dphi_rot -= 2*TMath::Pi();
            if(dphi_rot < -TMath::Pi()) dphi_rot += 2*TMath::Pi();            
            for (unsigned int iconst = 0; iconst < pythia_jet_constituents.size(); iconst++)
            {

                // get constituent info
                double phi_org = pythia_jet_constituents.at(iconst).phi();
                double eta_const = pythia_jet_constituents.at(iconst).eta();
                double pt_const = pythia_jet_constituents.at(iconst).pt();
                double pz_const = pythia_jet_constituents.at(iconst).pz();
                double e_const = pythia_jet_constituents.at(iconst).e();

                double delta_phi_old = phi_org - jet_truth_phi_orig;
                if(delta_phi_old > TMath::Pi()) delta_phi_old -= 2*TMath::Pi();
                if(delta_phi_old < -TMath::Pi()) delta_phi_old += 2*TMath::Pi();
                double delta_eta_old = eta_const - jet_eta_truth;

                // rotate constituent momentum
                double new_phi = phi_org + dphi_rot;
                double new_px = pt_const*TMath::Cos(new_phi);
                double new_py = pt_const*TMath::Sin(new_phi);
                fastjet::PseudoJet rotated_const = fastjet::PseudoJet(new_px, new_py, pz_const, e_const);

                double delta_phi = rotated_const.phi() - jet_phi_truth;
                if(delta_phi > TMath::Pi()) delta_phi -= 2*TMath::Pi();
                if(delta_phi < -TMath::Pi()) delta_phi += 2*TMath::Pi();
                double delta_eta = rotated_const.eta() - jet_eta_truth;

                double dR_old = TMath::Sqrt(delta_phi_old*delta_phi_old + delta_eta_old*delta_eta_old);
                double dR_new = TMath::Sqrt(delta_phi*delta_phi + delta_eta*delta_eta);

                h1_pre_rotated_jet_const_deltaR->Fill(dR_old);
                h1_post_rotated_jet_const_deltaR->Fill(dR_new);


                // push back rotated constituent
                particles_all.push_back(rotated_const);

            }
        }
        else 
        {
            jet_phi_truth = jet_truth_phi_orig;
            jet_v2_truth = 0 ;
            jet_v3_truth = 0 ;
            jet_v4_truth = 0 ;
            for (unsigned int iconst = 0; iconst < pythia_jet_constituents.size(); iconst++)
            {
                particles_all.push_back(pythia_jet_constituents.at(iconst));
            }
        }

        // ================================================================
        // ===================== BACKGROUND EVENT =========================
        // ================================================================
       
        // loop over particles
        int M = 0;
        double rho = 0.0;
        for( unsigned int ispecies = 0; ispecies < toymodel::n_species; ispecies++)
        {
            double mass = toymodel::particle_masses.at(ispecies);
            int KF = toymodel::particle_ids.at(ispecies);
            TF1 *pT_distro = pT_distos.at(ispecies);
            unsigned int vn_idx = static_cast<unsigned int>(ispecies/2);
            for (unsigned int ipart = 0; ipart < toymodel::particle_yeilds.at(ispecies); ipart++)
            {
                double pT = pT_distro->GetRandom(pt_rand); // get pT
                double eta = eta_rand->Uniform(- toymodel_settings.Main.Particle.MaxEta, toymodel_settings.Main.Particle.MaxEta); // get eta

                // vns
                double v2 = v2_distos.at(vn_idx)->Eval(pT);
                double v3 = v3_distos.at(vn_idx)->Eval(pT);
                double v4 = v4_distos.at(vn_idx)->Eval(pT);
                double v1 = v2 - 0.02;

                

                
                double tmp_dndpi = 1.0 + 2.0*(TMath::Abs(v1) + TMath::Abs(v2) + TMath::Abs(v3) + TMath::Abs(v4));
                double phi = 0.0;
                while(true)
                {
                    double y1 = phi_rand->Uniform(0.0,tmp_dndpi);
                    double x = phi_rand->Uniform(0.0,2*TMath::Pi());
                    double y2 = 1.0+2.0*v1*TMath::Cos(x-psi1)+2.0*v2*TMath::Cos(2.0*(x-psi2))+2.0*v3*TMath::Cos(3.0*(x-psi3))+2.0*v4*TMath::Cos(4.0*(x-psi4));
                    
                    if(y1 < y2){ phi = x; break; }
                }

                // transform to cartesian
                double px = pT*TMath::Cos(phi);
                double py = pT*TMath::Sin(phi);
                double pz = pT*TMath::SinH(eta);
                double e = TMath::Sqrt(px*px + py*py + pz*pz + mass*mass);
                particles_all.push_back(fastjet::PseudoJet(px, py, pz, e));
               
                M++;
                rho+=pT;
                h1_tenngen_dNdpT->Fill(pT);
                h1_tenngen_dNdeta->Fill(eta);
                h1_tenngen_dNdphi->Fill(phi);

                // do forward rapidity particles
                pT = pT_distro->GetRandom(pt_rand); // get pT
                v2 = v2_distos.at(vn_idx)->Eval(pT);
                v3 = v3_distos.at(vn_idx)->Eval(pT);
                v4 = v4_distos.at(vn_idx)->Eval(pT);
                v1 = v2 - 0.02;

                v2_sum += v2;
                v3_sum += v3;
                v4_sum += v4;
                n_tenngen_particles++;

                tmp_dndpi = 1.0 + 2.0*(TMath::Abs(v1) + TMath::Abs(v2) + TMath::Abs(v3) + TMath::Abs(v4));
                phi = 0.0;
                while(true)
                {
                    double y1 = phi_rand->Uniform(0.0,tmp_dndpi);
                    double x = phi_rand->Uniform(0.0,2*TMath::Pi());
                    double y2 = 1.0+2.0*v1*TMath::Cos(x-psi1)+2.0*v2*TMath::Cos(2.0*(x-psi2))+2.0*v3*TMath::Cos(3.0*(x-psi3))+2.0*v4*TMath::Cos(4.0*(x-psi4));
                    
                    if(y1 < y2){ phi = x; break; }
                }

                Q2_i += TMath::Sin(2.0*phi);
                Q2_r += TMath::Cos(2.0*phi);
                Q3_i += TMath::Sin(3.0*phi);
                Q3_r += TMath::Cos(3.0*phi);
                Q4_i += TMath::Sin(4.0*phi);
                Q4_r += TMath::Cos(4.0*phi);
                Q6_i += TMath::Sin(6.0*phi);
                Q6_r += TMath::Cos(6.0*phi);
                Q8_i += TMath::Sin(8.0*phi);
                Q8_r += TMath::Cos(8.0*phi);

                n_forward_particles++;

            } // end of particle loop   
        
        } 

        // // get average pt density
        rho /= M;
        average_pt = rho;



        // ================================================================
        // ===================== RECO JETS ===============================
        // ================================================================

        // cluster particles
        fastjet::ClusterSequenceArea clust_seq_mixed(particles_all, jet_def_antikt, area_def);
        std::vector<fastjet::PseudoJet>  jets = sorted_by_pt(reco_jet_select(clust_seq_mixed.inclusive_jets()));
        
        //================================================================
        //==================== JET MATCHING ==============================
        //================================================================
  
        // match jets
        n_unmatched_reco_jets = 0;
        unmatched_reco_phi.clear();
        unmatched_reco_nconst.clear();
        unmatched_reco_pt.clear();

        double dRMin = dRMax;
        for (unsigned int jjet = 0; jjet < jets.size(); jjet++)
        {


            // mulitplicity subtraction
            h1_mixed_event_dNjetdpT->Fill(jets.at(jjet).pt());
            int nconst = comp_select(jets.at(jjet).constituents()).size();
            if (nconst == 0)
            {
                continue;
            }            
            if (jets.at(jjet).pt() < toymodel_settings.Main.Jet.MinPt) 
            {
                continue;
            }


            // cut reco jets
            if(TMath::Abs(jets.at(jjet).eta()) > toymodel_settings.Main.Jet.MaxEta)
            {
                continue;
            }
            if(jets.at(jjet).area() < toymodel_settings.Main.Jet.MinArea*TMath::Pi()*toymodel_settings.Main.Jet.R*toymodel_settings.Main.Jet.R)
            {
                continue;
            }
            std::vector<fastjet::PseudoJet> reco_jet_constituents = comp_select(jets.at(jjet).constituents());
            n_merged_jets++;
            for (unsigned int iconst = 0; iconst < reco_jet_constituents.size(); iconst++)
            {
                h1_mixed_jet_const_dNdpT->Fill(reco_jet_constituents.at(iconst).pt());
            }

            double deta_jet = jet_eta_truth - jets.at(jjet).eta();
            double dphi_jet = jet_phi_truth - jets.at(jjet).phi();
            if (dphi_jet > TMath::Pi()) dphi_jet -= 2.0*TMath::Pi();
            if (dphi_jet < -TMath::Pi()) dphi_jet += 2.0*TMath::Pi();
            
            double dR = TMath::Sqrt(deta_jet*deta_jet + dphi_jet*dphi_jet);
            if (dR <= dRMin)
            {
                dRMin = dR;
                jet_pt_reco = jets.at(jjet).pt();
                jet_phi_reco = jets.at(jjet).phi();
                jet_nconst_reco = nconst;
                if(jet_pt_reco > max_pt_reco) max_pt_reco = jet_pt_reco;

            }
            else 
            {

                unmatched_reco_pt.push_back(jets.at(jjet).pt());
                unmatched_reco_phi.push_back(jets.at(jjet).phi());
                unmatched_reco_nconst.push_back(nconst);
                n_unmatched_reco_jets++;
            
            }
        }
        if(jet_pt_reco > 0)
        {
            percentage_matched+=1;
        }

        // ================================================================
        // ===================== FILL TREE =================================
        // ================================================================

       // fill tree
        event_tree->Fill();


        iEvent++;
        // print progress
        if (iEvent % progess_update == 0)
        {
            std::cout << "Progress: " << progress << "%" << std::endl;
            progress += 10;
        }

    } // end of event loop
    std::cout << "Progress: 100%" << std::endl;
    std::cout << "Event loop finished" << std::endl;
    percentage_matched/=toymodel_settings.Main.NumEvents;
    std::cout << "Percentage of matched jets: " << percentage_matched << std::endl;
    n_accepted_events = toymodel_settings.Main.NumEvents;
    n_processed_events = info.nAccepted();
    sum_of_weights = info.weightSum();
    integrated_luminosity = (info.nAccepted()) /(info.sigmaGen()*1e9);
    xsec_over_eventweight = (info.sigmaGen() / info.weightSum());

    // fill event info tree
    config_tree->Fill();
    
    //============================================================
    //==================== WRITE OUTPUT ==========================
    //============================================================

    // write output file
    outFile->cd();
    config_tree->Write();
    event_tree->Write();
    h1_pre_rotated_jet_const_deltaR->Write();
    h1_post_rotated_jet_const_deltaR->Write();
    h1_tenngen_dNdpT->Write();
    h1_pythia_dNdpT->Write();
    h1_tenngen_dNdeta->Write();
    h1_pythia_dNdeta->Write();
    h1_tenngen_dNdphi->Write();
    h1_pythia_dNdphi->Write();
    h1_mixed_event_dNjetdpT->Write();
    h1_truth_jet_const_dNdpT->Write();
    h1_mixed_jet_const_dNdpT->Write();
    
    // Convert TF1 to TH1D
    for (auto f : pT_distos)
    {
        TH1D * h = (TH1D*)f->GetHistogram();
        h->SetName(Form("%s_hist", f->GetName()));
        h->Write();
    }
    for (auto f : v2_distos)
    {
        TH1D * h = (TH1D*)f->GetHistogram();
        h->SetName(Form("%s_hist", f->GetName()));
        h->Write();
    }
    for (auto f : v3_distos)
    {
        TH1D * h = (TH1D*)f->GetHistogram();
        h->SetName(Form("%s_hist", f->GetName()));
        h->Write();
    }
    for (auto f : v4_distos)
    {
        TH1D * h = (TH1D*)f->GetHistogram();
        h->SetName(Form("%s_hist", f->GetName()));
        h->Write();
    }

    std::vector<TF1 *> jet_vn_funcs = {jet_v2_func, jet_v3_func, jet_v4_func};
    for (auto f : jet_vn_funcs)
    {
        TH1D * h = (TH1D*)f->GetHistogram();
        h->SetName(Form("%s_hist", f->GetName()));
        h->Write();
    }
    
    outFile->Close();

    std::cout << "v2 ref man : " << v2_sum/n_tenngen_particles << std::endl;
    std::cout << "v3 ref man : " << v3_sum/n_tenngen_particles << std::endl;
    std::cout << "v4 ref man : " << v4_sum/n_tenngen_particles << std::endl;

    std::cout << "Done" << std::endl;
    return ;

}

}