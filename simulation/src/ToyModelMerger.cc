// ROOT includes
#include <TFile.h>
#include <TTree.h>
#include <TMath.h>
#include <TF1.h>
#include <TH1D.h>
#include <TRandom3.h>
#include <TTimeStamp.h>
#include <TKey.h>

// FastJet includes
#include <fastjet/ClusterSequence.hh>
#include <fastjet/PseudoJet.hh>
#include <fastjet/ClusterSequenceArea.hh>
#include <fastjet/config.h>

// Custom includes
#include "Settings.h"
#include "ToyModelMerger.h"

// C++ includes
#include <iostream>
#include <vector>

namespace toymodel  
{


void ToyModelMerger::Run()
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

  
    if(toymodel_settings.Main.Verbosity > 0)
    {
        toymodel_settings.print(std::cout);
    }


    // read input signal file
    TFile *signalFile = new TFile(m_input_signal_filename.c_str(), "READ");
    TTree *signalTree = (TTree*)signalFile->Get("tree");
    TTree *signalConfigTree = (TTree*)signalFile->Get("config");
    int n_generated_events = 0;
    int n_accepted_events = 0;
    double sum_of_weights = 0.0;
    double integrated_luminosity = 0.0;
    float pT_hard_max = 0.0;
    float pT_hard_min = 0.0;
    float min_pythia_jet_pt = 0.0;
    signalConfigTree->SetBranchAddress("n_generated_events", &n_generated_events);
    signalConfigTree->SetBranchAddress("n_accepted_events", &n_accepted_events);
    signalConfigTree->SetBranchAddress("sum_of_weights", &sum_of_weights);
    signalConfigTree->SetBranchAddress("integrated_luminosity", &integrated_luminosity);
    signalConfigTree->SetBranchAddress("pT_hard_max", &pT_hard_max);
    signalConfigTree->SetBranchAddress("pT_hard_min", &pT_hard_min);
    signalConfigTree->SetBranchAddress("min_pythia_jet_pt", &min_pythia_jet_pt);
    signalConfigTree->GetEntry(0);

    // input signal variables
    int event_id_sig = 0;
    double weight = 0.0;
    double jet_pt_truth = 0.0;
    double jet_truth_phi_orig = 0.0;
    double jet_eta_truth= 0.0;
    int jet_nconst_truth = 0;
    std::vector<double> * const_pt_sig = 0;
    std::vector<double> * const_eta_sig = 0;
    std::vector<double> * const_delta_phi_sig = 0;
    std::vector<double> * const_pz_sig = 0;
    std::vector<double> * const_e_sig = 0;
    signalTree->SetBranchAddress("event_id", &event_id_sig);
    signalTree->SetBranchAddress("weight", &weight);
    signalTree->SetBranchAddress("jet_pt", &jet_pt_truth);
    signalTree->SetBranchAddress("jet_phi", &jet_truth_phi_orig);
    signalTree->SetBranchAddress("jet_eta", &jet_eta_truth);
    signalTree->SetBranchAddress("jet_nconst", &jet_nconst_truth);
    signalTree->SetBranchAddress("const_pt", &const_pt_sig);
    signalTree->SetBranchAddress("const_eta", &const_eta_sig);
    signalTree->SetBranchAddress("const_delta_phi", &const_delta_phi_sig);
    signalTree->SetBranchAddress("const_pz", &const_pz_sig);
    signalTree->SetBranchAddress("const_e", &const_e_sig);
    int n_signal_events = signalTree->GetEntries();

    // read input background file
    TFile *bkgdFile = new TFile(m_input_background_filename.c_str(), "READ");
    TTree *bkgdTree = (TTree*)bkgdFile->Get("tree");
    TTree *bkgdConfigTree = (TTree*)bkgdFile->Get("config");
    double v2_ref_truth = 0.0;
    double v3_ref_truth = 0.0;
    double v4_ref_truth = 0.0;
    bkgdConfigTree->SetBranchAddress("v2_ref_truth", &v2_ref_truth);
    bkgdConfigTree->SetBranchAddress("v3_ref_truth", &v3_ref_truth);
    bkgdConfigTree->SetBranchAddress("v4_ref_truth", &v4_ref_truth);
    bkgdConfigTree->GetEntry(0);

    // input background variables
    int event_id_bkgd = 0;
    double psi1 = 0.0;
    double psi2 = 0.0;
    double psi3 = 0.0;
    double psi4 = 0.0;
    double average_pt = 0.0;
    int n_forward_particles = 0;
    double Q2_i = 0.0;
    double Q2_r = 0.0;
    double Q3_i = 0.0;
    double Q3_r = 0.0;
    double Q4_i = 0.0;
    double Q4_r = 0.0;
    double Q6_i = 0.0;
    double Q6_r = 0.0;
    double Q8_i = 0.0;
    double Q8_r = 0.0;
    int n_midrapidity_particles = 0;
    std::vector<double> * bkgd_part_px = 0;
    std::vector<double> * bkgd_part_py = 0;
    std::vector<double> * bkgd_part_pz = 0;
    std::vector<double> * bkgd_part_e = 0;
    bkgdTree->SetBranchAddress("event_id", &event_id_bkgd);
    bkgdTree->SetBranchAddress("psi1", &psi1);
    bkgdTree->SetBranchAddress("psi2", &psi2);
    bkgdTree->SetBranchAddress("psi3", &psi3);
    bkgdTree->SetBranchAddress("psi4", &psi4);
    bkgdTree->SetBranchAddress("average_pt", &average_pt);
    bkgdTree->SetBranchAddress("n_forward_particles", &n_forward_particles);
    bkgdTree->SetBranchAddress("Q2_i", &Q2_i);
    bkgdTree->SetBranchAddress("Q2_r", &Q2_r);
    bkgdTree->SetBranchAddress("Q3_i", &Q3_i);
    bkgdTree->SetBranchAddress("Q3_r", &Q3_r);
    bkgdTree->SetBranchAddress("Q4_i", &Q4_i);
    bkgdTree->SetBranchAddress("Q4_r", &Q4_r);
    bkgdTree->SetBranchAddress("Q6_i", &Q6_i);
    bkgdTree->SetBranchAddress("Q6_r", &Q6_r);
    bkgdTree->SetBranchAddress("Q8_i", &Q8_i);
    bkgdTree->SetBranchAddress("Q8_r", &Q8_r);
    bkgdTree->SetBranchAddress("n_midrapidity_particles", &n_midrapidity_particles);
    bkgdTree->SetBranchAddress("part_px", &bkgd_part_px);
    bkgdTree->SetBranchAddress("part_py", &bkgd_part_py);
    bkgdTree->SetBranchAddress("part_pz", &bkgd_part_pz);
    bkgdTree->SetBranchAddress("part_e", &bkgd_part_e);
    int n_bkgd_events = bkgdTree->GetEntries();

    //============================================================
    //==================== OUTPUT CONFIGURATION ===================
    //============================================================
    // create output file
    TFile *outFile = new TFile(m_output_filename.c_str(),"RECREATE");
    // create output trees
    TTree * config_tree = new TTree("config", "config");
    TTree * event_tree = new TTree("tree", "tree");
   
    // config branches
    int max_nconst_truth = 0;
    double max_pt_reco = 0;
    double percentage_matched = 0;
    int nEvents = std::min(n_signal_events, n_bkgd_events);
    config_tree->Branch("nEvents", &nEvents, "nEvents/I");
    config_tree->Branch("seed", &toymodel_settings.Main.RandomSeed , "seed/i");
    config_tree->Branch("max_nconst_truth", &max_nconst_truth, "max_nconst_truth/I");
    config_tree->Branch("max_pt_reco", &max_pt_reco, "max_pt_reco/D");
    config_tree->Branch("percentage_matched", &percentage_matched, "percentage_matched/D");
    config_tree->Branch("max_particle_eta", & toymodel_settings.Main.Particle.MaxEta, "max_particle_eta/F");
    config_tree->Branch("min_particle_pt", & toymodel_settings.Main.Particle.MinPt, "min_particle_pt/F");
    config_tree->Branch("R", &toymodel_settings.Main.Jet.R, "R/F");
    config_tree->Branch("max_jet_eta", &toymodel_settings.Main.Jet.MaxEta, "max_jet_eta/F");
    config_tree->Branch("min_jet_pt", &toymodel_settings.Main.Jet.MinPt, "min_jet_pt/F");
    config_tree->Branch("min_pythia_jet_pt", &min_pythia_jet_pt, "min_pythia_jet_pt/F");
    config_tree->Branch("min_jet_area", &toymodel_settings.Main.Jet.MinArea, "min_jet_area/F");
    config_tree->Branch("dRMax", &dRMax, "dRMax/F");
    config_tree->Branch("pT_hard_min", &pT_hard_min, "pT_hard_min/F");
    config_tree->Branch("pT_hard_max", &pT_hard_max, "pT_hard_max/F");
    config_tree->Branch("n_accepted_events", &n_accepted_events, "n_accepted_events/I");
    config_tree->Branch("n_generated_events", &n_generated_events, "n_generated_events/I");
    config_tree->Branch("sum_of_weights", &sum_of_weights, "sum_of_weights/D");
    config_tree->Branch("integrated_luminosity", &integrated_luminosity, "integrated_luminosity/D");
    config_tree->Branch("v2_ref_truth", &v2_ref_truth, "v2_ref_truth/D");
    config_tree->Branch("v3_ref_truth", &v3_ref_truth, "v3_ref_truth/D");
    config_tree->Branch("v4_ref_truth", &v4_ref_truth, "v4_ref_truth/D");
   
    // output variables
    int event_id;
    double weight_reweight = 1/ integrated_luminosity;

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
    event_tree->Branch("weight", &weight_reweight, "weight/D");

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

    TRandom3 *phi_rand = new TRandom3(static_cast<unsigned int>(toymodel_settings.Main.RandomSeed*7));

   
    // ================================================================
    // ===================== FASTJET CONFIG  =============================
    // ================================================================
    
    // area definition
    fastjet::GhostedAreaSpec ghost_area_spec( toymodel_settings.Main.Particle.MaxEta, 1, 0.01);
    fastjet::AreaDefinition area_def(fastjet::active_area, ghost_area_spec);
    // jet definition
    fastjet::JetDefinition jet_def_antikt(fastjet::antikt_algorithm, toymodel_settings.Main.Jet.R, fastjet::E_scheme, fastjet::Best);
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
    TH1D * h1_pythia_jet_const_deltaR_postrotation = new TH1D("h1_jet_const_deltaR_postrotation", "h1_jet_const_deltaR_postrotation", 100, 0, 1);
    TH1D * h1_pythia_jet_const_delta_phi_postrotation = new TH1D("h1_jet_const_delta_phi_postrotation", "h1_jet_const_delta_phi_postrotation", 100, -TMath::Pi(), TMath::Pi());
    TH1D * h1_reco_jet_pt = new TH1D("h1_reco_jet_pt", "h1_reco_jet_pt", 100, 0, 100);
    TH1D * h1_reco_jet_eta = new TH1D("h1_reco_jet_eta", "h1_reco_jet_eta", 100, -5, 5);
    TH1D * h1_reco_jet_phi = new TH1D("h1_reco_jet_phi", "h1_reco_jet_phi", 100, -TMath::Pi(), TMath::Pi());
    TH1D * h1_pythia_jet_phi_postrotation = new TH1D("h1_pythia_jet_phi_postrotation", "h1_pythia_jet_phi_postrotation", 100, -TMath::Pi(), TMath::Pi());


    // ================================================================
    // ===================== EVENT LOOP ===============================
    // ================================================================
    if(toymodel_settings.Main.Verbosity > 0) std::cout << "Starting event loop" << std::endl;
    int progess_update = int(nEvents/10);
    int progress = 0;

    // loop over events
    int iEvent = 0;
    int n_merged_jets = 0;
    while(iEvent < nEvents)
    {
   
        // reset event variables
        event_id = iEvent;
        
        // get signal event
        signalTree->GetEntry(iEvent);
        // get background event
        bkgdTree->GetEntry(iEvent);


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
     
        if(jet_nconst_truth > max_nconst_truth) max_nconst_truth = jet_nconst_truth;

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
                double x = phi_rand->Uniform(0.0, 2*TMath::Pi());
                double y2 = 1.0
                    +2.0*(jet_v2_truth-0.02)*TMath::Cos(x-psi1)
                    +2.0*jet_v2_truth*TMath::Cos(2.0*(x-psi2))
                    +2.0*jet_v3_truth*TMath::Cos(3.0*(x-psi3))
                    +2.0*jet_v4_truth*TMath::Cos(4.0*(x-psi4));
                
                if(y1 < y2){ phi_tmp = x; break; }
            }

            // normalize phi to -pi to pi
            if(phi_tmp > TMath::Pi()) phi_tmp -= 2*TMath::Pi();
            if(phi_tmp < -TMath::Pi()) phi_tmp += 2*TMath::Pi();

            jet_phi_truth = phi_tmp;
            h1_pythia_jet_phi_postrotation->Fill(jet_truth_phi_orig);

            // rotate jet constituents
            // double dphi_rot = jet_phi_truth - jet_truth_phi_orig;
            // if(dphi_rot > TMath::Pi()) dphi_rot -= 2*TMath::Pi();
            // if(dphi_rot < -TMath::Pi()) dphi_rot += 2*TMath::Pi();            
            for (unsigned int iconst = 0; iconst < jet_nconst_truth; iconst++)
            {

                // get constituent info
                double delta_phi_old = const_delta_phi_sig->at(iconst);
                double eta_const = const_eta_sig->at(iconst);
                double pt_const = const_pt_sig->at(iconst);
                double pz_const = const_pz_sig->at(iconst);
                double e_const =  const_e_sig->at(iconst);

                // rotate constituent momentum
                // double new_phi = delta_phi_old + dphi_rot;
                double new_phi = delta_phi_old + jet_phi_truth;
                if(new_phi > TMath::Pi()) new_phi -= 2*TMath::Pi();
                if(new_phi < -TMath::Pi()) new_phi += 2*TMath::Pi();
                double new_px = pt_const*TMath::Cos(new_phi);
                double new_py = pt_const*TMath::Sin(new_phi);
                fastjet::PseudoJet rotated_const = fastjet::PseudoJet(new_px, new_py, pz_const, e_const);

                double delta_phi = rotated_const.phi_std() - jet_phi_truth;
                if(delta_phi > TMath::Pi()) delta_phi -= 2*TMath::Pi();
                if(delta_phi < -TMath::Pi()) delta_phi += 2*TMath::Pi();
                double delta_eta = rotated_const.eta() - jet_eta_truth;

                double dR_new = TMath::Sqrt(delta_phi*delta_phi + delta_eta*delta_eta);

                h1_pythia_jet_const_deltaR_postrotation->Fill(dR_new);
                h1_pythia_jet_const_delta_phi_postrotation->Fill(delta_phi);


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
            for (unsigned int iconst = 0; iconst < jet_nconst_truth; iconst++)
            {
                double phi_i = const_delta_phi_sig->at(iconst) + jet_phi_truth;
                if(phi_i > TMath::Pi()) phi_i -= 2*TMath::Pi();
                if(phi_i < -TMath::Pi()) phi_i += 2*TMath::Pi();
                double px = const_pt_sig->at(iconst)*TMath::Cos(phi_i);
                double py = const_pt_sig->at(iconst)*TMath::Sin(phi_i);
                double pz = const_pz_sig->at(iconst);
                double e = const_e_sig->at(iconst);
                particles_all.push_back(fastjet::PseudoJet(px, py, pz, e));
            }
        }

        // ================================================================
        // ===================== BACKGROUND EVENT =========================
        // ================================================================
       
        // background particles
        for (unsigned int ipart = 0; ipart < n_midrapidity_particles; ipart++)
        {
            double px = bkgd_part_px->at(ipart);
            double py = bkgd_part_py->at(ipart);
            double pz = bkgd_part_pz->at(ipart);
            double e = bkgd_part_e->at(ipart);
            particles_all.push_back(fastjet::PseudoJet(px, py, pz, e));
        }



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
            h1_reco_jet_pt->Fill(jets.at(jjet).pt());
            h1_reco_jet_eta->Fill(jets.at(jjet).eta());
            h1_reco_jet_phi->Fill(jets.at(jjet).phi_std());

            double deta_jet = jet_eta_truth - jets.at(jjet).eta();
            double dphi_jet = jet_phi_truth - jets.at(jjet).phi_std();
            if (dphi_jet > TMath::Pi()) dphi_jet -= 2.0*TMath::Pi();
            if (dphi_jet < -TMath::Pi()) dphi_jet += 2.0*TMath::Pi();
            
            double dR = TMath::Sqrt(deta_jet*deta_jet + dphi_jet*dphi_jet);
            if (dR <= dRMin)
            {
                dRMin = dR;
                jet_pt_reco = jets.at(jjet).pt();
                jet_phi_reco = jets.at(jjet).phi_std();
                jet_nconst_reco = nconst;
                if(jet_pt_reco > max_pt_reco) max_pt_reco = jet_pt_reco;

            }
            else 
            {

                unmatched_reco_pt.push_back(jets.at(jjet).pt());
                unmatched_reco_phi.push_back(jets.at(jjet).phi_std());
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
    percentage_matched/=nEvents;
    std::cout << "Percentage of matched jets: " << percentage_matched << std::endl;
    // fill event info tree
    config_tree->Fill();
    
    //============================================================
    //==================== WRITE OUTPUT ==========================
    //============================================================

    // write output file
    outFile->cd();
    config_tree->Write();
    event_tree->Write();
    h1_pythia_jet_const_deltaR_postrotation->Write();
    h1_pythia_jet_const_delta_phi_postrotation->Write();
    h1_reco_jet_pt->Write();
    h1_reco_jet_eta->Write();
    h1_reco_jet_phi->Write();
    h1_pythia_jet_phi_postrotation->Write();

    // write jet vn functions as histograms
    std::vector<TF1 *> jet_vn_funcs = {jet_v2_func, jet_v3_func, jet_v4_func};
    for (auto f : jet_vn_funcs)
    {
        TH1D * h = (TH1D*)f->GetHistogram();
        h->SetName(Form("%s_hist", f->GetName()));
        outFile->cd();
        h->Write();
    }

    // get key values from input files
    TList * list = signalFile->GetListOfKeys();
    TIter next(list);
    TKey * key;
    while(key = (TKey*)next())
    {
        if (key->GetClassName() == TString("TH1D"))
        {
            TH1D * h1 = (TH1D*)key->ReadObj();
            outFile->cd();
            h1->Write();
        }
    }

    list = bkgdFile->GetListOfKeys();
    next = TIter(list);
    while(key = (TKey*)next())
    {
        if (key->GetClassName() == TString("TH1D"))
        {
            TH1D * h1 = (TH1D*)key->ReadObj();
            outFile->cd();
            h1->Write();
        }
    }
    
    outFile->Close();
    signalFile->Close();
    bkgdFile->Close();

    std::cout << "Done" << std::endl;
    return ;

}

}