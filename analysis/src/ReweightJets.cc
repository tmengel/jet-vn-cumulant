#include "ReweightJets.h"
#include "Utils.h"

#include <iostream>
#include <string>
#include <vector>
#include <cstdlib>

#include <TFile.h>
#include <TTree.h>
#include <TSystem.h>
#include <TKey.h>
#include <TH1D.h>


int ReweightJets::Reweight()
{
    // make output directory
    if (m_output_dir.empty())
    {
        std::cerr << "Output directory not set" << std::endl;
        exit(1);
    }

    if(gSystem->AccessPathName(m_output_dir.c_str()))
    {
        std::cout << "Creating output directory " << m_output_dir << std::endl;
        gSystem->mkdir(m_output_dir.c_str(), kTRUE);
        if(gSystem->AccessPathName(m_output_dir.c_str()))
        {
            std::cout << "Error: output directory " << m_output_dir << " does not exist" << std::endl;
            exit(1);
        }
    }

    m_input_files.clear();
    m_input_files = Utils::GetFilesFromDir(m_input_dir);
    std::cout << "Number of input files: " << m_input_files.size() << std::endl;

    if (m_input_files.size() == 0)
    {
        std::cerr << "No input files found in " << m_input_dir << std::endl;
        exit(1);
    }

    for (auto &input_file : m_input_files)
    {

        std::string output_file = m_output_dir + "/" + Utils::GetBaseName(input_file) + "_reweighted.root";

        TFile * fin = new TFile(input_file.c_str(), "READ");
        if(!fin->IsOpen())
        {
            std::cout << "Error: could not open input file " << input_file << std::endl;
            exit(1);
        }

        TTree * event_tree = (TTree*)fin->Get("tree");
        if(!event_tree)
        {
            std::cout << "Error: could not find tree in input file " << input_file << std::endl;
            exit(1);
        }

        TTree * config_tree = (TTree*)fin->Get("config");
        if(!config_tree)
        {
            std::cout << "Error: could not find config tree in input file " << input_file << std::endl;
            exit(1);
        }


        // input tree variables
        int event_id = 0;
        double weight = 0.0;    
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
        // single truth jet variables
        double jet_pt_truth = 0.0;
        int jet_nconst_truth = 0;
        double jet_phi_truth = 0.0;
        // truth jet vn variables
        double jet_v2_truth = 0.0;
        double jet_v3_truth = 0.0;
        double jet_v4_truth = 0.0;
        // single reco jet variables
        double jet_pt_reco = 0.0;
        int jet_nconst_reco = 0;
        double jet_phi_reco = 0.0;
        // unmatched reco jet variables
        int n_unmatched_reco_jets = 0;
        std::vector<double> * unmatched_reco_pt = 0;
        std::vector<int> * unmatched_reco_nconst = 0;
        std::vector<double> * unmatched_reco_phi  = 0;
        // set input tree branches
        event_tree->SetBranchAddress("event_id", &event_id);
        event_tree->SetBranchAddress("weight", &weight);
        event_tree->SetBranchAddress("psi1", &psi1);
        event_tree->SetBranchAddress("psi2", &psi2);
        event_tree->SetBranchAddress("psi3", &psi3);
        event_tree->SetBranchAddress("psi4", &psi4);
        event_tree->SetBranchAddress("average_pt", &average_pt);
        event_tree->SetBranchAddress("n_forward_particles", &n_forward_particles);
        event_tree->SetBranchAddress("Q2_i", &Q2_i);
        event_tree->SetBranchAddress("Q2_r", &Q2_r);
        event_tree->SetBranchAddress("Q3_i", &Q3_i);
        event_tree->SetBranchAddress("Q3_r", &Q3_r);
        event_tree->SetBranchAddress("Q4_i", &Q4_i);
        event_tree->SetBranchAddress("Q4_r", &Q4_r);
        event_tree->SetBranchAddress("Q6_i", &Q6_i);
        event_tree->SetBranchAddress("Q6_r", &Q6_r);
        event_tree->SetBranchAddress("Q8_i", &Q8_i);
        event_tree->SetBranchAddress("Q8_r", &Q8_r);
        event_tree->SetBranchAddress("jet_pt_truth", &jet_pt_truth);
        event_tree->SetBranchAddress("jet_nconst_truth", &jet_nconst_truth);
        event_tree->SetBranchAddress("jet_phi_truth", &jet_phi_truth);
        event_tree->SetBranchAddress("jet_v2_truth", &jet_v2_truth);
        event_tree->SetBranchAddress("jet_v3_truth", &jet_v3_truth);
        event_tree->SetBranchAddress("jet_v4_truth", &jet_v4_truth);
        event_tree->SetBranchAddress("jet_pt_reco", &jet_pt_reco);
        event_tree->SetBranchAddress("jet_nconst_reco", &jet_nconst_reco);
        event_tree->SetBranchAddress("jet_phi_reco", &jet_phi_reco);
        event_tree->SetBranchAddress("n_unmatched_reco_jets", &n_unmatched_reco_jets);
        event_tree->SetBranchAddress("unmatched_reco_pt", &unmatched_reco_pt);
        event_tree->SetBranchAddress("unmatched_reco_nconst", &unmatched_reco_nconst);
        event_tree->SetBranchAddress("unmatched_reco_phi", &unmatched_reco_phi);


        // input config variables
        int nEvents = 0;
        unsigned int seed = 0;
        int max_nconst_truth = 0;
        double max_pt_reco = 0.0;
        double percentage_matched = 0.0;
        float max_particle_eta = 0.0;
        float min_particle_pt = 0.0;
        float R = 0.0;
        float max_jet_eta = 0.0;
        float min_jet_pt = 0.0;
        float min_pythia_jet_pt = 0.0;
        float min_jet_area = 0.0;
        float dRMax = 0.0;
        float pT_hard_min = 0.0;
        float pT_hard_max = 0.0;
        double xsec_over_eventweight = 0.0;
        int n_accepted_events = 0;
        int n_generated_events = 0;
        int n_processed_events = 0;
        double sum_of_weights = 0.0;
        double integrated_luminosity = 0.0;
        double v2_ref_truth = 0.0;
        double v3_ref_truth = 0.0;
        double v4_ref_truth = 0.0;
        double v2_sum = 0.0;
        double v3_sum = 0.0;
        double v4_sum = 0.0;
        int n_tenngen_particles = 0;
        config_tree->SetBranchAddress("nEvents", &nEvents);
        config_tree->SetBranchAddress("seed", &seed);
        config_tree->SetBranchAddress("max_nconst_truth", &max_nconst_truth);
        config_tree->SetBranchAddress("max_pt_reco", &max_pt_reco);
        config_tree->SetBranchAddress("percentage_matched", &percentage_matched);
        config_tree->SetBranchAddress("max_particle_eta", &max_particle_eta);
        config_tree->SetBranchAddress("min_particle_pt", &min_particle_pt);
        config_tree->SetBranchAddress("R", &R);
        config_tree->SetBranchAddress("max_jet_eta", &max_jet_eta);
        config_tree->SetBranchAddress("min_jet_pt", &min_jet_pt);
        config_tree->SetBranchAddress("min_pythia_jet_pt", &min_pythia_jet_pt);
        config_tree->SetBranchAddress("min_jet_area", &min_jet_area);
        config_tree->SetBranchAddress("dRMax", &dRMax);
        config_tree->SetBranchAddress("pT_hard_min", &pT_hard_min);
        config_tree->SetBranchAddress("pT_hard_max", &pT_hard_max);
        config_tree->SetBranchAddress("xsec_over_eventweight", &xsec_over_eventweight);
        config_tree->SetBranchAddress("n_accepted_events", &n_accepted_events);
        config_tree->SetBranchAddress("n_generated_events", &n_generated_events);
        config_tree->SetBranchAddress("n_processed_events", &n_processed_events);
        config_tree->SetBranchAddress("sum_of_weights", &sum_of_weights);
        config_tree->SetBranchAddress("integrated_luminosity", &integrated_luminosity);
        config_tree->SetBranchAddress("v2_ref_truth", &v2_ref_truth);
        config_tree->SetBranchAddress("v3_ref_truth", &v3_ref_truth);
        config_tree->SetBranchAddress("v4_ref_truth", &v4_ref_truth);
        config_tree->SetBranchAddress("v2_sum", &v2_sum);
        config_tree->SetBranchAddress("v3_sum", &v3_sum);
        config_tree->SetBranchAddress("v4_sum", &v4_sum);
        config_tree->SetBranchAddress("n_tenngen_particles", &n_tenngen_particles);


        // output file
        TFile * fout = new TFile(output_file.c_str(), "RECREATE");
        // create same trees in output file
        TTree * event_tree_out = new TTree("tree", "tree");
        TTree * config_tree_out = new TTree("config", "config");

        double weight_out = 0.0;

        // output tree variables
        config_tree_out->Branch("nEvents", &nEvents, "nEvents/I");
        config_tree_out->Branch("seed", &seed, "seed/i");
        config_tree_out->Branch("max_nconst_truth", &max_nconst_truth, "max_nconst_truth/I");
        config_tree_out->Branch("max_pt_reco", &max_pt_reco, "max_pt_reco/D");
        config_tree_out->Branch("percentage_matched", &percentage_matched, "percentage_matched/D");
        config_tree_out->Branch("max_particle_eta", &max_particle_eta, "max_particle_eta/F");
        config_tree_out->Branch("min_particle_pt", &min_particle_pt, "min_particle_pt/F");
        config_tree_out->Branch("R", &R, "R/F");
        config_tree_out->Branch("max_jet_eta", &max_jet_eta, "max_jet_eta/F");
        config_tree_out->Branch("min_jet_pt", &min_jet_pt, "min_jet_pt/F");
        config_tree_out->Branch("min_pythia_jet_pt", &min_pythia_jet_pt, "min_pythia_jet_pt/F");
        config_tree_out->Branch("min_jet_area", &min_jet_area, "min_jet_area/F");
        config_tree_out->Branch("dRMax", &dRMax, "dRMax/F");
        config_tree_out->Branch("pT_hard_min", &pT_hard_min, "pT_hard_min/F");
        config_tree_out->Branch("pT_hard_max", &pT_hard_max, "pT_hard_max/F");
        config_tree_out->Branch("xsec_over_eventweight", &xsec_over_eventweight, "xsec_over_eventweight/D");
        config_tree_out->Branch("n_accepted_events", &n_accepted_events, "n_accepted_events/I");
        config_tree_out->Branch("n_generated_events", &n_generated_events, "n_generated_events/I");
        config_tree_out->Branch("n_processed_events", &n_processed_events, "n_processed_events/I");
        config_tree_out->Branch("sum_of_weights", &sum_of_weights, "sum_of_weights/D");
        config_tree_out->Branch("integrated_luminosity", &integrated_luminosity, "integrated_luminosity/D");
        config_tree_out->Branch("v2_ref_truth", &v2_ref_truth, "v2_ref_truth/D");
        config_tree_out->Branch("v3_ref_truth", &v3_ref_truth, "v3_ref_truth/D");
        config_tree_out->Branch("v4_ref_truth", &v4_ref_truth, "v4_ref_truth/D");
        config_tree_out->Branch("v2_sum", &v2_sum, "v2_sum/D");
        config_tree_out->Branch("v3_sum", &v3_sum, "v3_sum/D");
        config_tree_out->Branch("v4_sum", &v4_sum, "v4_sum/D");
        config_tree_out->Branch("n_tenngen_particles", &n_tenngen_particles, "n_tenngen_particles/I");

        event_tree_out->Branch("event_id", &event_id, "event_id/I");
        event_tree_out->Branch("weight", &weight_out, "weight/D");
        event_tree_out->Branch("psi1", &psi1, "psi1/D");
        event_tree_out->Branch("psi2", &psi2, "psi2/D");
        event_tree_out->Branch("psi3", &psi3, "psi3/D");
        event_tree_out->Branch("psi4", &psi4, "psi4/D");
        event_tree_out->Branch("average_pt", &average_pt, "average_pt/D");
        event_tree_out->Branch("n_forward_particles", &n_forward_particles, "n_forward_particles/I");
        event_tree_out->Branch("Q2_i", &Q2_i, "Q2_i/D");
        event_tree_out->Branch("Q2_r", &Q2_r, "Q2_r/D");
        event_tree_out->Branch("Q3_i", &Q3_i, "Q3_i/D");
        event_tree_out->Branch("Q3_r", &Q3_r, "Q3_r/D");
        event_tree_out->Branch("Q4_i", &Q4_i, "Q4_i/D");
        event_tree_out->Branch("Q4_r", &Q4_r, "Q4_r/D");
        event_tree_out->Branch("Q6_i", &Q6_i, "Q6_i/D");
        event_tree_out->Branch("Q6_r", &Q6_r, "Q6_r/D");
        event_tree_out->Branch("Q8_i", &Q8_i, "Q8_i/D");
        event_tree_out->Branch("Q8_r", &Q8_r, "Q8_r/D");
        event_tree_out->Branch("jet_pt_truth", &jet_pt_truth, "jet_pt_truth/D");
        event_tree_out->Branch("jet_phi_truth", &jet_phi_truth, "jet_phi_truth/D");
        event_tree_out->Branch("jet_nconst_truth", &jet_nconst_truth, "jet_nconst_truth/I");
        event_tree_out->Branch("jet_v2_truth", &jet_v2_truth, "jet_v2_truth/D");
        event_tree_out->Branch("jet_v3_truth", &jet_v3_truth, "jet_v3_truth/D");
        event_tree_out->Branch("jet_v4_truth", &jet_v4_truth, "jet_v4_truth/D");
        event_tree_out->Branch("jet_pt_reco", &jet_pt_reco, "jet_pt_reco/D");
        event_tree_out->Branch("jet_nconst_reco", &jet_nconst_reco, "jet_nconst_reco/I");
        event_tree_out->Branch("jet_phi_reco", &jet_phi_reco, "jet_phi_reco/D");
        event_tree_out->Branch("n_unmatched_reco_jets", &n_unmatched_reco_jets, "n_unmatched_reco_jets/I");
        event_tree_out->Branch("unmatched_reco_pt", &unmatched_reco_pt);
        event_tree_out->Branch("unmatched_reco_nconst", &unmatched_reco_nconst);
        event_tree_out->Branch("unmatched_reco_phi", &unmatched_reco_phi);


        
        // loop over events
        int n_entries = event_tree->GetEntries();
        for (int i = 0; i < n_entries; i++)
        {
            event_tree->GetEntry(i);
            weight_out = 0.0;
            event_tree_out->Fill();
        }


        // get list of TH1Ds in file
        TList * list = fin->GetListOfKeys();
        TIter next(list);
        TKey * key;
        while(key = (TKey*)next())
        {
            if(key->GetClassName() == TString("TTree")) continue;
            else if (key->GetClassName() == TString("TH1D"))
            {
                TH1D * h1 = (TH1D*)key->ReadObj();
                fout->cd();
                h1->Write();
            }
        }


        fin->Close();


    }

    std::cout << "Reweighting complete" << std::endl;
   

    return 0;
}

