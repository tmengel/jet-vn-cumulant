#include "Utils.h"
#include "MultSubtractor.h"

#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <cstdlib>

#include <TSystem.h>
#include <TString.h>
#include <TFile.h>
#include <TTree.h>
#include <TKey.h>
#include <TProfile.h>
#include <TH1D.h>
#include <TH2D.h>


int MultSubtractor::Subtract(bool overwrite)
{

    // make output directory
    if (m_output_dir.empty())
    {
        std::cerr << "Output directory not set" << std::endl;
        exit(1);
    }

    if(gSystem->AccessPathName(m_output_dir.c_str()))
    {
        std::cout << "Error: output directory " << m_output_dir << " does not exist" << std::endl;
        std::cout << "Creating output directory " << m_output_dir << std::endl;
        
        gSystem->mkdir(m_output_dir.c_str(), kTRUE);
    }

    m_input_files.clear();
    m_input_files = Utils::GetFilesFromDir(m_input_dir);
    std::cout << "Number of input files: " << m_input_files.size() << std::endl;

    
    if(!overwrite)
    {
        // get files that already exist
        std::vector<std::string> existing_files = Utils::GetFilesFromDir(m_output_dir);
        if(existing_files.size() != 0)
        {
            std::cout << "Files already exist in output directory. Use --overwrite to overwrite" << std::endl;
            return 1;
        }
       
    }


    // read the calibration file
    TFile * fcal = new TFile(m_calib_file_name.c_str(), "READ");
    if(!fcal->IsOpen())
    {
        std::cout << "Error: could not open calibration file " << m_calib_file_name << std::endl;
        exit(1);
    }

    TH1D * h1_nconst_truth_vs_reco_pt = (TH1D*)fcal->Get("h1_nconst_truth_vs_reco_pt");
    if (!h1_nconst_truth_vs_reco_pt)
    {
        std::cout << "Error: could not find histogram h1_nconst_truth_vs_reco_pt in calibration file " << m_calib_file_name << std::endl;
        exit(1);
    }

    // m_input_list_name = m_output_dir + "/input_list.txt";
    // std::vector<std::string> input_files = Utils::GetListFromDir(m_input_dir, m_input_list_name);
    for (auto &input_file : m_input_files)
    {
        // subtract the multiplicity
        std::string subtracted_file = GetSubFileName(input_file);
       
        std::cout << "Subtracting multiplicity from " << input_file << " to " << subtracted_file << std::endl;
        // read the input file
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
        
            // input variables
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

        // set up input branches
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

        // create the output file
        TFile * fout = new TFile(subtracted_file.c_str(), "RECREATE");

        // create new output tree
        TTree * event_tree_subtracted = new TTree("tree", "tree");
        // output variables
        double jet_pt_reco_subtracted = 0.0;
        double jet_nsignal_reco = 0.0;
        std::vector<double> unmatched_reco_pt_subtracted;
        std::vector<double> unmatched_reco_nsingal;

        // set up output branches
        event_tree_subtracted->Branch("event_id", &event_id, "event_id/I");
        event_tree_subtracted->Branch("weight", &weight, "weight/D");
        event_tree_subtracted->Branch("average_pt", &average_pt, "average_pt/D");
        event_tree_subtracted->Branch("psi1", &psi1, "psi1/D");
        event_tree_subtracted->Branch("psi2", &psi2, "psi2/D");
        event_tree_subtracted->Branch("psi3", &psi3, "psi3/D");
        event_tree_subtracted->Branch("psi4", &psi4, "psi4/D");
        event_tree_subtracted->Branch("n_forward_particles", &n_forward_particles, "n_forward_particles/I");
        event_tree_subtracted->Branch("Q2_i", &Q2_i, "Q2_i/D");
        event_tree_subtracted->Branch("Q2_r", &Q2_r, "Q2_r/D");
        event_tree_subtracted->Branch("Q3_i", &Q3_i, "Q3_i/D");
        event_tree_subtracted->Branch("Q3_r", &Q3_r, "Q3_r/D");
        event_tree_subtracted->Branch("Q4_i", &Q4_i, "Q4_i/D");
        event_tree_subtracted->Branch("Q4_r", &Q4_r, "Q4_r/D");
        event_tree_subtracted->Branch("Q6_i", &Q6_i, "Q6_i/D");
        event_tree_subtracted->Branch("Q6_r", &Q6_r, "Q6_r/D");
        event_tree_subtracted->Branch("Q8_i", &Q8_i, "Q8_i/D");
        event_tree_subtracted->Branch("Q8_r", &Q8_r, "Q8_r/D");
        event_tree_subtracted->Branch("jet_pt_truth", &jet_pt_truth, "jet_pt_truth/D");
        event_tree_subtracted->Branch("jet_nconst_truth", &jet_nconst_truth, "jet_nconst_truth/I");
        event_tree_subtracted->Branch("jet_phi_truth", &jet_phi_truth, "jet_phi_truth/D");
        event_tree_subtracted->Branch("jet_v2_truth", &jet_v2_truth, "jet_v2_truth/D");
        event_tree_subtracted->Branch("jet_v3_truth", &jet_v3_truth, "jet_v3_truth/D");
        event_tree_subtracted->Branch("jet_v4_truth", &jet_v4_truth, "jet_v4_truth/D");
        event_tree_subtracted->Branch("jet_pt_reco_unsubtracted", &jet_pt_reco, "jet_pt_reco/D");
        event_tree_subtracted->Branch("jet_pt_reco", &jet_pt_reco_subtracted, "jet_pt_reco/D");
        event_tree_subtracted->Branch("jet_nconst_reco", &jet_nconst_reco, "jet_nconst_reco/I");
        event_tree_subtracted->Branch("jet_nsignal_reco", &jet_nsignal_reco, "jet_nsignal_reco/D");
        event_tree_subtracted->Branch("jet_phi_reco", &jet_phi_reco, "jet_phi_reco/D");
        event_tree_subtracted->Branch("n_unmatched_reco_jets", &n_unmatched_reco_jets, "n_unmatched_reco_jets/I");
        event_tree_subtracted->Branch("unmatched_reco_pt_unsubtracted", &unmatched_reco_pt);
        event_tree_subtracted->Branch("unmatched_reco_nconst", &unmatched_reco_nconst);
        event_tree_subtracted->Branch("unmatched_reco_nsignal", &unmatched_reco_nsingal);
        event_tree_subtracted->Branch("unmatched_reco_pt", &unmatched_reco_pt_subtracted);
        event_tree_subtracted->Branch("unmatched_reco_phi", &unmatched_reco_phi);


        int n_entries = event_tree->GetEntries();
        // std::cout << "Looping over " << n_entries << " entries" << std::endl;
        for (int i = 0; i < n_entries; i++)
        {
            event_tree->GetEntry(i);
            if(jet_pt_reco > 0)
            {
                jet_nsignal_reco = h1_nconst_truth_vs_reco_pt->GetBinContent(h1_nconst_truth_vs_reco_pt->FindBin(jet_pt_reco));
                double correction = average_pt*(jet_nconst_reco - jet_nsignal_reco);
                jet_pt_reco_subtracted = jet_pt_reco - correction;
            }
            else
            {
                jet_nsignal_reco = -1;
                jet_pt_reco_subtracted = -1;
            }
            unmatched_reco_pt_subtracted.clear();
            unmatched_reco_nsingal.clear();
            for (int j = 0; j < n_unmatched_reco_jets; j++)
            {
                double n_signal_unmatched = h1_nconst_truth_vs_reco_pt->GetBinContent(h1_nconst_truth_vs_reco_pt->FindBin(unmatched_reco_pt->at(j)));
                double correction_unmatched = average_pt*(unmatched_reco_nconst->at(j) - n_signal_unmatched);
                unmatched_reco_nsingal.push_back(n_signal_unmatched);
                unmatched_reco_pt_subtracted.push_back(unmatched_reco_pt->at(j) - correction_unmatched);
            }

            event_tree_subtracted->Fill();
        }

        // copy all contents of the input file to the output file except the tree
        TTree * config_tree = (TTree*)fin->Get("config");
        if(config_tree)
        {
            fout->cd();
            TTree * config_tree_copy = config_tree->CloneTree(-1, "fast");
            config_tree_copy->Write();
        }

        fout->cd();
        event_tree_subtracted->Write();

        // copy histograms from input file to output file
        TList * list = fin->GetListOfKeys();
        TIter next(list);
        TKey * key;
        while(key = (TKey*)next())
        {
            if(key->GetClassName() == TString("TTree")) continue;
            else if (key->GetClassName() == TString("TProfile"))
            {
                TProfile * prof = (TProfile*)key->ReadObj();
                fout->cd();
                prof->Write();
            }
            else if (key->GetClassName() == TString("TH1D"))
            {
                TH1D * h1 = (TH1D*)key->ReadObj();
                fout->cd();
                h1->Write();
            }
            else if (key->GetClassName() == TString("TH2D"))
            {
                TH2D * h2 = (TH2D*)key->ReadObj();
                fout->cd();
                h2->Write();
            }
        }

        // copy multiplicity calibration histograms
        list = fcal->GetListOfKeys();
        next = TIter(list);
        while(key = (TKey*)next())
        {
            if (key->GetClassName() == TString("TH1D"))
            {
                TH1D * h1 = (TH1D*)key->ReadObj();
                fout->cd();
                h1->Write();
            }
            else if (key->GetClassName() == TString("TH2D"))
            {
                TH2D * h2 = (TH2D*)key->ReadObj();
                fout->cd();
                h2->Write();
            }
            else if (key->GetClassName() == TString("TProfile"))
            {
                TProfile * prof = (TProfile*)key->ReadObj();
                fout->cd();
                prof->Write();
            }
        }

        fin->Close();
        fout->Close();
    }

    fcal->Close();

   
    std::cout << "Done subtracting multiplicity" << std::endl;
    return 0;
}

std::string MultSubtractor::GetSubFileName(std::string inputfile)
{
     // subtract the multiplicity
    TString inputfile_base = inputfile.c_str();
    inputfile_base.ReplaceAll(".root", "");
    // remove the path (everything before the last /)
    inputfile_base = inputfile_base(inputfile_base.Last('/')+1, inputfile_base.Length()-inputfile_base.Last('/')+1);
    TString subtracted_file = Form("%s/%s_subtracted.root", m_output_dir.c_str(), inputfile_base.Data());
    return subtracted_file.Data();
}