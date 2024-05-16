#include "MultSubCalibrator.h"
#include "Utils.h"

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

int MultSubCalibrator::Calibrate(bool overwrite)
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
    
    if (m_calib_file_name.empty())
    {
        std::cout << "Calibration file name not set. Using default name" << std::endl;
        m_calib_file_name = m_output_dir + "/mult_calib.root";
    }
    // m_calib_file_name = m_output_dir + "/mult_calib.root";

    if(!overwrite && !gSystem->AccessPathName(m_calib_file_name.c_str()))
    {
        std::cout << "Calibration file " << m_calib_file_name << " already exists. Use --overwrite to overwrite" << std::endl;
        return 1;
    }

    std::cout << "Calibrating multiplicity correction" << std::endl;
    std::cout << "Merging config trees to " << m_calib_file_name << std::endl;

    // create input list
    m_input_list_name = m_output_dir + "/input_list.txt";
    std::vector<std::string> input_files = Utils::GetListFromDir(m_input_dir, m_input_list_name);
    Utils::MergeTrees(m_input_list_name, m_calib_file_name, "config");

    // pT bins
    double pt_max = Utils::BranchMax(m_calib_file_name, "config", "max_pt_reco");
    const int n_bins_pt = (int)(pt_max + 1.0);
    const double delta_pT = (pt_max + 1.0)/n_bins_pt;
    double pt_bins[n_bins_pt+1];
    for (int i = 0; i < n_bins_pt+1; i++){  pt_bins[i] =  i*delta_pT; }
    std::vector<double> pt_bins_vec;

    if(m_user_pt_binning)
    {
        pt_bins_vec = m_pt_bins;
    }
    else
    {
        for (int i = 0; i < n_bins_pt+1; i++){  pt_bins_vec.push_back(pt_bins[i]); }
    }
    
    
    // nconst bins
    double max_nconst_truth = Utils::BranchMax(m_calib_file_name, "config", "max_nconst_truth");
    const int n_bins_nconst = (int)(max_nconst_truth + 1.0);
    const double delta_nconst = (max_nconst_truth + 1.0)/n_bins_nconst;
    double nconst_bins[n_bins_nconst+1];
    for (int i = 0; i < n_bins_nconst+1; i++){  nconst_bins[i] = i*delta_nconst; }

    // create the histogram
    TH2D * h2_nconst_truth_vs_reco_pt = new TH2D("h2_nconst_truth_vs_reco_pt", ";p_{T}^{reco} [GeV];n_{const}^{truth}"
    , pt_bins_vec.size()-1, pt_bins_vec.data(), n_bins_nconst, nconst_bins);

    TH1D * h1_jet_pt_reco = new TH1D("h1_jet_pt_reco", ";p_{T}^{reco} [GeV]", pt_bins_vec.size()-1, pt_bins_vec.data());
    // TH1D * h1_jet_nconst_truth = new TH1D("h1_jet_nconst_truth", ";n_{const}^{truth}", n_bins_nconst, nconst_bins);
    std::cout << "Filling calibration histograms" << std::endl;
    // loop over input files
    for (auto &input_file : input_files)
    {
        // read the input file
        TFile * fin = new TFile(input_file.c_str(), "READ");
        if(!fin->IsOpen())
        {
            std::cout << "Error: could not open input file " << input_file << std::endl;
            exit(-1);
        }

        TTree * event_tree = (TTree*)fin->Get("tree");
        if(!event_tree)
        {
            std::cout << "Error: could not find tree in input file " << input_file << std::endl;
            exit(-1);
        }

        // input variables
        double jet_pt_reco = 0.0;
        int jet_nconst_truth = 0;

        // set up input branches
        event_tree->SetBranchAddress("jet_pt_reco", &jet_pt_reco);
        event_tree->SetBranchAddress("jet_nconst_truth", &jet_nconst_truth);

        // loop over the events
        int n_entries = event_tree->GetEntries();
        for (int i = 0; i < n_entries; i++)
        {
            event_tree->GetEntry(i);
            h2_nconst_truth_vs_reco_pt->Fill(jet_pt_reco, jet_nconst_truth);
            h1_jet_pt_reco->Fill(jet_pt_reco);
        }

        fin->Close();

    }

    // create the output histograms
    TProfile * p1_nconst_truth_vs_reco_pt = h2_nconst_truth_vs_reco_pt->ProfileX();
    p1_nconst_truth_vs_reco_pt->SetName("p1_nconst_truth_vs_reco_pt");            
    TH1D * h1_nconst_truth_vs_reco_pt = (TH1D*)p1_nconst_truth_vs_reco_pt->ProjectionX();
    h1_nconst_truth_vs_reco_pt->SetName("h1_nconst_truth_vs_reco_pt");


    // save the histogram to config file
   
    std::cout << "Successfully created calibration histograms. Writing file." << std::endl;
    TFile * fout = new TFile(m_calib_file_name.c_str(), "UPDATE");
    h2_nconst_truth_vs_reco_pt->Write();
    p1_nconst_truth_vs_reco_pt->Write();
    h1_jet_pt_reco->Write();
    h1_nconst_truth_vs_reco_pt->Write();
    
    // copy histograms from input file to output file
    TFile * fin = new TFile(input_files[0].c_str(), "READ");
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
    
    fout->Close();
    fin->Close();
    
    std::cout << "Done calibrating multiplicity correction" << std::endl;
    return 0;
}