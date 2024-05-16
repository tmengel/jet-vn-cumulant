#include "Utils.h"
#include "FlowCalculator.h"

#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <cstdlib>
#include <cmath>
#include <algorithm>
#include <utility>

#include <TSystem.h>
#include <TString.h>
#include <TFile.h>
#include <TTree.h>
#include <TKey.h>
#include <TProfile.h>
#include <TH1D.h>
#include <TH2D.h>

int FlowCalculator::Calculate()
{

    // make output directory
    std::cout << "Starting FlowCalculator::Calculate()" << std::endl;
    // make output directory
    if (m_output_file_name.empty())
    {
        std::cerr << "Output filename not set" << std::endl;
        exit(1);
    }

    std::cout << "Calculating Flow" << std::endl;

  
    // read the reference flow values
    ReadRefFile();

    // read the calibration file
    // ReadCalibFile();
    std::cout << "Reading calib file " << m_calib_file_name << std::endl;
    TFile * f_calib = new TFile(m_calib_file_name.c_str(), "READ");
    if(!f_calib->IsOpen())
    {
        std::cerr << "Error: could not open calib file " << m_calib_file_name << std::endl;
        exit(1);
    }
    

    m_h1_jet_v2_func = (TH1D*)f_calib->Get("jet_v2_func_hist")->Clone("h1_jet_v2_func");
    m_h1_jet_v3_func = (TH1D*)f_calib->Get("jet_v3_func_hist")->Clone("h1_jet_v3_func");
    m_h1_jet_v4_func = (TH1D*)f_calib->Get("jet_v4_func_hist")->Clone("h1_jet_v4_func");

    // f_calib->Close();

    std::cout << "Read calib file " << m_calib_file_name << std::endl;

    std::cout << "Creating output file " << m_output_file_name << std::endl;
    TFile * fout = new TFile(m_output_file_name.c_str(), "RECREATE");

    // save the reference flow values
    TTree * tree = new TTree("tree", "tree");
    tree->Branch("v2_ref_truth", &m_v2_ref_truth);
    tree->Branch("v3_ref_truth", &m_v3_ref_truth);
    tree->Branch("v4_ref_truth", &m_v4_ref_truth);
    tree->Branch("v2_two_ref", &m_v2_two_ref);
    tree->Branch("v3_two_ref", &m_v3_two_ref);
    tree->Branch("v4_two_ref", &m_v4_two_ref);
    tree->Branch("v2_four_ref", &m_v2_four_ref);
    tree->Branch("v3_four_ref", &m_v3_four_ref);
    tree->Branch("v4_four_ref", &m_v4_four_ref);
    tree->Branch("v2_two_ref_err", &m_v2_two_ref_err);
    tree->Branch("v3_two_ref_err", &m_v3_two_ref_err);
    tree->Branch("v4_two_ref_err", &m_v4_two_ref_err);
    tree->Branch("v2_four_ref_err", &m_v2_four_ref_err);
    tree->Branch("v3_four_ref_err", &m_v3_four_ref_err);
    tree->Branch("v4_four_ref_err", &m_v4_four_ref_err);

    tree->Fill();

    fout->cd();
    tree->Write();
    m_h1_jet_v2_func->Write();
    m_h1_jet_v3_func->Write();
    m_h1_jet_v4_func->Write();

    std::cout << "Finished writing reference flow values" << std::endl;
    for (auto &input_file : m_unfold_input_files_harmonics)
    {
        std::string input_dir = input_file.first;
        int harmonic = input_file.second;
        std::cout << "Calculating Diff flow for harmonic " << harmonic << std::endl;
        m_input_files.clear();
        m_input_files = Utils::GetFilesFromDir(input_dir);
        std::cout << "Number of input files in unfolded directory " << input_dir << " is " << m_input_files.size() << std::endl;

        std::pair<TH1D*, TH1D*> vn_hists = GetVnHists(harmonic, m_input_files, "event_tree_vecs", "two_part_diff_measured", "four_part_diff_measured", "measured", 0);
        TH1D * h_vn2_measured = (TH1D*)vn_hists.first->Clone(Form("h1_v%d_two_diff_meas", harmonic));
        TH1D * h_vn4_measured = (TH1D*)vn_hists.second->Clone(Form("h1_v%d_four_diff_meas", harmonic));

        vn_hists = GetVnHists(harmonic, m_input_files, "event_tree_vecs", "two_part_diff_truth", "four_part_diff_truth", "truth", 0);
        TH1D * h_vn2_truth = (TH1D*)vn_hists.first->Clone(Form("h1_v%d_two_diff_truth", harmonic));
        TH1D * h_vn4_truth = (TH1D*)vn_hists.second->Clone(Form("h1_v%d_four_diff_truth", harmonic));

        // open first file to get number of unfolded vectors
        TFile * f = new TFile(m_input_files[0].c_str(), "READ");
        TTree * tree = (TTree*)f->Get("event_tree_unfolded");
        int n_iter_start = 0;
        tree->SetBranchAddress("iterunfolded", &n_iter_start);
        tree->GetEntry(0);
        int n_unfolded = tree->GetEntries();
        int niter = n_iter_start;
        f->Close();

        for (int i = 0; i < n_unfolded; i++)
        {
            vn_hists = GetVnHists(harmonic, m_input_files, "event_tree_unfolded", "two_part_diff_unfolded", "four_part_diff_unfolded", "unfolded", i);
            TH1D * h_vn2_unfolded = (TH1D*)vn_hists.first->Clone(Form("h1_v%d_two_diff_unfolded_%d", harmonic, i+n_iter_start));
            TH1D * h_vn4_unfolded = (TH1D*)vn_hists.second->Clone(Form("h1_v%d_four_diff_unfolded_%d", harmonic, i+n_iter_start));

            fout->cd();
            h_vn2_unfolded->Write();
            h_vn4_unfolded->Write();
        }

        fout->cd();
        h_vn2_measured->Write();
        h_vn4_measured->Write();
        h_vn2_truth->Write();
        h_vn4_truth->Write();
    }

    for (auto &input_file : m_input_files_harmonics)
    {
        std::string input_dir = input_file.first;
        int harmonic = input_file.second;
        m_input_files.clear();
        m_input_files = Utils::GetFilesFromDir(input_dir);
        std::cout << "Number of input files in directory " << input_dir << " is " << m_input_files.size() << std::endl;

        MakeTruthHists(harmonic, m_input_files);
    }

    fout->cd();
    m_h1_jet_v2_truth->Write();
    m_h1_jet_v3_truth->Write();
    m_h1_jet_v4_truth->Write();
    

    fout->Close();
    f_calib->Close();


    std::cout << "Finished calculating correlations" << std::endl;

    return 0;

}

void FlowCalculator::MakeTruthHists(const int harmonic, const std::vector<std::string> input_files)
{

    // get binning from first file
    TFile * f = new TFile(input_files[0].c_str(), "READ");
    if(!f->IsOpen())
    {
        std::cerr << "Error: could not open input file " << input_files[0] << std::endl;
        exit(1);
    }
    TTree * bin_tree = (TTree*)f->Get("bin_tree");
    if(!bin_tree)
    {
        std::cerr << "Error: could not find bin_tree in input file " << input_files[0] << std::endl;
        exit(1);
    }

    std::vector<double> * jet_vn_bins = 0;
    std::vector<double> * jet_pt_bins = 0;
    bin_tree->SetBranchAddress(Form("jet_v%d_truth_bins", harmonic), &jet_vn_bins);
    bin_tree->SetBranchAddress("pt_bins", &jet_pt_bins);
    bin_tree->GetEntry(0);
    std::vector<double> jet_vn_truth_bins = *jet_vn_bins;
    std::vector<double> jet_pt_truth_bins = *jet_pt_bins;

    f->Close();

    TH2D * h2_jet_vn_truth = new TH2D(Form("h2_jet_v%d_truth", harmonic), Form("h2_jet_v%d_truth", harmonic), jet_pt_truth_bins.size()-1, jet_pt_truth_bins.data(), jet_vn_truth_bins.size()-1, jet_vn_truth_bins.data());

    for (auto &input_file : input_files)
    {
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


        double jet_vn_truth = 0.0;
        double jet_pt_truth = 0.0;
        event_tree->SetBranchAddress(Form("jet_v%d_truth", harmonic), &jet_vn_truth);
        event_tree->SetBranchAddress("jet_pt_truth", &jet_pt_truth);

        for (int i = 0; i < event_tree->GetEntries(); i++)
        {
            event_tree->GetEntry(i);
            h2_jet_vn_truth->Fill(jet_pt_truth, jet_vn_truth);
        }

        fin->Close();
    }

    TProfile * p = h2_jet_vn_truth->ProfileX(Form("p_jet_v%d_truth", harmonic));
    TH1D * h1_jet_vn_truth = (TH1D*)p->ProjectionX(Form("h1_jet_v%d_truth_hist", harmonic));
    
    if(harmonic == 2)
    {
        m_h1_jet_v2_truth = (TH1D*)h1_jet_vn_truth->Clone("h1_jet_v2_truth");
    }
    else if(harmonic == 3)
    {
        m_h1_jet_v3_truth = (TH1D*)h1_jet_vn_truth->Clone("h1_jet_v3_truth");
    }
    else if(harmonic == 4)
    {
        m_h1_jet_v4_truth = (TH1D*)h1_jet_vn_truth->Clone("h1_jet_v4_truth");
    }

    return;
}

void FlowCalculator::ReadCalibFile()
{
    // open CalibFile

    std::cout << "Reading calib file " << m_calib_file_name << std::endl;
    TFile * f_calib = new TFile(m_calib_file_name.c_str(), "READ");
    if(!f_calib->IsOpen())
    {
        std::cerr << "Error: could not open calib file " << m_calib_file_name << std::endl;
        exit(1);
    }
    

    m_h1_jet_v2_func = (TH1D*)f_calib->Get("jet_v2_func_hist")->Clone("h1_jet_v2_func");
    m_h1_jet_v3_func = (TH1D*)f_calib->Get("jet_v3_func_hist")->Clone("h1_jet_v3_func");
    m_h1_jet_v4_func = (TH1D*)f_calib->Get("jet_v4_func_hist")->Clone("h1_jet_v4_func");

    f_calib->Close();

    std::cout << "Read calib file " << m_calib_file_name << std::endl;

    return;

    
}

void FlowCalculator::ReadRefFile()
{
    // open RefFile
    std::cout << "Reading ref file " << m_ref_flow_file_name << std::endl;
    TFile * f_ref = new TFile(m_ref_flow_file_name.c_str(), "READ");
    if(!f_ref->IsOpen())
    {
        std::cerr << "Error: could not open ref file " << m_ref_flow_file_name << std::endl;
        exit(1);
    }
    TTree * ref_tree = (TTree*)f_ref->Get("tree");
    if(!ref_tree)
    {
        std::cerr << "Error: could not find tree in ref file " << m_ref_flow_file_name << std::endl;
        exit(1);
    }

    double v2_ref_truth = 0.0;
    double v3_ref_truth = 0.0;
    double v4_ref_truth = 0.0;
    double v2_two_ref = 0.0;
    double v3_two_ref = 0.0;
    double v4_two_ref = 0.0;
    double v2_four_ref = 0.0;
    double v3_four_ref = 0.0;
    double v4_four_ref = 0.0;
    double v2_two_ref_err = 0.0;
    double v3_two_ref_err = 0.0;
    double v4_two_ref_err = 0.0;
    double v2_four_ref_err = 0.0;
    double v3_four_ref_err = 0.0;
    double v4_four_ref_err = 0.0;

    ref_tree->SetBranchAddress("v2_ref_truth", &v2_ref_truth);
    ref_tree->SetBranchAddress("v3_ref_truth", &v3_ref_truth);
    ref_tree->SetBranchAddress("v4_ref_truth", &v4_ref_truth);
    ref_tree->SetBranchAddress("v2_two_ref", &v2_two_ref);
    ref_tree->SetBranchAddress("v3_two_ref", &v3_two_ref);
    ref_tree->SetBranchAddress("v4_two_ref", &v4_two_ref);
    ref_tree->SetBranchAddress("v2_four_ref", &v2_four_ref);
    ref_tree->SetBranchAddress("v3_four_ref", &v3_four_ref);
    ref_tree->SetBranchAddress("v4_four_ref", &v4_four_ref);
    ref_tree->SetBranchAddress("v2_two_ref_err", &v2_two_ref_err);
    ref_tree->SetBranchAddress("v3_two_ref_err", &v3_two_ref_err);
    ref_tree->SetBranchAddress("v4_two_ref_err", &v4_two_ref_err);
    ref_tree->SetBranchAddress("v2_four_ref_err", &v2_four_ref_err);
    ref_tree->SetBranchAddress("v3_four_ref_err", &v3_four_ref_err);
    ref_tree->SetBranchAddress("v4_four_ref_err", &v4_four_ref_err);
    ref_tree->GetEntry(0);

    m_v2_ref_truth = v2_ref_truth;
    m_v3_ref_truth = v3_ref_truth;
    m_v4_ref_truth = v4_ref_truth;
    m_v2_two_ref = v2_two_ref;
    m_v3_two_ref = v3_two_ref;
    m_v4_two_ref = v4_two_ref;
    m_v2_four_ref = v2_four_ref;
    m_v3_four_ref = v3_four_ref;
    m_v4_four_ref = v4_four_ref;
    m_v2_two_ref_err = v2_two_ref_err;
    m_v3_two_ref_err = v3_two_ref_err;
    m_v4_two_ref_err = v4_two_ref_err;
    m_v2_four_ref_err = v2_four_ref_err;
    m_v3_four_ref_err = v3_four_ref_err;
    m_v4_four_ref_err = v4_four_ref_err;

    f_ref->Close();

    std::cout << "Read ref file " << m_ref_flow_file_name << std::endl;
    return;

}

std::pair<TH1D*, TH1D*> FlowCalculator::GetVnHists(const int harmonic, const std::vector<std::string> input_files, const std::string & treename, const std::string & v2branch, const std::string & v4branch, const std::string &name, const int entry)
{
 
    std::vector<std::vector<double>> vn_two_values{};
    std::vector<std::vector<double>> vn_four_values{};
    int n_subsamples = input_files.size();
    int nptbins =0;
    std::vector<double> ptbins{};

    for (auto &input_file : input_files)
    {
        // read the input file
        TFile * fin = new TFile(input_file.c_str(), "READ");
        if(!fin->IsOpen())
        {
            std::cout << "Error: could not open input file " << input_file << std::endl;
            exit(1);
        }

        TTree * event_tree_vecs = (TTree*)fin->Get(treename.c_str());
        if(!event_tree_vecs)
        {
            std::cout << "Error: could not find tree in input file " << input_file << std::endl;
            exit(1);
        }

        TTree * ref_tree = (TTree*)fin->Get("ref_tree");
        if(!ref_tree)
        {
            std::cout << "Error: could not find ref tree in input file " << input_file << std::endl;
            exit(1);
        }

        double incn_two_ref, incn_four_ref;
        ref_tree->SetBranchAddress("vn_two_ref", &incn_two_ref);
        ref_tree->SetBranchAddress("vn_four_ref", &incn_four_ref);
        ref_tree->GetEntry(0);
        
        double cn_two_ref = incn_two_ref*incn_two_ref;
        double cn_four_ref = -1.0*TMath::Power(incn_four_ref, 4.0);
        // input variables
        std::vector<double> * two_part_diff = 0;
        std::vector<double> * four_part_diff = 0;
        std::vector<double> * pt_bins = 0;

        int n_pt_bins = 0;
       
            event_tree_vecs->SetBranchAddress("n_pt_bins", &n_pt_bins);
            event_tree_vecs->SetBranchAddress("pt_bins", &pt_bins);
    
        event_tree_vecs->SetBranchAddress(v2branch.c_str(), &two_part_diff);
        event_tree_vecs->SetBranchAddress(v4branch.c_str(), &four_part_diff);
        
        
        event_tree_vecs->GetEntry(0);
        std::vector<double> two_part_vec = {};
        std::vector<double> four_part_vec = {};
        ptbins.clear();
        nptbins= n_pt_bins;
        for (int i = 0; i < n_pt_bins; i++)
        {
            two_part_vec.push_back(two_part_diff->at(i));
            four_part_vec.push_back(four_part_diff->at(i));
            ptbins.push_back(pt_bins->at(i));
        }

        
        std::vector<double> dn_two = D2Vec(two_part_vec);
        std::vector<double> dn_four = D4Vec(four_part_vec, two_part_vec, cn_two_ref);
        
        std::vector<double> vn_two = vn2Vec(dn_two, cn_two_ref);
        std::vector<double> vn_four = vn4Vec(dn_four, cn_four_ref);

        std::cout << "two_part_vec: ";
        for (auto &v : two_part_vec)
        {
            std::cout << v << " ";
        }
        std::cout << std::endl;
        std::cout << "four_part_vec: ";
        for (auto &v : four_part_vec)
        {
            std::cout << v << " ";
        }
        std::cout << std::endl;
        std::cout << "dn_two: ";
        for (auto &v : dn_two)
        {
            std::cout << v << " ";
        }
        std::cout << std::endl;
        std::cout << "dn_four: ";
        for (auto &v : dn_four)
        {
            std::cout << v << " ";
        }
        std::cout << std::endl;
        std::cout << "vn_two: ";
        for (auto &v : vn_two)
        {
            std::cout << v << " ";
        }
        std::cout << std::endl;
        std::cout << "vn_four: ";
        for (auto &v : vn_four)
        {
            std::cout << v << " ";
        }
        std::cout << std::endl;
        std::cout << "cn_two_ref = " << cn_two_ref << std::endl;
        std::cout << "cn_four_ref = " << cn_four_ref << std::endl;
        std::cout << std::endl;


        vn_two_values.push_back(vn_two);
        vn_four_values.push_back(vn_four);

        fin->Close();
    }

    std::vector<double> vn2{};
    std::vector<double> vn2_err{};
    std::vector<double> vn4{};
    std::vector<double> vn4_err{};

    for (int i = 0; i < nptbins; i++)
    {
        std::vector<double> vn_two_subsample{};
        std::vector<double> vn_four_subsample{};
        
        for (int j = 0; j < n_subsamples; j++)
        {
            vn_two_subsample.push_back(vn_two_values.at(j).at(i));
            vn_four_subsample.push_back(vn_four_values.at(j).at(i));
        }

        std::pair<double,double> vn_two_mean_err = MeanErrVec(vn_two_subsample);
        std::pair<double,double> vn_four_mean_err = MeanErrVec(vn_four_subsample);

        vn2.push_back(vn_two_mean_err.first);
        vn2_err.push_back(vn_two_mean_err.second);
        vn4.push_back(vn_four_mean_err.first);
        vn4_err.push_back(vn_four_mean_err.second);
    }

    TH1D * h_vn2 = GetHistFromVec(vn2, vn2_err, "h1_v" + std::to_string(harmonic) + "_two_diff_"+name, ptbins);
    TH1D * h_vn4 = GetHistFromVec(vn4, vn4_err, "h1_v" + std::to_string(harmonic) + "_four_diff_"+name, ptbins);

    return std::make_pair(h_vn2, h_vn4);
}


TH1D * FlowCalculator::GetHistFromVec(const std::vector<double> &vec, const std::vector<double> &err, const std::string &name, const std::vector<double> &bins)
{
    TH1D * h = new TH1D(name.c_str(), name.c_str(), bins.size()-1, bins.data());
    for (int i = 0; i < vec.size(); i++)
    {
        h->SetBinContent(i+1, vec.at(i));
        h->SetBinError(i+1, err.at(i));
    }
    return h;
}

std::vector<double> FlowCalculator::D2Vec(const std::vector<double> &two_avg_vec)
{
   
    std::vector<double> d2_vec{};
    for (auto &v : two_avg_vec)
    {
        d2_vec.push_back(v);
    }
    return d2_vec;
}

std::vector<double> FlowCalculator::vn2Vec(const std::vector<double> &d2_vec, const double &cn2)
{
    
    std::vector<double> vn2_vec{};
    for (auto &v : d2_vec)
    {
        double vn2 = v/TMath::Sqrt(cn2);
        vn2_vec.push_back(vn2);
    }
    return vn2_vec;
}

std::vector<double> FlowCalculator::D4Vec(const std::vector<double> &four_avg_vec, const std::vector<double> &two_avg_vec, const double &cn2)
{
    std::vector<double> d4_vec{};
    for (int i = 0; i < four_avg_vec.size(); i++)
    {
        double d4 = four_avg_vec.at(i) - (2.0*cn2*two_avg_vec.at(i));
        d4_vec.push_back(d4);
    }
    return d4_vec;
}

std::vector<double> FlowCalculator::vn4Vec(const std::vector<double> &d4_vec, const double &cn4)
{
    std::vector<double> vn4_vec{};
    for (auto &v : d4_vec)
    {
        double neg_c4 = -1.0*cn4;
        double denom = TMath::Power(neg_c4, 3.0/4.0);
        double neg_v = -1.0*v;
        double vn4 = neg_v/denom;
        vn4_vec.push_back(vn4);
    }
    return vn4_vec;
}

std::pair<double,double> FlowCalculator::MeanErrVec(const std::vector<double> &vec)
{
    double sum = 0.0;
    double sum_sq = 0.0;
    for (auto &v : vec)
    {
        sum += v;
        sum_sq += v*v;
    }
    double mean = sum / vec.size();
    double std = TMath::Sqrt(sum_sq / vec.size() - mean*mean);
    double std_err = std / TMath::Sqrt(vec.size());
    return std::make_pair(mean, std_err);

}


