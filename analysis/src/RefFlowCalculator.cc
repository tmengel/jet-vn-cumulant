#include "Utils.h"
#include "RefFlowCalculator.h"
#include "CorrFunctions.h"

#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <cstdlib>
#include <map>
#include <utility>

#include <TSystem.h>
#include <TString.h>
#include <TFile.h>
#include <TTree.h>
#include <TComplex.h>

void RefFlowCalculator::ReadTruthFile()
{
    TFile * f = new TFile(m_truth_file_name.c_str(), "READ");
    if(!f->IsOpen())
    {
        std::cout << "Error: could not open truh file " << m_truth_file_name << std::endl;
        exit(1);
    }

    TTree * t = (TTree*)f->Get("config");
    if(!t)
    {
        std::cout << "Error: could not find config tree in truth file " << m_truth_file_name << std::endl;
        exit(1);
    }
    double v2_ref_truth = 0.0;
    double v3_ref_truth = 0.0;
    double v4_ref_truth = 0.0;

    t->SetBranchAddress("v2_ref_truth", &v2_ref_truth);
    t->SetBranchAddress("v3_ref_truth", &v3_ref_truth);
    t->SetBranchAddress("v4_ref_truth", &v4_ref_truth);
    t->GetEntry(0);

    m_v2_ref_truth = v2_ref_truth;
    m_v3_ref_truth = v3_ref_truth;
    m_v4_ref_truth = v4_ref_truth;

    f->Close();

    return;
}

int RefFlowCalculator::Calculate()
{

    std::cout << "Starting RefFlowCalculator::Calculate()" << std::endl;
    // make output directory
    if (m_output_file_name.empty())
    {
        std::cerr << "Output filename not set" << std::endl;
        exit(1);
    }

    std::cout << "Calculating Ref flow" << std::endl;
    ReadTruthFile();

    // create input list
    // m_input_list_name = m_output_dir + "/input_list.txt";
    // std::vector<std::string> input_files = Utils::GetFilesFromDir(m_input_dir);
    m_input_files.clear();
    m_input_files = Utils::GetFilesFromDir(m_input_dir);
    std::cout << "Number of input files: " << m_input_files.size() << std::endl;

    std::map<int, std::vector<double>> two_part_event_avg {};
    std::map<int, std::vector<double>> four_part_event_avg {};
    std::map<int, std::vector<double>> two_part_cn {};
    std::map<int, std::vector<double>> four_part_cn {};
    std::map<int, std::vector<double>> two_part_vn {};
    std::map<int, std::vector<double>> four_part_vn {};
    std::map<int, std::pair<double, double>> vn2_mean_std {};
    std::map<int, std::pair<double, double>> vn4_mean_std {};

    for (int iharm =2; iharm < 5; iharm++)
    {
        std::vector<double> two_part_ref_vec{};
        std::vector<double> four_part_ref_vec{};

        CalculateHarm(iharm, two_part_ref_vec, four_part_ref_vec);

        two_part_event_avg[iharm] = two_part_ref_vec;
        four_part_event_avg[iharm] = four_part_ref_vec;

        std::vector<double> c2_vec = C2Vec(two_part_ref_vec);
        std::vector<double> c4_vec = C4Vec(two_part_ref_vec, four_part_ref_vec);
        std::vector<double> vn2_vec = vn2Vec(c2_vec);
        std::vector<double> vn4_vec = vn4Vec(c4_vec);

        two_part_cn[iharm] = c2_vec;
        four_part_cn[iharm] = c4_vec;
        two_part_vn[iharm] = vn2_vec;
        four_part_vn[iharm] = vn4_vec;

        std::pair<double, double> vn2_mean_std_pair = MeanErrVec(vn2_vec);
        std::pair<double, double> vn4_mean_std_pair = MeanErrVec(vn4_vec);

        vn2_mean_std[iharm] = vn2_mean_std_pair;
        vn4_mean_std[iharm] = vn4_mean_std_pair;
    }

    // write the output file
    TFile * fout = new TFile(m_output_file_name.c_str(), "RECREATE");
    TTree * event_tree_results = new TTree("tree", "tree");

    event_tree_results->Branch("v2_ref_truth", &m_v2_ref_truth, "v2_ref_truth/D");
    event_tree_results->Branch("v3_ref_truth", &m_v3_ref_truth, "v3_ref_truth/D");
    event_tree_results->Branch("v4_ref_truth", &m_v4_ref_truth, "v4_ref_truth/D");
    event_tree_results->Branch("v2_two_ref", &vn2_mean_std[2].first, "v2_ref/D");
    event_tree_results->Branch("v2_two_ref_err", &vn2_mean_std[2].second, "v2_ref_err/D");
    event_tree_results->Branch("v3_two_ref", &vn2_mean_std[3].first, "v3_ref/D");
    event_tree_results->Branch("v3_two_ref_err", &vn2_mean_std[3].second, "v3_ref_err/D");
    event_tree_results->Branch("v4_two_ref", &vn2_mean_std[4].first, "v4_ref/D");
    event_tree_results->Branch("v4_two_ref_err", &vn2_mean_std[4].second, "v4_ref_err/D");
    event_tree_results->Branch("v2_four_ref", &vn4_mean_std[2].first, "v2_ref/D");
    event_tree_results->Branch("v2_four_ref_err", &vn4_mean_std[2].second, "v2_ref_err/D");
    event_tree_results->Branch("v3_four_ref", &vn4_mean_std[3].first, "v3_ref/D");
    event_tree_results->Branch("v3_four_ref_err", &vn4_mean_std[3].second, "v3_ref_err/D");
    event_tree_results->Branch("v4_four_ref", &vn4_mean_std[4].first, "v4_ref/D");
    event_tree_results->Branch("v4_four_ref_err", &vn4_mean_std[4].second, "v4_ref_err/D");
    event_tree_results->Fill();

    TTree * event_tree_vecs = new TTree("vecs", "vecs");
    double two_part_event_avg_2_val, two_part_event_avg_3_val, two_part_event_avg_4_val;
    double four_part_event_avg_2_val, four_part_event_avg_3_val, four_part_event_avg_4_val;
    double two_part_cn_2_val, two_part_cn_3_val, two_part_cn_4_val;
    double four_part_cn_2_val, four_part_cn_3_val, four_part_cn_4_val;
    double two_part_vn_2_val, two_part_vn_3_val, two_part_vn_4_val;
    double four_part_vn_2_val, four_part_vn_3_val, four_part_vn_4_val;
    int n_subsamples = two_part_event_avg[2].size();
    event_tree_vecs->Branch("n_subsamples", &n_subsamples, "n_subsamples/I");
    event_tree_vecs->Branch("two_part_event_avg_2", &two_part_event_avg_2_val, "two_part_event_avg_2/D");
    event_tree_vecs->Branch("two_part_event_avg_3", &two_part_event_avg_3_val, "two_part_event_avg_3/D");
    event_tree_vecs->Branch("two_part_event_avg_4", &two_part_event_avg_4_val, "two_part_event_avg_4/D");
    event_tree_vecs->Branch("four_part_event_avg_2", &four_part_event_avg_2_val, "four_part_event_avg_2/D");
    event_tree_vecs->Branch("four_part_event_avg_3", &four_part_event_avg_3_val, "four_part_event_avg_3/D");
    event_tree_vecs->Branch("four_part_event_avg_4", &four_part_event_avg_4_val, "four_part_event_avg_4/D");
    event_tree_vecs->Branch("two_part_cn_2", &two_part_cn_2_val, "two_part_cn_2/D");
    event_tree_vecs->Branch("two_part_cn_3", &two_part_cn_3_val, "two_part_cn_3/D");
    event_tree_vecs->Branch("two_part_cn_4", &two_part_cn_4_val, "two_part_cn_4/D");
    event_tree_vecs->Branch("four_part_cn_2", &four_part_cn_2_val, "four_part_cn_2/D");
    event_tree_vecs->Branch("four_part_cn_3", &four_part_cn_3_val, "four_part_cn_3/D");
    event_tree_vecs->Branch("four_part_cn_4", &four_part_cn_4_val, "four_part_cn_4/D");
    event_tree_vecs->Branch("two_part_vn_2", &two_part_vn_2_val, "two_part_vn_2/D");
    event_tree_vecs->Branch("two_part_vn_3", &two_part_vn_3_val, "two_part_vn_3/D");
    event_tree_vecs->Branch("two_part_vn_4", &two_part_vn_4_val, "two_part_vn_4/D");
    event_tree_vecs->Branch("four_part_vn_2", &four_part_vn_2_val, "four_part_vn_2/D");
    event_tree_vecs->Branch("four_part_vn_3", &four_part_vn_3_val, "four_part_vn_3/D");
    event_tree_vecs->Branch("four_part_vn_4", &four_part_vn_4_val, "four_part_vn_4/D");

    for (int i = 0; i < n_subsamples; i++)
    {
        two_part_event_avg_2_val = two_part_event_avg[2][i];
        two_part_event_avg_3_val = two_part_event_avg[3][i];
        two_part_event_avg_4_val = two_part_event_avg[4][i];
        four_part_event_avg_2_val = four_part_event_avg[2][i];
        four_part_event_avg_3_val = four_part_event_avg[3][i];
        four_part_event_avg_4_val = four_part_event_avg[4][i];
        two_part_cn_2_val = two_part_cn[2][i];
        two_part_cn_3_val = two_part_cn[3][i];
        two_part_cn_4_val = two_part_cn[4][i];
        four_part_cn_2_val = four_part_cn[2][i];
        four_part_cn_3_val = four_part_cn[3][i];
        four_part_cn_4_val = four_part_cn[4][i];
        two_part_vn_2_val = two_part_vn[2][i];
        two_part_vn_3_val = two_part_vn[3][i];
        two_part_vn_4_val = two_part_vn[4][i];
        four_part_vn_2_val = four_part_vn[2][i];
        four_part_vn_3_val = four_part_vn[3][i];
        four_part_vn_4_val = four_part_vn[4][i];
        event_tree_vecs->Fill();
    }

    fout->cd();
    event_tree_results->Write();
    event_tree_vecs->Write();
    fout->Close();   
    // fcal->Close();

    return 0;

}

int RefFlowCalculator::CalculateHarm(int harmonic, std::vector<double> &two_part_ref_vec, std::vector<double> &four_part_ref_vec)
{

    two_part_ref_vec.clear();
    four_part_ref_vec.clear();

    for (auto &input_file : m_input_files)
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
        
        // input variables
        int event_id = 0;
        int n_forward_particles = 0;
        double Qn_i = 0.0;
        double Qn_r = 0.0;
        double Q2n_i = 0.0;
        double Q2n_r = 0.0;

        // set up input branches
        event_tree->SetBranchAddress("event_id", &event_id);
        event_tree->SetBranchAddress("n_forward_particles", &n_forward_particles);
        event_tree->SetBranchAddress(Form("Q%d_i", harmonic), &Qn_i);
        event_tree->SetBranchAddress(Form("Q%d_r", harmonic), &Qn_r);
        event_tree->SetBranchAddress(Form("Q%d_i", int(2*harmonic)), &Q2n_i);
        event_tree->SetBranchAddress(Form("Q%d_r", int(2*harmonic)), &Q2n_r);

        double two_part_ref_sum = 0.0;
        double four_part_ref_sum = 0.0;
        double two_part_ref_weight_sum = 0.0;
        double four_part_ref_weight_sum = 0.0;

        
        // loop over the events
        int n_entries = event_tree->GetEntries();
        for (int i = 0; i < n_entries; i++)
        {
            event_tree->GetEntry(i);
            
            TComplex Qn(Qn_r, Qn_i);
            TComplex Q2n(Q2n_r, Q2n_i);

            double two_part_ref = CorrFunctions::TwoPartRef(Qn, n_forward_particles);
            double four_part_ref = CorrFunctions::FourPartRef(Qn, Q2n, n_forward_particles);
            double two_part_ref_weight = CorrFunctions::TwoPartRefWeight(n_forward_particles);
            double four_part_ref_weight = CorrFunctions::FourPartRefWeight(n_forward_particles);

            two_part_ref_sum += two_part_ref*two_part_ref_weight;
            four_part_ref_sum += four_part_ref*four_part_ref_weight;
            two_part_ref_weight_sum += two_part_ref_weight;
            four_part_ref_weight_sum += four_part_ref_weight;
        }

        double two_part_avg_ref = two_part_ref_sum / two_part_ref_weight_sum;
        double four_part_avg_ref = four_part_ref_sum / four_part_ref_weight_sum;

        two_part_ref_vec.push_back(two_part_avg_ref);
        four_part_ref_vec.push_back(four_part_avg_ref);
        fin->Close();
    }

    return 0;


}


std::pair<double,double> RefFlowCalculator::MeanErrVec(const std::vector<double> &vec)
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

std::vector<double> RefFlowCalculator::C2Vec(const std::vector<double> &two_avg_vec)
{
    std::vector<double> c2_vec{};
    for (auto &v : two_avg_vec)
    {
        c2_vec.push_back(v);
    }
    return c2_vec;
}

std::vector<double> RefFlowCalculator::vn2Vec(const std::vector<double> &c2_vec)
{
    std::vector<double> vn2_vec{};
    for (auto &v : c2_vec)
    {
        double vn2 = TMath::Sqrt(v);
        vn2_vec.push_back(vn2);
    }
    return vn2_vec;
}

std::vector<double> RefFlowCalculator::C4Vec(const std::vector<double> &two_avg_vec, const std::vector<double> &four_avg_vec)
{
    std::vector<double> c4_vec{};
    for (int i = 0; i < two_avg_vec.size(); i++)
    {
        double c4 = four_avg_vec[i] - 2.0*two_avg_vec[i]*two_avg_vec[i];
        c4_vec.push_back(c4);
    }
    return c4_vec;
}

std::vector<double> RefFlowCalculator::vn4Vec( const std::vector<double> &c4_vec)
{
    std::vector<double> vn4_vec{};
    for (auto &v : c4_vec)
    {
        double x = -1.0*v;
        double vn4 = TMath::Power(x, 1.0/4.0);
        vn4_vec.push_back(vn4);
    }
    return vn4_vec;
}

