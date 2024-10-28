#include "SimpleDiffFlow.h"

#include "CorrFunctions.h"
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
#include <TComplex.h>

int SimpleDiffFlow::Calculate()
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
    // m_input_files = GetFilesFromDir(m_input);
    if(m_input.find(".root") != std::string::npos)
    {
        std::cout << "Input is a file" << std::endl;
        m_input_files.push_back(m_input);
    }
    else
    {
        std::cout << "Input is a directory" << std::endl;
        m_input_files = Utils::GetFilesFromDir(m_input);
    }
    std::cout << "Found " << m_input_files.size() << " input files" << std::endl;

    for (auto &input_file : m_input_files)
    {
        // calculate correlations
        std::string cum_file = GetOutputFileName(input_file);
        std::cout << "Creating corr file " << cum_file << std::endl;
  
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

        TTree * bin_tree = (TTree*)fin->Get("bin_tree");
        if(!bin_tree)
        {
            std::cout << "Error: could not find bin_tree in input file " << input_file << std::endl;
            exit(1);
        }

        TTree * ref_tree = (TTree*)fin->Get("ref_tree");
        if(!ref_tree)
        {
            std::cout << "Error: could not find ref_tree in input file " << input_file << std::endl;
            exit(1);
        }
        
        // input variables
        int event_id = 0;
        double weight = 1.0;
        double jet_v2_truth = 0.0;
        double jet_v3_truth = 0.0;
        double jet_v4_truth = 0.0;
        int n_forward_particles = 0;        
        double jet_pt_truth = 0.0;
        double two_part_diff_truth = 0.0;
        double four_part_diff_truth = 0.0;
        double jet_pt_reco = 0.0;
        double two_part_diff_reco = 0.0;
        double four_part_diff_reco = 0.0;
        std::vector<double> * jet_pt_unmatched = 0;
        std::vector<double> * two_part_diff_unmatched = 0;
        std::vector<double> * four_part_diff_unmatched = 0;
       
        event_tree->SetBranchAddress("event_id", &event_id);
        event_tree->SetBranchAddress("weight", &weight);
        event_tree->SetBranchAddress("jet_v2_truth", &jet_v2_truth);
        event_tree->SetBranchAddress("jet_v3_truth", &jet_v3_truth);
        event_tree->SetBranchAddress("jet_v4_truth", &jet_v4_truth);
        event_tree->SetBranchAddress("n_forward_particles", &n_forward_particles);
        event_tree->SetBranchAddress("jet_pt_truth", &jet_pt_truth);
        event_tree->SetBranchAddress("two_part_diff_truth", &two_part_diff_truth);
        event_tree->SetBranchAddress("four_part_diff_truth", &four_part_diff_truth);
        event_tree->SetBranchAddress("jet_pt_reco", &jet_pt_reco);
        event_tree->SetBranchAddress("two_part_diff_reco", &two_part_diff_reco);
        event_tree->SetBranchAddress("four_part_diff_reco", &four_part_diff_reco);
        event_tree->SetBranchAddress("jet_pt_unmatched", &jet_pt_unmatched);
        event_tree->SetBranchAddress("two_part_diff_unmatched", &two_part_diff_unmatched);
        event_tree->SetBranchAddress("four_part_diff_unmatched", &four_part_diff_unmatched);


      
        std::vector<double> * in_pt_bins = 0;
      
        bin_tree->SetBranchAddress("pt_bins", &in_pt_bins);
        bin_tree->GetEntry(0);
        std::vector<double> pt_bins = *in_pt_bins;

        double two_part_event_avg_ref, four_part_event_avg_ref;
        double cn_two_ref, cn_four_ref;
        double vn_two_ref, vn_four_ref;
        ref_tree->SetBranchAddress("two_part_event_avg_ref", &two_part_event_avg_ref);
        ref_tree->SetBranchAddress("four_part_event_avg_ref", &four_part_event_avg_ref);
        ref_tree->SetBranchAddress("cn_two_ref", &cn_two_ref);
        ref_tree->SetBranchAddress("cn_four_ref", &cn_four_ref);
        ref_tree->SetBranchAddress("vn_two_ref", &vn_two_ref);
        ref_tree->SetBranchAddress("vn_four_ref", &vn_four_ref);
        ref_tree->GetEntry(0);
      

        //=======================================================================================================


        // create the output file
        TFile * fout = new TFile(cum_file.c_str(), "RECREATE");
      
        const int n_pt_bins = pt_bins.size()-1;
        std::vector<double> two_part_diff_truth_vec(n_pt_bins, 0.0);
        std::vector<double> four_part_diff_truth_vec(n_pt_bins, 0.0);
        std::vector<double> two_part_diff_meas_vec(n_pt_bins, 0.0);
        std::vector<double> four_part_diff_meas_vec(n_pt_bins, 0.0);

        std::vector<double> two_part_diff_w_2_sum(n_pt_bins, 0.0);
        std::vector<double> four_part_diff_w_4_sum(n_pt_bins, 0.0);
        std::vector<double> two_part_diff_w_2_meas_sum(n_pt_bins, 0.0);
        std::vector<double> four_part_diff_w_4_meas_sum(n_pt_bins, 0.0);

        int n_entries = event_tree->GetEntries();
        double w_2 = 0;
        double w_4 = 0;
        // fill the histograms
        for (int i = 0; i < n_entries; i++)
        {
            
            
            event_tree->GetEntry(i);
            w_2 = CorrFunctions::TwoPartDiffWeight(n_forward_particles, 1);
            w_4 = CorrFunctions::FourPartDiffWeight(n_forward_particles, 1);
           // fill vectors
            for (int j = 0; j < n_pt_bins; j++)
            {
                if (jet_pt_truth >= pt_bins[j] && jet_pt_truth < pt_bins[j+1])
                {
                    two_part_diff_truth_vec[j] += two_part_diff_truth * w_2;
                    four_part_diff_truth_vec[j] += four_part_diff_truth * w_4;
                    // two_part_diff_truth_vec[j] += two_part_diff_truth;
                    // four_part_diff_truth_vec[j] += four_part_diff_truth;
                    // if(two_part_diff_w_2_sum[j] == 0)
                    // {
                    //     two_part_diff_w_2_sum[j] = w_2;
                    // }
                    // if(four_part_diff_w_4_sum[j] == 0)
                    // {
                    //     four_part_diff_w_4_sum[j] = w_4;
                    // }
                    two_part_diff_w_2_sum[j] += w_2;
                    four_part_diff_w_4_sum[j] += w_4;
                }
                if (jet_pt_reco >= pt_bins[j] && jet_pt_reco < pt_bins[j+1])
                {
                    two_part_diff_meas_vec[j] += two_part_diff_reco * w_2;
                    four_part_diff_meas_vec[j] += four_part_diff_reco * w_4;
                    // two_part_diff_meas_vec[j] += two_part_diff_reco;
                    // four_part_diff_meas_vec[j] += four_part_diff_reco;
                    // if(two_part_diff_w_2_meas_sum[j] == 0)
                    // {
                    //     two_part_diff_w_2_meas_sum[j] = w_2;
                    // }
                    // if(four_part_diff_w_4_meas_sum[j] == 0)
                    // {
                    //     four_part_diff_w_4_meas_sum[j] = w_4;
                    // }
                    two_part_diff_w_2_meas_sum[j] += w_2;
                    four_part_diff_w_4_meas_sum[j] += w_4;
                }
            }
        
        }

        // normalize
        for (int i = 0; i < n_pt_bins; i++)
        {
            if (two_part_diff_w_2_sum[i] > 0)
            {
                two_part_diff_truth_vec[i] /= two_part_diff_w_2_sum[i];
            }
            if (four_part_diff_w_4_sum[i] > 0)
            {
                four_part_diff_truth_vec[i] /= four_part_diff_w_4_sum[i];
            }

            if (two_part_diff_w_2_meas_sum[i] > 0)
            {
                two_part_diff_meas_vec[i] /= two_part_diff_w_2_meas_sum[i];
            }

            if (four_part_diff_w_4_meas_sum[i] > 0)
            {
                four_part_diff_meas_vec[i] /= four_part_diff_w_4_meas_sum[i];
            }
        }

        // calc d2 and d4
        std::vector<double> d2_vec = D2Vec(two_part_diff_truth_vec);
        std::vector<double> d4_vec = D4Vec(four_part_diff_truth_vec, two_part_diff_truth_vec, cn_two_ref);
        std::vector<double> vn2_vec = vn2Vec(d2_vec, cn_two_ref);
        std::vector<double> vn4_vec = vn4Vec(d4_vec, cn_four_ref);

        std::vector<double> d2_vec_meas = D2Vec(two_part_diff_meas_vec);
        std::vector<double> d4_vec_meas = D4Vec(four_part_diff_meas_vec, two_part_diff_meas_vec, cn_two_ref);
        std::vector<double> vn2_vec_meas = vn2Vec(d2_vec_meas, cn_two_ref);
        std::vector<double> vn4_vec_meas = vn4Vec(d4_vec_meas, cn_four_ref);


        // create histograms
        TH1D * h_t2 = new TH1D("t2", "t2", n_pt_bins, pt_bins.data());
        TH1D * h_t4 = new TH1D("t4", "t4", n_pt_bins, pt_bins.data());
        TH1D * h_d2 = new TH1D("d2", "d2", n_pt_bins, pt_bins.data());
        TH1D * h_d4 = new TH1D("d4", "d4", n_pt_bins, pt_bins.data());
        TH1D * h_vn2 = new TH1D("vn2", "vn2", n_pt_bins, pt_bins.data());
        TH1D * h_vn4 = new TH1D("vn4", "vn4", n_pt_bins, pt_bins.data());

        TH1D * h_t2_meas = new TH1D("t2_meas", "t2_meas", n_pt_bins, pt_bins.data());
        TH1D * h_t4_meas = new TH1D("t4_meas", "t4_meas", n_pt_bins, pt_bins.data());
        TH1D * h_d2_meas = new TH1D("d2_meas", "d2_meas", n_pt_bins, pt_bins.data());
        TH1D * h_d4_meas = new TH1D("d4_meas", "d4_meas", n_pt_bins, pt_bins.data());
        TH1D * h_vn2_meas = new TH1D("vn2_meas", "vn2_meas", n_pt_bins, pt_bins.data());
        TH1D * h_vn4_meas = new TH1D("vn4_meas", "vn4_meas", n_pt_bins, pt_bins.data());


        for (int i = 0; i < n_pt_bins; i++)
        {
            h_t2->SetBinContent(i+1, two_part_diff_truth_vec[i]);
            h_t4->SetBinContent(i+1, four_part_diff_truth_vec[i]);
            h_d2->SetBinContent(i+1, d2_vec[i]);
            h_d4->SetBinContent(i+1, d4_vec[i]);
            h_vn2->SetBinContent(i+1, vn2_vec[i]);
            h_vn4->SetBinContent(i+1, vn4_vec[i]);

            h_t2_meas->SetBinContent(i+1, two_part_diff_meas_vec[i]);
            h_t4_meas->SetBinContent(i+1, four_part_diff_meas_vec[i]);
            h_d2_meas->SetBinContent(i+1, d2_vec_meas[i]);
            h_d4_meas->SetBinContent(i+1, d4_vec_meas[i]);
            h_vn2_meas->SetBinContent(i+1, vn2_vec_meas[i]);
            h_vn4_meas->SetBinContent(i+1, vn4_vec_meas[i]);

        }

        fout->cd();
        h_t2->Write();
        h_t4->Write();
        h_d2->Write();
        h_d4->Write();
        h_vn2->Write();
        h_vn4->Write();
        
        h_t2_meas->Write();
        h_t4_meas->Write();
        h_d2_meas->Write();
        h_d4_meas->Write();
        h_vn2_meas->Write();
        h_vn4_meas->Write();

        
        fout->Close();

        fin->Close();






        
    }

    std::cout << "Done" << std::endl;
    return 0;

}

std::vector<double> SimpleDiffFlow::D2Vec(const std::vector<double> &two_avg_vec)
{
   
    std::vector<double> d2_vec{};
    for (auto &v : two_avg_vec)
    {
        d2_vec.push_back(v);
    }
    return d2_vec;
}

std::vector<double> SimpleDiffFlow::vn2Vec(const std::vector<double> &d2_vec, const double &cn2)
{
    
    std::vector<double> vn2_vec{};
    for (auto &v : d2_vec)
    {
        double vn2 = v/TMath::Sqrt(cn2);
        vn2_vec.push_back(vn2);
    }
    return vn2_vec;
}

std::vector<double> SimpleDiffFlow::D4Vec(const std::vector<double> &four_avg_vec, const std::vector<double> &two_avg_vec, const double &cn2)
{
    std::vector<double> d4_vec{};
    for (int i = 0; i < four_avg_vec.size(); i++)
    {
        double d4 = four_avg_vec.at(i) - (2.0*cn2*two_avg_vec.at(i));
        d4_vec.push_back(d4);
    }
    return d4_vec;
}

std::vector<double> SimpleDiffFlow::vn4Vec(const std::vector<double> &d4_vec, const double &cn4)
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



std::vector<double> SimpleDiffFlow::ConvertToVec(TH2D * h2,  double w)
{
 
    TProfile * p = h2->ProfileX();
    TH1D * h1 = (TH1D*)p->ProjectionX();
    if (w != 1.0)
    {
        h1->Scale(w);
    }
    std::vector<double> vec;
    for (int i = 1; i <= h1->GetNbinsX(); i++)
    {
        vec.push_back(h1->GetBinContent(i));
    }
    return vec;
}

TH1D * SimpleDiffFlow::Make1D(TH2D * h2, std::string name, double w)
{
    TProfile * p = h2->ProfileX();
    TH1D * h1 = (TH1D*)p->ProjectionX(name.c_str());
    if (w != 1.0)
    {
        h1->Scale(w);
    }   
    return h1;
}

std::string SimpleDiffFlow::GetOutputFileName(std::string inputfile)
{
    // subtract the multiplicity
    TString inputfile_base = inputfile.c_str();
    inputfile_base.ReplaceAll(".root", "");
    // remove the path (everything before the last /)
    inputfile_base = inputfile_base(inputfile_base.Last('/')+1, inputfile_base.Length()-inputfile_base.Last('/')+1);
    TString outputfile = Form("%s/%s_simple.root", m_output_dir.c_str(), inputfile_base.Data());
    return outputfile.Data();
}