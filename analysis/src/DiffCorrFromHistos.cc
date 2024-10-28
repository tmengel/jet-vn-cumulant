#include "CorrFunctions.h"
#include "Utils.h"
#include "DiffCorrFromHistos.h"

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


int DiffCorrFromHistos::Calculate()
{
    
    std::cout << "Calculating differential correlations from histograms" << std::endl;
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


        std::vector<double> * in_jet_v2_truth_bins = 0;
        std::vector<double> * in_jet_v3_truth_bins = 0;
        std::vector<double> * in_jet_v4_truth_bins = 0;
        std::vector<double> * in_pt_bins = 0;
        std::vector<double> * in_two_part_diff_bins = 0;
        std::vector<double> * in_four_part_diff_bins = 0;

        bin_tree->SetBranchAddress("pt_bins", &in_pt_bins);
        // bin_tree->SetBranchAddress("truth_pt_bins", &in_truth_pt_bins);
        bin_tree->SetBranchAddress("jet_v2_truth_bins", &in_jet_v2_truth_bins);
        bin_tree->SetBranchAddress("jet_v3_truth_bins", &in_jet_v3_truth_bins);
        bin_tree->SetBranchAddress("jet_v4_truth_bins", &in_jet_v4_truth_bins);
        bin_tree->SetBranchAddress("two_diff_bins", &in_two_part_diff_bins);
        bin_tree->SetBranchAddress("four_diff_bins", &in_four_part_diff_bins);

        bin_tree->GetEntry(0);
        std::vector<double> jet_v2_truth_bins = *in_jet_v2_truth_bins;
        std::vector<double> jet_v3_truth_bins = *in_jet_v3_truth_bins;
        std::vector<double> jet_v4_truth_bins = *in_jet_v4_truth_bins;
        std::cout << "jet_v2_truth_bins: " << std::endl;
        // for (auto &bin : jet_v2_truth_bins)
        // {
        //     std::cout << "bin: " << bin << std::endl;
        // }
        std::vector<double> pt_bins = *in_pt_bins;
        std::vector<double> two_part_diff_bins = *in_two_part_diff_bins;
        std::vector<double> four_part_diff_bins = *in_four_part_diff_bins;

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
        TTree * event_tree_out = new TTree("ref_tree", "ref_tree");
        event_tree_out->Branch("two_part_event_avg_ref", &two_part_event_avg_ref);
        event_tree_out->Branch("four_part_event_avg_ref", &four_part_event_avg_ref);
        event_tree_out->Branch("cn_two_ref", &cn_two_ref);
        event_tree_out->Branch("cn_four_ref", &cn_four_ref);
        event_tree_out->Branch("vn_two_ref", &vn_two_ref);
        event_tree_out->Branch("vn_four_ref", &vn_four_ref);
        event_tree_out->Fill();

        int n_pt_bins = pt_bins.size()-1;
        std::vector<double> two_part_diff_truth_vec{};
        std::vector<double> two_part_diff_measured_vec{};
        std::vector<double> four_part_diff_truth_vec{};
        std::vector<double> four_part_diff_measured_vec{};
        std::vector<double> two_part_diff_unfolded_vec{};
        std::vector<double> four_part_diff_unfolded_vec{};
        TTree * event_tree_vecs = new TTree("event_tree_vecs", "event_tree_vecs");
        event_tree_vecs->Branch("n_pt_bins", &n_pt_bins, "n_pt_bins/I");
        event_tree_vecs->Branch("pt_bins", &pt_bins);
        event_tree_vecs->Branch("two_part_diff_truth", &two_part_diff_truth_vec);
        event_tree_vecs->Branch("two_part_diff_measured", &two_part_diff_measured_vec);
        event_tree_vecs->Branch("four_part_diff_truth", &four_part_diff_truth_vec);
        event_tree_vecs->Branch("four_part_diff_measured", &four_part_diff_measured_vec);

        TTree * event_tree_unfolded = new TTree("event_tree_unfolded", "event_tree_unfolded");
        int iterunfolded = 0;
        event_tree_unfolded->Branch("n_pt_bins", &n_pt_bins, "n_pt_bins/I");
        event_tree_unfolded->Branch("pt_bins", &pt_bins);
        event_tree_unfolded->Branch("iterunfolded", &iterunfolded, "iterunfolded/I");
        event_tree_unfolded->Branch("two_part_diff_unfolded", &two_part_diff_unfolded_vec);
        event_tree_unfolded->Branch("four_part_diff_unfolded", &four_part_diff_unfolded_vec);


        // declare histograms
        TH1D * h1_jet_pt_truth_all = new TH1D("h1_jet_pt_truth_all", "h1_jet_pt_truth_all", pt_bins.size()-1, pt_bins.data());
        TH1D * h1_jet_pt_fake = new TH1D("h1_jet_pt_fake", "h1_jet_pt_fake", pt_bins.size()-1, pt_bins.data());
        TH1D * h1_jet_pt_measured = new TH1D("h1_jet_pt_measured", "h1_jet_pt_measured", pt_bins.size()-1, pt_bins.data());

        // used for unfolding
        TH2D * h2_two_part_diff_truth = new TH2D("h2_two_part_diff_truth", "h2_two_part_diff_truth", pt_bins.size()-1, pt_bins.data() , two_part_diff_bins.size()-1, two_part_diff_bins.data());
        TH2D * h2_two_part_diff_measured = new TH2D("h2_two_part_diff_measured", "h2_two_part_diff_measured", pt_bins.size()-1, pt_bins.data() , two_part_diff_bins.size()-1, two_part_diff_bins.data());
        TH2D * h2_two_part_diff_reco_matched = new TH2D("h2_two_part_diff_reco_matched", "h2_two_part_diff_reco_matched", pt_bins.size()-1, pt_bins.data() , two_part_diff_bins.size()-1, two_part_diff_bins.data());
        TH2D * h2_two_part_diff_fake = new TH2D("h2_two_part_diff_fake", "h2_two_part_diff_fake", pt_bins.size()-1, pt_bins.data() , two_part_diff_bins.size()-1, two_part_diff_bins.data());
        
        TH2D * h2_four_part_diff_truth = new TH2D("h2_four_part_diff_truth", "h2_four_part_diff_truth", pt_bins.size()-1, pt_bins.data() , four_part_diff_bins.size()-1, four_part_diff_bins.data());
        TH2D * h2_four_part_diff_measured = new TH2D("h2_four_part_diff_measured", "h2_four_part_diff_measured", pt_bins.size()-1, pt_bins.data() , four_part_diff_bins.size()-1, four_part_diff_bins.data());
        TH2D * h2_four_part_diff_reco_matched = new TH2D("h2_four_part_diff_reco_matched", "h2_four_part_diff_reco_matched", pt_bins.size()-1, pt_bins.data() , four_part_diff_bins.size()-1, four_part_diff_bins.data());
        TH2D * h2_four_part_diff_fake = new TH2D("h2_four_part_diff_fake", "h2_four_part_diff_fake", pt_bins.size()-1, pt_bins.data() , four_part_diff_bins.size()-1, four_part_diff_bins.data());

        TH2D * h2_jet_v2_truth = new TH2D("h2_jet_v2_truth", "h2_jet_v2_truth", pt_bins.size()-1, pt_bins.data(), jet_v2_truth_bins.size()-1, jet_v2_truth_bins.data());
        TH2D * h2_jet_v3_truth = new TH2D("h2_jet_v3_truth", "h2_jet_v3_truth", pt_bins.size()-1, pt_bins.data(), jet_v3_truth_bins.size()-1, jet_v3_truth_bins.data());
        TH2D * h2_jet_v4_truth = new TH2D("h2_jet_v4_truth", "h2_jet_v4_truth", pt_bins.size()-1, pt_bins.data(), jet_v4_truth_bins.size()-1, jet_v4_truth_bins.data());
        
       

        int n_entries = event_tree->GetEntries();
        int n_entries_for_response = int(n_entries/3);
        double w_2 = 0;
        double w_4 = 0;
        // fill the histograms
        for (int i = 0; i < n_entries; i++)
        {
            
            
            event_tree->GetEntry(i);
            if(i ==0 )
            {
                w_2 = CorrFunctions::TwoPartDiffWeight(n_forward_particles, 1);
                w_4 = CorrFunctions::FourPartDiffWeight(n_forward_particles, 1);
            }

            h2_jet_v2_truth->Fill(jet_pt_truth, jet_v2_truth);
            h2_jet_v3_truth->Fill(jet_pt_truth, jet_v3_truth);
            h2_jet_v4_truth->Fill(jet_pt_truth, jet_v4_truth);


            // fill histograms
            if(jet_pt_truth > 0 && jet_pt_reco > -998)
            { // matched jets
                h1_jet_pt_truth_all->Fill(jet_pt_truth, weight);
                // h1_jet_pt_reco_matched->Fill(jet_pt_reco, weight);
                h1_jet_pt_measured->Fill(jet_pt_reco, weight);

                if(i > n_entries_for_response)
                {
                    // h2_two_part_diff_truth->Fill(jet_pt_truth, two_part_diff_truth, weight);
                    // h2_four_part_diff_truth->Fill(jet_pt_truth, four_part_diff_truth, weight);
                    // h2_two_part_diff_reco_matched->Fill(jet_pt_reco, two_part_diff_reco, weight);
                    // h2_four_part_diff_reco_matched->Fill(jet_pt_reco, four_part_diff_reco, weight);
                    // h2_two_part_diff_measured->Fill(jet_pt_reco, two_part_diff_reco, weight);
                    // h2_four_part_diff_measured->Fill(jet_pt_reco, four_part_diff_reco, weight);   
                    h2_two_part_diff_truth->Fill(jet_pt_truth, two_part_diff_truth, w_2);
                    h2_four_part_diff_truth->Fill(jet_pt_truth, four_part_diff_truth, w_4);
                    h2_two_part_diff_measured->Fill(jet_pt_reco, two_part_diff_reco, w_2);
                    h2_four_part_diff_measured->Fill(jet_pt_reco, four_part_diff_reco, w_4);
                }

            }
           
        
        }

        w_2 = 1.0;
        w_4 = 1.0;

 
        // write the histograms
        fout->cd();
        h1_jet_pt_truth_all->Write();
        // h1_jet_pt_reco_matched->Write();
        h1_jet_pt_fake->Write();
        h1_jet_pt_measured->Write();
        h2_jet_v2_truth->Write();
        h2_jet_v3_truth->Write();
        h2_jet_v4_truth->Write();

        TH1D * h1_jet_v2_truth = Make1D(h2_jet_v2_truth, "h1_jet_v2_truth");
        TH1D * h1_jet_v3_truth = Make1D(h2_jet_v3_truth, "h1_jet_v3_truth");
        TH1D * h1_jet_v4_truth = Make1D(h2_jet_v4_truth, "h1_jet_v4_truth");

        TH1D * h1_two_part_diff_truth = Make1D(h2_two_part_diff_truth, "h1_two_part_diff_truth", 1.0/w_2);
        TH1D * h1_two_part_diff_measured = Make1D(h2_two_part_diff_measured, "h1_two_part_diff_measured", 1.0/w_2);
        // TH1D * h1_two_part_diff_reco_matched = Make1D(h2_two_part_diff_reco_matched, "h1_two_part_diff_reco_matched", 1.0/w_2);
        TH1D * h1_two_part_diff_fake = Make1D(h2_two_part_diff_fake, "h1_two_part_diff_fake", 1.0/w_2);

        TH1D * h1_four_part_diff_truth = Make1D(h2_four_part_diff_truth, "h1_four_part_diff_truth", 1.0/w_4);
        TH1D * h1_four_part_diff_measured = Make1D(h2_four_part_diff_measured, "h1_four_part_diff_measured", 1.0/w_4);
        // TH1D * h1_four_part_diff_reco_matched = Make1D(h2_four_part_diff_reco_matched, "h1_four_part_diff_reco_matched", 1.0/w_4);
        TH1D * h1_four_part_diff_fake = Make1D(h2_four_part_diff_fake, "h1_four_part_diff_fake", 1.0/w_4);

        // convert to vectors
        two_part_diff_truth_vec = ConvertToVec(h2_two_part_diff_truth, 1.0/w_2);
        two_part_diff_measured_vec = ConvertToVec(h2_two_part_diff_measured, 1.0/w_2);
        four_part_diff_truth_vec = ConvertToVec(h2_four_part_diff_truth, 1.0/w_4);
        four_part_diff_measured_vec = ConvertToVec(h2_four_part_diff_measured, 1.0/w_4);
        event_tree_vecs->Fill();

        h2_two_part_diff_truth->Write();
        h2_two_part_diff_measured->Write();
        // h2_two_part_diff_reco_matched->Write();
        // h2_two_part_diff_fake->Write();
        h2_four_part_diff_truth->Write();
        h2_four_part_diff_measured->Write();
        // h2_four_part_diff_reco_matched->Write();
        // h2_four_part_diff_fake->Write();
        

        h1_jet_v2_truth->Write();
        h1_jet_v3_truth->Write();
        h1_jet_v4_truth->Write();
        h1_two_part_diff_truth->Write();
        h1_two_part_diff_measured->Write(); 
        // h1_two_part_diff_reco_matched->Write();
        // h1_two_part_diff_fake->Write();
        h1_four_part_diff_truth->Write();
        h1_four_part_diff_measured->Write();
        // h1_four_part_diff_reco_matched->Write();
        // h1_four_part_diff_fake->Write();

        event_tree_vecs->Write();
        event_tree_unfolded->Write();
        event_tree_out->Write();
        fout->Close();        
        fin->Close();
        
    }

    std::cout << "Done" << std::endl;
    return 0;

}


std::vector<double> DiffCorrFromHistos::ConvertToVec(TH2D * h2,  double w)
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

TH1D * DiffCorrFromHistos::Make1D(TH2D * h2, std::string name, double w)
{
    TProfile * p = h2->ProfileX();
    TH1D * h1 = (TH1D*)p->ProjectionX(name.c_str());
    if (w != 1.0)
    {
        h1->Scale(w);
    }   
    return h1;
}

std::string DiffCorrFromHistos::GetOutputFileName(std::string inputfile)
{
    // subtract the multiplicity
    TString inputfile_base = inputfile.c_str();
    inputfile_base.ReplaceAll(".root", "");
    // remove the path (everything before the last /)
    inputfile_base = inputfile_base(inputfile_base.Last('/')+1, inputfile_base.Length()-inputfile_base.Last('/')+1);
    TString outputfile = Form("%s/%s_fromhist.root", m_output_dir.c_str(), inputfile_base.Data());
    return outputfile.Data();
}

