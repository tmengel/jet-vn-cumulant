#ifndef DIFFCORRUNFOLDER_H
#define DIFFCORRUNFOLDER_H

#include "CorrFunctions.h"

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

#include <RooUnfoldResponse.h>
#include <RooUnfoldBayes.h>
#include <RooUnfoldErrors.h>

class DiffCorrUnfolder
{
    public:

        DiffCorrUnfolder(const std::string &input_dir) 
            : m_input(input_dir)
        {
        }
        ~DiffCorrUnfolder() {}

        void OutputDir(const std::string &name) { m_output_dir = name; }
        std::string OutputDir() { return m_output_dir; }

        int Run() { return Calculate(); }
        
        void MaxIter(int max_iter) { m_max_iter = max_iter; }
        int MaxIter() { return m_max_iter; }

        void MinIter(int min_iter) { m_min_iter = min_iter; }
        int MinIter() { return m_min_iter; }

        void SetCrossRotation(bool cross_rotation) { m_cross_rotation = cross_rotation; }
        bool GetCrossRotation() { return m_cross_rotation; }

        void SetInputFiles(const std::vector<std::string> &input_files) { m_input_files = input_files; }

    private:

        std::string m_input{""};
        std::vector<std::string> m_input_files{};
        std::string m_output_dir {""};
        

        int m_max_iter{6};
        int m_min_iter{3};

        bool m_cross_rotation{false};

        int Calculate();
        int CalculateCrossRotation();
        
        std::string GetOutputFileName(std::string inputfile);
        std::string GetOutputFileNameAlt(std::string inputfile, std::string res_file);
        static TH1D * Make1D(TH2D * h2, std::string name, double w=1.0);
        static std::vector<double> ConvertToVec(TH2D * h2, double w=1.0);

        static std::vector<std::string> GetFilesFromDir(const std::string &dir)
        {
            std::vector<std::string> input_files{};

            if(gSystem->AccessPathName(dir.c_str()))
            {
                std::cerr << "Error: directory " << dir << " does not exist" << std::endl;
                exit(1);
            }

            bool has_trailing_slash = dir.back() == '/';
            std::string ls_dir = dir;
            if(!has_trailing_slash) ls_dir += "/";

            std::string cmd = "ls " + ls_dir+ "*.root";
            FILE *fp = popen(cmd.c_str(), "r");
            if(!fp)
            {
                std::cerr << "Cannot open directory " << dir << std::endl;
                exit(1);
            }

            // read the output of the command
            char path[1024];
            while(fgets(path, sizeof(path), fp) != NULL)
            {
                std::string file_path = path;
                file_path.erase(std::remove(file_path.begin(), file_path.end(), '\n'), file_path.end());
                input_files.push_back(file_path);
            }

            pclose(fp);

            if(input_files.size() == 0)
            {
                std::cerr << "No files found in directory " << dir << std::endl;
                exit(1);
            }

            return input_files;
        }


};

int DiffCorrUnfolder::Calculate()
{
    if(m_cross_rotation)
    {
        return CalculateCrossRotation();
    }

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
        m_input_files = GetFilesFromDir(m_input);
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
        TH1D * h1_jet_pt_reco_matched = new TH1D("h1_jet_pt_reco_matched", "h1_jet_pt_reco_matched", pt_bins.size()-1, pt_bins.data());
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
        
        const int n_iterations = m_max_iter - m_min_iter + 1;
        TH2D * h2_unfolded_two_part_unfolded[n_iterations];
        TH2D * h2_unfolded_four_part_unfolded[n_iterations];

        RooUnfold::ErrorTreatment errorTreatment = RooUnfold::kCovariance;
        RooUnfoldResponse *response_two_part = new RooUnfoldResponse("response_two_part", "");
        RooUnfoldResponse *response_four_part = new RooUnfoldResponse("response_four_part", "");
        response_two_part->Setup(h2_two_part_diff_measured, h2_two_part_diff_truth);
        response_four_part->Setup(h2_four_part_diff_measured, h2_four_part_diff_truth);

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
                h1_jet_pt_truth_all->Fill(jet_pt_truth);
                h1_jet_pt_reco_matched->Fill(jet_pt_reco);
                h1_jet_pt_measured->Fill(jet_pt_reco);

                if(i < n_entries_for_response)
                {
                    // sample theorectical values
                    
                    response_two_part->Fill(jet_pt_reco, two_part_diff_reco, jet_pt_truth, two_part_diff_truth, 1.0);
                    response_four_part->Fill(jet_pt_reco, four_part_diff_reco, jet_pt_truth, four_part_diff_truth, 1.0);
                }
                else 
                {
                    h2_two_part_diff_truth->Fill(jet_pt_truth, two_part_diff_truth);
                    h2_four_part_diff_truth->Fill(jet_pt_truth, four_part_diff_truth);
                    h2_two_part_diff_reco_matched->Fill(jet_pt_reco, two_part_diff_reco);
                    h2_four_part_diff_reco_matched->Fill(jet_pt_reco, four_part_diff_reco);
                    h2_two_part_diff_measured->Fill(jet_pt_reco, two_part_diff_reco);
                    h2_four_part_diff_measured->Fill(jet_pt_reco, four_part_diff_reco);   
                }

            }
            // else if(jet_pt_truth > 0 && jet_pt_reco < -50)
            // { // missed jets

            //     h1_jet_pt_truth_all->Fill(jet_pt_truth);
            //     if(i < n_entries_for_response)
            //     {
            //         response_two_part->Miss(jet_pt_truth, two_part_diff_truth, 1.0);
            //         response_four_part->Miss(jet_pt_truth, four_part_diff_truth, 1.0);
            //     }
            //     else 
            //     {
            //         h2_two_part_diff_truth->Fill(jet_pt_truth, two_part_diff_truth);
            //         h2_four_part_diff_truth->Fill(jet_pt_truth, four_part_diff_truth);
            //     }
            // }
           
            // // unmatched jets       
            // for (int j = 0; j < jet_pt_unmatched->size(); j++)
            // {
            //     h1_jet_pt_fake->Fill(jet_pt_unmatched->at(j));
            //     h1_jet_pt_measured->Fill(jet_pt_unmatched->at(j));

            //     if(i < n_entries_for_response)
            //     {
            //         response_two_part->Fake(jet_pt_unmatched->at(j), two_part_diff_unmatched->at(j), 1.0);
            //         response_four_part->Fake(jet_pt_unmatched->at(j), four_part_diff_unmatched->at(j), 1.0);
            //     }
            //     else 
            //     {
            //         h2_two_part_diff_fake->Fill(jet_pt_unmatched->at(j), two_part_diff_unmatched->at(j));
            //         h2_four_part_diff_fake->Fill(jet_pt_unmatched->at(j), four_part_diff_unmatched->at(j));
            //         h2_two_part_diff_measured->Fill(jet_pt_unmatched->at(j), two_part_diff_unmatched->at(j));
            //         h2_four_part_diff_measured->Fill(jet_pt_unmatched->at(j), four_part_diff_unmatched->at(j));
            //     }
               
            // } 
        
        
        }

        // w_2 = 1.0;
        // w_4 = 1.0;
        // unfold the measured cumulants
        for (int i = 0; i < n_iterations; i++)
        {
            //clear vectors
            two_part_diff_unfolded_vec.clear();
            four_part_diff_unfolded_vec.clear();


            RooUnfoldBayes unfold_two_part(response_two_part, h2_two_part_diff_measured, i+m_min_iter);
            RooUnfoldBayes unfold_four_part(response_four_part, h2_four_part_diff_measured, i+m_min_iter);


            TH2D * h2_unfolded_two_part_tmp = (TH2D*)unfold_two_part.Hunfold(errorTreatment);
            TH2D * h2_unfolded_four_part_tmp = (TH2D*)unfold_four_part.Hunfold(errorTreatment);

            // convert to vectors
            two_part_diff_unfolded_vec = ConvertToVec(h2_unfolded_two_part_tmp, 1.0/w_2);
            four_part_diff_unfolded_vec = ConvertToVec(h2_unfolded_four_part_tmp, 1.0/w_4);
            iterunfolded = i+m_min_iter;
            event_tree_unfolded->Fill();

            h2_unfolded_two_part_unfolded[i] = (TH2D*)h2_unfolded_two_part_tmp->Clone(Form("h2_unfolded_two_part_unfolded_%d", i+m_min_iter));
            h2_unfolded_four_part_unfolded[i] = (TH2D*)h2_unfolded_four_part_tmp->Clone(Form("h2_unfolded_four_part_unfolded_%d", i+m_min_iter));
            

            TH1D * h1_unfolded_two_part_tmp = Make1D(h2_unfolded_two_part_tmp, Form("h1_unfolded_two_part_%d", i+m_min_iter), 1.0/w_2);
            TH1D * h1_unfolded_four_part_tmp = Make1D(h2_unfolded_four_part_tmp, Form("h1_unfolded_four_part_%d", i+m_min_iter), 1.0/w_4);

            fout->cd();
            h1_unfolded_two_part_tmp->Write();
            h1_unfolded_four_part_tmp->Write(); 
        }

        // write the histograms
        fout->cd();
        h1_jet_pt_truth_all->Write();
        h1_jet_pt_reco_matched->Write();
        h1_jet_pt_fake->Write();
        h1_jet_pt_measured->Write();

        TH1D * h1_jet_v2_truth = Make1D(h2_jet_v2_truth, "h1_jet_v2_truth");
        TH1D * h1_jet_v3_truth = Make1D(h2_jet_v3_truth, "h1_jet_v3_truth");
        TH1D * h1_jet_v4_truth = Make1D(h2_jet_v4_truth, "h1_jet_v4_truth");

        TH1D * h1_two_part_diff_truth = Make1D(h2_two_part_diff_truth, "h1_two_part_diff_truth", 1.0/w_2);
        TH1D * h1_two_part_diff_measured = Make1D(h2_two_part_diff_measured, "h1_two_part_diff_measured", 1.0/w_2);
        TH1D * h1_two_part_diff_reco_matched = Make1D(h2_two_part_diff_reco_matched, "h1_two_part_diff_reco_matched", 1.0/w_2);
        TH1D * h1_two_part_diff_fake = Make1D(h2_two_part_diff_fake, "h1_two_part_diff_fake", 1.0/w_2);

       

        TH1D * h1_four_part_diff_truth = Make1D(h2_four_part_diff_truth, "h1_four_part_diff_truth", 1.0/w_4);
        TH1D * h1_four_part_diff_measured = Make1D(h2_four_part_diff_measured, "h1_four_part_diff_measured", 1.0/w_4);
        TH1D * h1_four_part_diff_reco_matched = Make1D(h2_four_part_diff_reco_matched, "h1_four_part_diff_reco_matched", 1.0/w_4);
        TH1D * h1_four_part_diff_fake = Make1D(h2_four_part_diff_fake, "h1_four_part_diff_fake", 1.0/w_4);

        // convert to vectors
        two_part_diff_truth_vec = ConvertToVec(h2_two_part_diff_truth, 1.0/w_2);
        two_part_diff_measured_vec = ConvertToVec(h2_two_part_diff_measured, 1.0/w_2);
        four_part_diff_truth_vec = ConvertToVec(h2_four_part_diff_truth, 1.0/w_4);
        four_part_diff_measured_vec = ConvertToVec(h2_four_part_diff_measured, 1.0/w_4);
        event_tree_vecs->Fill();

        h2_two_part_diff_truth->Write();
        h2_two_part_diff_measured->Write();
        h2_two_part_diff_reco_matched->Write();
        h2_two_part_diff_fake->Write();
        h2_four_part_diff_truth->Write();
        h2_four_part_diff_measured->Write();
        h2_four_part_diff_reco_matched->Write();
        h2_four_part_diff_fake->Write();
        

        h1_jet_v2_truth->Write();
        h1_jet_v3_truth->Write();
        h1_jet_v4_truth->Write();
        h1_two_part_diff_truth->Write();
        h1_two_part_diff_measured->Write(); 
        h1_two_part_diff_reco_matched->Write();
        h1_two_part_diff_fake->Write();
        h1_four_part_diff_truth->Write();
        h1_four_part_diff_measured->Write();
        h1_four_part_diff_reco_matched->Write();
        h1_four_part_diff_fake->Write();

        event_tree_vecs->Write();
        event_tree_unfolded->Write();
        event_tree_out->Write();
        fout->Close();        
        fin->Close();
        
    }

    std::cout << "Done" << std::endl;
    return 0;

}

int DiffCorrUnfolder::CalculateCrossRotation()
{

    std::cout << "Calculating cross rotation" << std::endl;

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

    for (auto &input_file : m_input_files)
    { 

        for (auto &response_file : m_input_files)
        {
            if(response_file == input_file) continue;
            
            
            // calculate correlations
            std::string cum_file = GetOutputFileNameAlt(input_file, response_file);
            std::cout << "Creating corr file " << cum_file << std::endl;
    

            // this is the file that will be unfolded with all the other files
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
            bin_tree->SetBranchAddress("jet_v2_truth_bins", &in_jet_v2_truth_bins);
            bin_tree->SetBranchAddress("jet_v3_truth_bins", &in_jet_v3_truth_bins);
            bin_tree->SetBranchAddress("jet_v4_truth_bins", &in_jet_v4_truth_bins);
            bin_tree->SetBranchAddress("two_diff_bins", &in_two_part_diff_bins);
            bin_tree->SetBranchAddress("four_diff_bins", &in_four_part_diff_bins);

            bin_tree->GetEntry(0);
            std::vector<double> jet_v2_truth_bins = *in_jet_v2_truth_bins;
            std::vector<double> jet_v3_truth_bins = *in_jet_v3_truth_bins;
            std::vector<double> jet_v4_truth_bins = *in_jet_v4_truth_bins;
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
            TH1D * h1_jet_pt_reco_matched = new TH1D("h1_jet_pt_reco_matched", "h1_jet_pt_reco_matched", pt_bins.size()-1, pt_bins.data());
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
            
            const int n_iterations = m_max_iter - m_min_iter + 1;
            TH2D * h2_unfolded_two_part_unfolded[n_iterations];
            TH2D * h2_unfolded_four_part_unfolded[n_iterations];

            
            RooUnfold::ErrorTreatment errorTreatment = RooUnfold::kCovariance;
            RooUnfoldResponse *response_two_part = new RooUnfoldResponse("response_two_part", "");
            RooUnfoldResponse *response_four_part = new RooUnfoldResponse("response_four_part", "");
            response_two_part->Setup(h2_two_part_diff_measured, h2_two_part_diff_truth);
            response_four_part->Setup(h2_four_part_diff_measured, h2_four_part_diff_truth);
            
            double w_2 = 0;
            double w_4 = 0;
            int n_entries = event_tree->GetEntries();
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

                if(jet_pt_truth > 0 && jet_pt_reco > -998)
                { // matched jets
                    h1_jet_pt_truth_all->Fill(jet_pt_truth);
                    h1_jet_pt_reco_matched->Fill(jet_pt_reco);
                    h1_jet_pt_measured->Fill(jet_pt_reco);
                    h2_two_part_diff_truth->Fill(jet_pt_truth, two_part_diff_truth);
                    h2_four_part_diff_truth->Fill(jet_pt_truth, four_part_diff_truth);
                    h2_two_part_diff_reco_matched->Fill(jet_pt_reco, two_part_diff_reco);
                    h2_four_part_diff_reco_matched->Fill(jet_pt_reco, four_part_diff_reco);
                    h2_two_part_diff_measured->Fill(jet_pt_reco, two_part_diff_reco);
                    h2_four_part_diff_measured->Fill(jet_pt_reco, four_part_diff_reco);   
                }
            
            }

            // this is the file that will be unfolded with all the other files
            TFile * resfin = new TFile(response_file.c_str(), "READ");
            if(!resfin->IsOpen())
            {
                std::cout << "Error: could not open input file " << response_file << std::endl;
                exit(1);
            }

            TTree * res_event_tree = (TTree*)resfin->Get("tree");
            if(!res_event_tree)
            {
                std::cout << "Error: could not find tree in input file " << response_file << std::endl;
                exit(1);
            }


            // input variables
            int res_event_id = 0;
            int res_n_forward_particles = 0;        
            double res_jet_pt_truth = 0.0;
            double res_two_part_diff_truth = 0.0;
            double res_four_part_diff_truth = 0.0;
            double res_jet_pt_reco = 0.0;
            double res_two_part_diff_reco = 0.0;
            double res_four_part_diff_reco = 0.0;
            std::vector<double> * res_jet_pt_unmatched = 0;
            std::vector<double> * res_two_part_diff_unmatched = 0;
            std::vector<double> * res_four_part_diff_unmatched = 0;

            res_event_tree->SetBranchAddress("event_id", &res_event_id);
            // res_event_tree->SetBranchAddress("jet_v2_truth", &res_jet_v2_truth);
            // res_event_tree->SetBranchAddress("jet_v3_truth", &res_jet_v3_truth);
            // res_event_tree->SetBranchAddress("jet_v4_truth", &res_jet_v4_truth);
            res_event_tree->SetBranchAddress("n_forward_particles", &res_n_forward_particles);
            res_event_tree->SetBranchAddress("jet_pt_truth", &res_jet_pt_truth);
            res_event_tree->SetBranchAddress("two_part_diff_truth", &res_two_part_diff_truth);
            res_event_tree->SetBranchAddress("four_part_diff_truth", &res_four_part_diff_truth);
            res_event_tree->SetBranchAddress("jet_pt_reco", &res_jet_pt_reco);
            res_event_tree->SetBranchAddress("two_part_diff_reco", &res_two_part_diff_reco);
            res_event_tree->SetBranchAddress("four_part_diff_reco", &res_four_part_diff_reco);
            res_event_tree->SetBranchAddress("jet_pt_unmatched", &res_jet_pt_unmatched);
            res_event_tree->SetBranchAddress("two_part_diff_unmatched", &res_two_part_diff_unmatched);
            res_event_tree->SetBranchAddress("four_part_diff_unmatched", &res_four_part_diff_unmatched);

            int n_entries_res = res_event_tree->GetEntries();
            int n_entries_for_response = int(n_entries_res/3);

            w_2 = 0;
            w_4 = 0;
            for (int i = 0; i < n_entries_for_response; i++)
            {


                res_event_tree->GetEntry(i);
                if(i ==0 )
                {
                    w_2 = CorrFunctions::TwoPartDiffWeight(n_forward_particles, 1);
                    w_4 = CorrFunctions::FourPartDiffWeight(n_forward_particles, 1);
                }

                // fill histograms
                if(res_jet_pt_truth > 0 && res_jet_pt_reco > -998)
                { 

                    response_two_part->Fill(res_jet_pt_reco, res_two_part_diff_reco, res_jet_pt_truth, res_two_part_diff_truth, 1.0);
                    response_four_part->Fill(res_jet_pt_reco, res_four_part_diff_reco, res_jet_pt_truth, res_four_part_diff_truth, 1.0);
                }
                // else if(jet_pt_truth > 0 && jet_pt_reco < -50)
                // { // missed jets

                //     h1_jet_pt_truth_all->Fill(jet_pt_truth);
                //     if(i < n_entries_for_response)
                //     {
                //         response_two_part->Miss(jet_pt_truth, two_part_diff_truth, 1.0);
                //         response_four_part->Miss(jet_pt_truth, four_part_diff_truth, 1.0);
                //     }
                //     else 
                //     {
                //         h2_two_part_diff_truth->Fill(jet_pt_truth, two_part_diff_truth);
                //         h2_four_part_diff_truth->Fill(jet_pt_truth, four_part_diff_truth);
                //     }
                // }

                // unmatched jets       
                // for (int j = 0; j < jet_pt_unmatched->size(); j++)
                // {
                //     h1_jet_pt_fake->Fill(jet_pt_unmatched->at(j));
                //     h1_jet_pt_measured->Fill(jet_pt_unmatched->at(j));

                //     if(i < n_entries_for_response)
                //     {
                //         response_two_part->Fake(jet_pt_unmatched->at(j), two_part_diff_unmatched->at(j), 1.0);
                //         response_four_part->Fake(jet_pt_unmatched->at(j), four_part_diff_unmatched->at(j), 1.0);
                //     }
                //     else 
                //     {
                //         h2_two_part_diff_fake->Fill(jet_pt_unmatched->at(j), two_part_diff_unmatched->at(j));
                //         h2_four_part_diff_fake->Fill(jet_pt_unmatched->at(j), four_part_diff_unmatched->at(j));
                //         h2_two_part_diff_measured->Fill(jet_pt_unmatched->at(j), two_part_diff_unmatched->at(j));
                //         h2_four_part_diff_measured->Fill(jet_pt_unmatched->at(j), four_part_diff_unmatched->at(j));
                //     }
                
                // } 
            }

            // unfold the measured cumulants
            for (int i = 0; i < n_iterations; i++)
            {
                //clear vectors
                two_part_diff_unfolded_vec.clear();
                four_part_diff_unfolded_vec.clear();


                RooUnfoldBayes unfold_two_part(response_two_part, h2_two_part_diff_measured, i+m_min_iter);
                RooUnfoldBayes unfold_four_part(response_four_part, h2_four_part_diff_measured, i+m_min_iter);


                TH2D * h2_unfolded_two_part_tmp = (TH2D*)unfold_two_part.Hunfold(errorTreatment);
                TH2D * h2_unfolded_four_part_tmp = (TH2D*)unfold_four_part.Hunfold(errorTreatment);

                // convert to vectors
                two_part_diff_unfolded_vec = ConvertToVec(h2_unfolded_two_part_tmp, 1.0/w_2);
                four_part_diff_unfolded_vec = ConvertToVec(h2_unfolded_four_part_tmp, 1.0/w_4);
                iterunfolded = i+m_min_iter;
                event_tree_unfolded->Fill();

                h2_unfolded_two_part_unfolded[i] = (TH2D*)h2_unfolded_two_part_tmp->Clone(Form("h2_unfolded_two_part_unfolded_%d", i+m_min_iter));
                h2_unfolded_four_part_unfolded[i] = (TH2D*)h2_unfolded_four_part_tmp->Clone(Form("h2_unfolded_four_part_unfolded_%d", i+m_min_iter));
                

                TH1D * h1_unfolded_two_part_tmp = Make1D(h2_unfolded_two_part_tmp, Form("h1_unfolded_two_part_%d", i+m_min_iter), 1.0/w_2);
                TH1D * h1_unfolded_four_part_tmp = Make1D(h2_unfolded_four_part_tmp, Form("h1_unfolded_four_part_%d", i+m_min_iter), 1.0/w_4);

                fout->cd();
                h1_unfolded_two_part_tmp->Write();
                h1_unfolded_four_part_tmp->Write(); 
            }

            // write the histograms
            fout->cd();
            h1_jet_pt_truth_all->Write();
            h1_jet_pt_reco_matched->Write();
            h1_jet_pt_fake->Write();
            h1_jet_pt_measured->Write();

            TH1D * h1_two_part_diff_truth = Make1D(h2_two_part_diff_truth, "h1_two_part_diff_truth", 1.0/w_2);
            TH1D * h1_two_part_diff_measured = Make1D(h2_two_part_diff_measured, "h1_two_part_diff_measured", 1.0/w_2);
            TH1D * h1_two_part_diff_reco_matched = Make1D(h2_two_part_diff_reco_matched, "h1_two_part_diff_reco_matched", 1.0/w_2);
            TH1D * h1_two_part_diff_fake = Make1D(h2_two_part_diff_fake, "h1_two_part_diff_fake", 1.0/w_2);

            TH1D * h1_four_part_diff_truth = Make1D(h2_four_part_diff_truth, "h1_four_part_diff_truth", 1.0/w_4);
            TH1D * h1_four_part_diff_measured = Make1D(h2_four_part_diff_measured, "h1_four_part_diff_measured", 1.0/w_4);
            TH1D * h1_four_part_diff_reco_matched = Make1D(h2_four_part_diff_reco_matched, "h1_four_part_diff_reco_matched", 1.0/w_4);
            TH1D * h1_four_part_diff_fake = Make1D(h2_four_part_diff_fake, "h1_four_part_diff_fake", 1.0/w_4);

            // convert to vectors
            two_part_diff_truth_vec = ConvertToVec(h2_two_part_diff_truth, 1.0/w_2);
            two_part_diff_measured_vec = ConvertToVec(h2_two_part_diff_measured, 1.0/w_2);
            four_part_diff_truth_vec = ConvertToVec(h2_four_part_diff_truth, 1.0/w_4);
            four_part_diff_measured_vec = ConvertToVec(h2_four_part_diff_measured, 1.0/w_4);
            event_tree_vecs->Fill();

            h2_two_part_diff_truth->Write();
            h2_two_part_diff_measured->Write();
            h2_two_part_diff_reco_matched->Write();
            h2_two_part_diff_fake->Write();
            h2_four_part_diff_truth->Write();
            h2_four_part_diff_measured->Write();
            h2_four_part_diff_reco_matched->Write();
            h2_four_part_diff_fake->Write();
            

            // h1_jet_v2_truth->Write();
            // h1_jet_v3_truth->Write();
            // h1_jet_v4_truth->Write();
            h1_two_part_diff_truth->Write();
            h1_two_part_diff_measured->Write(); 
            h1_two_part_diff_reco_matched->Write();
            h1_two_part_diff_fake->Write();
            h1_four_part_diff_truth->Write();
            h1_four_part_diff_measured->Write();
            h1_four_part_diff_reco_matched->Write();
            h1_four_part_diff_fake->Write();

            event_tree_vecs->Write();
            event_tree_unfolded->Write();
            event_tree_out->Write();
            fout->Close();        
            fin->Close();
            resfin->Close();
        }
    }
    return 0;
}

std::vector<double> DiffCorrUnfolder::ConvertToVec(TH2D * h2,  double w)
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

TH1D * DiffCorrUnfolder::Make1D(TH2D * h2, std::string name, double w)
{
    TProfile * p = h2->ProfileX();
    TH1D * h1 = (TH1D*)p->ProjectionX(name.c_str());
    if (w != 1.0)
    {
        h1->Scale(w);
    }   
    return h1;
}

std::string DiffCorrUnfolder::GetOutputFileName(std::string inputfile)
{
    // subtract the multiplicity
    TString inputfile_base = inputfile.c_str();
    inputfile_base.ReplaceAll(".root", "");
    // remove the path (everything before the last /)
    inputfile_base = inputfile_base(inputfile_base.Last('/')+1, inputfile_base.Length()-inputfile_base.Last('/')+1);
    TString outputfile = Form("%s/%s_unfolded.root", m_output_dir.c_str(), inputfile_base.Data());
    return outputfile.Data();
}

std::string DiffCorrUnfolder::GetOutputFileNameAlt(std::string inputfile, std::string res_file)
{
    // subtract the multiplicity
    TString inputfile_base = inputfile.c_str();
    inputfile_base.ReplaceAll(".root", "");
    // remove the path (everything before the last /)
    inputfile_base = inputfile_base(inputfile_base.Last('/')+1, inputfile_base.Length()-inputfile_base.Last('/')+1);

    // get sub dir after the last /rootfiles/
    TString rotation_type_input = inputfile.c_str();
    rotation_type_input.ReplaceAll("/lustre/isaac/scratch/tmengel/jet-vn-cumulant/ana/rootfiles/", "");
    rotation_type_input = rotation_type_input(0, rotation_type_input.First('/'));
    TString rotation_type_res = res_file.c_str();
    rotation_type_res.ReplaceAll("/lustre/isaac/scratch/tmengel/jet-vn-cumulant/ana/rootfiles/", "");
    rotation_type_res = rotation_type_res(0, rotation_type_res.First('/'));
    TString outputfile = Form("%s/%s_%s_unfolded_w_%s.root", m_output_dir.c_str(), inputfile_base.Data(), rotation_type_res.Data(), rotation_type_input.Data());
    return outputfile.Data();
}

#endif // DIFFCORRUNFOLDER_H