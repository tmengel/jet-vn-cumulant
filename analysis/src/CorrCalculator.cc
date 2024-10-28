#include "Utils.h"
#include "CorrCalculator.h"
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

int CorrCalculator::Calculate()
{

    if(m_pt_bins.size() == 0)
    {
        std::cerr << "Pt binning not set" << std::endl;
        exit(1);
    }

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
    }

    m_input_files.clear();
    // see if the input is a directory or a file
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

    std::cout << "Calculating correlations" << std::endl;
    for (int iharm =2; iharm < 5; iharm++)
    {
        CalculateHarm(iharm);
    }

    std::cout << "Done calculating correlations" << std::endl;


    return 0;

}

int CorrCalculator::CalculateHarm(int harmonic)
{   
    std::cout << "Calculating " << harmonic << "-order correlations" << std::endl;

    std::string output_dir = m_output_dir + "/Order" + std::to_string(harmonic);
    if(gSystem->AccessPathName(output_dir.c_str()))
    {
        std::cout << "Creating output directory " << output_dir << std::endl;
        gSystem->mkdir(output_dir.c_str(), kTRUE);
    }

    for (auto &input_file : m_input_files)
    {
        // calculate correlations
        std::string corr_file = GetOutputFileName(input_file, output_dir, harmonic);
        std::cout << "Creating corr file " << corr_file << std::endl;

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
    

        int event_id = 0;
        double weight = 1.0;
        int n_forward_particles = 0;
        double Qn_i = 0.0;
        double Qn_r = 0.0;
        double Q2n_i = 0.0;
        double Q2n_r = 0.0;

        // single truth jet variables
        double jet_pt_truth = 0.0;
        double jet_phi_truth = 0.0;

        double jet_v2_truth = 0.0;
        double jet_v3_truth = 0.0;
        double jet_v4_truth = 0.0;

        // single reco jet variables
        double jet_pt_reco = 0.0;
        double jet_phi_reco = 0.0;

        // unmatched reco jet variables
        int n_unmatched_reco_jets = 0;
        std::vector<double> * unmatched_reco_pt = 0;
        std::vector<double> * unmatched_reco_phi  = 0;

        // set up input branches
        event_tree->SetBranchAddress("event_id", &event_id);
        event_tree->SetBranchAddress("weight", &weight);
        event_tree->SetBranchAddress("n_forward_particles", &n_forward_particles);
        event_tree->SetBranchAddress(Form("Q%d_i", harmonic), &Qn_i);
        event_tree->SetBranchAddress(Form("Q%d_r", harmonic), &Qn_r);
        event_tree->SetBranchAddress(Form("Q%d_i", int(2*harmonic)), &Q2n_i);
        event_tree->SetBranchAddress(Form("Q%d_r", int(2*harmonic)), &Q2n_r);
        event_tree->SetBranchAddress("jet_pt_truth", &jet_pt_truth);
        event_tree->SetBranchAddress("jet_phi_truth", &jet_phi_truth);
        event_tree->SetBranchAddress("jet_v2_truth", &jet_v2_truth);
        event_tree->SetBranchAddress("jet_v3_truth", &jet_v3_truth);
        event_tree->SetBranchAddress("jet_v4_truth", &jet_v4_truth);
        event_tree->SetBranchAddress("jet_pt_reco", &jet_pt_reco);
        event_tree->SetBranchAddress("jet_phi_reco", &jet_phi_reco);
        event_tree->SetBranchAddress("n_unmatched_reco_jets", &n_unmatched_reco_jets);
        event_tree->SetBranchAddress("unmatched_reco_pt", &unmatched_reco_pt);
        event_tree->SetBranchAddress("unmatched_reco_phi", &unmatched_reco_phi);



        // create the output file
        TFile * fout = new TFile(corr_file.c_str(), "RECREATE");
        TTree * event_tree_corr = new TTree("tree", "tree");

        // output variables
        double two_part_diff_truth;
        double four_part_diff_truth;
        double two_part_diff_reco;
        double four_part_diff_reco;
        std::vector<double> two_part_diff_unmatched;
        std::vector<double> four_part_diff_unmatched;
        event_tree_corr->Branch("event_id", &event_id, "event_id/I");
        event_tree_corr->Branch("weight", &weight, "weight/D");
        event_tree_corr->Branch("jet_v2_truth", &jet_v2_truth, "jet_v2_truth/D");
        event_tree_corr->Branch("jet_v3_truth", &jet_v3_truth, "jet_v3_truth/D");
        event_tree_corr->Branch("jet_v4_truth", &jet_v4_truth, "jet_v4_truth/D");
        event_tree_corr->Branch("n_forward_particles", &n_forward_particles, "n_forward_particles/I");
        event_tree_corr->Branch("jet_pt_truth", &jet_pt_truth, "jet_pt_truth/D");
        event_tree_corr->Branch("two_part_diff_truth", &two_part_diff_truth, "two_part_diff_truth/D");
        event_tree_corr->Branch("four_part_diff_truth", &four_part_diff_truth, "four_part_diff_truth/D");
        event_tree_corr->Branch("jet_pt_reco", &jet_pt_reco, "jet_pt_reco/D");
        event_tree_corr->Branch("two_part_diff_reco", &two_part_diff_reco, "two_part_diff_reco/D");
        event_tree_corr->Branch("four_part_diff_reco", &four_part_diff_reco, "four_part_diff_reco/D");
        event_tree_corr->Branch("jet_pt_unmatched", &unmatched_reco_pt);
        event_tree_corr->Branch("two_part_diff_unmatched", &two_part_diff_unmatched);
        event_tree_corr->Branch("four_part_diff_unmatched", &four_part_diff_unmatched);

        TTree * ref_tree = new TTree("ref_tree", "ref_tree");
        double two_part_event_avg_ref, four_part_event_avg_ref;
        double cn_two_ref, cn_four_ref;
        double vn_two_ref, vn_four_ref;
        ref_tree->Branch("two_part_event_avg_ref", &two_part_event_avg_ref, "two_part_event_avg_ref/D");
        ref_tree->Branch("four_part_event_avg_ref", &four_part_event_avg_ref, "four_part_event_avg_ref/D");
        ref_tree->Branch("cn_two_ref", &cn_two_ref, "cn_two_ref/D");
        ref_tree->Branch("cn_four_ref", &cn_four_ref, "cn_four_ref/D");
        ref_tree->Branch("vn_two_ref", &vn_two_ref, "vn_two_ref/D");
        ref_tree->Branch("vn_four_ref", &vn_four_ref, "vn_four_ref/D");

   

        // keep track of min max for binning
        std::pair<double, double> jet_v2_truth_minmax = {9999.0, -9999.0};
        std::pair<double, double> jet_v3_truth_minmax = {9999.0, -9999.0};
        std::pair<double, double> jet_v4_truth_minmax = {9999.0, -9999.0};

        std::pair<double,double> two_truth_minmax = {9999.0, -9999.0};
        std::pair<double,double> four_truth_minmax = {9999.0, -9999.0};
        std::pair<double,double> two_reco_minmax = {9999.0, -9999.0};
        std::pair<double,double> four_reco_minmax = {9999.0, -9999.0};
        std::pair<double,double> two_unmatched_minmax = {9999.0, -9999.0};
        std::pair<double,double> four_unmatched_minmax = {9999.0, -9999.0};

        std::pair<double,double> two_minmax = {9999.0, -9999.0};
        std::pair<double,double> four_minmax = {9999.0, -9999.0};
        
        // loop over the events
        int n_entries = event_tree->GetEntries();

        double two_part_ref_sum = 0.0;
        double four_part_ref_sum = 0.0;
        double two_part_ref_weight_sum = 0.0;
        double four_part_ref_weight_sum = 0.0;


       

        for (int i = 0; i < n_entries; i++)
        {
            event_tree->GetEntry(i);
            
            // clear the vectors
            two_part_diff_unmatched.clear();
            four_part_diff_unmatched.clear();
            two_part_diff_truth = -999;
            four_part_diff_truth = -999;
            two_part_diff_reco = -999;
            four_part_diff_reco = -999;

            if (jet_v2_truth < jet_v2_truth_minmax.first) jet_v2_truth_minmax.first = jet_v2_truth;
            if (jet_v2_truth > jet_v2_truth_minmax.second) jet_v2_truth_minmax.second = jet_v2_truth;

            if (jet_v3_truth < jet_v3_truth_minmax.first) jet_v3_truth_minmax.first = jet_v3_truth;
            if (jet_v3_truth > jet_v3_truth_minmax.second) jet_v3_truth_minmax.second = jet_v3_truth;

            if (jet_v4_truth < jet_v4_truth_minmax.first) jet_v4_truth_minmax.first = jet_v4_truth;
            if (jet_v4_truth > jet_v4_truth_minmax.second) jet_v4_truth_minmax.second = jet_v4_truth;
            

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

            if(jet_pt_truth > 0)
            {
                TComplex pn = CorrFunctions::pvec(jet_phi_truth, harmonic);
                double w_two = CorrFunctions::TwoPartDiffWeight(n_forward_particles, 1);
                double w_four = CorrFunctions::FourPartDiffWeight(n_forward_particles, 1, 0);
                double two_diff = CorrFunctions::TwoPartDiff(Qn, pn, 1, n_forward_particles);
                double four_diff = CorrFunctions::FourPartDiff(Qn, Q2n, pn, 1, n_forward_particles);
                // two_part_diff_truth = two_diff*w_two;
                // four_part_diff_truth = four_diff*w_four;
                two_part_diff_truth= two_diff;
                four_part_diff_truth = four_diff;

                if(two_part_diff_truth < two_truth_minmax.first) two_truth_minmax.first = two_part_diff_truth;
                if(two_part_diff_truth > two_truth_minmax.second) two_truth_minmax.second = two_part_diff_truth;

                if(four_part_diff_truth < four_truth_minmax.first) four_truth_minmax.first = four_part_diff_truth;
                if(four_part_diff_truth > four_truth_minmax.second) four_truth_minmax.second = four_part_diff_truth;

                if(two_part_diff_truth < two_minmax.first) two_minmax.first = two_part_diff_truth;
                if(two_part_diff_truth > two_minmax.second) two_minmax.second = two_part_diff_truth;
                
                if(four_part_diff_truth < four_minmax.first) four_minmax.first = four_part_diff_truth;
                if(four_part_diff_truth > four_minmax.second) four_minmax.second = four_part_diff_truth;
            
            }
            if(jet_pt_reco > 0)
            {
                TComplex pn = CorrFunctions::pvec(jet_phi_reco, harmonic);
                double w_two = CorrFunctions::TwoPartDiffWeight(n_forward_particles, 1);
                double w_four = CorrFunctions::FourPartDiffWeight(n_forward_particles, 1, 0);
                double two_diff = CorrFunctions::TwoPartDiff(Qn, pn, 1, n_forward_particles);
                double four_diff = CorrFunctions::FourPartDiff(Qn, Q2n, pn, 1, n_forward_particles);
                // two_part_diff_reco = two_diff*w_two;
                // four_part_diff_reco = four_diff*w_four;
                two_part_diff_reco = two_diff;
                four_part_diff_reco = four_diff;


                if(two_part_diff_reco < two_reco_minmax.first) two_reco_minmax.first = two_part_diff_reco;
                if(two_part_diff_reco > two_reco_minmax.second) two_reco_minmax.second = two_part_diff_reco;

                if(four_part_diff_reco < four_reco_minmax.first) four_reco_minmax.first = four_part_diff_reco;
                if(four_part_diff_reco > four_reco_minmax.second) four_reco_minmax.second = four_part_diff_reco;

                if(two_part_diff_reco < two_minmax.first) two_minmax.first = two_part_diff_reco;
                if(two_part_diff_reco > two_minmax.second) two_minmax.second = two_part_diff_reco;

                if(four_part_diff_reco < four_minmax.first) four_minmax.first = four_part_diff_reco;
                if(four_part_diff_reco > four_minmax.second) four_minmax.second = four_part_diff_reco;

            }

            for (int j = 0; j < unmatched_reco_phi->size(); j++)
            {
                TComplex pn = CorrFunctions::pvec(unmatched_reco_phi->at(j), harmonic);
                
                double w_two = CorrFunctions::TwoPartDiffWeight(n_forward_particles, 1);
                double w_four = CorrFunctions::FourPartDiffWeight(n_forward_particles, 1, 0);
                double two_part_diff = CorrFunctions::TwoPartDiff(Qn, pn, 1, n_forward_particles);
                double four_part_diff = CorrFunctions::FourPartDiff(Qn, Q2n, pn, 1, n_forward_particles);

                // two_part_diff_unmatched.push_back(two_part_diff*w_two);
                // four_part_diff_unmatched.push_back(four_part_diff*w_four);
                two_part_diff_unmatched.push_back(two_part_diff);
                four_part_diff_unmatched.push_back(four_part_diff);

                if(two_part_diff_unmatched.back() < two_unmatched_minmax.first) two_unmatched_minmax.first = two_part_diff_unmatched.back();
                if(two_part_diff_unmatched.back() > two_unmatched_minmax.second) two_unmatched_minmax.second = two_part_diff_unmatched.back();

                if(four_part_diff_unmatched.back() < four_unmatched_minmax.first) four_unmatched_minmax.first = four_part_diff_unmatched.back();
                if(four_part_diff_unmatched.back() > four_unmatched_minmax.second) four_unmatched_minmax.second = four_part_diff_unmatched.back();
           
            }
           
           
            event_tree_corr->Fill();
        }

        two_part_event_avg_ref = two_part_ref_sum/two_part_ref_weight_sum;
        four_part_event_avg_ref = four_part_ref_sum/four_part_ref_weight_sum;
        cn_two_ref = two_part_event_avg_ref;
        cn_four_ref = four_part_event_avg_ref - 2*two_part_event_avg_ref*two_part_event_avg_ref;
        vn_two_ref = TMath::Sqrt(cn_two_ref);
        double neg_fou_ref = -1.0*cn_four_ref;
        vn_four_ref = TMath::Power(neg_fou_ref, 1.0/4.0);

        ref_tree->Fill();

        // bin tree for unfolding
        TTree * bin_tree = new TTree("bin_tree", "bin_tree");
        jet_v2_truth_minmax.first*=0.5;
        jet_v2_truth_minmax.second*=1.5;
        jet_v3_truth_minmax.first*=0.5;
        jet_v3_truth_minmax.second*=1.5;
        jet_v4_truth_minmax.first*=0.5;
        jet_v4_truth_minmax.second*=1.5;
        

         std::vector<double> jet_v2_truth_bins = GetBinning(jet_v2_truth_minmax, 20);
        std::vector<double> jet_v3_truth_bins = GetBinning(jet_v3_truth_minmax, 20);
        std::vector<double> jet_v4_truth_bins = GetBinning(jet_v4_truth_minmax, 20);
        std::vector<double> two_diff_truth_bins = GetCorrBinning(two_truth_minmax, 100);
        std::vector<double> four_diff_truth_bins = GetCorrBinning(four_truth_minmax, 100);
        std::vector<double> two_diff_reco_bins = GetCorrBinning(two_reco_minmax, 100);
        std::vector<double> four_diff_reco_bins = GetCorrBinning(four_reco_minmax, 100);
        std::vector<double> two_diff_unmatched_bins = GetBinning(two_unmatched_minmax, 100);
        std::vector<double> four_diff_unmatched_bins = GetBinning(four_unmatched_minmax, 100);
        std::vector<double> two_diff_bins = GetCorrBinning(two_minmax, 150);
        std::vector<double> four_diff_bins = GetCorrBinning(four_minmax, 250);
    
        bin_tree->Branch("pt_bins", &m_pt_bins);
        bin_tree->Branch("jet_v2_truth_bins", &jet_v2_truth_bins);
        bin_tree->Branch("jet_v3_truth_bins", &jet_v3_truth_bins);
        bin_tree->Branch("jet_v4_truth_bins", &jet_v4_truth_bins);
        bin_tree->Branch("two_diff_truth_bins", &two_diff_truth_bins);
        bin_tree->Branch("four_diff_truth_bins", &four_diff_truth_bins);
        bin_tree->Branch("two_diff_reco_bins", &two_diff_reco_bins);
        bin_tree->Branch("four_diff_reco_bins", &four_diff_reco_bins);
        bin_tree->Branch("two_diff_unmatched_bins", &two_diff_unmatched_bins);
        bin_tree->Branch("four_diff_unmatched_bins", &four_diff_unmatched_bins);
        bin_tree->Branch("two_diff_bins", &two_diff_bins);
        bin_tree->Branch("four_diff_bins", &four_diff_bins);
        bin_tree->Fill();

        fout->cd();
        event_tree_corr->Write();
        bin_tree->Write();
        ref_tree->Write();
        fin->Close();
        fout->Close();
    }
 
    std::cout << "Done calculating "<< harmonic << "-particle correlations" << std::endl;

    return 0;
}


std::vector<double> CorrCalculator::GetBinning(std::pair<double,double> minmax, const int n_bins)
{
    std::vector<double> bins;
    const double delta = (1.05*minmax.second - 0.95*minmax.first)/n_bins;
    for (int i = 0; i < n_bins+1; i++){  bins.push_back(minmax.first + i*delta); }
    return bins;
}

// std::vector<double> CorrCalculator::GetCorrBinning(std::pair<double,double> minmax, const int n_bins)
// {
//     double max_abs = std::max(std::abs(minmax.first), std::abs(minmax.second));
//     max_abs*=0.8;
//     // // make sure the binning is symmetric around 0
//     double delta = 2.0*max_abs/n_bins;
//     std::vector<double> bins;
//     for (int i = 0; i < n_bins+1; i++){  bins.push_back(-max_abs + i*delta); }
//     // return bins;
//     // std::vector<double> bins;
//     // const double delta = (minmax.second - minmax.first)/n_bins;
//     // for (int i = 0; i < n_bins+1; i++){  bins.push_back(minmax.first + i*delta); }
//     return bins;
// }
std::vector<double> CorrCalculator::GetCorrBinning(std::pair<double,double> minmax, const int n_bins)
{
    double max_abs = std::max(std::abs(minmax.first), std::abs(minmax.second));
    max_abs*=0.5;
    // make sure the binning is symmetric around 0
    double delta = 2.0*max_abs/n_bins;
    std::vector<double> bins;
    for (int i = 0; i < n_bins+1; i++){  bins.push_back(-max_abs + i*delta); }
    return bins;
}

std::string CorrCalculator::GetOutputFileName(std::string inputfile, std::string output_dir, int harmonic)
{
    // subtract the multiplicity
    TString inputfile_base = inputfile.c_str();
    inputfile_base.ReplaceAll(".root", "");
    // remove the path (everything before the last /)
    inputfile_base = inputfile_base(inputfile_base.Last('/')+1, inputfile_base.Length()-inputfile_base.Last('/')+1);
    TString outputfile = Form("%s/%s_v%d_correlations.root", output_dir.c_str(), inputfile_base.Data(), harmonic);
    return outputfile.Data();
}