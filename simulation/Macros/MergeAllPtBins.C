#include "TFile.h"
#include "TTree.h"
#include "TH2.h"
#include "TString.h"
#include "TSystem.h"
#include "TChain.h"

#include <iostream>
#include <vector>

using namespace std;

void MergeAllPtBins(TString input_dir, Int_t n_pt_bins, Float_t jetparam, TString output_file) {
  
    // Create a chain of root trees
    TChain *chain = new TChain("tree");
    TString jet_res = Form("R0%.0f", jetparam*10);
    if(input_dir.EndsWith("/")) input_dir.Remove(input_dir.Length()-1);
    std::vector<TString> input_files;
    for (Int_t i = 0; i < n_pt_bins; i++) {
        TString input_file = Form("%s/200GeV_MixedEvents_ptbin%d_RealignJets_%s.root", input_dir.Data(), i, jet_res.Data());
        input_files.push_back(input_file);
        chain->Add(input_file);
    }
    cout << "Merging " << n_pt_bins << " pt bins into " << output_file << endl;
    chain->Merge(output_file.Data());

    // get histograms from each file and add them to histogram in output file
    cout << "Adding histograms to output file" << endl;
    TFile *output = new TFile(output_file.Data(), "UPDATE");
    
    TH2D* v2_pi = new TH2D("v2_pi_tenngen", "v2_pi_tenngen", 50,0,5,50,0,1);
    TH2D* v2_k = new TH2D("v2_k_tenngen", "v2_k_tenngen", 50,0,5,50,0,1);
    TH2D* v2_p = new TH2D("v2_p_tenngen", "v2_p_tenngen", 50,0,5,50,0,1);
    TH2D* v3_pi = new TH2D("v3_pi_tenngen", "v3_pi_tenngen", 50,0,5,50,0,1);
    TH2D* v3_k = new TH2D("v3_k_tenngen", "v3_k_tenngen", 50,0,5,50,0,1);
    TH2D* v3_p = new TH2D("v3_p_tenngen", "v3_p_tenngen", 50,0,5,50,0,1);
    TH2D* v4_pi = new TH2D("v4_pi_tenngen", "v4_pi_tenngen", 50,0,5,50,0,1);
    TH2D* v4_k = new TH2D("v4_k_tenngen", "v4_k_tenngen", 50,0,5,50,0,1);
    TH2D* v4_p = new TH2D("v4_p_tenngen", "v4_p_tenngen", 50,0,5,50,0,1);
    
    for (Int_t i = 0; i < n_pt_bins; i++) {
        TFile *input = new TFile(input_files[i].Data());
        TH2D* v2_pi_temp = (TH2D*)input->Get("v2_pi");
        TH2D* v2_k_temp = (TH2D*)input->Get("v2_k");
        TH2D* v2_p_temp = (TH2D*)input->Get("v2_p");
        TH2D* v3_pi_temp = (TH2D*)input->Get("v3_pi");
        TH2D* v3_k_temp = (TH2D*)input->Get("v3_k");
        TH2D* v3_p_temp = (TH2D*)input->Get("v3_p");
        TH2D* v4_pi_temp = (TH2D*)input->Get("v4_pi");
        TH2D* v4_k_temp = (TH2D*)input->Get("v4_k");
        TH2D* v4_p_temp = (TH2D*)input->Get("v4_p");
        v2_pi->Add(v2_pi_temp);
        v2_k->Add(v2_k_temp);
        v2_p->Add(v2_p_temp);
        v3_pi->Add(v3_pi_temp);
        v3_k->Add(v3_k_temp);
        v3_p->Add(v3_p_temp);
        v4_pi->Add(v4_pi_temp);
        v4_k->Add(v4_k_temp);
        v4_p->Add(v4_p_temp);
        input->Close();
    }
    cout << "Done" << endl;
    output->cd();
    v2_pi->Write();
    v2_k->Write();
    v2_p->Write();
    v3_pi->Write();
    v3_k->Write();
    v3_p->Write();
    v4_pi->Write();
    v4_k->Write();
    v4_p->Write();
    output->Close();


 
}

Int_t main(Int_t argc, Char_t** argv) {
    if (argc != 5) {
        cout << "Usage: ./MergeAllPtBins <input_dir> <n_pt_bins> <jetparam> <output_file>" << endl;
        return 1;
    }
    TString input_dir = argv[1];
    Int_t n_pt_bins = atoi(argv[2]);
    Float_t jetparam = atof(argv[3]);
    TString output_file = argv[4];
    MergeAllPtBins(input_dir, n_pt_bins, jetparam, output_file);
    return 0;
}