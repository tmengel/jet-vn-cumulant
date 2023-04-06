
#include "TFile.h"
#include "TTree.h"
#include "TH2.h"
#include "TString.h"
#include "TSystem.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"

#include "TRandom3.h"
#include "TTimeStamp.h"
#include "TMath.h"

#include "fastjet/ClusterSequence.hh"
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/config.h"

#include <iostream>
#include <vector>

using namespace std;
using namespace fastjet;


const Double_t JetV2 = 0.04;
const Double_t JetV3 = 0.01;
const Double_t JetV4 = 0.005;

const Int_t MAX_PARTICLES = 5000;
const Int_t MAX_JETS = 100;

// struct contains all information from TTree for a single event
struct EventReader_t{

    TFile *input_file;
    TTree *input_tree;
    // TTreeReader *input_tree_reader;
    // TTreeReaderValue<Int_t> *event_id_number;
    // TTreeReaderValue<Int_t> *event_cent_bin;
    // TTreeReaderValue<Int_t> *event_pt_hard_bin;
    // TTreeReaderValue<Int_t> *event_multiplicity;
    // TTreeReaderValue<Int_t> *event_num_pythia_jets;
    // TTreeReaderValue<Double_t> *event_weight;
    // TTreeReaderValue<Double_t> *event_psi1;
    // TTreeReaderValue<Double_t> *event_psi2;
    // TTreeReaderValue<Double_t> *event_psi3;
    // TTreeReaderValue<Double_t> *event_psi4;
    // TTreeReaderValue<Double_t> *event_truth_jet_v2;
    // TTreeReaderValue<Double_t> *event_truth_jet_v3;
    // TTreeReaderValue<Double_t> *event_truth_jet_v4;
    // TTreeReaderValue<Double_t> *event_pythia_phi_shift;

    // TTreeReaderArray<Double_t> *particle_px;
    // TTreeReaderArray<Double_t> *particle_py;
    // TTreeReaderArray<Double_t> *particle_pz;
    // TTreeReaderArray<Double_t> *particle_E;
    // TTreeReaderArray<Double_t> *pythia_jets_original_phi;
    // TTreeReaderArray<Double_t> *pythia_jets_pt;
    // TTreeReaderArray<Double_t> *pythia_jets_realigned_phi;
    Int_t current_entry;
    Int_t tree_entries;
    Int_t event_id_number;
    Int_t event_cent_bin;
    Int_t event_pt_hard_bin;
    Int_t event_multiplicity;
    Int_t event_num_pythia_jets;
    Double_t event_weight;
    Double_t event_psi1;
    Double_t event_psi2;
    Double_t event_psi3;
    Double_t event_psi4;
    Double_t event_truth_jet_v2;
    Double_t event_truth_jet_v3;
    Double_t event_truth_jet_v4;
    Double_t event_pythia_phi_shift;

    Double_t particle_px[MAX_PARTICLES];
    Double_t particle_py[MAX_PARTICLES];
    Double_t particle_pz[MAX_PARTICLES];
    Double_t particle_E[MAX_PARTICLES];
    Double_t pythia_jets_original_phi[MAX_JETS];
    Double_t pythia_jets_pt[MAX_JETS];
    Double_t pythia_jets_realigned_phi[MAX_JETS];

    // TH2D *v2_pi_tenngen;
    // TH2D *v2_p_tenngen;
    // TH2D *v2_K_tenngen;
    // TH2D *v3_pi_tenngen;
    // TH2D *v3_p_tenngen;
    // TH2D *v3_K_tenngen;
    // TH2D *v4_pi_tenngen;
    // TH2D *v4_p_tenngen;
    // TH2D *v4_K_tenngen;

    std::vector<PseudoJet> particles;
    std::vector<PseudoJet> jets;

};

void InitEventReader(struct EventReader_t *reader, TString inputfile){

    reader->input_file = new TFile(inputfile.Data(), "READ");
    reader->input_tree = (TTree*)reader->input_file->Get("tree");
    reader->input_tree->SetBranchAddress("event_id_number", &reader->event_id_number);
    reader->input_tree->SetBranchAddress("event_cent_bin", &reader->event_cent_bin);
    reader->input_tree->SetBranchAddress("event_pt_hard_bin", &reader->event_pt_hard_bin);
    reader->input_tree->SetBranchAddress("event_multiplicity", &reader->event_multiplicity);
    reader->input_tree->SetBranchAddress("event_num_pythia_jets", &reader->event_num_pythia_jets);
    reader->input_tree->SetBranchAddress("event_weight", &reader->event_weight);
    reader->input_tree->SetBranchAddress("event_psi1", &reader->event_psi1);
    reader->input_tree->SetBranchAddress("event_psi2", &reader->event_psi2);
    reader->input_tree->SetBranchAddress("event_psi3", &reader->event_psi3);
    reader->input_tree->SetBranchAddress("event_psi4", &reader->event_psi4);
    reader->input_tree->SetBranchAddress("event_truth_jet_v2", &reader->event_truth_jet_v2);
    reader->input_tree->SetBranchAddress("event_truth_jet_v3", &reader->event_truth_jet_v3);
    reader->input_tree->SetBranchAddress("event_truth_jet_v4", &reader->event_truth_jet_v4);
    reader->input_tree->SetBranchAddress("event_pythia_phi_shift", &reader->event_pythia_phi_shift);
    reader->input_tree->SetBranchAddress("particle_px", reader->particle_px);
    reader->input_tree->SetBranchAddress("particle_py", reader->particle_py);
    reader->input_tree->SetBranchAddress("particle_pz", reader->particle_pz);
    reader->input_tree->SetBranchAddress("particle_E", reader->particle_E);
    reader->input_tree->SetBranchAddress("pythia_jets_original_phi", reader->pythia_jets_original_phi);
    reader->input_tree->SetBranchAddress("pythia_jets_pt", reader->pythia_jets_pt);
    reader->input_tree->SetBranchAddress("pythia_jets_realigned_phi", reader->pythia_jets_realigned_phi);

    reader->current_entry = 0;
    reader->tree_entries = reader->input_tree->GetEntries();

    // reader->input_tree_reader = new TTreeReader(reader->input_tree);
    // reader->event_id_number = new TTreeReaderValue<Int_t>(*reader->input_tree_reader, "event_id_number");
    // reader->event_cent_bin = new TTreeReaderValue<Int_t>(*reader->input_tree_reader, "event_cent_bin");
    // reader->event_pt_hard_bin = new TTreeReaderValue<Int_t>(*reader->input_tree_reader, "event_pt_hard_bin");
    // reader->event_multiplicity = new TTreeReaderValue<Int_t>(*reader->input_tree_reader, "event_multiplicity");
    // reader->event_num_pythia_jets = new TTreeReaderValue<Int_t>(*reader->input_tree_reader, "event_num_pythia_jets");
    // reader->event_weight = new TTreeReaderValue<Double_t>(*reader->input_tree_reader, "event_weight");
    // reader->event_psi1 = new TTreeReaderValue<Double_t>(*reader->input_tree_reader, "event_psi1");
    // reader->event_psi2 = new TTreeReaderValue<Double_t>(*reader->input_tree_reader, "event_psi2");
    // reader->event_psi3 = new TTreeReaderValue<Double_t>(*reader->input_tree_reader, "event_psi3");
    // reader->event_psi4 = new TTreeReaderValue<Double_t>(*reader->input_tree_reader, "event_psi4");
    // reader->event_truth_jet_v2 = new TTreeReaderValue<Double_t>(*reader->input_tree_reader, "event_truth_jet_v2");
    // reader->event_truth_jet_v3 = new TTreeReaderValue<Double_t>(*reader->input_tree_reader, "event_truth_jet_v3");
    // reader->event_truth_jet_v4 = new TTreeReaderValue<Double_t>(*reader->input_tree_reader, "event_truth_jet_v4");
    // reader->event_pythia_phi_shift = new TTreeReaderValue<Double_t>(*reader->input_tree_reader, "event_pythia_phi_shift");
    // reader->particle_px = new TTreeReaderArray<Double_t>(*reader->input_tree_reader, "particle_px");
    // reader->particle_py = new TTreeReaderArray<Double_t>(*reader->input_tree_reader, "particle_py");
    // reader->particle_pz = new TTreeReaderArray<Double_t>(*reader->input_tree_reader, "particle_pz");
    // reader->particle_E = new TTreeReaderArray<Double_t>(*reader->input_tree_reader, "particle_E");
    // reader->pythia_jets_original_phi = new TTreeReaderArray<Double_t>(*reader->input_tree_reader, "pythia_jets_original_phi");
    // reader->pythia_jets_pt = new TTreeReaderArray<Double_t>(*reader->input_tree_reader, "pythia_jets_pt");
    // reader->pythia_jets_realigned_phi = new TTreeReaderArray<Double_t>(*reader->input_tree_reader, "pythia_jets_realigned_phi");

    // reader->v2_pi_tenngen = (TH2D*)reader->input_file->Get("v2_pi");
    // reader->v2_p_tenngen = (TH2D*)reader->input_file->Get("v2_p");
    // reader->v2_K_tenngen = (TH2D*)reader->input_file->Get("v2_K");
    // reader->v3_pi_tenngen = (TH2D*)reader->input_file->Get("v3_pi");
    // reader->v3_p_tenngen = (TH2D*)reader->input_file->Get("v3_p");
    // reader->v3_K_tenngen = (TH2D*)reader->input_file->Get("v3_K");
    // reader->v4_pi_tenngen = (TH2D*)reader->input_file->Get("v4_pi");
    // reader->v4_p_tenngen = (TH2D*)reader->input_file->Get("v4_p");
    // reader->v4_K_tenngen = (TH2D*)reader->input_file->Get("v4_K");

    // reader->v2_pi_tenngen->SetName("v2_pi_tenngen");
    // reader->v2_p_tenngen->SetName("v2_p_tenngen");
    // reader->v2_K_tenngen->SetName("v2_K_tenngen");
    // reader->v3_pi_tenngen->SetName("v3_pi_tenngen");
    // reader->v3_p_tenngen->SetName("v3_p_tenngen");
    // reader->v3_K_tenngen->SetName("v3_K_tenngen");
    // reader->v4_pi_tenngen->SetName("v4_pi_tenngen");
    // reader->v4_p_tenngen->SetName("v4_p_tenngen");
    // reader->v4_K_tenngen->SetName("v4_K_tenngen");

    reader->particles.clear();
    reader->jets.clear();

    if(!reader->input_file->IsOpen()){ cout << "Input file not open" << endl; exit(-1);}

}

void FindJets(struct EventReader_t *reader, Double_t jetparam, Double_t min_jet_pt){

    for (Int_t i = 0; i < (reader->event_multiplicity); i++){
        PseudoJet particle_temp((reader->particle_px)[i], (reader->particle_py)[i], (reader->particle_pz)[i], (reader->particle_E)[i]);
        reader->particles.push_back(particle_temp);
    }

    Double_t particle_max_eta = 0.9;
    Double_t ghost_max_eta = 0.9;
    Double_t jet_min_pt_antikt = min_jet_pt;
    Double_t antikt_jet_abs_eta_max = particle_max_eta - jetparam;
    fastjet::GhostedAreaSpec ghost_area_spec(ghost_max_eta);
    fastjet::JetDefinition antikt_jet_def(fastjet::antikt_algorithm, jetparam);
    fastjet::AreaDefinition area_def(fastjet::active_area, ghost_area_spec);  
    fastjet::ClusterSequenceArea antikt_cs(reader->particles, antikt_jet_def, area_def);

    reader->jets = fastjet::sorted_by_pt(antikt_cs.inclusive_jets(jet_min_pt_antikt));
    
}

Int_t NextEvent(struct EventReader_t *reader){
    
    if(reader->current_entry >= reader->tree_entries){ cout << "No more events" << endl; return 0;}
    reader->input_tree->GetEntry(reader->current_entry);
    reader->current_entry++;
    reader->particles.clear();
    reader->jets.clear();
    return 1;
}

void CloseEventReader(struct EventReader_t *reader){
    reader->input_file->Close();
}

Int_t main(Int_t argc, Char_t** argv){

    // check if correct number of arguments
    if(argc != 4){ cout << "Incorrect number of arguments" << endl; exit(-1);}

    // get arguments
    TString input_file = argv[1];
    Double_t jetparam = atof(argv[2]);
    Double_t min_jet_pt = 10.0;
    TString output_file = argv[3];

    // initialize event reader
    struct EventReader_t reader;
    InitEventReader(&reader, input_file);

    TFile *fout = new TFile(output_file.Data(), "RECREATE");
    TTree *outTree = new TTree("tree", "tree");

    // output variables
    Int_t event_number;
    Int_t event_multiplicity;
    Int_t event_num_jets;
    Int_t event_num_pythia_jets;
    Double_t event_truth_psi1, event_truth_psi2, event_truth_psi3, event_truth_psi4;
    Double_t event_truth_jetv2, event_truth_jetv3, event_truth_jetv4;
    // Double_t jet_phi[MAX_JETS], jet_pt[MAX_JETS];
    // Double_t particle_phi[MAX_PARTICLES], particle_pt[MAX_PARTICLES];

    // output branches
    outTree->Branch("event_number", &event_number, "event_number/I");
    outTree->Branch("event_multiplicity", &event_multiplicity, "event_multiplicity/I");
    outTree->Branch("event_num_jets", &event_num_jets, "event_num_jets/I");
    outTree->Branch("event_num_pythia_jets", &event_num_pythia_jets, "event_num_pythia_jets/I");
    outTree->Branch("event_truth_psi1", &event_truth_psi1, "event_truth_psi1/D");
    outTree->Branch("event_truth_psi2", &event_truth_psi2, "event_truth_psi2/D");
    outTree->Branch("event_truth_psi3", &event_truth_psi3, "event_truth_psi3/D");
    outTree->Branch("event_truth_psi4", &event_truth_psi4, "event_truth_psi4/D");
    outTree->Branch("event_truth_jetv2", &event_truth_jetv2, "event_truth_jetv2/D");
    outTree->Branch("event_truth_jetv3", &event_truth_jetv3, "event_truth_jetv3/D");
    outTree->Branch("event_truth_jetv4", &event_truth_jetv4, "event_truth_jetv4/D");
    // outTree->Branch("jet_phi", jet_phi, "jet_phi[event_num_jets]/D");
    // outTree->Branch("jet_pt", jet_pt, "jet_pt[event_num_jets]/D");
    // outTree->Branch("particle_phi", particle_phi, "particle_phi[event_multiplicity]/D");
    // outTree->Branch("particle_pt", particle_pt, "particle_pt[event_multiplicity]/D");


    TH2D *v2_pi_truth = (TH2D*)reader.input_file->Get("v2_pi");
    TH2D *v2_p_truth = (TH2D*)reader.input_file->Get("v2_p");
    TH2D *v2_K_truth = (TH2D*)reader.input_file->Get("v2_k");
    TH2D *v3_pi_truth = (TH2D*)reader.input_file->Get("v3_pi");
    TH2D *v3_p_truth = (TH2D*)reader.input_file->Get("v3_p");
    TH2D *v3_K_truth = (TH2D*)reader.input_file->Get("v3_k");
    TH2D *v4_pi_truth = (TH2D*)reader.input_file->Get("v4_pi");
    TH2D *v4_p_truth = (TH2D*)reader.input_file->Get("v4_p");
    TH2D *v4_K_truth = (TH2D*)reader.input_file->Get("v4_k");
    
    TH2D *v2_jet = new TH2D("v2_jet", "v2_jet", 100, 0, 100, 100, 0, 1);
    TH2D *v3_jet = new TH2D("v3_jet", "v3_jet", 100, 0, 100, 100, 0, 1);
    TH2D *v4_jet = new TH2D("v4_jet", "v4_jet", 100, 0, 100, 100, 0, 1);

    TH2D *v2_jet_truth = new TH2D("v2_jet_truth", "v2_jet_truth", 100, 0, 100, 100, 0, 1);
    TH2D *v3_jet_truth = new TH2D("v3_jet_truth", "v3_jet_truth", 100, 0, 100, 100, 0, 1);
    TH2D *v4_jet_truth = new TH2D("v4_jet_truth", "v4_jet_truth", 100, 0, 100, 100, 0, 1);

    TH1D *dNdPhi_Particles = new TH1D("dNdPhi_Particles", "dNdPhi_Particles", 100, 0, 2*TMath::Pi());
    TH1D *dNdPhi_Jets = new TH1D("dNdPhi_Jets", "dNdPhi_Jets", 100, 0, 2*TMath::Pi());

    TH1D *Jet_Pt_Spectrum = new TH1D("Jet_Pt_Spectrum", "Jet_Pt_Spectrum", 100, 0, 100);

    // loop over events

    cout << "Starting event loop" << endl;
    while(NextEvent(&reader)){

        // find jets
        FindJets(&reader, jetparam, min_jet_pt);
        // loop over jets
        
        std::vector<Double_t> particle_phi;
        std::vector<Double_t> particle_pt;

        std::vector<Double_t> jet_phi;
        std::vector<Double_t> jet_pt;

        for (Int_t ijet = 0; ijet< reader.jets.size(); ijet++)
        {
            jet_phi.push_back((reader.jets.at(ijet)).phi());
            jet_pt.push_back((reader.jets.at(ijet)).pt());

            Jet_Pt_Spectrum->Fill((reader.jets.at(ijet)).pt());
            dNdPhi_Jets->Fill((reader.jets.at(ijet)).phi());

            v2_jet_truth->Fill((reader.jets.at(ijet)).pt(), reader.event_truth_jet_v2);
            v3_jet_truth->Fill((reader.jets.at(ijet)).pt(), reader.event_truth_jet_v3);
            v4_jet_truth->Fill((reader.jets.at(ijet)).pt(), reader.event_truth_jet_v4);

        }

        for (Int_t ipart = 0; ipart < reader.event_multiplicity; ipart++)
        {
            particle_phi.push_back((reader.particles.at(ipart)).phi());
            particle_pt.push_back((reader.particles.at(ipart)).pt());

            dNdPhi_Particles->Fill((reader.particles.at(ipart)).phi());
        }

        event_number= reader.event_id_number;
        event_num_jets = reader.jets.size();
        event_num_pythia_jets = reader.event_num_pythia_jets;
        event_multiplicity = reader.event_multiplicity;
        event_truth_psi1 = reader.event_psi1;
        event_truth_psi2 = reader.event_psi2;
        event_truth_psi3 = reader.event_psi3;
        event_truth_psi4 = reader.event_psi4;
        event_truth_jetv2 = reader.event_truth_jet_v2;
        event_truth_jetv3 = reader.event_truth_jet_v3;
        event_truth_jetv4 = reader.event_truth_jet_v4;

        outTree->Fill();      
    }

    // write output file

    fout->cd();

    v2_pi_truth->Write();
    v2_p_truth->Write();
    v2_K_truth->Write();
    v3_pi_truth->Write();
    v3_p_truth->Write();
    v3_K_truth->Write();
    v4_pi_truth->Write();
    v4_p_truth->Write();
    v4_K_truth->Write();

    v2_jet->Write();
    v3_jet->Write();
    v4_jet->Write();

    v2_jet_truth->Write();
    v3_jet_truth->Write();
    v4_jet_truth->Write();

    dNdPhi_Particles->Write();
    dNdPhi_Jets->Write();

    Jet_Pt_Spectrum->Write();

    fout->Write();
    fout->Close();

    CloseEventReader(&reader);
    return 0;

}

