
#include "TFile.h"
#include "TTree.h"
#include "TH2.h"
#include "TH1D.h"
#include "TString.h"
#include "TSystem.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"
#include "TClass.h"

#include "TRandom3.h"
#include "TTimeStamp.h"
#include "TMath.h"

#include "fastjet/ClusterSequence.hh"
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/config.h"

#include "TComplex.h"

#include <iostream>
#include <vector>

using namespace std;
using namespace fastjet;
using namespace ROOT;

const Double_t t_JetV2 = 0.04;
const Double_t t_JetV3 = 0.01;
const Double_t t_JetV4 = 0.005;

const Int_t MAX_PARTICLES = 10000;
const Int_t MAX_JETS = 100;
const Double_t MIN_JET_PT = 10.0;

const Int_t NUM_JET_PT_BINS = 7;
const Double_t JET_PT_BINS[NUM_JET_PT_BINS+1] = {10, 20, 30, 40, 50, 60, 70, 80};



// TComplex Q_vector_second_order_mag(std::vector<Double_t> phi_vec, Int_t n){
//     Int_t M = phi_vec.size();
//     TComplex Qn(M, 0);
//     for (Int_t i = 0; i < phi_vec.size(); i++){
//         for (Int_t j = 0; j < phi_vec.size(); j++){
//             if(i != j){
//                 Qn += TComplex(TMath::Cos(n*(phi_vec[i] - phi_vec[j])), TMath::Sin(n*(phi_vec[i] - phi_vec[j])));
//             }
//         }
//     }
//     return Qn;   
// }

// TComplex Q_vector_fourth_order_mag(std::vector<Double_t> phi_vec, Int_t n){
//     Int_t M = phi_vec.size();
//     TComplex Qn(0, 0);
//     for (Int_t i = 0; i < phi_vec.size(); i++)
//     {
//         for (Int_t j = 0; j < phi_vec.size(); j++)
//         {
//             for (Int_t k = 0; k < phi_vec.size(); k++)
//             {
//                 for (Int_t l = 0; l < phi_vec.size(); l++)
//                 {
//                     Qn += TComplex(TMath::Cos(n*(phi_vec[i] + phi_vec[j] - phi_vec[k] - phi_vec[l])), TMath::Sin(n*(phi_vec[i] + phi_vec[j] - phi_vec[k] - phi_vec[l])));
//                 }
//             }
//         }
//     }
//     return Qn;   
// }

// TComplex flow_vector(std::vector<Double_t> phi_vec, Int_t n){
//     TComplex Qn(0, 0);
//     for (Int_t i = 0; i < phi_vec.size(); i++){
//         Qn += TComplex(TMath::Cos(n*phi_vec[i]), TMath::Sin(n*phi_vec[i]));
//     }
//     return Qn;
// }

// TComplex two_particle_correlations(std::vector<Double_t> phi_vec, Int_t n){
//     TComplex Qn(0, 0);
//     for (Int_t i = 0; i < phi_vec.size(); i++){
//         for (Int_t j = 0; j < phi_vec.size(); j++){
//             if(i != j){
//                 Qn += TComplex(TMath::Cos(n*(phi_vec[i] - phi_vec[j])), TMath::Sin(n*(phi_vec[i] - phi_vec[j])));
//             }
//         }
//     }
//     Int_t M = phi_vec.size();
//     Double_t P_M2 = TMath::Factorial(M)/(TMath::Factorial(M-2));
//     return Qn/P_M2;
// }

// TComplex four_particle_correlations(std::vector<Double_t> phi_vec, Int_t n){
//     TComplex Qn(0, 0);
//     for (Int_t i = 0; i < phi_vec.size(); i++)
//     {
//         for (Int_t j = 0; j < phi_vec.size(); j++)
//         {
//             for (Int_t k = 0; k < phi_vec.size(); k++)
//             {
//                 for (Int_t l = 0; l < phi_vec.size(); l++)
//                 {
//                     if(i != j && i != k && i != l && j != k && j != l && k != l)
//                     {
//                         Qn += TComplex(TMath::Cos(n*(phi_vec[i] + phi_vec[j] - phi_vec[k] - phi_vec[l])), TMath::Sin(n*(phi_vec[i] + phi_vec[j] - phi_vec[k] - phi_vec[l])));
//                     }
//                 }
//             }
//         }
//     }
//     Int_t M = phi_vec.size();
//     Double_t P_M4 = TMath::Factorial(M)/(TMath::Factorial(M-4));
//     return Qn/P_M4;
// }

// TComplex two_particle_single_event_reduced_correlation(std::vector<Double_t> POI_vec, std::vector<Double_t> RFP_vec, std::vector<Double_t> POI_and_RFP_vec, Int_t n){
    
//     std::vector<Double_t> RFP_total_vec = RFP_vec;
//     std::vector<Double_t> POI_total_vec = POI_vec;
//     for (Int_t i = 0; i < POI_and_RFP_vec.size(); i++){
//         RFP_total_vec.push_back(POI_and_RFP_vec[i]);
//         POI_total_vec.push_back(POI_and_RFP_vec[i]);
//     }

//     Double_t mp = 1.0*POI_vec.size();
//     Double_t mq = 1.0*POI_and_RFP_vec.size();
//     Double_t M = 1.0*(RFP_total_vec.size() + POI_total_vec.size());

//     TComplex pn = flow_vector(POI_total_vec, n);
//     TComplex qn = flow_vector(POI_and_RFP_vec, n);
//     TComplex Qn = flow_vector(RFP_total_vec, n);
//     TComplex Qn_conj = TComplex::Conjugate(Qn);

//     return (pn*Qn_conj - mq)/(mp*M - mq);

// }

// TComplex four_particle_single_event_reduced_correlation(std::vector<Double_t> POI_vec, std::vector<Double_t> RFP_vec, std::vector<Double_t> POI_and_RFP_vec, Int_t n){
    
//     std::vector<Double_t> RFP_total_vec = RFP_vec;
//     std::vector<Double_t> POI_total_vec = POI_vec;
//     for (Int_t i = 0; i < POI_and_RFP_vec.size(); i++){
//         RFP_total_vec.push_back(POI_and_RFP_vec[i]);
//         POI_total_vec.push_back(POI_and_RFP_vec[i]);
//     }

//     Double_t mp = 1.0*POI_vec.size();
//     Double_t mq = 1.0*POI_and_RFP_vec.size();
//     Double_t M = 1.0*(RFP_total_vec.size() + POI_total_vec.size());

//     TComplex pn = flow_vector(POI_total_vec, n);
//     TComplex qn = flow_vector(POI_and_RFP_vec, n);
//     TComplex qn_conj = TComplex::Conjugate(qn);
//     TComplex q2n = flow_vector(POI_and_RFP_vec, 2*n);
//     TComplex Qn = flow_vector(RFP_total_vec, n);
//     TComplex Q2n = flow_vector(RFP_total_vec, 2*n);
//     TComplex Qn_conj = TComplex::Conjugate(Qn);
//     TComplex Q2n_conj = TComplex::Conjugate(Q2n);

//     return (pn*Qn*Qn_conj*Q2n_conj 
//     - q2n*Qn_conj*Qn_conj 
//     - pn*Qn*Q2n_conj 
//     - 2.0*M*pn*Qn_conj 
//     - 2.0*mq*Qn*Qn_conj 
//     + 7.0*qn*Qn_conj
//     - Qn*qn_conj
//     + q2n*Q2n_conj
//     + 2.0*pn*Qn_conj
//     + 2.0*mq*M 
//     - 6.0*mq ) /((mp*M-3*mq)*(M-1)*(M-2));

// }

TComplex FlowVector(Int_t n, std::vector<Double_t> phi_vec){
    TComplex Qn = TComplex(0,0);
    for (Int_t i = 0; i < phi_vec.size(); i++){
        Qn += TComplex(TMath::Cos(n*phi_vec[i]), TMath::Sin(n*phi_vec[i]));
    }
    return Qn;
}

Double_t SingleEventTwoParticleCorrelation(Int_t n, std::vector<Double_t> phi_vec){
    TComplex Qn = TComplex(0,0);
    Double_t Pmn = TMath::Factorial(phi_vec.size())/(TMath::Factorial(phi_vec.size()-n));
    for (Int_t i = 0; i < phi_vec.size(); i++){
      for (Int_t j = 0; j < phi_vec.size(); j++){
        if (i == j) continue;
        Qn += TComplex(TMath::Cos(n*(phi_vec[i] - phi_vec[j])), TMath::Sin(n*(phi_vec[i] - phi_vec[j])));
      }
    }
    Double_t realQn = Qn.Re();
    //Double_t imagQn = Qn.Im();
    return realQn;
}

Int_t main(Int_t argc, Char_t** argv){

    // check if correct number of arguments
    if(argc != 4){ cout << "Incorrect number of arguments" << endl; exit(-1);}

    // get arguments
    TString input_file = argv[1];
    Double_t jetparam = atof(argv[2]);
    Double_t min_jet_pt = MIN_JET_PT;
    TString output_file = argv[3];

    cout << "Input file: " << input_file << endl;
    cout << "Jet param: " << jetparam << endl;
    cout << "Min jet pt: " << min_jet_pt << endl;
    cout << "Output file: " << output_file << endl;

    // Open input file and create tree branch addresses
    cout << "Opening input file" << endl;
    TFile fin(input_file, "READ");
    if(!fin.IsOpen()){ cout << "Error opening input file" << endl; exit(-1);}
    TTreeReader reader("tree", &fin);

    // branches
    TTreeReaderValue<Int_t> event_number(reader, "event_id_number");
    // TTreeReaderValue<Int_t> event_num_jets(reader, "event_num_jets");
    TTreeReaderValue<Int_t> event_num_pythia_jets(reader, "event_num_pythia_jets");
    TTreeReaderValue<Int_t> event_cent_bin(reader, "event_cent_bin");
    // TTreeReaderValue<Int_t> event_pt_hard_bin(reader, "event_pt_hard_bin");
    TTreeReaderValue<Int_t> event_multiplicity(reader, "event_multiplicity");
    TTreeReaderValue<Double_t> event_weight(reader, "event_weight");
    TTreeReaderValue<Double_t> event_psi1(reader, "event_psi1");
    TTreeReaderValue<Double_t> event_psi2(reader, "event_psi2");
    TTreeReaderValue<Double_t> event_psi3(reader, "event_psi3");
    TTreeReaderValue<Double_t> event_psi4(reader, "event_psi4");
    TTreeReaderValue<Double_t> event_truth_jet_v2(reader, "event_truth_jet_v2");
    TTreeReaderValue<Double_t> event_truth_jet_v3(reader, "event_truth_jet_v3");
    TTreeReaderValue<Double_t> event_truth_jet_v4(reader, "event_truth_jet_v4");
    TTreeReaderValue<Double_t> event_pythia_phi_shift(reader, "event_pythia_phi_shift");
    TTreeReaderArray<Double_t> particle_px_array(reader, "particle_px");
    TTreeReaderArray<Double_t> particle_py_array(reader, "particle_py");
    TTreeReaderArray<Double_t> particle_pz_array(reader, "particle_pz");
    TTreeReaderArray<Double_t> particle_E_array(reader, "particle_E");
    TTreeReaderArray<Double_t> pythia_jets_original_phi_array(reader, "pythia_jets_original_phi");
    TTreeReaderArray<Double_t> pythia_jets_pt_array(reader, "pythia_jets_pt");
    TTreeReaderArray<Double_t> pythia_jets_realigned_phi_array(reader, "pythia_jets_realigned_phi");

    cout << "Getting tenngen vn histograms" << endl;
    // get tenngen vn histograms from input file
    TH2D *v2_pi_tenngen = (TH2D*)fin.Get("v2_pi_tenngen");
    TH2D *v3_pi_tenngen = (TH2D*)fin.Get("v3_pi_tenngen");
    TH2D *v4_pi_tenngen = (TH2D*)fin.Get("v4_pi_tenngen");
    TH2D *v2_k_tenngen = (TH2D*)fin.Get("v2_k_tenngen");
    TH2D *v3_k_tenngen = (TH2D*)fin.Get("v3_k_tenngen");
    TH2D *v4_k_tenngen = (TH2D*)fin.Get("v4_k_tenngen");
    TH2D *v2_p_tenngen = (TH2D*)fin.Get("v2_p_tenngen");
    TH2D *v3_p_tenngen = (TH2D*)fin.Get("v3_p_tenngen");
    TH2D *v4_p_tenngen = (TH2D*)fin.Get("v4_p_tenngen");

    cout << "Initializing histograms" << endl;
    // create histograms
    TH1D *event_psi1_hist = new TH1D("event_psi1_hist", "event_psi1_hist", 100,0,2*TMath::Pi());
    TH1D *event_psi2_hist = new TH1D("event_psi2_hist", "event_psi2_hist", 100,0,2*TMath::Pi());
    TH1D *event_psi3_hist = new TH1D("event_psi3_hist", "event_psi3_hist", 100,0,2*TMath::Pi());
    TH1D *event_psi4_hist = new TH1D("event_psi4_hist", "event_psi4_hist", 100,0,2*TMath::Pi());

    TH1D *pythia_event_phi_shift_hist = new TH1D("pythia_event_phi_shift_hist", "pythia_event_phi_shift_hist", 100,-TMath::Pi(),TMath::Pi());

    TH2D *v2_pi_calulated = new TH2D("v2_pi_calulated", "v2_pi_calulated",10,0,5,50,0,1);
    TH2D *v3_pi_calulated = new TH2D("v3_pi_calulated", "v3_pi_calulated",10,0,5,50,0,1);
    TH2D *v4_pi_calulated = new TH2D("v4_pi_calulated", "v4_pi_calulated",10,0,5,50,0,1);
    TH2D *v2_k_calulated = new TH2D("v2_k_calulated", "v2_k_calulated",10,0,5,50,0,1);
    TH2D *v3_k_calulated = new TH2D("v3_k_calulated", "v3_k_calulated",10,0,5,50,0,1);
    TH2D *v4_k_calulated = new TH2D("v4_k_calulated", "v4_k_calulated",10,0,5,50,0,1);
    TH2D *v2_p_calulated = new TH2D("v2_p_calulated", "v2_p_calulated",10,0,5,50,0,1);
    TH2D *v3_p_calulated = new TH2D("v3_p_calulated", "v3_p_calulated",10,0,5,50,0,1);
    TH2D *v4_p_calulated = new TH2D("v4_p_calulated", "v4_p_calulated",10,0,5,50,0,1);

    TH2D *v2_jet = new TH2D("v2_jet", "v2_jet", NUM_JET_PT_BINS, JET_PT_BINS, 100, -0.5, 1);
    TH2D *v3_jet = new TH2D("v3_jet", "v3_jet", NUM_JET_PT_BINS, JET_PT_BINS, 100, -0.5, 1);
    TH2D *v4_jet = new TH2D("v4_jet", "v4_jet", NUM_JET_PT_BINS, JET_PT_BINS, 100, -0.5, 1);

    TH2D *v2_jet_truth = new TH2D("v2_jet_truth", "v2_jet_truth", NUM_JET_PT_BINS, JET_PT_BINS, 100, -0.5, 1);
    TH2D *v3_jet_truth = new TH2D("v3_jet_truth", "v3_jet_truth", NUM_JET_PT_BINS, JET_PT_BINS, 100, -0.5, 1);
    TH2D *v4_jet_truth = new TH2D("v4_jet_truth", "v4_jet_truth", NUM_JET_PT_BINS, JET_PT_BINS, 100, -0.5, 1);

    TH1D *dNdPhi_Particles = new TH1D("dNdPhi_Particles", "dNdPhi_Particles", 100, 0, 2*TMath::Pi());

    TH1D *dNdPhi_Jets = new TH1D("dNdPhi_Jets", "dNdPhi_Jets", 100, 0, 2*TMath::Pi());
    TH1D *dNdPhi_pythia_jets_original = new TH1D("dNdPhi_pythia_jets_original", "dNdPhi_pythia_jets_original", 100,0.0,2.0*TMath::Pi());
    TH1D *dNdPhi_pythia_jets = new TH1D("dNdPhi_pythia_jets", "dNdPhi_pythia_jets", 100, 0.0,2.0*TMath::Pi());

    TH1D *Jet_Pt_Spectrum = new TH1D("Jet_Pt_Spectrum", "Jet_Pt_Spectrum", NUM_JET_PT_BINS, JET_PT_BINS);
    TH1D *Jet_Pythia_Pt_Spectrum = new TH1D("Jet_Pythia_Pt_Spectrum", "Jet_Pythia_Pt_Spectrum", NUM_JET_PT_BINS, JET_PT_BINS);

    TH1D* particle_mass_hist = new TH1D("particle_mass_hist", "particle_mass_hist", 100, 0, 5);

    // fast jet stuff   
    Double_t particle_max_eta_midrapidity = 0.9;
    Double_t ghost_max_eta = 1.0;
    Double_t jet_min_pt_antikt = MIN_JET_PT;
    Double_t antikt_jet_abs_eta_max = particle_max_eta_midrapidity - jetparam;
    fastjet::GhostedAreaSpec ghost_area_spec(ghost_max_eta);
    fastjet::JetDefinition antikt_jet_def(fastjet::antikt_algorithm, jetparam);
    fastjet::AreaDefinition area_def(fastjet::active_area, ghost_area_spec);  
    

    std::vector<fastjet::PseudoJet> mid_rapidity_particles;
    std::vector<fastjet::PseudoJet> forward_particles;
    std::vector<fastjet::PseudoJet> mid_rapidity_jets;
    // loop over events
  
    Int_t num_particle_pt_bins = 10;
    Double_t particle_pt_bins[11] = {0.0};
    for (Int_t ix = 1; ix<v2_pi_calulated->GetNbinsX()+1; ix++){
        particle_pt_bins[ix-1] = v2_pi_calulated->GetXaxis()->GetBinLowEdge(ix);
    }
    particle_pt_bins[10] = v2_pi_calulated->GetXaxis()->GetBinUpEdge(v2_pi_calulated->GetNbinsX());


    Double_t particle_pt_bins_center[10];
    for (Int_t ibin=1;ibin<v2_pi_calulated->GetNbinsX()+1; ibin++){
        particle_pt_bins_center[ibin-1] = v2_pi_calulated->GetXaxis()->GetBinCenter(ibin);
    }


    Double_t W_two_particle_pi[10]={0};
    Double_t omega_two_particle_differential_pi[10]={0};
    Double_t second_order_differential_cumulant_pi_v2[10]={0};
    Double_t second_order_differential_cumulant_pi_v3[10]={0};
    Double_t second_order_differential_cumulant_pi_v4[10]={0};

    Double_t second_order_cumulant_pi_v2[10]={0};
    Double_t second_order_cumulant_pi_v3[10]={0};
    Double_t second_order_cumulant_pi_v4[10]={0};

    cout << "Starting event loop" << endl;

    Int_t event_counter = 0;
    while(reader.Next()){
        if(event_counter >1000) break;
        
        //Clear vectors
        mid_rapidity_particles.clear();
        forward_particles.clear();
        mid_rapidity_jets.clear();

        // read in midrapidty and forward particles
        for (Int_t ipart = 0; ipart < *event_multiplicity; ipart++){
            PseudoJet particle_temp(particle_px_array[ipart], particle_py_array[ipart], particle_pz_array[ipart], particle_E_array[ipart]);
            if (particle_temp.pt() < 0.15) continue;
            //if (particle_temp.m() == 0) continue;

            if (TMath::Abs(particle_temp.eta()) < 0.9){
                mid_rapidity_particles.push_back(particle_temp);
                dNdPhi_Particles->Fill(particle_temp.phi());
            }
            
            if(TMath::Abs(particle_temp.eta()) < 3.5 && TMath::Abs(particle_temp.eta()) > 3.0){ 
                forward_particles.push_back(particle_temp);
                dNdPhi_Particles->Fill(particle_temp.phi());
            }
        }

        // get pythia jet pt and phi
        for (Int_t ijet = 0; ijet < *event_num_pythia_jets; ijet++){
            dNdPhi_pythia_jets_original->Fill(pythia_jets_original_phi_array[ijet]);
            dNdPhi_pythia_jets->Fill(pythia_jets_realigned_phi_array[ijet]);
            Jet_Pythia_Pt_Spectrum->Fill(pythia_jets_pt_array[ijet]);
        }

        // get misc event info
        event_psi1_hist->Fill(*event_psi1);
        event_psi2_hist->Fill(*event_psi2);
        event_psi3_hist->Fill(*event_psi3);
        event_psi4_hist->Fill(*event_psi4);
        pythia_event_phi_shift_hist->Fill(*event_pythia_phi_shift);

        // cluster midrapidity particles
        fastjet::ClusterSequenceArea antikt_cs(mid_rapidity_particles, antikt_jet_def, area_def);
        mid_rapidity_jets = sorted_by_pt(antikt_cs.inclusive_jets(jet_min_pt_antikt));

        // fill jet pt spectrum and dNdPhi
        for (Int_t ijet = 0; ijet < mid_rapidity_jets.size(); ijet++){
            Jet_Pt_Spectrum->Fill(mid_rapidity_jets[ijet].pt());
            dNdPhi_Jets->Fill(mid_rapidity_jets[ijet].phi());
            v2_jet_truth->Fill(mid_rapidity_jets[ijet].pt(), t_JetV2);
            v3_jet_truth->Fill(mid_rapidity_jets[ijet].pt(), t_JetV3);
            v4_jet_truth->Fill(mid_rapidity_jets[ijet].pt(), t_JetV4);
        }

        for (Int_t ipart = 0; ipart < mid_rapidity_particles.size(); ipart++){
            particle_mass_hist->Fill(mid_rapidity_particles[ipart].m());

        }


        
        std::vector<Double_t> particle_poi_phi;
        std::vector<Double_t> particle_rfp_phi;
        // loop over all particle pt bins
        for (Int_t ibin=0;ibin<num_particle_pt_bins+1; ibin++){
            
            // pion = poi, kaon = rfp
            // clear vectors
            particle_poi_phi.clear();
            particle_rfp_phi.clear();

            // loop over midrapidity particles 
            for (Int_t ipart = 0; ipart < mid_rapidity_particles.size(); ipart++){
                // check if particle is in pt bin
                if (mid_rapidity_particles[ipart].pt() > particle_pt_bins[ibin] && mid_rapidity_particles[ipart].pt() < particle_pt_bins[ibin+1]){
                    // check if particle is a pion
                    if(mid_rapidity_particles[ipart].m() < 0.14 && mid_rapidity_particles[ipart].m() > 0.13){
                        particle_poi_phi.push_back(mid_rapidity_particles[ipart].phi());
                    }
                    // check if particle is a kaon
                    if(mid_rapidity_particles[ipart].m() < 0.5 && mid_rapidity_particles[ipart].m() > 0.4){
                        particle_rfp_phi.push_back(mid_rapidity_particles[ipart].phi());
                    }

                    TComplex ref_flow_Q_2 = FlowVector(2,particle_rfp_phi);
                    TComplex ref_flow_Q_3 = FlowVector(3,particle_rfp_phi);
                    TComplex ref_flow_Q_4 = FlowVector(4,particle_rfp_phi);

                    TComplex poi_flow_Q_2 = FlowVector(2,particle_poi_phi);
                    TComplex poi_flow_Q_3 = FlowVector(3,particle_poi_phi);
                    TComplex poi_flow_Q_4 = FlowVector(4,particle_poi_phi);

                    Double_t second_order_omega = 1.0*particle_poi_phi.size()*particle_rfp_phi.size(); // mq =0
                    Double_t second_order_W = 1.0*particle_rfp_phi.size()*(particle_rfp_phi.size()-1); // mq = 1

                    Double_t single_event_two_particle_correlation_v2 = SingleEventTwoParticleCorrelation(2,particle_rfp_phi);
                    Double_t single_event_two_particle_correlation_v3 = SingleEventTwoParticleCorrelation(3,particle_rfp_phi);
                    Double_t single_event_two_particle_correlation_v4 = SingleEventTwoParticleCorrelation(4,particle_rfp_phi);

                    Double_t single_event_differential_two_particle_correlation_v2 = (poi_flow_Q_2*TComplex::Conjugate(ref_flow_Q_2)).Re();
                    Double_t single_event_differential_two_particle_correlation_v3 = (poi_flow_Q_3*TComplex::Conjugate(ref_flow_Q_3)).Re();
                    Double_t single_event_differential_two_particle_correlation_v4 = (poi_flow_Q_4*TComplex::Conjugate(ref_flow_Q_4)).Re();

                    W_two_particle_pi[ibin]+=second_order_W;
                    omega_two_particle_differential_pi[ibin]+=second_order_omega;
                    
                    second_order_differential_cumulant_pi_v2[ibin]+=(single_event_differential_two_particle_correlation_v2);
                    second_order_differential_cumulant_pi_v3[ibin]+=(single_event_differential_two_particle_correlation_v3);
                    second_order_differential_cumulant_pi_v4[ibin]+=(single_event_differential_two_particle_correlation_v4);

                    second_order_cumulant_pi_v2[ibin]+=single_event_two_particle_correlation_v2;
                    second_order_cumulant_pi_v3[ibin]+=single_event_two_particle_correlation_v3;
                    second_order_cumulant_pi_v4[ibin]+=single_event_two_particle_correlation_v4;
                }
            }
        
        }


        event_counter++;

    }

    //renormalize
    for (Int_t ibin=0;ibin<num_particle_pt_bins+1; ibin++){
        cout << "W_two_particle_pi: " << W_two_particle_pi[ibin] << endl;
        cout << "omega_two_particle_differential_pi: " << omega_two_particle_differential_pi[ibin] << endl;
        cout << "second_order_cumulant_pi_v2: " << second_order_cumulant_pi_v2[ibin] << endl;
        cout << "second_order_cumulant_pi_v3: " << second_order_cumulant_pi_v3[ibin] << endl;
        second_order_cumulant_pi_v2[ibin]/=W_two_particle_pi[ibin];
        second_order_cumulant_pi_v3[ibin]/=W_two_particle_pi[ibin];
        second_order_cumulant_pi_v4[ibin]/=W_two_particle_pi[ibin];

        second_order_differential_cumulant_pi_v2[ibin]/=omega_two_particle_differential_pi[ibin];
        second_order_differential_cumulant_pi_v3[ibin]/=omega_two_particle_differential_pi[ibin];
        second_order_differential_cumulant_pi_v4[ibin]/=omega_two_particle_differential_pi[ibin];
    }

    // calculate v2, v3, v4
    for (Int_t ibin=0;ibin<num_particle_pt_bins; ibin++){
        cout << "Normalization: " << endl;
        cout << "second_order_cumulant_pi_v2: " << second_order_cumulant_pi_v2[ibin] << endl;
        cout << "second_order_cumulant_pi_v3: " << second_order_cumulant_pi_v3[ibin] << endl;

        Double_t pi_v2 = second_order_differential_cumulant_pi_v2[ibin]/TMath::Sqrt(second_order_cumulant_pi_v2[ibin]);
        Double_t pi_v3 = second_order_differential_cumulant_pi_v3[ibin]/TMath::Sqrt(second_order_cumulant_pi_v3[ibin]);
        Double_t pi_v4 = second_order_differential_cumulant_pi_v4[ibin]/TMath::Sqrt(second_order_cumulant_pi_v4[ibin]);
        cout<<"pt bin: "<<ibin<<" v2: "<<pi_v2<<" v3: "<<pi_v3<<" v4: "<<pi_v4<<endl;
    
        v2_pi_calulated->Fill(particle_pt_bins[ibin],pi_v2);
        v3_pi_calulated->Fill(particle_pt_bins[ibin],pi_v3);
        v4_pi_calulated->Fill(particle_pt_bins[ibin],pi_v4);
    }

    // write histograms to file
    cout << "Writing histograms to file" << endl;
    // create output file
    TFile *fout = new TFile(output_file.Data(), "RECREATE");
    fout->cd();
    v2_pi_tenngen->Write();
    v3_pi_tenngen->Write();
    v4_pi_tenngen->Write();
    v2_k_tenngen->Write();
    v3_k_tenngen->Write();
    v4_k_tenngen->Write();
    v2_p_tenngen->Write();
    v3_p_tenngen->Write();
    v4_p_tenngen->Write();

    v2_pi_calulated->Write();
    v3_pi_calulated->Write();
    v4_pi_calulated->Write();
    v2_k_calulated->Write();
    v3_k_calulated->Write();
    v4_k_calulated->Write();
    v2_p_calulated->Write();
    v3_p_calulated->Write();
    v4_p_calulated->Write();

    event_psi1_hist->Write();
    event_psi2_hist->Write();
    event_psi3_hist->Write();
    event_psi4_hist->Write();

    pythia_event_phi_shift_hist->Write();
    
    v2_jet_truth->Write();
    v3_jet_truth->Write();
    v4_jet_truth->Write();

    v2_jet->Write();
    v3_jet->Write();
    v4_jet->Write();

    dNdPhi_Particles->Write();

    dNdPhi_Jets->Write();
    dNdPhi_pythia_jets->Write();
    dNdPhi_pythia_jets_original->Write();

    Jet_Pt_Spectrum->Write();
    Jet_Pythia_Pt_Spectrum->Write();

    particle_mass_hist->Write();

    cout << "Writing to file complete" << endl;
    fout->Close();
    fin.Close();

    return 0;



    // // initialize weights and cn's
    // Double_t W2_sum = 0;
    // Double_t W4_sum = 0;
    // TComplex event_averaged_two_particle_correlation[4];
    // TComplex event_averaged_four_particle_correlation[4];

    // // jet bins
    // Double_t w2_sum[NUM_JET_PT_BINS];
    // Double_t w4_sum[NUM_JET_PT_BINS];

    // TComplex event_averaged_two_particle_differential_correlation[NUM_JET_PT_BINS][4];
    // TComplex event_averaged_four_particle_differential_correlation[NUM_JET_PT_BINS][4];

    // TComplex C2[4];
    // TComplex C4[4];

    // TComplex d2[NUM_JET_PT_BINS][4];
    // TComplex d4[NUM_JET_PT_BINS][4];

    // TComplex vn_2[4];
    // TComplex vn_4[4];

    // TComplex vn_jet_2[NUM_JET_PT_BINS][4];
    // TComplex vn_jet_4[NUM_JET_PT_BINS][4];
    
    // // zero out the arrays
    // for (Int_t vn = 0; vn < 4; vn++){
    //     event_averaged_two_particle_correlation[vn] = TComplex(0,0);
    //     event_averaged_four_particle_correlation[vn] = TComplex(0,0);
    //     C2[vn] = TComplex(0,0);
    //     C4[vn] = TComplex(0,0);
    //     vn_2[vn] = TComplex(0,0);
    //     vn_4[vn] = TComplex(0,0);
    // }

    // for (Int_t jetbin; jetbin < NUM_JET_PT_BINS; jetbin++)
    // {
    //     w2_sum[jetbin] = 0;
    //     w4_sum[jetbin] = 0;
    //     for (Int_t i = 0; i < 4; i++)
    //     {
    //         event_averaged_two_particle_differential_correlation[jetbin][i] = TComplex(0,0);
    //         event_averaged_four_particle_differential_correlation[jetbin][i] = TComplex(0,0);
    //         d2[jetbin][i] = TComplex(0,0);
    //         d4[jetbin][i] = TComplex(0,0);
    //         vn_jet_2[jetbin][i] = TComplex(0,0);
    //         vn_jet_4[jetbin][i] = TComplex(0,0);
    //     }
    // }


    // cout << "Starting event loop" << endl;
    // Int_t event_counter = 0;
    // while(NextEvent(&reader)){
    //     if(event_counter == 1000) break;
    //     cout << "Event " << event_counter << endl;
    //     // find jets
    //     FindJets(&reader, jetparam, min_jet_pt);
    //     // loop over jets
        
    //     std::vector<Double_t> particle_phi;
    //     std::vector<Double_t> particle_pt;

    //     std::vector<Double_t> particle_phi_forward_rapidity;
    //     std::vector<Double_t> particle_pt_forward_rapidity;
        
    //     for (Int_t ijet = 0; ijet< reader.jets.size(); ijet++)
    //     {
    //         Jet_Pt_Spectrum->Fill((reader.jets.at(ijet)).pt());
    //         dNdPhi_Jets->Fill((reader.jets.at(ijet)).phi());
    //         v2_jet_truth->Fill((reader.jets.at(ijet)).pt(), reader.event_truth_jet_v2);
    //         v3_jet_truth->Fill((reader.jets.at(ijet)).pt(), reader.event_truth_jet_v3);
    //         v4_jet_truth->Fill((reader.jets.at(ijet)).pt(), reader.event_truth_jet_v4);

    //     }

    //     for (Int_t ipart = 0; ipart < reader.event_multiplicity; ipart++)
    //     {
    //         if(reader.particles.at(ipart).eta() > 3.0 && reader.particles.at(ipart).eta() < 3.5)
    //         {
    //             particle_phi_forward_rapidity.push_back((reader.particles.at(ipart)).phi());
    //             particle_pt_forward_rapidity.push_back((reader.particles.at(ipart)).pt());
    //         }
    //         if(reader.particles.at(ipart).eta() > 0.5 && reader.particles.at(ipart).eta() < 0.9)
    //         {
    //             particle_phi.push_back((reader.particles.at(ipart)).phi());
    //             particle_pt.push_back((reader.particles.at(ipart)).pt());
    //         }
    //         // particle_phi.push_back((reader.particles.at(ipart)).phi());
    //         // particle_pt.push_back((reader.particles.at(ipart)).pt());
    //         dNdPhi_Particles->Fill((reader.particles.at(ipart)).phi());
    //     }

    //     // cout << "Number of particles: " << particle_phi.size() << endl;
    //     // cout << "Number of jets: " << reader.jets.size() << endl;

    //     // calculate weights
    //     Double_t W_2 = particle_phi.size()*(particle_phi.size()-1);
    //     Double_t W_4 = W_2*(particle_phi.size()-2)*(particle_phi.size()-3);

    //     W2_sum += W_2;
    //     W4_sum += W_4;
        
    //     for (Int_t i = 0; i < 4; i++)
    //     {
    //         event_averaged_two_particle_correlation[i] += W_2*two_particle_correlations(particle_phi,i+1);
    //        // event_averaged_four_particle_correlation[i] += W_4*four_particle_correlations(particle_phi,i+1);
    //        event_averaged_four_particle_correlation[i]+= TComplex(0,0);
    //     }

    

    //     for (Int_t jetbin = 0; jetbin < NUM_JET_PT_BINS; jetbin ++)
    //     {
    //         std::vector<Double_t> jet_phi;
    //         std::vector<Double_t> jet_pt;

    //         for (Int_t ijet = 0; ijet< reader.jets.size(); ijet++)
    //         {
    //             if((reader.jets.at(ijet)).pt() > JET_PT_BINS[jetbin] && (reader.jets.at(ijet)).pt() < JET_PT_BINS[jetbin+1])
    //             {
    //                 jet_phi.push_back((reader.jets.at(ijet)).phi());
    //                 jet_pt.push_back((reader.jets.at(ijet)).pt());
    //             }

    //         }

    //         // differential flow
    //         Double_t w_2 = jet_phi.size()*particle_phi_forward_rapidity.size();
    //         Double_t w_4 = w_2*(particle_phi_forward_rapidity.size()-1)*(particle_phi_forward_rapidity.size()-2);

    //         w2_sum[jetbin] += w_2;
    //         w4_sum[jetbin] += w_4;

            
    //         for(Int_t i = 0; i < 4; i++)
    //         {
    //             TComplex Qn = flow_vector(particle_phi_forward_rapidity, i+1);
    //             //TComplex Q2n = flow_vector(particle_phi_forward_rapidity, 2*(i+1));
    //             TComplex pn = flow_vector(jet_phi, i+1);
    //             TComplex Qn_conj = TComplex::Conjugate(Qn);
    //             //TComplex Q2n_conj = TComplex::Conjugate(Q2n);

    //             event_averaged_two_particle_differential_correlation[jetbin][i] += w_2*(pn*Qn_conj)/w_2;
    //             // event_averaged_four_particle_differential_correlation[jetbin][i] += w_4*(  pn*Qn*Qn_conj*Qn_conj 
    //             //                                                             - pn*Qn*Q2n_conj
    //             //                                                             - 2.0*particle_phi_forward_rapidity.size()*pn*Qn_conj
    //             //                                                             - 2.0*pn*Qn_conj    )/w_4;
    //             event_averaged_four_particle_differential_correlation[jetbin][i] = TComplex(0,0);
    //         }

             
    //     }




    //     event_number = reader.event_id_number;
    //     event_num_jets = reader.jets.size();
    //     event_num_pythia_jets = reader.event_num_pythia_jets;
    //     event_multiplicity = reader.event_multiplicity;
    //     event_truth_psi1 = reader.event_psi1;
    //     event_truth_psi2 = reader.event_psi2;
    //     event_truth_psi3 = reader.event_psi3;
    //     event_truth_psi4 = reader.event_psi4;
    //     event_truth_jetv2 = reader.event_truth_jet_v2;
    //     event_truth_jetv3 = reader.event_truth_jet_v3;
    //     event_truth_jetv4 = reader.event_truth_jet_v4;

    //     outTree->Fill();      

    //     event_counter++;
    // }
    // cout << "Number of events: " << event_counter << endl;
    // // normalize event averaged quantities
    // for (Int_t i = 0; i < 4; i++)
    // {
    //     event_averaged_two_particle_correlation[i] /= W2_sum;
    //     //event_averaged_four_particle_correlation[i] /= W4_sum;
    //     for (Int_t jetbin = 0; jetbin < NUM_JET_PT_BINS; jetbin ++)
    //     {
    //         event_averaged_two_particle_differential_correlation[jetbin][i] /= w2_sum[jetbin];
    //         //event_averaged_four_particle_differential_correlation[jetbin][i] /= w4_sum[jetbin];
    //     }
    // }

    // // calculate particle cumulants
    // for (Int_t i = 0; i < 4; i++)
    // {
    //     C2[i] = event_averaged_two_particle_correlation[i];
    //     //C4[i] = event_averaged_four_particle_correlation[i] - 2.0*event_averaged_two_particle_correlation[i]*event_averaged_two_particle_correlation[i];
    //     for (Int_t jetbin = 0; jetbin < NUM_JET_PT_BINS; jetbin ++)
    //     {
    //         d2[jetbin][i] = event_averaged_two_particle_differential_correlation[jetbin][i];
    //         //d4[jetbin][i] = event_averaged_four_particle_differential_correlation[jetbin][i] - 2.0*event_averaged_two_particle_differential_correlation[jetbin][i]*event_averaged_two_particle_correlation[i];
    //     }
    // }

    // // calculate flow 
    // for (Int_t i = 0; i < 4; i++)
    // {
    //     vn_2[i] = TMath::Sqrt(C2[i]);
    //     //vn_4[i] = TMath::Power(-1.0*C4[i], 1.0/4.0);
    //     for (Int_t jetbin = 0; jetbin < NUM_JET_PT_BINS; jetbin ++)
    //     {
    //         vn_jet_2[jetbin][i] = d2[jetbin][i]/vn_2[i];
    //         //vn_jet_4[jetbin][i] = d4[jetbin][i]/(TMath::Power(-1.0*C4[i], 3.0/4.0));
    //     }
    // }
    // // cout << "finished calculating flow" << endl;
    // // print debug info
   
    // for(Int_t i = 0; i < NUM_JET_PT_BINS; i++)
    // {
    //     //cout << "jet pt bin: " << i << endl;
    //     //cout << "vn_jet_2: " << TComplex::Abs(vn_jet_2[i][0]) << endl;
    // //    cout << "jet pt bin: " << i << " vn_jet_2: " << TComplex::Abs(vn_jet_2[i][1]) << " vn_jet_3: " << vn_jet_2[i][1].Re() << " vn_jet_4: " << vn_jet_2[i][2].Re() << endl;
    //     v2_jet->SetBinContent(i+1, vn_jet_2[i][1].Re());
    //     v3_jet->SetBinContent(i+1, vn_jet_2[i][2].Re());
    //     v4_jet->SetBinContent(i+1, vn_jet_2[i][3].Re());
    // }
    // // cout << "finished filling jet v2" << endl;

    // // write output file

    // fout->cd();
    // cout << "writing output file" << endl;

    // // v2_pi_truth->Write();
    // // v2_p_truth->Write();
    // // v2_K_truth->Write();
    // // v3_pi_truth->Write();
    // // v3_p_truth->Write();
    // // v3_K_truth->Write();
    // // v4_pi_truth->Write();
    // // v4_p_truth->Write();
    // // v4_K_truth->Write();
    // // cout << "finished writing truth flow" << endl;
    // v2_jet->Write();
    // v3_jet->Write();
    // v4_jet->Write();
    // // // cout << "finished writing jet flow" << endl;
    // v2_jet_truth->Write();
    // v3_jet_truth->Write();
    // v4_jet_truth->Write();
    // // cout << "finished writing jet truth flow" << endl;
    // dNdPhi_Particles->Write();
    // dNdPhi_Jets->Write();
    // // cout << "finished writing dNdPhi" << endl;
    // Jet_Pt_Spectrum->Write();
    // // cout << "finished writing jet pt spectrum" << endl;
    // fout->Write();
    // cout << "finished writing output file" << endl;
    // fout->Close();
    // cout << "finished closing output file" << endl;

    // CloseEventReader(&reader);
    // cout << "finished closing event reader" << endl;
    // return 0;

}

