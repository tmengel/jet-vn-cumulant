
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
const Double_t MAX_DNDPHI = (1.0+2.0*(TMath::Abs(JetV2)+TMath::Abs(JetV3)+TMath::Abs(JetV4)));

const Int_t MAX_PARTICLES = 10000;
const Int_t MAX_JETS = 100;


Double_t dNdPhi(Double_t phiPart, Double_t Psi2, Double_t Psi3, Double_t Psi4, Double_t v2,Double_t v3,Double_t v4){
    return (1.0+2.0*(v2*TMath::Cos(2.0*(phiPart-Psi2)) + v3*TMath::Cos(3.0*(phiPart-Psi3))+ v4*TMath::Cos(4.0*(phiPart-Psi4)) ) );
}


Int_t RealignJets(TString input_file_pythia, TString input_file_tenngen, TString output_file, Double_t jetparam, Int_t nEvents = -1){
    
    // Initialize random number generator
    TTimeStamp *timestamp = new TTimeStamp();
    UInt_t Seed = UInt_t(timestamp->GetSec());
    UInt_t seed_uniform_random_phi = (Seed*3)*4;
    TRandom3 *Harmonics_Phi_Dist_Rand = new TRandom3(seed_uniform_random_phi);

    // print out starting message
    cout << "Starting Jet Alignment" << endl;
    // open input files
    TFile pythia_file(input_file_pythia.Data());
    TFile tenngen_file(input_file_tenngen.Data());
    // check if files are open
    if(!pythia_file.IsOpen()){ cout << "PYTHIA file not open" << endl; exit(-1);}
    if(!tenngen_file.IsOpen()){ cout << "TennGen file not open" << endl; exit(-1);}

    // tree names
    TString partTree = "tree";
    TString eventTree = "eventInfo";
    // create tree readers
    
    // PYTHIA event info
    TTreeReader pythia_event(eventTree.Data(), &pythia_file);
    TTreeReaderValue<Int_t> pythia_nevent(pythia_event, "nevent");
    // TTreeReaderValue<Double_t> pythia_ptmax(pythia_event, "ptmin");
    // TTreeReaderValue<Double_t> pythia_ptmin(pythia_event, "ptmax");
    TTreeReaderValue<Int_t> pythia_ptbin(pythia_event, "ptbinID");
    TTreeReaderValue<Double_t> pythia_xsec(pythia_event, "xsec_over_eventweight");
    
    // PYTHIA particles
    TTreeReader pythia_parts(partTree.Data(), &pythia_file);
    TTreeReaderValue<Double_t> pythia_weight(pythia_parts, "weight");
    TTreeReaderValue<Int_t> pythia_numparts(pythia_parts, "nparts");
    TTreeReaderArray<Double_t> pythia_px(pythia_parts, "particle_px");
    TTreeReaderArray<Double_t> pythia_py(pythia_parts, "particle_py");
    TTreeReaderArray<Double_t> pythia_pz(pythia_parts, "particle_pz");
    TTreeReaderArray<Double_t> pythia_E(pythia_parts, "particle_e");
    TTreeReaderArray<Int_t> pythia_id(pythia_parts, "particle_id");
    
    // TennGen event info
    TTreeReader tenngen_event(eventTree.Data(), &tenngen_file);
    TTreeReaderValue<Int_t> tenngen_nevent(tenngen_event, "nevent");
    TTreeReaderValue<Int_t> tenngen_centbin(tenngen_event, "centbin");
    TTreeReaderValue<Double_t> tenngen_etarange(tenngen_event, "etarange");
    
    // TennGen particles
    TTreeReader tenngen_parts(partTree.Data(), &tenngen_file);
    TTreeReaderValue<Int_t> tenngen_numparts(tenngen_parts, "nparts");
    TTreeReaderArray<Double_t> tenngen_px(tenngen_parts, "particle_px");
    TTreeReaderArray<Double_t> tenngen_py(tenngen_parts, "particle_py");
    TTreeReaderArray<Double_t> tenngen_pz(tenngen_parts, "particle_pz");
    TTreeReaderArray<Double_t> tenngen_E(tenngen_parts, "particle_e");
    TTreeReaderArray<Int_t> tenngen_id(tenngen_parts, "particle_id");
    TTreeReaderValue<Double_t> tenngen_psi1(tenngen_parts, "psi1");
    TTreeReaderValue<Double_t> tenngen_psi2(tenngen_parts, "psi2");
    TTreeReaderValue<Double_t> tenngen_psi3(tenngen_parts, "psi3");
    TTreeReaderValue<Double_t> tenngen_psi4(tenngen_parts, "psi4");
    
    // Read TH2D of particle vn's from TennGen file
    TH2D *v2_pi = (TH2D*)tenngen_file.Get("v2_pi");
    TH2D *v2_k = (TH2D*)tenngen_file.Get("v2_k");
    TH2D *v2_p = (TH2D*)tenngen_file.Get("v2_p");
    TH2D *v3_pi = (TH2D*)tenngen_file.Get("v3_pi");
    TH2D *v3_k = (TH2D*)tenngen_file.Get("v3_k");
    TH2D *v3_p = (TH2D*)tenngen_file.Get("v3_p");
    TH2D *v4_pi = (TH2D*)tenngen_file.Get("v4_pi");
    TH2D *v4_k = (TH2D*)tenngen_file.Get("v4_k");
    TH2D *v4_p = (TH2D*)tenngen_file.Get("v4_p");

    //////////////////////// Create output file //////////////////////////
    TFile *fout = new TFile(output_file.Data(), "RECREATE");
    
    // output variables
    Int_t event_id_number, event_cent_bin, event_pt_hard_bin, event_multiplicity, event_num_pythia_jets;
    Double_t event_weight, event_psi1, event_psi2, event_psi3, event_psi4, event_truth_jet_v2, event_truth_jet_v3, event_truth_jet_v4, event_pythia_phi_shift;
    Double_t particle_px[MAX_PARTICLES], particle_py[MAX_PARTICLES], particle_pz[MAX_PARTICLES], particle_E[MAX_PARTICLES];
    Double_t pythia_jets_original_phi[MAX_JETS], pythia_jets_pt[MAX_JETS], pythia_jets_realigned_phi[MAX_JETS];

    // Output tree
    TTree *outTree = new TTree("tree", "tree");
    outTree->Branch("event_id_number", &event_id_number, "event_id_number/I");
    outTree->Branch("event_cent_bin", &event_cent_bin, "event_cent_bin/I");
    outTree->Branch("event_pt_hard_bin", &event_pt_hard_bin, "event_pt_hard_bin/I");
    outTree->Branch("event_multiplicity", &event_multiplicity, "event_multiplicity/I");
    outTree->Branch("event_num_pythia_jets", &event_num_pythia_jets, "event_num_pythia_jets/I");
    outTree->Branch("event_weight", &event_weight, "event_weight/D");
    outTree->Branch("event_psi1", &event_psi1, "event_psi1/D");
    outTree->Branch("event_psi2", &event_psi2, "event_psi2/D");
    outTree->Branch("event_psi3", &event_psi3, "event_psi3/D");
    outTree->Branch("event_psi4", &event_psi4, "event_psi4/D");
    outTree->Branch("event_truth_jet_v2", &event_truth_jet_v2, "event_truth_jet_v2/D");
    outTree->Branch("event_truth_jet_v3", &event_truth_jet_v3, "event_truth_jet_v3/D");
    outTree->Branch("event_truth_jet_v4", &event_truth_jet_v4, "event_truth_jet_v4/D");
    outTree->Branch("event_pythia_phi_shift", &event_pythia_phi_shift, "event_pythia_phi_shift/D");
    outTree->Branch("particle_px", particle_px, "particle_px[event_multiplicity]/D");
    outTree->Branch("particle_py", particle_py, "particle_py[event_multiplicity]/D");
    outTree->Branch("particle_pz", particle_pz, "particle_pz[event_multiplicity]/D");
    outTree->Branch("particle_E", particle_E, "particle_E[event_multiplicity]/D");
    outTree->Branch("pythia_jets_original_phi", pythia_jets_original_phi, "pythia_jets_original_phi[event_num_pythia_jets]/D");
    outTree->Branch("pythia_jets_pt", pythia_jets_pt, "pythia_jets_pt[event_num_pythia_jets]/D");
    outTree->Branch("pythia_jets_realigned_phi", pythia_jets_realigned_phi, "pythia_jets_realigned_phi[event_num_pythia_jets]/D");
    
    ////////////////////////////////////////////////////////////////////////////////
    // check if the trees are empty
    if(!pythia_event.Next()){  cout << "PYTHIA tree is empty" << endl; exit(-1);}
    if(!tenngen_event.Next()){ cout << "TennGen tree is empty" << endl; exit(-1);}
    
    event_cent_bin = *tenngen_centbin;
    event_pt_hard_bin = *pythia_ptbin;
    Double_t pythia_xsec_over_eventweight = *pythia_xsec;
    Int_t nEvents_in_file = *pythia_nevent;
    if(*tenngen_nevent < *pythia_nevent) nEvents_in_file = *tenngen_nevent;
    Double_t  etaRange = *tenngen_etarange;

    // initialize fastjet
    Double_t particle_max_eta = 0.9;
    Double_t ghost_max_eta = 0.9;
    Double_t jet_min_pt_antikt = 10.0;
    Double_t antikt_jet_abs_eta_max = particle_max_eta - jetparam;
    fastjet::GhostedAreaSpec ghost_area_spec(ghost_max_eta);
    fastjet::JetDefinition antikt_jet_def(fastjet::antikt_algorithm, jetparam);


    ////////////////////////////////////////////////////////////////////////////////
    Int_t event_total;  
    if(nEvents < 0 || nEvents > nEvents_in_file) event_total = nEvents_in_file;
    else event_total = nEvents;
    cout << "Number of events to process: " << event_total << " for cent bin " << event_cent_bin << " and jet radius " << jetparam << endl;

    for(Int_t ievent = 0; ievent < event_total; ievent++){

            if(!tenngen_parts.Next()){ cout << "Reached end of tenngen file" << endl; break;}
            if(!pythia_parts.Next()){ cout << "Reached end of pythia file" << endl; break;}
            ////////////////////////////////////////////////////////////////////////////////
            // Fill the tree
            event_id_number = ievent;
            Int_t nParts_pythia = *pythia_numparts;
            Int_t nParts_tenngen = *tenngen_numparts;
            // Get event weight
            event_weight = *pythia_weight;
            event_weight *= pythia_xsec_over_eventweight;
            // Get event plane
            event_psi1 = *tenngen_psi1;
            event_psi2 = *tenngen_psi2;
            event_psi3 = *tenngen_psi3;
            event_psi4 = *tenngen_psi4;

            event_truth_jet_v2 = JetV2;
            event_truth_jet_v3 = JetV3;
            event_truth_jet_v4 = JetV4;
         
            // Get pythia particles
            event_multiplicity=0;
            std::vector<fastjet::PseudoJet> particles;
            for(Int_t ipart = 0; ipart < nParts_pythia; ipart++){   
                PseudoJet particle_temp(pythia_px[ipart], pythia_py[ipart], pythia_pz[ipart], pythia_E[ipart]);
                if(particle_temp.pt() < 0.15) continue; // remove particles with pt < 0.15 GeV
                if(( TMath::Abs(particle_temp.eta())  >  1.1 ) && !(particle_temp.eta() < 4.9 && particle_temp.eta() > 2.0)) continue;
                if(particle_temp.px() == 0 && particle_temp.py() == 0 && particle_temp.pz() == 0 && particle_temp.E() == 0) continue; // remove particles with 0 momentum
                
                particles.push_back(particle_temp);
                //particles.at(event_multiplicity).set_user_index(pythia_id[ipart]);    // set the user index to the particle id (not used in this analysis)
                event_multiplicity++; // increment the number of particles
            }

            // find jets in pythia event
            fastjet::AreaDefinition area_def(fastjet::active_area, ghost_area_spec);  
            fastjet::ClusterSequenceArea antikt_cs(particles, antikt_jet_def, area_def);
            std::vector<fastjet::PseudoJet> jets = antikt_cs.inclusive_jets(jet_min_pt_antikt);
            // cout << "Number of jets: " << jets.size() << endl;
            event_num_pythia_jets = 0;
            Double_t leading_jet_pt = 0.0;
            Double_t leading_jet_phi = 0.0;
            
            for(Int_t ijet = 0; ijet < jets.size(); ijet++){ // loop over antikt jets
            
                // cout << "Jet " << ijet << endl;
                // cout << "Jet " << ijet << " has constituents " << jets[ijet].constituents().size() << endl;
                // cout << "Jet " << ijet << " has pt " << jets[ijet].pt() << endl;
                // cout << "Jet " << ijet << " has area " << jets[ijet].area() << endl;
                // if(!jets[ijet].has_area()) continue; // check if the jet has an area
                //  cout << "Jet " << ijet << " has area " << jets[ijet].area() << endl;
                if(jets[ijet].constituents().size() == 0) continue; // check if the jet has constituents
                // cout << "Jet " << ijet << " has constituents " << jets[ijet].constituents().size() << endl;
                if(jets[ijet].pt() < jet_min_pt_antikt) continue; // check if the jet has pt > 20
                // cout << "Jet " << ijet << " has pt " << jets[ijet].pt() << endl;
                if(TMath::Abs(jets[ijet].eta()) > antikt_jet_abs_eta_max) continue; // check if the jet is in the eta range
                // fill the jet info
                pythia_jets_original_phi[event_num_pythia_jets] = jets[ijet].phi();
                pythia_jets_pt[event_num_pythia_jets] = jets[ijet].pt();

                // find leading jet
                if(pythia_jets_pt[event_num_pythia_jets] > leading_jet_pt){
                    leading_jet_pt = jets[ijet].pt();
                    leading_jet_phi = jets[ijet].phi();
                }
                // increment the number of jets
                event_num_pythia_jets++;

            }
            // clear the vector
            jets.clear();
            // skip events with no jets
            if(event_num_pythia_jets == 0) {continue;}

            Double_t jet_phi_new = -1.0;
            // find new phi for leading jet
            while(1){ // von Neumann rejection method
                Double_t dndphi_max = Harmonics_Phi_Dist_Rand->Uniform(0.0,MAX_DNDPHI);
                Double_t test_value = Harmonics_Phi_Dist_Rand->Uniform(0.0,2.0*TMath::Pi());
                Double_t test_phi = dNdPhi(test_value, event_psi2, event_psi3, event_psi4, JetV2, JetV3, JetV4);
                if(dndphi_max > test_value) continue;
                jet_phi_new = test_phi;
                break;
            }
            // check for errors
            if(jet_phi_new == -1.0) {cout << "ERROR: Could not find new phi for jet in event " << ievent << endl; return 0;}

            // calculate difference between old and new phi
            event_pythia_phi_shift = jet_phi_new - leading_jet_phi;
            // cout << "Jet phi shift " << event_pythia_phi_shift << endl;
            // cout << "Jet phi new " << jet_phi_new << endl;
            // cout << "Jet phi old " << leading_jet_phi << endl;
            // realign jet phis
            for (Int_t ijet = 0; ijet < event_num_pythia_jets; ijet++){
                Double_t new_phi = pythia_jets_original_phi[ijet] + event_pythia_phi_shift;
                if (new_phi < 0.0) new_phi += 2.0*TMath::Pi();
                if (new_phi > 2.0*TMath::Pi()) new_phi -= 2.0*TMath::Pi();
                pythia_jets_realigned_phi[ijet] = new_phi;
                // cout << "New phi " << new_phi << endl;
                // cout << "Jet " << ijet << " has phi " << pythia_jets_original_phi[ijet] << " and new phi " << pythia_jets_realigned_phi[ijet] << endl;
            }

            // rotate pythia particles
            for (Int_t ipart = 0; ipart < event_multiplicity; ipart++){

                Double_t phi_temp = particles.at(ipart).phi() + event_pythia_phi_shift;
                if(phi_temp < 0.0) phi_temp += 2.0*TMath::Pi();
                if(phi_temp > 2.0*TMath::Pi()) phi_temp -= 2.0*TMath::Pi(); 
                
                particle_px[ipart] = (particles.at(ipart).pt())*TMath::Cos(phi_temp);
                particle_py[ipart] = (particles.at(ipart).pt())*TMath::Sin(phi_temp);
                particle_pz[ipart] = (particles.at(ipart).pt())*TMath::SinH(particles.at(ipart).eta());
                particle_E[ipart] = particles.at(ipart).E();
            }
            // clear particles vector
            particles.clear();

            // Get TennGen Particles
            for(Int_t ipart = 0; ipart < nParts_tenngen; ipart++){
                PseudoJet particle_temp(tenngen_px[ipart], tenngen_py[ipart], tenngen_pz[ipart], tenngen_E[ipart]);
                if(particle_temp.pt() < 0.15) continue; // remove particles with pt < 0.15 GeV
                if(( TMath::Abs(particle_temp.eta())  >  1.1 ) && !(particle_temp.eta() < 4.9 && particle_temp.eta() > 2.0)) continue;
                // if(TMath::Abs(particle_temp.eta()) > etaRange) continue; // remove particles outside of eta range
                if(particle_temp.px() == 0 && particle_temp.py() == 0 && particle_temp.pz() == 0 && particle_temp.E() == 0) continue; // remove particles with 0 momentum
               
                // fill the particle info
                particle_px[event_multiplicity] = particle_temp.px();
                particle_py[event_multiplicity] = particle_temp.py();
                particle_pz[event_multiplicity] = particle_temp.pz();
                particle_E[event_multiplicity] = particle_temp.E();
                event_multiplicity++;

            }
            // fill the tree
            outTree->Fill();

            // clear the arrays
            for(Int_t i = 0; i < event_multiplicity; i++){
                particle_px[i] = 0.0;
                particle_py[i] = 0.0;
                particle_pz[i] = 0.0;
                particle_E[i] = 0.0;
            }

            for(Int_t i =0; i < event_num_pythia_jets; i++){
                pythia_jets_pt[i] = 0.0;
                pythia_jets_original_phi[i] = 0.0;
                pythia_jets_realigned_phi[i] = 0.0;
            }
    } //end of event/cent loop

    fout->cd();
    v2_pi->Write();
    v3_pi->Write();
    v4_pi->Write();
    v2_k->Write();
    v3_k->Write();
    v4_k->Write();
    v2_p->Write();
    v3_p->Write();
    v4_p->Write();
    fout->Write();
    fout->Close();
    // Close the input file
    tenngen_file.Close();
    pythia_file.Close();
    // Write the tree to the output file  
    return 0;
}


Int_t main(Int_t argc, Char_t** argv){


    if(argc < 5){
        cout << "Usage: ./RealignJets <input pythia file> <input tenngen file> <prefix> <jet resolution parameter> <nEvents>" << endl;
        return 0;
    }
    cout << "Starting" << endl;
    TString pythia_file_name = argv[1];
    TString tenngen_file_name = argv[2];
    TString prefix = argv[3];
    Double_t jetparam = atof(argv[4]);
    Double_t nEvents = argc > 5 ? atoi(argv[5]) : -1;

    // get Pwd
    TString pwd = gSystem->pwd();
    cout << "Current directory: " << pwd << endl;
    // if(!pythia_file_name.Contains(pwd.Data()))  pythia_file_name = pwd + "/" + pythia_file_name;
    // if(!tenngen_file_name.Contains(pwd.Data()))  tenngen_file_name = pwd + "/" + tenngen_file_name;
    // if(prefix.EndsWith("/") && prefix.Length() > 1) prefix.Remove(prefix.Length()-1);
    // if(!prefix.Contains(pwd.Data()))  prefix = pwd + "/" + prefix;
    // check that the input files exist
    // if(!gSystem->AccessPathName(pythia_file_name.Data())){
    //     cout << "Pythia file " << pythia_file_name << " does not exist" << endl;
    //     return 0;
    // }
    // if(!gSystem->AccessPathName(tenngen_file_name.Data())){
    //     cout << "TennGen file " << tenngen_file_name << " does not exist" << endl;
    //     return 0;
    // }

    // check that prefix is valid file location, if not try to create it
    if(gSystem->AccessPathName(prefix)){
        cout << "Prefix " << prefix << " does not exist, trying to create it" << endl;
        gSystem->mkdir(prefix);
        if(gSystem->AccessPathName(prefix)){
            cout << "Could not create prefix " << prefix << endl;
            return 0;
        }
    }

    TString ptbin = tenngen_file_name(tenngen_file_name.Index("ptbin") + 5, 2);
    if (ptbin.Contains("_")) ptbin.Remove(ptbin.Index("_"));
    if (ptbin.Contains(".")) ptbin.Remove(ptbin.Index("."));

    TString output_file_name = Form("%s/200GeV_MixedEvents_ptbin%s_RealignJets_R0%.0f.root", prefix.Data(), ptbin.Data(), jetparam*10.0);

    // print out the input parameters
    cout << "Realigning Pythia file " << pythia_file_name << " with TennGen file " << tenngen_file_name << endl;
    cout << "Output file name: " << output_file_name << endl;
    cout << "Jet resolution parameter: " << jetparam << endl;
    if (nEvents > 0) cout << "Number of events: " << nEvents << endl;
    else cout << "Number of events: all" << endl;

    cout << "Starting to process events" << endl;
    // process the events
    

    RealignJets(pythia_file_name, tenngen_file_name, output_file_name, jetparam, nEvents);

    cout << "Finished processing events" << endl;
    return 0;

}
