// ROOT includes
#include <TFile.h>
#include <TTree.h>
#include <TMath.h>
#include <TF1.h>
#include <TH1D.h>
#include <TRandom3.h>
#include <TTimeStamp.h>
 
// Pythia includes
#include <Pythia8/Pythia.h>

// FastJet includes
#include <fastjet/ClusterSequence.hh>
#include <fastjet/PseudoJet.hh>
#include <fastjet/ClusterSequenceArea.hh>
#include <fastjet/config.h>

// Custom includes
#include "Settings.h"
#include "BkgdFunctions.h"
#include "BkgdGen.h"

// C++ includes
#include <iostream>
#include <vector>

namespace toymodel  
{


void BkgdGen::Run()
{

    // ================================================================
    // ===================== CONFIGURATION ============================
    // ================================================================
    ToyModelSettings toymodel_settings;
    toymodel_settings.read(m_config_filename);

    if(toymodel_settings.Main.RandomSeed == 0 )
    {
        TTimeStamp *time = new TTimeStamp();
        toymodel_settings.Main.RandomSeed = static_cast<unsigned int>(time->GetNanoSec());
        delete time;
    }
    if (toymodel_settings.Main.RandomSeed > 900000000) toymodel_settings.Main.RandomSeed = toymodel_settings.Main.RandomSeed % 900000000;

    
    if(toymodel_settings.Main.Verbosity > 0)
    {
        toymodel_settings.print(std::cout);
    }



    double v2_ref_truth = toymodel::TruthRefVn(0, toymodel_settings.Main.Particle.MinPt);
    double v3_ref_truth = toymodel::TruthRefVn(1, toymodel_settings.Main.Particle.MinPt);
    double v4_ref_truth = toymodel::TruthRefVn(2, toymodel_settings.Main.Particle.MinPt);
    //============================================================
    //==================== OUTPUT CONFIGURATION ===================
    //============================================================
    // create output file
    TFile *outFile = new TFile(m_output_filename.c_str(),"RECREATE");
    
    // create output trees
    TTree * config_tree = new TTree("config", "config");
    TTree * event_tree = new TTree("tree", "tree");
   
    // config branches
    config_tree->Branch("nEvents", &toymodel_settings.Main.NumEvents, "nEvents/I");
    config_tree->Branch("seed", &toymodel_settings.Main.RandomSeed , "seed/i");
    config_tree->Branch("max_particle_eta", & toymodel_settings.Main.Particle.MaxEta, "max_particle_eta/F");
    config_tree->Branch("min_particle_pt", & toymodel_settings.Main.Particle.MinPt, "min_particle_pt/F");
    config_tree->Branch("v2_ref_truth", &v2_ref_truth, "v2_ref_truth/D");
    config_tree->Branch("v3_ref_truth", &v3_ref_truth, "v3_ref_truth/D");
    config_tree->Branch("v4_ref_truth", &v4_ref_truth, "v4_ref_truth/D");

    // output variables (bkdg)
    int event_id;
    double weight;
    double psi1, psi2, psi3, psi4; // event plane angles    
    double average_pt;

    int n_forward_particles; // forward rapidity multiplicity
    double Q2_i, Q2_r; // q-vector real and imaginary parts
    double Q3_i, Q3_r; // q-vector real and imaginary parts
    double Q4_i, Q4_r; // q-vector real and imaginary parts
    double Q6_i, Q6_r; // q-vector real and imaginary parts
    double Q8_i, Q8_r; // q-vector real and imaginary parts
    
    int n_midrapidity_particles; // midrapidity multiplicity
    std::vector<double> part_px{};
    std::vector<double> part_py{};
    std::vector<double> part_pz{};
    std::vector<double> part_e{};

    // set up output branches
    event_tree->Branch("event_id", &event_id, "event_id/I");
    event_tree->Branch("weight", &weight, "weight/D");

    event_tree->Branch("psi1", &psi1, "psi1/D");
    event_tree->Branch("psi2", &psi2, "psi2/D");
    event_tree->Branch("psi3", &psi3, "psi3/D");
    event_tree->Branch("psi4", &psi4, "psi4/D");
    event_tree->Branch("average_pt", &average_pt, "average_pt/D");

    event_tree->Branch("n_forward_particles", &n_forward_particles, "n_forward_particles/I");
    event_tree->Branch("Q2_i", &Q2_i, "Q2_i/D");
    event_tree->Branch("Q2_r", &Q2_r, "Q2_r/D");
    event_tree->Branch("Q3_i", &Q3_i, "Q3_i/D");
    event_tree->Branch("Q3_r", &Q3_r, "Q3_r/D");
    event_tree->Branch("Q4_i", &Q4_i, "Q4_i/D");
    event_tree->Branch("Q4_r", &Q4_r, "Q4_r/D");
    event_tree->Branch("Q6_i", &Q6_i, "Q6_i/D");
    event_tree->Branch("Q6_r", &Q6_r, "Q6_r/D");
    event_tree->Branch("Q8_i", &Q8_i, "Q8_i/D");
    event_tree->Branch("Q8_r", &Q8_r, "Q8_r/D");
    
    event_tree->Branch("n_midrapidity_particles", &n_midrapidity_particles, "n_midrapidity_particles/I");
    event_tree->Branch("part_px", &part_px);
    event_tree->Branch("part_py", &part_py);
    event_tree->Branch("part_pz", &part_pz);
    event_tree->Branch("part_e", &part_e);


    
    //================================================================
    //==================== BKGD CONFIGURATION =======================
    //================================================================
    // initialize pT distributions
    if(toymodel_settings.Main.Verbosity > 0) std::cout << "Configuring Bkgd distributions" << std::endl;
    
    TF1 * f1_piPluspT = new TF1("f1_piPluspT", toymodel::BGBW,  toymodel_settings.Main.Particle.MinPt, 50, 5);
    TF1 * f1_piMinuspT = new TF1("f1_piMinuspT", toymodel::BGBW,  toymodel_settings.Main.Particle.MinPt, 50, 5);
    TF1 * f1_KPluspT = new TF1("f1_KPluspT", toymodel::BGBW,  toymodel_settings.Main.Particle.MinPt, 50, 5);
    TF1 * f1_KMinuspT = new TF1("f1_KMinuspT", toymodel::BGBW,  toymodel_settings.Main.Particle.MinPt, 50, 5);
    TF1 * f1_protonpT = new TF1("f1_protonpT", toymodel::BGBW,  toymodel_settings.Main.Particle.MinPt, 50, 5);
    TF1 * f1_pbarpT = new TF1("f1_pbarpT", toymodel::BGBW,  toymodel_settings.Main.Particle.MinPt, 50, 5);

    // initialize vn functions
    std::cout << "Configuring vn functions" << std::endl;
    TF1 * v2_pi =  new TF1("v2_pi", toymodel::VnFunction,  toymodel_settings.Main.Particle.MinPt, 50, 7);
    TF1 * v3_pi =  new TF1("v3_pi", toymodel::VnFunction,  toymodel_settings.Main.Particle.MinPt, 50, 7);
    TF1 * v4_pi =  new TF1("v4_pi", toymodel::VnFunction,  toymodel_settings.Main.Particle.MinPt, 50, 7);
    TF1 * v2_k =  new TF1("v2_k", toymodel::VnFunction,  toymodel_settings.Main.Particle.MinPt, 50, 7);
    TF1 * v3_k =  new TF1("v3_k", toymodel::VnFunction,  toymodel_settings.Main.Particle.MinPt, 50, 7);
    TF1 * v4_k =  new TF1("v4_k", toymodel::VnFunction,  toymodel_settings.Main.Particle.MinPt, 50, 7);
    TF1 * v2_pro =  new TF1("v2_pro", toymodel::VnFunction,  toymodel_settings.Main.Particle.MinPt, 50, 7);
    TF1 * v3_pro =  new TF1("v3_pro", toymodel::VnFunction,  toymodel_settings.Main.Particle.MinPt, 50, 7);
    TF1 * v4_pro =  new TF1("v4_pro", toymodel::VnFunction,  toymodel_settings.Main.Particle.MinPt, 50, 7);
    
    std::vector<TF1*> pT_distos = { f1_piPluspT,
                                    f1_piMinuspT, 
                                    f1_KPluspT, 
                                    f1_KMinuspT, 
                                    f1_protonpT, 
                                    f1_pbarpT};
    for (auto f : pT_distos) f->SetNpx(1000);

    for (int i=0;i<6;i++)
    {
        for (int j=0;j<5;j++)
        {
            pT_distos.at(i)->SetParameter(j,toymodel::blastwave_params[j][i]);
        }
    }
  
    std::vector<TF1*> v2_distos = {v2_pi, 
                                   v2_k, 
                                   v2_pro};
    std::vector<TF1*> v3_distos = {v3_pi, 
                                   v3_k, 
                                   v3_pro};
    std::vector<TF1*> v4_distos = {v4_pi, 
                                   v4_k, 
                                   v4_pro};
    for (auto f : v2_distos) f->SetNpx(1000);
    for (auto f : v3_distos) f->SetNpx(1000);
    for (auto f : v4_distos) f->SetNpx(1000);
    for (int i=0;i<3;i++)
    {
        if(i == 0)
        {
            for (int j=0;j<7;j++)
            {
                v2_distos.at(i)->SetParameter(j,toymodel::pi_vn_params[0][j]);
                v3_distos.at(i)->SetParameter(j,toymodel::pi_vn_params[1][j]);
                v4_distos.at(i)->SetParameter(j,toymodel::pi_vn_params[2][j]);
            }
        }
        else if(i == 1)
        {
            for (int j=0;j<7;j++)
            {
                v2_distos.at(i)->SetParameter(j,toymodel::k_vn_params[0][j]);
                v3_distos.at(i)->SetParameter(j,toymodel::k_vn_params[1][j]);
                v4_distos.at(i)->SetParameter(j,toymodel::k_vn_params[2][j]);
            }
        }
        else if(i == 2)
        {
            for (int j=0;j<7;j++)
            {
                v2_distos.at(i)->SetParameter(j,toymodel::pro_vn_params[0][j]);
                v3_distos.at(i)->SetParameter(j,toymodel::pro_vn_params[1][j]);
                v4_distos.at(i)->SetParameter(j,toymodel::pro_vn_params[2][j]);
            }
        }
       
    }

    // initialize random number generators
    TRandom3 *eta_rand = new TRandom3(static_cast<unsigned int>(toymodel_settings.Main.RandomSeed));
    TRandom3 *phi_rand = new TRandom3(static_cast<unsigned int>(toymodel_settings.Main.RandomSeed*7));
    TRandom3 *pt_rand = new TRandom3(static_cast<unsigned int>(toymodel_settings.Main.RandomSeed*3));
    TRandom3 *psi_rand = new TRandom3(static_cast<unsigned int>(toymodel_settings.Main.RandomSeed*5));

    
    TH1D * h1_tenngen_part_pt = new TH1D("h1_tenngen_part_pt", "h1_tenngen_part_pt", 100, 0, 100);
    TH1D * h1_tenngen_part_eta = new TH1D("h1_tenngen_part_eta", "h1_tenngen_part_eta", 100, -5, 5);
    TH1D * h1_tenngen_part_phi = new TH1D("h1_tenngen_part_phi", "h1_tenngen_part_phi", 100, -TMath::Pi(), TMath::Pi());
    TH1D * h1_tenngen_forward_part_phi = new TH1D("h1_tenngen_forward_part_phi", "h1_tenngen_forward_part_phi", 100, -TMath::Pi(), TMath::Pi());
    TH1D * h1_tenngen_forward_part_eta = new TH1D("h1_tenngen_forward_part_eta", "h1_tenngen_forward_part_eta", 100, -5, 5);


    // ================================================================
    // ===================== EVENT LOOP ===============================
    // ================================================================
    if(toymodel_settings.Main.Verbosity > 0) std::cout << "Starting event loop" << std::endl;
    int progess_update = int(toymodel_settings.Main.NumEvents / 10);
    int progress = 0;

    // loop over events
    int iEvent = 0;
    while(iEvent < toymodel_settings.Main.NumEvents)
    {
   
        // reset event variables
        event_id = iEvent;
        weight = 0.0;

        // get psi
        if(toymodel_settings.Bkgd.ConstEventPlane1 < 0){ psi1 = psi_rand->Uniform(0,2*TMath::Pi()); }
        else { psi1 = toymodel_settings.Bkgd.ConstEventPlane1; }
        if(toymodel_settings.Bkgd.ConstEventPlane2 < 0){ psi2 = psi_rand->Uniform(0,2*TMath::Pi()); }
        else { psi2 = toymodel_settings.Bkgd.ConstEventPlane2; }
        if(toymodel_settings.Bkgd.ConstEventPlane3 < 0){ psi3 = psi_rand->Uniform(0,2*TMath::Pi()); }
        else { psi3 = toymodel_settings.Bkgd.ConstEventPlane3; }
        if(toymodel_settings.Bkgd.ConstEventPlane4 < 0){ psi4 = psi_rand->Uniform(0,2*TMath::Pi()); }
        else { psi4 = toymodel_settings.Bkgd.ConstEventPlane4; }


        // clear variables
        average_pt = 0.0;
        n_forward_particles = 0;
        Q2_i = 0.0;
        Q2_r = 0.0;
        Q3_i = 0.0;
        Q3_r = 0.0;
        Q4_i = 0.0;
        Q4_r = 0.0;
        Q6_i = 0.0;
        Q6_r = 0.0;
        Q8_i = 0.0;
        Q8_r = 0.0;

        n_midrapidity_particles = 0;
        part_px.clear();
        part_py.clear();
        part_pz.clear();
        part_e.clear();

        // loop over particles
        int M = 0;
        double rho = 0.0;
        for( unsigned int ispecies = 0; ispecies < toymodel::n_species; ispecies++)
        {
            double mass = toymodel::particle_masses.at(ispecies);
            int KF = toymodel::particle_ids.at(ispecies);
            TF1 *pT_distro = pT_distos.at(ispecies);
            unsigned int vn_idx = static_cast<unsigned int>(ispecies/2);
            for (unsigned int ipart = 0; ipart < toymodel::particle_yeilds.at(ispecies); ipart++)
            {
                double pT = pT_distro->GetRandom(pt_rand); // get pT
                double eta = eta_rand->Uniform(- toymodel_settings.Main.Particle.MaxEta, toymodel_settings.Main.Particle.MaxEta); // get eta

                // vns
                double v2 = v2_distos.at(vn_idx)->Eval(pT);
                double v3 = v3_distos.at(vn_idx)->Eval(pT);
                double v4 = v4_distos.at(vn_idx)->Eval(pT);
                double v1 = v2 - 0.02;

                

                
                double tmp_dndpi = 1.0 + 2.0*(TMath::Abs(v1) + TMath::Abs(v2) + TMath::Abs(v3) + TMath::Abs(v4));
                double phi = 0.0;
                while(true)
                {
                    double y1 = phi_rand->Uniform(0.0,tmp_dndpi);
                    double x = phi_rand->Uniform(0.0, 2*TMath::Pi());
                    double y2 = 1.0+2.0*v1*TMath::Cos(x-psi1)+2.0*v2*TMath::Cos(2.0*(x-psi2))+2.0*v3*TMath::Cos(3.0*(x-psi3))+2.0*v4*TMath::Cos(4.0*(x-psi4));
                    
                    if(y1 < y2){ phi = x; break; }
                }
                // transform phi to -pi,pi
                phi-=TMath::Pi();
                if(phi < -TMath::Pi()) phi += 2.0*TMath::Pi();
                if(phi > TMath::Pi()) phi -= 2.0*TMath::Pi();


                // transform to cartesian
                double px = pT*TMath::Cos(phi);
                double py = pT*TMath::Sin(phi);
                double pz = pT*TMath::SinH(eta);
                double e = TMath::Sqrt(px*px + py*py + pz*pz + mass*mass);
                
                n_midrapidity_particles++;
                part_px.push_back(px);
                part_py.push_back(py);
                part_pz.push_back(pz);
                part_e.push_back(e);

               
                M++;
                rho+=pT;
                h1_tenngen_part_pt->Fill(pT);
                h1_tenngen_part_eta->Fill(eta);
                h1_tenngen_part_phi->Fill(phi);



                // do forward rapidity particles
                pT = pT_distro->GetRandom(pt_rand); // get pT
                v2 = v2_distos.at(vn_idx)->Eval(pT);
                v3 = v3_distos.at(vn_idx)->Eval(pT);
                v4 = v4_distos.at(vn_idx)->Eval(pT);
                v1 = v2 - 0.02;

                tmp_dndpi = 1.0 + 2.0*(TMath::Abs(v1) + TMath::Abs(v2) + TMath::Abs(v3) + TMath::Abs(v4));
                phi = 0.0;
                while(true)
                {
                    double y1 = phi_rand->Uniform(0.0,tmp_dndpi);
                    double x = phi_rand->Uniform(0.0,2*TMath::Pi());
                    double y2 = 1.0+2.0*v1*TMath::Cos(x-psi1)+2.0*v2*TMath::Cos(2.0*(x-psi2))+2.0*v3*TMath::Cos(3.0*(x-psi3))+2.0*v4*TMath::Cos(4.0*(x-psi4));
                    
                    if(y1 < y2){ phi = x; break; }
                }

                // transform phi to -pi,pi
                phi-=TMath::Pi();
                if(phi < -TMath::Pi()) phi += 2.0*TMath::Pi();
                if(phi > TMath::Pi()) phi -= 2.0*TMath::Pi();


                Q2_i += TMath::Sin(2.0*phi);
                Q2_r += TMath::Cos(2.0*phi);
                Q3_i += TMath::Sin(3.0*phi);
                Q3_r += TMath::Cos(3.0*phi);
                Q4_i += TMath::Sin(4.0*phi);
                Q4_r += TMath::Cos(4.0*phi);
                Q6_i += TMath::Sin(6.0*phi);
                Q6_r += TMath::Cos(6.0*phi);
                Q8_i += TMath::Sin(8.0*phi);
                Q8_r += TMath::Cos(8.0*phi);
                n_forward_particles++;

                h1_tenngen_forward_part_phi->Fill(phi);
                h1_tenngen_forward_part_eta->Fill(eta);


            } // end of particle loop   
        
        } 

        // // get average pt density
        rho /= M;
        average_pt = rho;

        // ================================================================
        // ===================== FILL TREE =================================
        // ================================================================

       // fill tree
        event_tree->Fill();


        iEvent++;
        // print progress
        if (iEvent % progess_update == 0)
        {
            std::cout << "Progress: " << progress << "%" << std::endl;
            progress += 10;
        }

    } // end of event loop
    std::cout << "Progress: 100%" << std::endl;
    std::cout << "Event loop finished" << std::endl;
    // fill event info tree
    config_tree->Fill();
    
    //============================================================
    //==================== WRITE OUTPUT ==========================
    //============================================================

    // write output file
    outFile->cd();
    config_tree->Write();
    event_tree->Write();
    h1_tenngen_part_pt->Write();
    h1_tenngen_part_eta->Write();
    h1_tenngen_part_phi->Write();
    h1_tenngen_forward_part_phi->Write();
    h1_tenngen_forward_part_eta->Write();

    // Convert TF1 to TH1D
    for (auto f : pT_distos)
    {
        TH1D * h = (TH1D*)f->GetHistogram();
        h->SetName(Form("%s_hist", f->GetName()));
        h->Write();
    }
    for (auto f : v2_distos)
    {
        TH1D * h = (TH1D*)f->GetHistogram();
        h->SetName(Form("%s_hist", f->GetName()));
        h->Write();
    }
    for (auto f : v3_distos)
    {
        TH1D * h = (TH1D*)f->GetHistogram();
        h->SetName(Form("%s_hist", f->GetName()));
        h->Write();
    }
    for (auto f : v4_distos)
    {
        TH1D * h = (TH1D*)f->GetHistogram();
        h->SetName(Form("%s_hist", f->GetName()));
        h->Write();
    }

    outFile->Close();
    std::cout << "Done" << std::endl;
    return ;

}

}