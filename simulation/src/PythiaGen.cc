// ROOT includes
#include <TFile.h>
#include <TTree.h>
#include <TMath.h>
#include <TRandom3.h>
#include <TTimeStamp.h>
#include <TH1D.h>
 
// Pythia includes
#include <Pythia8/Pythia.h>

// FastJet includes
#include <fastjet/ClusterSequence.hh>
#include <fastjet/PseudoJet.hh>
#include <fastjet/config.h>

// Custom includes
#include "Settings.h"
#include "PythiaGen.h"

// C++ includes
#include <iostream>
#include <vector>
#include <string>
#include <cstdio>

namespace toymodel  
{
   
    bool PythiaGen::accept_event(std::vector<fastjet::PseudoJet> &jets, double pTmin, double pTmax)
    {

            if(pTmin  < 0 ) return true;
            if(pTmax < 0) pTmax = 1000000; // if pTmax is negative, set to large number

            fastjet::JetDefinition jetdef(fastjet::antikt_algorithm, 0.4, fastjet::E_scheme, fastjet::Best);
            fastjet::ClusterSequence clust_seq(jets, jetdef);
            std::vector<fastjet::PseudoJet> inclusive_jets = fastjet::sorted_by_pt(clust_seq.inclusive_jets(0));
            
            double leading_jet_pt = 0;
            for (auto jet : jets)
            {
                if (jet.pt() > leading_jet_pt) leading_jet_pt = jet.pt();
            }
            
            if(leading_jet_pt < pTmin || leading_jet_pt > pTmax)
            {
                return false;
            }


            return true;
    }



    void PythiaGen::Run()
    {

        // ================================================================
        // ===================== CONFIGURATION ============================
        // ================================================================
        if (m_output_filename == "")
        {
            std::cerr << "No output file provided" << std::endl;
            exit(1);
        }

        if (m_config_filename == "")
        {
            std::cerr << "No config file provided" << std::endl;
            exit(1);
        }

        ToyModelSettings toymodel_settings;
        toymodel_settings.read(m_config_filename);

        if(toymodel_settings.Main.RandomSeed == 0 )
        {
            TTimeStamp *time = new TTimeStamp();
            toymodel_settings.Main.RandomSeed = static_cast<unsigned int>(time->GetNanoSec());
            delete time;
        }
        if (toymodel_settings.Main.RandomSeed > 900000000) toymodel_settings.Main.RandomSeed = toymodel_settings.Main.RandomSeed % 900000000;

        
        if(toymodel_settings.Main.Jet.MaxEta  < 0) toymodel_settings.Main.Jet.MaxEta  = toymodel_settings.Main.Particle.MaxEta - toymodel_settings.Main.Jet.R;
        if (toymodel_settings.Signal.Jet.MinPt < 0) toymodel_settings.Signal.Jet.MinPt = toymodel_settings.Main.Jet.MinPt;

        std::string pythia_config_file = toymodel_settings.Signal.Pythia.Config;
        std::vector<std::string> pythia_commands = toymodel_settings.Signal.Pythia.Commands;

        if(toymodel_settings.Main.Verbosity > 0)
        {
            toymodel_settings.print(std::cout);
        }


        double xsec_over_eventweight = 0.0;
        int n_accepted_events = 0;
        int n_generated_events = 0;
        int n_processed_events = 0;
        double sum_of_weights = 0.0;
        double integrated_luminosity = 0.0;

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
        config_tree->Branch("R", &toymodel_settings.Main.Jet.R, "R/F");
        config_tree->Branch("max_jet_eta", &toymodel_settings.Main.Jet.MaxEta, "max_jet_eta/F");
        config_tree->Branch("min_pythia_jet_pt", &toymodel_settings.Signal.Jet.MinPt, "min_pythia_jet_pt/F");
        config_tree->Branch("pT_hard_min", &toymodel_settings.Signal.Pythia.PtHardMin, "pT_hard_min/F");
        config_tree->Branch("pT_hard_max", &toymodel_settings.Signal.Pythia.PtHardMax, "pT_hard_max/F");
        config_tree->Branch("xsec_over_eventweight", &xsec_over_eventweight, "xsec_over_eventweight/D");
        config_tree->Branch("n_accepted_events", &n_accepted_events, "n_accepted_events/I");
        config_tree->Branch("n_generated_events", &n_generated_events, "n_generated_events/I");
        config_tree->Branch("sum_of_weights", &sum_of_weights, "sum_of_weights/D");
        config_tree->Branch("integrated_luminosity", &integrated_luminosity, "integrated_luminosity/D");

        // output variables (bkdg)
        int event_id;
        double weight;

        // single truth jet variables
        double jet_pt;
        double jet_phi;
        double jet_eta;
        int jet_nconst;
        std::vector<double> const_pt{};
        std::vector<double> const_eta{};
        std::vector<double> const_delta_phi{};
        std::vector<double> const_pz{};
        std::vector<double> const_e{};

        // set up output branches
        event_tree->Branch("event_id", &event_id, "event_id/I");
        event_tree->Branch("weight", &weight, "weight/D");
        event_tree->Branch("jet_pt", &jet_pt, "jet_pt/D");
        event_tree->Branch("jet_phi", &jet_phi, "jet_phi/D");
        event_tree->Branch("jet_eta", &jet_eta, "jet_eta/D");
        event_tree->Branch("jet_nconst", &jet_nconst, "jet_nconst/I");
        event_tree->Branch("const_pt", &const_pt);
        event_tree->Branch("const_eta", &const_eta);
        event_tree->Branch("const_delta_phi", &const_delta_phi);
        event_tree->Branch("const_pz", &const_pz);
        event_tree->Branch("const_e", &const_e);

        TH1D * h1_jet_const_delta_phi  = new TH1D("h1_jet_const_delta_phi", "h1_jet_const_delta_phi", 100, -TMath::Pi(), TMath::Pi());
        TH1D * h1_jet_const_deltaR  = new TH1D("h1_jet_const_deltaR", "h1_jet_const_deltaR", 100, 0, 1);
        TH1D * h1_pythia_part_pt = new TH1D("h1_pythia_part_pt", "h1_pythia_part_pt", 100, 0, 100);
        TH1D * h1_pythia_part_eta = new TH1D("h1_pythia_part_eta", "h1_pythia_part_eta", 100, -5, 5);   
        TH1D * h1_pythia_part_phi = new TH1D("h1_pythia_part_phi", "h1_pythia_part_phi", 100, -TMath::Pi(), TMath::Pi());
        TH1D * h1_pythia_jet_pt = new TH1D("h1_pythia_jet_pt", "h1_pythia_jet_pt", 100, 0, 100);
        TH1D * h1_pythia_jet_eta = new TH1D("h1_pythia_jet_eta", "h1_pythia_jet_eta", 100, -5, 5);
        TH1D * h1_pythia_jet_phi = new TH1D("h1_pythia_jet_phi", "h1_pythia_jet_phi", 100, -TMath::Pi(), TMath::Pi());

        //============================================================
        //==================== PYTHIA CONFIGURATION ==================
        //============================================================
        if(toymodel_settings.Main.Verbosity > 0) std::cout << "Configuring Pythia" << std::endl;

        Pythia8::Pythia pythia;
        Pythia8::Settings& settings = pythia.settings;
        const Pythia8::Info& info = pythia.info;
        Pythia8::Event& event = pythia.event;
        if(pythia_config_file != "")
        {
            pythia.readFile(pythia_config_file);
        }
        for (auto command : pythia_commands)
        {
            pythia.readString(command);
        }
        if(toymodel_settings.Signal.Pythia.PtHardMin > 0)
        {
            settings.parm("PhaseSpace:pTHatMin", toymodel_settings.Signal.Pythia.PtHardMin);
            if(toymodel_settings.Signal.Pythia.PtHardMax > 0)
            {
                settings.parm("PhaseSpace:pTHatMax", toymodel_settings.Signal.Pythia.PtHardMax);
            }
        }
        settings.readString("Random:setSeed = on");
        settings.parm("Random:seed", toymodel_settings.Main.RandomSeed);
        pythia.init();
    
        // ================================================================
        // ===================== FASTJET CONFIG  =============================
        // ================================================================
        
        // jet definition
        fastjet::JetDefinition jet_def_antikt(fastjet::antikt_algorithm, toymodel_settings.Main.Jet.R, fastjet::E_scheme, fastjet::Best);

        // pythia jet selector
        fastjet::Selector pythia_jet_select = fastjet::SelectorNHardest(1) * (fastjet::SelectorAbsEtaMax(toymodel_settings.Main.Jet.MaxEta) && fastjet::SelectorPtMin(toymodel_settings.Signal.Jet.MinPt));

        // ================================================================
        // ===================== EVENT LOOP ===============================
        // ================================================================
        if(toymodel_settings.Main.Verbosity > 0) std::cout << "Starting event loop" << std::endl;
        int progess_update = int(toymodel_settings.Main.NumEvents / 10);
        int progress = 0;

        // loop over events
        int iEvent = 0;
        std::vector<fastjet::PseudoJet> pythia_particles {};
        while(iEvent < toymodel_settings.Main.NumEvents)
        {
    
            // reset event variables
            event_id = iEvent;
            weight = 0.0;
            jet_pt = 0.0;
            jet_phi = 0.0;
            jet_eta = 0.0;
            jet_nconst = 0;
            const_pt.clear();
            const_eta.clear();
            const_delta_phi.clear();
            const_pz.clear();
            const_e.clear();

            // ================================================================
            // ===================== PYTHIA EVENT =============================
            // ================================================================
            bool acceptEvent = false;
            while(!acceptEvent)
            {   
                
                
                while (!pythia.next()) 
                {
                    n_generated_events++;
                }
                n_generated_events++;
                
                // get event info
                weight = info.weight();
                pythia_particles.clear();
                for (unsigned int ipart = 0; ipart < event.size(); ipart++)
                {

                    
                    // only consider final state charged particles
                    if(!event[ipart].isFinal())
                    {
                        continue;
                    }
                    // only consider charged particles
                    if(!event[ipart].isCharged())
                    {
                        continue;
                    }

                    // min pT cut
                    if(event[ipart].pT() <  toymodel_settings.Main.Particle.MinPt)
                    {
                        continue;
                    }

                    // abs(eta) cut
                    if(TMath::Abs(event[ipart].eta()) >  toymodel_settings.Main.Particle.MaxEta)
                    {
                        continue;
                    }


                    // push back particle
                    pythia_particles.push_back(fastjet::PseudoJet(event[ipart].px(), event[ipart].py(), event[ipart].pz(), event[ipart].e()));
                }

                
                acceptEvent = accept_event(pythia_particles, toymodel_settings.Signal.Pythia.PtHardMin, toymodel_settings.Signal.Pythia.PtHardMax);
            }


            // ================================================================
            // ===================== PYTHIA JET ===============================
            // ================================================================

            fastjet::ClusterSequence clust_seq_pythia(pythia_particles, jet_def_antikt); 
            std::vector<fastjet::PseudoJet> pythia_jets = sorted_by_pt(pythia_jet_select(clust_seq_pythia.inclusive_jets()));
            if(pythia_jets.size() == 0){ continue; }

            // get leading jet 
            fastjet::PseudoJet pythia_jet = pythia_jets.at(0);
            std::vector<fastjet::PseudoJet> pythia_jet_constituents = pythia_jet.constituents();

            jet_pt = pythia_jet.pt();
            jet_nconst = pythia_jet_constituents.size();
            jet_phi = pythia_jet.phi_std();
            jet_eta = pythia_jet.eta();

            h1_pythia_jet_pt->Fill(jet_pt);
            h1_pythia_jet_eta->Fill(jet_eta);
            h1_pythia_jet_phi->Fill(jet_phi);

            for (unsigned int iconst = 0; iconst < pythia_jet_constituents.size(); iconst++)
            {
                double pt = pythia_jet_constituents.at(iconst).pt();
                double eta = pythia_jet_constituents.at(iconst).eta();
                double phi = pythia_jet_constituents.at(iconst).phi_std();
                double pz = pythia_jet_constituents.at(iconst).pz();
                double e = pythia_jet_constituents.at(iconst).e();

                double dphi = phi - jet_phi;
                if(dphi > TMath::Pi()) dphi -= 2*TMath::Pi();
                if(dphi < -TMath::Pi()) dphi += 2*TMath::Pi();

                double deltaR = TMath::Sqrt(dphi*dphi + (eta - jet_eta)*(eta - jet_eta));

                const_pt.push_back(pt);
                const_eta.push_back(eta);
                const_delta_phi.push_back(dphi);
                const_pz.push_back(pz);
                const_e.push_back(e);

                h1_jet_const_delta_phi->Fill(dphi);
                h1_jet_const_deltaR->Fill(deltaR);
                h1_pythia_part_pt->Fill(pt);
                h1_pythia_part_eta->Fill(eta);
                h1_pythia_part_phi->Fill(phi);
            }
            
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
       
        // ============================================================
        // ==================== Fill Config Tree ======================
        // ============================================================

        n_accepted_events = iEvent;
        n_generated_events = info.nAccepted();
        sum_of_weights = info.weightSum();
        integrated_luminosity = (info.nAccepted()) /(info.sigmaGen()*1e9);
        xsec_over_eventweight = (info.sigmaGen() / info.weightSum());
        // fill event info tree
        config_tree->Fill();
        
        // ============================================================
        // ==================== WRITE OUTPUT ==========================
        // ============================================================

        // write output file
        outFile->cd();
        h1_jet_const_delta_phi->Write();
        h1_jet_const_deltaR->Write();
        h1_pythia_part_pt->Write();
        h1_pythia_part_eta->Write();
        h1_pythia_part_phi->Write();
        h1_pythia_jet_pt->Write();
        h1_pythia_jet_eta->Write();
        h1_pythia_jet_phi->Write();
        
        config_tree->Write();
        event_tree->Write();
        outFile->Close();
        return ;

    }

}