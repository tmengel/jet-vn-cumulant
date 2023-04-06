// #include "TennGen.h"
#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TSystem.h"

#include <iostream>
#include "Pythia8/Pythia.h"

using namespace Pythia8;
// using namespace tenngen;
using namespace std;

const Int_t MAXPARTS = 3000;
const Int_t nPtBins = 26;
const Double_t pthardbin[nPtBins]= {5.,7.,8.,10.,11.,12.,15.,16.,17.,20.,21.,23.,25.,27.,30.,35.,40.,45.,50.,55.0,60.0,65.0,70.0,75.0,80.0,-1};

Int_t GeneratePythiaPP(Int_t collen, Int_t nevent, Int_t ptbin){

    Int_t nparts;
    Double_t particle_px[MAXPARTS], particle_py[MAXPARTS], particle_pz[MAXPARTS], particle_e[MAXPARTS];
    Int_t particle_id[MAXPARTS];
    Double_t weight, ptmin, ptmax, xsec_over_eventweight;

    TString datadir = "../root-files/PP/";
    if (gSystem->AccessPathName(datadir.Data())) gSystem->mkdir(datadir.Data());
    
    TString filename = Form("%s%dGeV_PP_ptbin%d",datadir.Data(),collen,ptbin);
    filename= filename + ".root";
    TFile *outFile = new TFile(filename.Data(),"RECREATE");
    TTree* outTree = new TTree("tree", "tree");
    TTree* eventInfo = new TTree("eventInfo","eventInfo");

    eventInfo->Branch("nevent",&nevent,"nevent/I");
    eventInfo->Branch("ptmin",&ptmin,"ptmin/D");
    eventInfo->Branch("ptmax",&ptmax,"ptmax/D");
    eventInfo->Branch("ptbinID", &ptbin, "ptbinID/I");
    eventInfo->Branch("xsec_over_eventweight",&xsec_over_eventweight,"xsec_over_eventweight/D");

    outTree->Branch("weight", &weight, "weight/D");
    outTree->Branch("nparts",&nparts,"nparts/I");
    outTree->Branch("particle_px", particle_px, "particle_px[nparts]/D");
    outTree->Branch("particle_py", particle_py, "particle_py[nparts]/D");
    outTree->Branch("particle_pz", particle_pz, "particle_pz[nparts]/D");
    outTree->Branch("particle_e", particle_e, "particle_e[nparts]/D");
    outTree->Branch("particle_id", particle_id, "particle_id[nparts]/I");

    Pythia pythia;
    Pythia8::Settings& settings = pythia.settings;
    const Pythia8::Info& info = pythia.info;
    Event& event = pythia.event;
    if(collen == 200 )pythia.readFile("pythia-configs/pythiaSettings_RHIC.cmnd");
    else if(collen == 2760) pythia.readFile("pythia-configs/pythiaSettings_LHC.cmnd");
    else{
        cout << "Invalid collision energy" << endl;
        return 1;
    }
    settings.parm("Main:numberOfEvents", nevent);
    settings.parm("PhaseSpace:pTHatMin", pthardbin[ptbin]);
    settings.parm("PhaseSpace:pTHatMax", pthardbin[ptbin+1]);
    pythia.init();

    for (Int_t ievent=0; ievent<nevent; ievent++){
        if (!pythia.next()) continue;
        weight = info.weight();
        nparts = 0;
        for (Int_t i = 0; i < event.size(); i++){
            //if (event[i].isFinal() && event[i].isCharged()){
            if (event[i].isFinal()){
               particle_px[nparts] = event[i].px();
                particle_py[nparts] = event[i].py();
                particle_pz[nparts] = event[i].pz();
                particle_e[nparts] = event[i].e();
                particle_id[nparts] = event[i].id();
                nparts++;
            }
        }

        outTree->Fill();
        for (Int_t k =0; k<nparts; k++){
            particle_px[k] = 0.0;
            particle_py[k] = 0.0;
            particle_pz[k] = 0.0;
            particle_e[k] = 0.0;
        }



    }

    
    ptmin = pthardbin[ptbin];
    ptmax = pthardbin[ptbin+1];
    xsec_over_eventweight = (info.sigmaGen() / info.weightSum());
    eventInfo->Fill();
    outFile->Write();
    outFile->Close();
    delete outFile;   
    cout << "Done" << endl;
    return 0;
}

Int_t main(Int_t argc, char** argv){
    
    Int_t nevent = 100;
    Int_t ptbin = 0;
    if(argc != 4){
        cout << "Usage: ./GeneratePythiaPP <collen> <nevents> <ptbin>" << endl;
        return 1;
    }
    Int_t collen = atoi(argv[1]);
    nevent = atoi(argv[2]);
    ptbin = atoi(argv[3]);
    cout << "Generating " << nevent << " events in pt bin " << ptbin << endl;
    cout << "Collision energy: " << collen << endl;
    return GeneratePythiaPP(collen,nevent, ptbin);
}
