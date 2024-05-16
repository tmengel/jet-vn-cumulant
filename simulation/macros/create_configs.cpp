#include <Settings.h>

#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <utility>



// std::map<std::string,std::pair<float,float>> pTbin_map =
// {
//     {"MB", std::make_pair(-1.0,-1.0)},
//     {"5to30_GeV", std::make_pair(5.0,30.0)},
//     {"30toInf_GeV", std::make_pair(30.0,-1.0)}
// };
std::map<std::string,std::pair<float,float>> pTbin_map =
{
    {"15toInf_GeV", std::make_pair(15.0,-1.0)}
};

std::map<std::string,std::vector<std::string>> jet_vn_map =
{
    // {"vn_0", {"0.04", "0.01", "0.005"}},
    {"2vn_0", {"0.08", "0.02", "0.01"}},
    // {"5vn_0", {"0.2", "0.05", "0.025"}},
    // {"linear", {"0.04 + 0.001 * x", "0.01 + 0.001 * x", "0.005 + 0.001 * x"}},
    // {"logarithmic", {"0.04 + 0.02*TMath::Log(x)", "0.01 + 0.02*TMath::Log(x)", "0.005 + 0.02*TMath::Log(x)"}},
    // {"cosine", {"0.08 + 0.04*TMath::Cos(0.2*x)", "0.02 + 0.01*TMath::Cos(0.2*x)", "0.01 + 0.005*TMath::Cos(0.2*x)"}}
};

void MakeConfig(
    std::string filename,
    int numEvents,
    float areaCut, 
    std::pair<float,float> ptHardBins,
    std::vector<std::string> jet_vn_strs)
{

    toymodel::ToyModelSettings settings;
    
    settings.Main.NumEvents = numEvents;

    // main settings (standard for all configurations)
    settings.Main.RandomSeed = 0; // 0 means random seed
    settings.Main.Verbosity = 1;
    settings.Main.Particle.MinPt = 0.5;
    settings.Main.Particle.MaxEta = 1.1;
    settings.Main.Jet.MinPt = 5.0;
    settings.Main.Jet.MinArea = areaCut;
    settings.Main.Jet.R = 0.4;

    // bkdg settings
    settings.Bkgd.Multiplicity = 363;
    settings.Bkgd.ConstEventPlane1 = -1; // -1 means random
    settings.Bkgd.ConstEventPlane2 = 0; // fixed to 0
    settings.Bkgd.ConstEventPlane3 = -1; // -1 means random
    settings.Bkgd.ConstEventPlane4 = 0; // fixed to 0

    settings.Signal.NumJets = 1;
    settings.Signal.Jet.MinPt = 5;

    // rotation settings
    if(jet_vn_strs.size() == 0)
    {
        settings.Signal.Jet.Rotate = false;
    }
    else
    {
        settings.Signal.Jet.Rotate = true;
        settings.Signal.Jet.v2Function = jet_vn_strs.at(0);
        settings.Signal.Jet.v3Function = jet_vn_strs.at(1);
        settings.Signal.Jet.v4Function = jet_vn_strs.at(2);
    }
   
    settings.Signal.Pythia.Config = "";

    settings.Signal.Pythia.Commands =  {
        "Beams:idA = 2212",
        "Beams:idB = 2212",
        "Beams:eCM = 200.0",
        "PDF:pSet = 13",
        "HardQCD:all = on",
        "Init:showProcesses = off",
        "Init:showMultipartonInteractions = off",
        "Init:showChangedSettings = on",
        "Next:numberCount = 1000000000",
        "Next:numberShowInfo = 0",
        "Next:numberShowProcess = 0",
        "Next:numberShowEvent = 0",
    };

    settings.Signal.Pythia.PtHardMin = ptHardBins.first;
    settings.Signal.Pythia.PtHardMax = ptHardBins.second;

    settings.write(filename);

    return;

}

int main()
{
    std::string output_dir = "/lustre/isaac/scratch/tmengel/jet-vn-cumulant/simulation/configs/";
    std::string output_filename = "toy_model_config_run2_";
    std::string output_extension = ".txt";

    int numEvents = 400000;
    float areaCut = 0.6;
    for (auto const &pTbin : pTbin_map)
    {
        for (auto const &jet_vn : jet_vn_map)
        {
            std::string pTbin_str = pTbin.first;
            std::pair<float,float> pTbin_vals = pTbin.second;
            std::string jet_vn_str = jet_vn.first;
            std::vector<std::string> jet_vn_strs = jet_vn.second;

            std::string output_filename_full = output_dir + output_filename + pTbin_str + "_" + jet_vn_str + output_extension;
            MakeConfig(output_filename_full, numEvents, areaCut, pTbin_vals, jet_vn_strs);
        }
    }

    return 0; 
}