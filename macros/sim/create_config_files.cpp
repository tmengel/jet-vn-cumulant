#include <Settings.h>

#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <utility>
#include <filesystem>

std::map<std::string, std::pair<float, float>> get_pTbin_map(std::vector<float> pt_binning)
{
    std::map<std::string, std::pair<float, float>> pTbin_map{};
    for (int i = 0; i < pt_binning.size() - 1; i++)
    {
        std::string pTbin_str = std::to_string(int(pt_binning.at(i))) + "to" + std::to_string(int(pt_binning.at(i + 1))) + "_GeV";
        std::cout << "pTbin_str: " << pTbin_str << std::endl;
        std::pair<float, float> pTbin_vals = std::make_pair(pt_binning.at(i), pt_binning.at(i + 1));
    
        pTbin_map[pTbin_str] = pTbin_vals;
    }
    return pTbin_map;
};

std::map<std::string,std::vector<std::string>> jet_vn_map =
{
    {"vn_0", {"0.04", "0.01", "0.005"}},
    {"2vn_0", {"0.08", "0.02", "0.01"}},
    {"5vn_0", {"0.2", "0.05", "0.025"}},
    {"linear", {"0.04 + 0.001 * x", "0.01 + 0.001 * x", "0.005 + 0.001 * x"}},
    {"logarithmic", {"0.04 + 0.02*TMath::Log(x)", "0.01 + 0.02*TMath::Log(x)", "0.005 + 0.02*TMath::Log(x)"}},
    {"cosine", {"0.08 + 0.04*TMath::Cos(0.2*x)", "0.02 + 0.01*TMath::Cos(0.2*x)", "0.01 + 0.005*TMath::Cos(0.2*x)"}}
};

void MakeConfig(
    std::string filename,
    int numEvents,
    float areaCut, 
    std::pair<float,float> ptHardBins,
    std::vector<std::string> jet_vn_strs = {})
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
    settings.Signal.Jet.MinPt = ptHardBins.first > 0 ? 1.0*ptHardBins.first : 5.0;

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
    std::string output_dir = "/lustre/isaac/scratch/tmengel/JetVnCumulantMethod/configs/";
    std::string output_filename = "toy_model_config_";
    std::string output_extension = ".txt";

    
    std::vector<float> pt_binning = {5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 35.0, 40.0, 50.0, 60.0, 80.0};
    std::vector<int> numJobs = {100,100,100,100,100,100,100,100,100,100,100};
    int numEvents = 1000000;
    float areaCut = 0.6;
    int iconfig = 0;
    for (int iptbin = 0; iptbin < pt_binning.size() - 1; iptbin++)
    {
        std::string pTbin_str = std::to_string(int(pt_binning.at(iptbin))) + "to" + std::to_string(int(pt_binning.at(iptbin + 1))) + "_GeV";
        std::pair<float,float> pTbin_vals = std::make_pair(pt_binning.at(iptbin), pt_binning.at(iptbin + 1));
        std::vector<std::string> jet_vn_strs = {};
        std::string output_filename_full = output_dir + output_filename + pTbin_str + output_extension;
        std::cout << "Creating config file: " << output_filename_full << std::endl;
        int Nevents = int(numEvents/numJobs.at(iconfig));
        MakeConfig(output_filename_full, Nevents, areaCut, pTbin_vals, jet_vn_strs);
        iconfig++;
    }

    numEvents = -1;
    std::pair<float,float> ptHardBins = std::make_pair(5.0, 80.0);
    output_dir = "/lustre/isaac/scratch/tmengel/JetVnCumulantMethod/configs/rotation/";
    if (!std::filesystem::exists(output_dir))
    {
        std::filesystem::create_directories(output_dir);
    }
    for (auto const &jet_vn : jet_vn_map)
    {
        std::string jet_vn_str = jet_vn.first;
        std::vector<std::string> jet_vn_strs = jet_vn.second;
        std::string output_filename_full = output_dir + output_filename + jet_vn_str + output_extension;
        std::cout << "Creating config file: " << output_filename_full << std::endl;
        int Nevents = int(numEvents/numJobs.at(iconfig));
        MakeConfig(output_filename_full, Nevents, areaCut, ptHardBins, jet_vn_strs);
        iconfig++;
    }

    return 0; 
}