#include <string>
#include <iostream>
#include <cstdlib>
#include <vector>
#include <algorithm>
#include <filesystem>
#include <chrono>

#include <SubSampler.h>
#include <ToyModelMerger.h>
#include <Utils.h>

static std::string register_directory(const std::string & dir)
{
    // if(clear && std::filesystem::exists(dir))
    // {
    //     // std::filesystem::remove_all(dir, std::filesystem
    // }

    if(!std::filesystem::exists(dir))
    {
        std::filesystem::create_directory(dir);
        if (!std::filesystem::exists(dir))
        {
            std::cerr << "Error: could not create directory " << dir << std::endl;
            exit(1);
        }
    }
    return dir;
}

int main(int argc, char* argv[])
{
    if(argc != 3)
    {
        std::cerr << "Usage: " << argv[0] << " <config_file> <file idx>" << std::endl;
        return 1;
    }
    
    // start timer
    auto start = std::chrono::high_resolution_clock::now();

    std::string config_dir = argv[1];
    int file_idx = std::stoi(argv[2]);
    std::cout << "Config file: " << config_dir << std::endl;
    std::string rootfile_path="/lustre/isaac/scratch/tmengel/JetVnCumulantMethod/rootfiles/sim";

    std::vector<std::string> config_files{};
    for (const auto & entry : std::filesystem::directory_iterator(config_dir))
    {
        std::string config_file = entry.path().string();
        config_files.push_back(config_file);
    }

    for (const auto & config_file : config_files)
    {
        std::string config_base = Utils::GetBaseName(config_file);
        std::string::size_type pos = config_base.find("toy_model_config_");
        if(pos != std::string::npos)
        {
            config_base.erase(pos, 17);
        }

        std::cout << "=============================================" << std::endl;
        std::cout << "Config file: " << config_file << std::endl;
        std::cout << "Config base: " << config_base << std::endl;
        std::cout << "=============================================" << std::endl;

        for (const auto & entry : std::filesystem::directory_iterator(rootfile_path))
        {
            std::string bkgd_dir = entry.path().string() + "/bkgd";
            std::string sig_dir = entry.path().string() + "/sig";
            std::string rotated_dir = register_directory(entry.path().string() + "/rotated");
            std::vector<std::string> bkgd_files = Utils::GetFilesFromDir(bkgd_dir);
            std::vector<std::string> sig_files = Utils::GetFilesFromDir(sig_dir);
            
            if(file_idx >= bkgd_files.size() || file_idx >= sig_files.size())
            {
                std::cerr << "Error: file index " << file_idx << " out of range" << std::endl;
                exit(1);
            }
            std::string output_dir = register_directory(rotated_dir + "/" + config_base);

            std::string bkgd_file = bkgd_files.at(file_idx);
            std::string sig_file = sig_files.at(file_idx);
            std::string bkgd_base = Utils::GetBaseName(bkgd_file);
            std::string sig_base = Utils::GetBaseName(sig_file);
            std::string output_file = output_dir + "/" + bkgd_base + "_" + config_base + ".root";
            std::cout << "Background file: " << bkgd_file << std::endl;
            std::cout << "Signal file: " << sig_file << std::endl;
            std::cout << "Output file: " << output_file << std::endl;

            toymodel::ToyModelMerger * merger = new toymodel::ToyModelMerger();
            merger->SetConfig(config_file);
            merger->SetInputSignal(sig_file);
            merger->SetInputBackground(bkgd_file);
            merger->SetOutput(output_file);
            merger->Run();
        }
    }

    // end timer
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    std::cout << "Elapsed time: " << elapsed.count() << " s\n";

    
    return 0;
}