#include <string>
#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <vector>
#include <filesystem>
#include <chrono>

#include <SubSampler.h>
#include <Utils.h>

int main(int argc, char* argv[])
{
    if(argc != 2)
    {
        std::cerr << "Usage: " << argv[0] << "<input dir>" << std::endl;
        return 1;
    }
    
    // start timer
    auto start = std::chrono::high_resolution_clock::now();

    // get rotation directories
    std::string input_dir = argv[1];
    if(input_dir.back() == '/')
    {
        input_dir.pop_back();
    }
    std::vector<std::vector<std::string>> file_merging_vec{};
    int n_files_per_pt_bin = 10;
    std::vector<std::string> sub_dirs{};
    for (const auto & entry : std::filesystem::directory_iterator(input_dir))
    {
       std::vector<std::string> files = Utils::GetFilesFromDir(entry.path().string());
       if(files.size() != n_files_per_pt_bin)
       {
           std::cerr << "Error: " << files.size() << " files in " << entry.path().string() << std::endl;
           exit(1);
       }

        file_merging_vec.push_back(files);
        sub_dirs.push_back(entry.path().string());
    }

    std::string output_dir = input_dir;
    std::string rotation = input_dir.substr(input_dir.find_last_of("/")+1);
    std::cout << "Rotation: " << rotation << std::endl;
    for (int i = 0; i < n_files_per_pt_bin; i++)
    {
        std::vector<std::string> files_to_merge{};
        for (const auto & files : file_merging_vec)
        {
            files_to_merge.push_back(files.at(i));
        }

        std::string output_file = output_dir + "/" + "toy_model_" + rotation + "_allpTbins_subsample_" + std::to_string(i) + ".root";
        SubSampler::MergeFileVec(files_to_merge, output_file);
        files_to_merge.clear();
    }

    std::vector<std::string> merged_files = Utils::GetFilesFromDir(output_dir);
    if(merged_files.size() != n_files_per_pt_bin)
    {
        std::cerr << "Error: " << merged_files.size() << " files in " << output_dir << std::endl;
        exit(1);
    }

    // rm subdirectories
    for (const auto & sub_dir : sub_dirs)
    {
        // std::string cmd = "rm -rf " + sub_dir;
        std::cout << "Removing " << sub_dir << std::endl;
        std::string cmd = "rm -rf " + sub_dir;
        system(cmd.c_str());
    }

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    std::cout << "Elapsed time: " << elapsed.count() << " s\n";

    
    return 0;
}