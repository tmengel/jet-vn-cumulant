#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <cstdlib>
#include <algorithm>
#include <filesystem>

#include <SubSampler.h>
#include <Utils.h>


int main(int argc, char** argv)
{
    if(argc !=4)
    {
        std::cout << "Usage: " << argv[0] << " <input_dir_A> <input_dir_B> <output_dir>" << std::endl;
        return 1;
    }

    std::string input_dir_A = argv[1];
    std::string input_dir_B = argv[2];
    std::string output_dir = argv[3];

    // see if input directories exist
    if (!std::filesystem::exists(input_dir_A) || !std::filesystem::exists(input_dir_B))
    {
        std::cout << "Input directory A or B does not exist" << std::endl;
        return 1;
    }
    // see if output directory exists
    if (!std::filesystem::exists(output_dir))
    {
        std::cout << "Creating output directory " << output_dir << std::endl;
        std::filesystem::create_directory(output_dir);
    }

    std::vector<std::string> input_files_A_vec = Utils::GetFilesFromDir(input_dir_A);
    std::vector<std::string> input_files_B_vec = Utils::GetFilesFromDir(input_dir_B);

    int n_files_A = input_files_A_vec.size();
    int n_files_B = input_files_B_vec.size();

    if(n_files_A == 0 || n_files_B == 0)
    {
        std::cout << "No input files found in input directories" << std::endl;
        return 1;
    }

    bool need_to_clean_tmp = false;
    if(n_files_A != n_files_B)
    {
        need_to_clean_tmp = true;
        std::string larger_dir = "";
        int n_sub_samples = 0;
        
        if (n_files_A > n_files_B)
        {
            larger_dir = input_dir_A;
            n_sub_samples = n_files_B;
        }
        else
        {
            larger_dir = input_dir_B;
            n_sub_samples = n_files_A;
        }

        std::cout << "Number of files in input directories do not match" << std::endl;
        std::cout << "Larger directory: " << larger_dir << std::endl;
        std::cout << "Subsampling to match number of files" << std::endl;

        std::string tmp_output_dir = output_dir + "/tmp";
        std::filesystem::create_directory(tmp_output_dir);
        
        std::cout << "Starting SubSampler" << std::endl;
        SubSampler * ss = new SubSampler(larger_dir);
        ss->OutputDir(tmp_output_dir);
        ss->nSubSamples(n_sub_samples);
        ss->Run();

        std::cout << "Subsampling complete" << std::endl;
        if (n_files_A > n_files_B)
        {
            input_dir_A = tmp_output_dir;
            input_files_A_vec.clear();
            input_files_A_vec = Utils::GetFilesFromDir(input_dir_A);
        }
        else
        {
            input_dir_B = tmp_output_dir;
            input_files_B_vec.clear();
            input_files_B_vec = Utils::GetFilesFromDir(input_dir_B);
        }

        // check that they are the same size
        n_files_A = input_files_A_vec.size();
        n_files_B = input_files_B_vec.size();
        if (n_files_A != n_files_B)
        {
            std::cout << "Error: Subsampling did not produce same number of files" << std::endl;
            return 1;
        }

    }


    // merge files
    int n_files = n_files_A;
    if(n_files != n_files_B || n_files == 0)
    {
        std::cout << "Error: number of files in input directories do not match" << std::endl;
        return 1;
    }

    std::vector<std::string> input_files_to_merge {};
    for (int i = 0; i < n_files; i++)
    {
        input_files_to_merge.clear();
        std::string output_file = output_dir + "/merged_" + std::to_string(i) + ".root";
        
        input_files_to_merge.push_back(input_files_A_vec.at(i));
        input_files_to_merge.push_back(input_files_B_vec.at(i));
        
        SubSampler::MergeFileVec(input_files_to_merge, output_file);
    }

    if (need_to_clean_tmp)
    {
        std::cout << "Cleaning up temporary files" << std::endl;
        std::string cmd = "rm -rf " + output_dir + "/tmp";
        system(cmd.c_str());
    }

    std::cout << "Merging complete" << std::endl;
    // check that output files exist
    std::vector<std::string> output_files = Utils::GetFilesFromDir(output_dir);
    if (output_files.size() != n_files)
    {
        std::cout << "Error: number of output files does not match number of input files" << std::endl;
        return 1;
    }

    return 0;
}