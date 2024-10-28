#include <string>
#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <vector>
#include <filesystem>
#include <chrono>

#include <SubSampler.h>
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
        std::cerr << "Usage: " << argv[0] << " <input dir> <output dir>" << std::endl;
        return 1;
    }
    
    // start timer
    auto start = std::chrono::high_resolution_clock::now();
    std::string input_dir = argv[1];
    std::string output_dir = argv[2];
    output_dir = register_directory(output_dir);
    int n_files = 10;
    SubSampler * ss = new SubSampler(input_dir);
    ss->OutputDir(output_dir);
    ss->nSubSamples(n_files);
    ss->Run();
   
    // end timer
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    std::cout << "Elapsed time: " << elapsed.count() << " s\n";

    
    return 0;
}
