#include <string>
#include <iostream>
#include <chrono>

#include <PythiaGen.h>
#include <BkgdGen.h>

int main(int argc, char* argv[])
{
    if(argc != 4)
    {
        std::cerr << "Usage: " << argv[0] << " <config_file> <signal_output_file> <background_output_file>" << std::endl;
        return 1;
    }
    
    // start timer
    auto start = std::chrono::high_resolution_clock::now();

    std::string config_file = argv[1];
    std::string signal_output_file = argv[2];
    std::string background_output_file = argv[3];

    toymodel::PythiaGen * model = new toymodel::PythiaGen();
    model->SetConfigFile(config_file);
    model->SetOutputFile(signal_output_file);
    model->Run();

    toymodel::BkgdGen * model_bkgd = new toymodel::BkgdGen();
    model_bkgd->SetConfigFile(config_file);
    model_bkgd->SetOutputFile(background_output_file);
    model_bkgd->Run();

    // end timer
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    std::cout << "Elapsed time: " << elapsed.count() << " s\n";

    
    return 0;
}