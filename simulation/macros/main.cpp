#include <ToyModel.h>
#include <string>
#include <iostream>
#include <chrono>


int main(int argc, char* argv[])
{
    if(argc != 3)
    {
        std::cerr << "Usage: " << argv[0] << " <config file> <output file>" << std::endl;
        return 1;
    }
    
    // start timer
    auto start = std::chrono::high_resolution_clock::now();

    std::string config_file = argv[1];
    std::string output_file = argv[2];

    toymodel::ToyModel model;
    model.set_config(config_file);
    model.set_output(output_file);
    model.run();

    // end timer
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    std::cout << "Elapsed time: " << elapsed.count() << " s\n";

    
    return 0;
}