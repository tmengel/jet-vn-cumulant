#include <iostream>
#include <string>
#include <vector>
#include <fstream>

#include <CorrCalculator.h>
#include <Unfolder.h>

int main(int argc, char** argv)
{
    if(argc !=3)
    {
        std::cout << "Usage: " << argv[0] << " <input_dir> <output_dir>" << std::endl;
        return 1;
    }

    bool correlate = true;
    bool unfold = true;

    std::string sub_inputfile = argv[1];
    std::string output_dir = argv[2];
    std::string corr_dir = output_dir+"/DiffCorrs";
   
    // std::vector<double> PtBins =  {5.0, 10.0, 15.0, 20.0, 30.0, 40.0, 60.0};
    std::vector<double> PtBins = {5.0, 10.0, 15.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 100.0};
    // std::vector<double> RecoPtBins = {5.0, 10.0, 15.0, 20.0,30.0, 40.0, 60.0};
    // std::vector<double> TruthPtBins = {5.0, 10.0, 15.0, 20.0, 30.0, 40.0, 60.0};

    // start CorrCalculator
    std::cout << "Starting CorrCalculator" << std::endl;
    CorrCalculator * cc = new CorrCalculator(sub_inputfile);
    cc->OutputDir(corr_dir);
    cc->SetPtBinning(PtBins);
    // cc->SetTruthPtBinning(TruthPtBins);
    if(correlate)
    {
        cc->Run();
    }
    std::cout << "Finished CorrCalculator" << std::endl;
    if(unfold)
    {
        Unfolder * uf = new Unfolder();
        for (int i = 2; i < 5; i++)
        {
            std::string corr_file =  cc->GetOutputFileName(sub_inputfile, corr_dir+"/Order"+std::to_string(i), i);
            uf->AddPathPair(corr_file, corr_dir + "/UnfoldedOrder"+std::to_string(i));
            std::cout << "Adding path pair: " << corr_file << " " << corr_dir + "/UnfoldedOrder"+std::to_string(i) << std::endl;
        }

        uf->SetMaxIter(9);
        uf->SetMinIter(5);
        uf->Run();
    }

  


    
    return 0;
}
