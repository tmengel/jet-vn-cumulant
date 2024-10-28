#include <iostream>
#include <string>
#include <vector>
#include <fstream>

#include <CorrCalculator.h>
#include <Unfolder.h>
#include <SimpleDiffFlow.h>
#include <DiffCorrFromHistos.h>

int main(int argc, char** argv)
{
    if(argc !=3)
    {
        std::cout << "Usage: " << argv[0] << " <input_dir> <output_dir>" << std::endl;
        return 1;
    }

    bool correlate = false;
    bool unfold = true;
    bool simple = true;
    bool diffcorrfromhists = true;

    std::string sub_inputfile = argv[1];
    std::string output_dir = argv[2];
    std::string corr_dir = output_dir+"/DiffCorrs";
   
    // std::vector<double> PtBins = {5.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 80.0};
    std::vector<double> PtBins = {5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 35.0, 40.0, 50.0, 60.0, 80.0};
    // std::vector<double> PtBins = {5.0, 7.5, 10.0, 12.5, 15.0, 17.5, 20.0, 25.0, 30.0, 35.0, 40.0, 45.0, 50.0};

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
        
        for (int i = 2; i < 5; i++)
        {
            std::string corr_file =  cc->GetOutputFileName(sub_inputfile, corr_dir+"/Order"+std::to_string(i), i);
            Unfolder * uf = new Unfolder(corr_file);
            uf->OutputDir(corr_dir + "/UnfoldedOrder"+std::to_string(i));
            uf->SetMaxIter(5);
            uf->SetMinIter(3);
            uf->Run();
        }
    }

    if(simple)
    {
        for (int i = 2; i < 5; i++)
        {
            std::string corr_file =  cc->GetOutputFileName(sub_inputfile, corr_dir+"/Order"+std::to_string(i), i);
            SimpleDiffFlow * sdf = new SimpleDiffFlow(corr_file);
            sdf->OutputDir(corr_dir + "/SimpleOrder"+std::to_string(i));
            sdf->Run();
        }
    }

    if(diffcorrfromhists)
    {
        for (int i = 2; i < 5; i++)
        {
            std::string corr_file =  cc->GetOutputFileName(sub_inputfile, corr_dir+"/Order"+std::to_string(i), i);
            DiffCorrFromHistos * dcfh = new DiffCorrFromHistos(corr_file);
            dcfh->OutputDir(corr_dir + "/DiffCorrHistosOrder"+std::to_string(i));
            dcfh->Run();
        }
    }
  


    
    return 0;
}
