#include <iostream>
#include <string>
#include <vector>
#include <fstream>

#include <SubSampler.h>
#include <Utils.h>

#include <MultSubCalibrator.h>
#include <MultSubtractor.h>

#include <CorrCalculator.h>
#include <RefFlowCalculator.h>

#include <Unfolder.h>

#include <FlowCalculator.h>


namespace Ana 
{
    bool Calibrate = false;
    bool Subtract = false;
    bool RefFlow = false;
    bool Correlate = false;
    bool Unfold = false;
    bool Flow = true;
    bool DiffCorrFromHistos =true;

    std::string InputDir = "";
    std::string OutputDir = "";
    std::string CalibFile = "";
    std::string SubtractedDir = "";
    std::string RefFlowFile = "";
    std::string CorrDir = "";
}

int main(int argc, char** argv)
{
    if(argc !=3)
    {
        std::cout << "Usage: " << argv[0] << " <input_dir> <output_dir>" << std::endl;
        return 1;
    }

    Ana::InputDir = argv[1];
    Ana::OutputDir = argv[2];
    Ana::CalibFile = Ana::OutputDir+"/CalibrationFile.root";
    Ana::SubtractedDir = Ana::OutputDir+"/MultSubtracted";
    Ana::CorrDir = Ana::OutputDir+"/DiffCorrs";
    Ana::RefFlowFile = Ana::OutputDir+"/ref_flow_results.root";
   
   
    std::vector<double> PtBins = {5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 35.0, 40.0, 50.0, 60.0, 80.0};

   
    // start FlowCalculator
    if(Ana::Flow)
    {
        std::cout << "Starting FlowCalculator" << std::endl;
        FlowCalculator * fc = new FlowCalculator();
        fc->OutputFile(Ana::OutputDir+"/flow_results.root");
        fc->CalibFile(Ana::CalibFile);
        fc->RefFile(Ana::RefFlowFile);
        for (int i = 2; i < 5; i++)
        {
            fc->AddUnfoldedDir(Ana::CorrDir + "/UnfoldedOrder"+std::to_string(i), i);
            fc->AddInputDir(Ana::CorrDir + "/Order"+std::to_string(i), i);
        }
        fc->Run();

        delete fc;
    }

    // start FlowCalculator
    if(Ana::DiffCorrFromHistos)
    {
        std::cout << "Starting FlowCalculator" << std::endl;
        FlowCalculator * fc = new FlowCalculator();
        fc->OutputFile(Ana::OutputDir+"/flow_results_nounfolding.root");
        fc->CalibFile(Ana::CalibFile);
        fc->RefFile(Ana::RefFlowFile);
        for (int i = 2; i < 5; i++)
        {
            fc->AddUnfoldedDir(Ana::CorrDir + "/DiffCorrHistosOrder"+std::to_string(i), i);
            fc->AddInputDir(Ana::CorrDir + "/Order"+std::to_string(i), i);
        }
        fc->Run();
    }


    
    return 0;
}
