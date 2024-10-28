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
    bool Calibrate = true;
    bool Subtract = true;
    bool RefFlow = true;
    bool Correlate = false;
    bool Unfold = false;
    bool Flow = false;

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

    // start MultSubCalibrator
    if(Ana::Calibrate)
    {
        std::cout << "Starting MultSubCalibrator" << std::endl;
        MultSubCalibrator * msc = new MultSubCalibrator(Ana::InputDir);
        msc->SetPtBinning(PtBins);
        msc->SetCalibFileName(Ana::CalibFile);
        msc->OutputDir(Ana::OutputDir);
        msc->Run(true); // overwrite = false
    }
    
    // start MultSubtractor
    if(Ana::Subtract)
    {
        std::cout << "Starting MultSubtractor" << std::endl;
        MultSubtractor * ms = new MultSubtractor(Ana::InputDir, Ana::CalibFile);
        ms->OutputDir(Ana::SubtractedDir);
        ms->Run(true); // overwrite = false
    }

    // start RefFlowCalculator
    if(Ana::RefFlow)
    {
        std::cout << "Starting RefFlowCalculator" << std::endl;
        RefFlowCalculator * rfc = new RefFlowCalculator(Ana::InputDir);
        rfc->TruthFile(Ana::CalibFile);
        rfc->OutputFileName(Ana::RefFlowFile);
        rfc->Run();
    }


    
    return 0;
}
