#if !(defined(__CINT__) || defined(__CLING__)) || defined(__ACLIC__)

#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <cstdlib>
#include <cstdio>

#include <TFile.h>
#include <TString.h>
#include <TChain.h>
#include <TTree.h>
#include <TBranch.h>
#include <TObjArray.h>
#include <TObject.h>
#include <TKey.h>
#include <TProfile.h>
#include <TSystem.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TComplex.h>

#include <RooUnfoldResponse.h>
#include <RooUnfoldBayes.h>
#include <RooUnfoldErrors.h>

#endif

#include "DiffCorrUnfolder.h"
// #include "SpecUnfolder.h"

void UnfoldCorrs(std::string input_dir, std::string output_dir, int min_iter, int max_iter)
{
    // start EventAverageCorrs
    cout << "Starting DiffCorrUnfolder" << endl;
    DiffCorrUnfolder * cc = new DiffCorrUnfolder(input_dir);
    cc->OutputDir(output_dir);
    cc->MaxIter(max_iter);
    cc->MinIter(min_iter);
    cc->Run();
    // SpecUnfolder * cc = new SpecUnfolder(input_dir);
    // cc->OutputDir(output_dir);
    // std::cout << "MaxIter: " << max_iter << std::endl;
    // std::cout << "MinIter: " << min_iter << std::endl;
    // cc->Run();


    return;
}

int main(int argc, char** argv)
{
    if(argc != 5)
    {
        cout << "Usage: UnfoldCorrs <input_dir> <output_dir> <miniter> <maxiter>" << endl;
        return 1;
    }
    UnfoldCorrs(argv[1], argv[2], atoi(argv[3]), atoi(argv[4]));

    return 0;
}