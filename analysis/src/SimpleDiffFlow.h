#ifndef SimpleDiffFlow_H
#define SimpleDiffFlow_H

#include "CorrFunctions.h"

#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <cstdlib>

#include <TSystem.h>
#include <TString.h>
#include <TFile.h>
#include <TTree.h>
#include <TKey.h>
#include <TProfile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TComplex.h>


class SimpleDiffFlow
{
    public:

        SimpleDiffFlow(const std::string &input_dir) 
            : m_input(input_dir)
        {
        }
        ~SimpleDiffFlow() {}

        void OutputDir(const std::string &name) { m_output_dir = name; }
        std::string OutputDir() { return m_output_dir; }

        int Run() { return Calculate(); }
        
        void SetInputFiles(const std::vector<std::string> &input_files) { m_input_files = input_files; }

    private:

        std::string m_input{""};
        std::vector<std::string> m_input_files{};
        std::string m_output_dir {""};
        

        int Calculate();

        
        std::string GetOutputFileName(std::string inputfile);
        static TH1D * Make1D(TH2D * h2, std::string name, double w=1.0);
        static std::vector<double> ConvertToVec(TH2D * h2, double w=1.0);

        static std::vector<double> D2Vec(const std::vector<double> &two_avg_vec);
        static std::vector<double> vn2Vec(const std::vector<double> &d2_vec, const double &cn2);
        static std::vector<double> D4Vec(const std::vector<double> &four_avg_vec, const std::vector<double> &two_avg_vec, const double &cn2);
        static std::vector<double> vn4Vec(const std::vector<double> &d4_vec, const double &cn4);
        static std::pair<double,double> MeanErrVec(const std::vector<double> &vec);


};




#endif // SimpleDiffFlow_H