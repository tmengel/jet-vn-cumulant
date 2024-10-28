#ifndef DIFFCORRFROMHISTOS_H
#define DIFFCORRFROMHISTOS_H

#include "CorrFunctions.h"
#include "Utils.h"

#include <iostream>
#include <string>
#include <vector>

#include <TH1D.h>
#include <TH2D.h>


class DiffCorrFromHistos
{
    public:

        DiffCorrFromHistos(const std::string &input_dir) 
            : m_input(input_dir)
        {
        }
        ~DiffCorrFromHistos() {}

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

       


};


#endif // DIFFCORRFROMHISTOS_H