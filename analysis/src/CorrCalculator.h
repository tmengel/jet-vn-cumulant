#ifndef CORRCALCULATOR_H
#define CORRCALCULATOR_H

#include <string>

#include <vector>
#include <iostream>



class CorrCalculator
{
    public:

        CorrCalculator(const std::string &input) 
            : m_input(input)
        {
        }
        ~CorrCalculator() {}

        void OutputDir(const std::string &name) { m_output_dir = name; }
        std::string OutputDir() { return m_output_dir; }

        void SetPtBinning(const std::vector<double> &pt_bins) { m_pt_bins = pt_bins; }
        std::vector<double> PtBinning() { return m_pt_bins; }

        // void SetTruthPtBinning(const std::vector<double> &pt_bins) { m_truth_pt_bins = pt_bins; }
        // std::vector<double> TruthPtBinning() { return m_truth_pt_bins; }

        int Run() { return Calculate(); }

        std::string GetOutputFileName(std::string inputfile,std::string output_dir, int harmonic);

    private:

        std::string m_input{""};
        std::string m_output_dir {""};
        std::vector<double> m_pt_bins{}; 
        // std::vector<double> m_truth_pt_bins{};
        std::vector<std::string> m_input_files{};


        int Calculate();
        int CalculateHarm(int harmonic);
        static std::vector<double> GetBinning(std::pair<double,double> minmax, const int n_bins);
        static std::vector<double> GetCorrBinning(std::pair<double,double> minmax, const int n_bins);
};


#endif // CORRCALCULATOR_H