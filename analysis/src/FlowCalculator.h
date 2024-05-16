#ifndef FLOWCALCULATOR
#define FLOWCALCULATOR


#include <string>
#include <vector>
#include <utility>
#include <iostream>

#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>

class FlowCalculator
{
    public:

        FlowCalculator() 
        {
        }
        ~FlowCalculator() {}

        void OutputFile(const std::string &name) { m_output_file_name = name; }
        std::string OutputFile() { return m_output_file_name; }

        void CalibFile(const std::string &name) { m_calib_file_name = name; }
        std::string CalibFile() { return m_calib_file_name; }

        void RefFile(const std::string &name) { m_ref_flow_file_name = name; }
        std::string RefFile() { return m_ref_flow_file_name; }

        void AddUnfoldedDir(const std::string &dir, const int harmonic){ m_unfold_input_files_harmonics.push_back(std::make_pair(dir, harmonic)); }
        void AddInputDir(const std::string &dir, const int harmonic){ m_input_files_harmonics.push_back(std::make_pair(dir, harmonic)); }

        int Run() { return Calculate(); }


    private:

        std::vector<std::pair<std::string, int>> m_input_files_harmonics{};
        std::vector<std::pair<std::string, int>> m_unfold_input_files_harmonics{};

        std::vector<std::string> m_input_files{};

        std::string m_calib_file_name {""};
        std::string m_ref_flow_file_name {""};
        std::string m_output_file_name {""};

        int Calculate();

        double m_v2_ref_truth = 0.0;
        double m_v3_ref_truth = 0.0;
        double m_v4_ref_truth = 0.0;
        double m_v2_two_ref = 0.0;
        double m_v3_two_ref = 0.0;
        double m_v4_two_ref = 0.0;
        double m_v2_four_ref = 0.0;
        double m_v3_four_ref = 0.0;
        double m_v4_four_ref = 0.0;
        double m_v2_two_ref_err = 0.0;
        double m_v3_two_ref_err = 0.0;
        double m_v4_two_ref_err = 0.0;
        double m_v2_four_ref_err = 0.0;
        double m_v3_four_ref_err = 0.0;
        double m_v4_four_ref_err = 0.0;

        TH1D * m_h1_jet_v2_func{nullptr};
        TH1D * m_h1_jet_v3_func{nullptr};
        TH1D * m_h1_jet_v4_func{nullptr};

        TH1D * m_h1_jet_v2_truth{nullptr};
        TH1D * m_h1_jet_v3_truth{nullptr};
        TH1D * m_h1_jet_v4_truth{nullptr};

        void MakeTruthHists(const int harmonic, const std::vector<std::string> input_files);
        void ReadCalibFile();
        void ReadRefFile();
        std::pair<TH1D*, TH1D*> GetVnHists(const int harmonic, const std::vector<std::string> input_files, const std::string & treename, 
                const std::string & v2branch, const std::string & v4branch, 
                const std::string &name, const int entry=0);
        static TH1D * GetHistFromVec(const std::vector<double> &vec, const std::vector<double> &err,
             const std::string &name, const std::vector<double> &bins);
        static std::vector<double> D2Vec(const std::vector<double> &two_avg_vec);
        static std::vector<double> vn2Vec(const std::vector<double> &d2_vec, const double &cn2);
        static std::vector<double> D4Vec(const std::vector<double> &four_avg_vec, const std::vector<double> &two_avg_vec, const double &cn2);
        static std::vector<double> vn4Vec(const std::vector<double> &d4_vec, const double &cn4);
        static std::pair<double,double> MeanErrVec(const std::vector<double> &vec);

}; 


#endif // FLOWCALCULATOR