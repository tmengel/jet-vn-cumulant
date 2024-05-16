#ifndef REFFLOWCCALCULATOR_H
#define REFFLOWCCALCULATOR_H

#include <string>

class RefFlowCalculator
{
    public:

        RefFlowCalculator(const std::string &input_dir) 
            : m_input_dir(input_dir)
        {
        }
        ~RefFlowCalculator() {}

        // void OutputDir(const std::string &name) { m_output_dir = name; }
        // std::string OutputDir() { return m_output_dir; }

        void TruthFile(const std::string &name) { m_truth_file_name = name; }
        std::string TruthFile() { return m_truth_file_name; }

        void OutputFileName(const std::string &name) { m_output_file_name = name; }
        std::string OutputFileName() { return m_output_file_name; }

        int Run() { return Calculate(); }

    private:

        std::string m_input_dir{""};
        // std::string m_input_list_name {""};
        std::vector<std::string> m_input_files{};
        // std::string m_output_dir {""};

        std::string m_truth_file_name {""};
        std::string m_output_file_name {""};
       

        int Calculate();
        int CalculateHarm(int harmonic, std::vector<double> &two_part_ref_vec, std::vector<double> &four_part_ref_vec);
 

        std::pair<double,double> MeanErrVec(const std::vector<double> &vec);
        std::vector<double> C2Vec(const std::vector<double> &two_avg_vec);
        std::vector<double> vn2Vec(const std::vector<double> &c2_vec);
        std::vector<double> C4Vec(const std::vector<double> &two_avg_vec, const std::vector<double> &four_avg_vec);
        std::vector<double> vn4Vec( const std::vector<double> &c4_vec);

        double m_v2_ref_truth{0.0};
        double m_v3_ref_truth{0.0};
        double m_v4_ref_truth{0.0};
        void ReadTruthFile();


};


#endif // REFFLOWCCALCULATOR_H