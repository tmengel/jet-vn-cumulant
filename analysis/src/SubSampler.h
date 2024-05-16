#ifndef SUBSAMPLE_H
#define SUBSAMPLE_H

#include <string>
#include <vector>

class SubSampler
{
    public:

        SubSampler(const std::string &input_dir ) 
            : m_input_dir(input_dir)
        {
        }
        ~SubSampler() {}

        void OutputDir(const std::string &name) { m_output_dir = name; }
        std::string OutputDir() { return m_output_dir; }

        void nSubSamples(int n) { m_n_sub_samples = n; }
        int Run() { return SubSample(); }

        static void MergeFileVec(const std::vector<std::string> &vec, const std::string &output_file);
    private:

        std::string m_input_dir{""};
        
        int m_n_sub_samples{0};
        std::vector<std::string> m_input_files{};
        std::string m_output_dir {""};


        int SubSample();
        void MergeFiles();

        std::string m_outfile_base{""};
        void SetOutFileBase();
        std::string GetOutputFileName(int i);

        
};


#endif // SUBSAMPLE_H