#ifndef UNFOLDER_H
#define UNFOLDER_H

#include <string>
#include <vector>
#include <iostream>
#include <utility>

class Unfolder
{
    public:

        Unfolder(const std::string &input_dir) 
            : m_input(input_dir)
        {
        }
        ~Unfolder() {}
    
        void OutputDir(const std::string &name) { m_output_dir = name; }
        std::string OutputDir() { return m_output_dir; }

        // void AddPathPair(const std::string &input_dir, const std::string &output_dir = "");
        void SetMaxIter(int max_iter) { m_max_iter = max_iter; }
        void SetMinIter(int min_iter) { m_min_iter = min_iter; }


        int Run() { return Calculate(); }


    private:

        std::string m_input{""};
        std::string m_output_dir {""};
        int m_max_iter {10};
        int m_min_iter {6};

        int Calculate();

        // std::vector<std::pair<std::string, std::string>> m_path_pairs{};
        const std::string root_exe{"/lustre/isaac/scratch/tmengel/JetVnCumulantMethod/jet-vn-cumulant/analysis/src/UnfoldCorrs.C"};
};


#endif // UNFOLDER_H