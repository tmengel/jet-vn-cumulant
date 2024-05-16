#ifndef UNFOLDER_H
#define UNFOLDER_H

#include <string>
#include <vector>
#include <iostream>
#include <utility>

class Unfolder
{
    public:

        Unfolder(){}
        ~Unfolder() {}
    
        void AddPathPair(const std::string &input_dir, const std::string &output_dir = "");
        void SetMaxIter(int max_iter) { m_max_iter = max_iter; }
        void SetMinIter(int min_iter) { m_min_iter = min_iter; }


        int Run();

        typedef std::pair<std::string, std::string> PathPair;

    private:

        
        int m_max_iter {10};
        int m_min_iter {6};

        std::vector<std::pair<std::string, std::string>> m_path_pairs{};
        const std::string root_exe{"/lustre/isaac/scratch/tmengel/jet-vn-cumulant/ana/src/UnfoldCorrs.C"};
};


#endif // UNFOLDER_H