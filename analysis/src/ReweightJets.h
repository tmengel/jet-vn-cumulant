#ifndef REWEIGHTJETS_H
#define REWEIGHTJETS_H

#include <string>
#include <vector>

class ReweightJets
{
    public:

        ReweightJets(const std::string &input_dir ) 
            : m_input_dir(input_dir)
        {
        }
        ~ReweightJets() {}

        void OutputDir(const std::string &name) { m_output_dir = name; }
        std::string OutputDir() { return m_output_dir; }

        int Run() { return Reweight(); }

    private:

        std::string m_input_dir{""};
        std::vector<std::string> m_input_files{};
        std::string m_output_dir {""};
        int Reweight();        
};


#endif // REWEIGHTJETS_H