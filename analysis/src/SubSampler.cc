#include "SubSampler.h"
#include "Utils.h"

#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <cstdlib>
#include <algorithm>

#include <TSystem.h>
#include <TString.h>
#include <TFile.h>
#include <TTree.h>
#include <TKey.h>
#include <TProfile.h>
#include <TH1D.h>
#include <TH2D.h>

#include <ROOT/RDataFrame.hxx>
#include <ROOT/RVec.hxx>

int SubSampler::SubSample()
{
    // make output directory
    if (m_output_dir.empty())
    {
        std::cerr << "Output directory not set" << std::endl;
        exit(1);
    }

    if(gSystem->AccessPathName(m_output_dir.c_str()))
    {
        std::cout << "Error: output directory " << m_output_dir << " does not exist" << std::endl;
        std::cout << "Creating output directory " << m_output_dir << std::endl;
        gSystem->mkdir(m_output_dir.c_str(), kTRUE);
    }

    m_input_files.clear();
    m_input_files = Utils::GetFilesFromDir(m_input_dir);
    std::cout << "Number of input files: " << m_input_files.size() << std::endl;

    if (m_input_files.size() == 0)
    {
        std::cerr << "No input files found in " << m_input_dir << std::endl;
        exit(1);
    }

    // suffle input files if more than one
    if (m_input_files.size() > 1)
    {
        std::random_shuffle(m_input_files.begin(), m_input_files.end());
    }

    
    if(m_n_sub_samples > m_input_files.size())
    {
        std::cerr << "Number of sub samples greater than number of input files" << std::endl;
        exit(1);
    }
    else
    {
        MergeFiles();
    }

    return 0;
}

void SubSampler::MergeFiles()
{
    std::cout << "Merging files" << std::endl;


    int n_files = m_input_files.size();
    int n_files_per_sample = int(n_files / m_n_sub_samples);
    SetOutFileBase();


    std::vector<std::string> input_files_to_merge {};
    for (int ifile =0; ifile < m_n_sub_samples; ifile++)
    {

        std::string output_file = GetOutputFileName(ifile);

    
        // get first n_files_per_sample files in m_input_files vector
        input_files_to_merge.clear();
        for (int jfile = 0; jfile < n_files_per_sample; jfile++)
        {
            input_files_to_merge.push_back(m_input_files.at(jfile));

        }   

        //remove first n_files_per_sample files from m_input_files vector
        for (int jfile = 0; jfile < n_files_per_sample; jfile++)
        {
            m_input_files.erase(m_input_files.begin());
        }

        // if not divisible by n_files_per_sample, add the remaining files to this sample
        if (m_input_files.size() < n_files_per_sample)
        {
            for (int jfile = 0; jfile < m_input_files.size(); jfile++)
            {
                input_files_to_merge.push_back(m_input_files.at(jfile));
            }
            m_input_files.clear();
        }

        // print input files
        std::cout << "Input files to merge to " << output_file << std::endl;
        for (auto file : input_files_to_merge)
        {
            std::cout << "\t " << file << std::endl;
        }

        // merge files using hadd
        std::string hadd_cmd = "hadd -f " + output_file;
        for (auto file : input_files_to_merge)
        {
            hadd_cmd += " " + file;
        }

        std::cout << "Merging files to " << output_file << std::endl;
        std::cout << hadd_cmd << std::endl;
        std::system(hadd_cmd.c_str());

        input_files_to_merge.clear();

    }

    if (m_input_files.size() > 0)
    {
        std::cout << "ERROR: Remaining files: " << m_input_files.size() << std::endl;
        for (auto file : m_input_files)
        {
            std::cout << "\t " << file << std::endl;
        }
    }

    std::cout << "Done" << std::endl;
    return;
}


std::string SubSampler::GetOutputFileName(int i)
{
    std::string output_file = m_outfile_base + std::to_string(i) + ".root";
    return output_file;
}

void SubSampler::SetOutFileBase()
{
    m_outfile_base = m_output_dir;
    if (m_outfile_base.back() != '/') m_outfile_base += "/";

    // get first file name without path
    std::string first_file = m_input_files[0];
    size_t pos = first_file.find_last_of("/");
    if (pos != std::string::npos)
    {
        first_file = first_file.substr(pos+1);
    }

    // remove extension
    pos = first_file.find_last_of(".");
    if (pos != std::string::npos)
    {
        first_file = first_file.substr(0, pos);
    }

    // find "job_#"
    
    pos = first_file.find("job_");
    if (pos != std::string::npos)
    {
        first_file = first_file.substr(0, pos);
    }

    m_outfile_base += first_file + "_sub_";

    std::cout << "Output file base: " << m_outfile_base << std::endl;

    return;
}

void SubSampler::MergeFileVec(const std::vector<std::string> &vec, const std::string &output_file)
{
    if(vec.size() == 0)
    {
        std::cerr << "No files to merge" << std::endl;
        return;
    }

    if(output_file.empty())
    {
        std::cerr << "No output file specified" << std::endl;
        return;
    }

    // if(gSystem->AccessPathName(output_file.c_str()))
    // {
    //     std::cout << "Error: output file " << output_file << " does not exist" << std::endl;
    //     return;
    // }

    std::string hadd_cmd = "hadd -f " + output_file;
    for(auto &f : vec)
    {
        hadd_cmd += " " + f;
    }
    
    std::cout << "Merging files to " << output_file << std::endl;
    std::cout << hadd_cmd << std::endl;
    std::system(hadd_cmd.c_str());


    return;
}