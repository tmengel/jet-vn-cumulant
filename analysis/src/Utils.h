#ifndef UTILS_H
#define UTILS_H

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>

#include <cstdlib>
#include <cstdio>
#include <regex>
#include <filesystem>

#include <TFile.h>
#include <TString.h>
#include <TChain.h>
#include <TTree.h>
#include <TBranch.h>
#include <TObjArray.h>
#include <TObject.h>
#include <TSystem.h>

#include <ROOT/RDataFrame.hxx>
#include <ROOT/RDF/RInterface.hxx>

namespace Utils 
{

static std::string GetBaseName(const std::string &path)
{
    std::string base_name = path;
    size_t pos = base_name.find_last_of("/\\");
    if(pos != std::string::npos)
    {
        base_name = base_name.substr(pos+1);
    }

    pos = base_name.find_last_of(".");
    if(pos != std::string::npos)
    {
        base_name = base_name.substr(0, pos);
    }
    
    return base_name;
}

static std::vector<std::string> ReadListFromFile(const std::string & inputlist)
{
    std::vector<std::string> input_files;
    std::ifstream file(inputlist);
    std::string line;
    
    while(std::getline(file, line))
    {
        TString input_file = Form("%s", line.c_str());
        if(!input_file.EndsWith(".root")) continue;
        
        TFile * fin = new TFile(input_file, "READ");
        if(!fin->IsOpen())
        {
            std::cerr << "Cannot open file " << input_file << std::endl;
            continue;
        }
        fin->Close();
        std::string input_file_str = input_file.Data();
        input_files.push_back(input_file_str);
    }

    if (input_files.size() == 0)
    {
        std::cerr << "No input files found in " << inputlist << std::endl;
        exit(1);
    }

    std::cout << "Read " << input_files.size() << " files from " << inputlist << std::endl;

    return input_files;

}

static std::vector<std::string> GetListFromDir(const std::string &dir, const std::string &output_file, const std::string &ext=".root")
{
    std::vector<std::string> input_files;

    // check if inputlist is a file
    if(!std::ifstream(dir))
    {
        std::cerr << "Input directory " << dir << " does not exist" << std::endl;
        exit(1);
    }

    // get all files in the directory
    std::string cmd = "ls " + dir + "/*" + ext + " > " + output_file;
    std::system(cmd.c_str());


    input_files = ReadListFromFile(output_file);

    // add the directory to the file names
    // for(auto &f : input_files)
    // {
    //     f = dir + "/" + f;
    // }


    return input_files;
}

static std::vector<std::string> GetFilesFromDir(const std::string &dir)
{
    std::vector<std::string> input_files{};

    if(!std::filesystem::exists(dir))
    {
        std::cerr << "Input directory " << dir << " does not exist" << std::endl;
        exit(1);
    }

    bool has_trailing_slash = dir.back() == '/';
    std::string ls_dir = dir;
    if(!has_trailing_slash) ls_dir += "/";

    std::string cmd = "ls " + ls_dir+ "*.root";
    FILE *fp = popen(cmd.c_str(), "r");
    if(!fp)
    {
        std::cerr << "Cannot open directory " << dir << std::endl;
        exit(1);
    }

    // read the output of the command
    char path[1024];
    while(fgets(path, sizeof(path), fp) != NULL)
    {
        std::string file_path = path;
        file_path.erase(std::remove(file_path.begin(), file_path.end(), '\n'), file_path.end());
        input_files.push_back(file_path);
    }

    pclose(fp);

    if(input_files.size() == 0)
    {
        std::cerr << "No files found in directory " << dir << std::endl;
        exit(1);
    }

    return input_files;
}
  
static std::string WriteListToFile(const std::vector<std::string> &files, const std::string &output_file)
{
    std::ofstream file(output_file);
    if(!file.is_open())
    {
        std::cerr << "Cannot open output file " << output_file << std::endl;
        exit(1);
    }

    for(auto &f : files)
    {
        file << f<< std::endl;
    }
    file.close();
    return output_file;
}

static void MergeTrees(const std::string &filelist_name, const std::string &output_file, const std::string &tree_name)
{

    std::vector<std::string> input_files = ReadListFromFile(filelist_name);

    TChain *chain = new TChain(tree_name.c_str());
    for(auto &f : input_files)
    {
        chain->Add(f.c_str());
    }

    // merge the trees
    chain->Merge(output_file.c_str());

    delete chain;

    std::cout << "Done merging " << input_files.size() << " files into " << output_file << std::endl;

    return;


}

static double BranchMax(const std::string &file_name, const std::string &tree_name, const std::string &branch_name)
{

    double max = 0.0;

    // enable multithreading
    ROOT::EnableImplicitMT();

    ROOT::RDataFrame df(tree_name, file_name);

    max = 1.0*df.Max(branch_name).GetValue(); // convert to double in case of integer branches

    ROOT::DisableImplicitMT();

    std::cout << "Max value of branch " << branch_name << " in tree " << tree_name << " in file " << file_name << " is " << max << std::endl;
    
    return max;
}

static double BranchMin(const std::string &file_name, const std::string &tree_name, const std::string &branch_name)
{

    double min = 0.0;

    // enable multithreading
    ROOT::EnableImplicitMT();

    ROOT::RDataFrame df(tree_name, file_name);

    min = 1.0*df.Max(branch_name).GetValue(); // convert to double in case of integer branches

    ROOT::DisableImplicitMT();

    std::cout << "Min value of branch " << branch_name << " in tree " << tree_name << " in file " << file_name << " is " << min << std::endl;
    
    return min;
}

static void AddMinMaxTree(const std::string &file_name,const std::string &tree_name)
{
    // add min and max tree to the file for all branches in given tree
    // enable multithreading
    // ROOT::EnableImplicitMT();

    TFile *fin = new TFile(file_name.c_str(), "UPDATE");
    if(!fin->IsOpen())
    {
        std::cerr << "Cannot open file " << file_name << std::endl;
        exit(1);
    }

    TTree *tree = (TTree*)fin->Get(tree_name.c_str());
    if(!tree)
    {
        std::cerr << "Cannot find tree " << tree_name << " in file " << file_name << std::endl;
        exit(1);
    }

    TTree *min_max_tree = new TTree("min_max_tree", "Tree with min and max values of branches");

    // get the list of branches
    TObjArray *branches = tree->GetListOfBranches();
    for(int i = 0; i < branches->GetEntries(); i++)
    {
        TBranch *branch = (TBranch*)branches->At(i);
        std::string branch_name = branch->GetName();
        double min = 0.0;
        double max = 0.0;

        // get the min and max values of the branch
        min = 1.0*BranchMin(file_name, tree_name, branch_name);
        max = 1.0*BranchMax(file_name, tree_name, branch_name);

        min_max_tree->Branch((branch_name + "_min").c_str(), &min);
        min_max_tree->Branch((branch_name + "_max").c_str(), &max);
    }

    min_max_tree->Fill();

    fin->cd();
    min_max_tree->Write();
    fin->Close();

    // ROOT::DisableImplicitMT();

    std::cout << "Added min and max tree to " << file_name << std::endl;

    return;


}

}

#endif