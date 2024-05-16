#include "Unfolder.h"

void Unfolder::AddPathPair(const std::string &input_dir, const std::string &output_dir)
{
    if (output_dir.empty())
    {
        m_path_pairs.push_back(std::make_pair(input_dir, input_dir));
    }
    else
    {
        m_path_pairs.push_back(std::make_pair(input_dir, output_dir));
    }

    return;
}

int Unfolder::Run()
{
    for (auto &path_pair : m_path_pairs)
    {
        std::string input_dir = path_pair.first;
        std::string output_dir = path_pair.second;

        std::string command = "root -l -b -q '" + root_exe + "(\"" + input_dir + "\", \"" + output_dir + "\", " + std::to_string(m_min_iter) + ", " + std::to_string(m_max_iter) + ")'";
        std::cout << "Running command: " << command << std::endl;
        system(command.c_str());
    }

    return 0;
}