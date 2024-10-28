#include "Unfolder.h"

#include <iostream>
#include <string>
#include <cstdlib>
#include <vector>
#include <cstdio>

int Unfolder::Calculate()
{
  
    std::string command = "root -l -b -q '" + root_exe + "(\"" + m_input + "\", \"" + m_output_dir + "\", " + std::to_string(m_min_iter) + ", " + std::to_string(m_max_iter) + ")'";
    std::cout << "Running command: " << command << std::endl;
    system(command.c_str());

    return 0;
}