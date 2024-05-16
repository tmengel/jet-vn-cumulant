#ifndef MULTSUBTRACTOR_H
#define MULTSUBTRACTOR_H

#include <string>
#include <vector>

class MultSubtractor
{
    public:

        MultSubtractor(const std::string &input_dir, const std::string & calib_file_name) 
            : m_input_dir(input_dir), m_calib_file_name(calib_file_name)
        {
        }
        ~MultSubtractor() {}

        void OutputDir(const std::string &name) { m_output_dir = name; }
        std::string OutputDir() { return m_output_dir; }

        void CalibFile(const std::string &name) { m_calib_file_name = name; }
        std::string CalibFile() { return m_calib_file_name; }

        int Run(bool overwrite = false) { return Subtract(overwrite); }

    private:

        std::string m_input_dir{""};
        // std::string m_input_list_name {""};
        std::vector<std::string> m_input_files{};
        std::string m_output_dir {""};
        std::string m_calib_file_name {""};

        int Subtract(bool overwrite);
        std::string GetSubFileName(std::string inputfile);
    

};

#endif // MULTSUBTRACTOR_H