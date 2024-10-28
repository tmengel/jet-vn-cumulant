#ifndef BKDGGEN_H
#define BKDGGEN_H

#include <string>
#include <vector>

#include "Settings.h"

namespace toymodel  
{
    class BkgdGen
    {
        public:

            BkgdGen(const std::string &config="", const std::string &outfile="")
                : m_config_filename(config), m_output_filename(outfile) {};
            ~BkgdGen() {};

            void SetConfigFile(const std::string& filename) { m_config_filename = filename; }
            std::string SetConfigFile() const { return m_config_filename; }
     
            void SetOutputFile(const std::string& filename) { m_output_filename = filename; }
            std::string GetOutputFile() const { return m_output_filename; }

            void Run();

        private:

            std::string m_config_filename{""};
            std::string m_output_filename{""};
    };

}

#endif
