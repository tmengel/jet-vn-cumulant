#ifndef TOYMODEL_H
#define TOYMODEL_H

#include <string>
#include <vector>
#include <fastjet/PseudoJet.hh>

namespace toymodel  
{
    class ToyModel
    {
        public:
            ToyModel() {};
            ~ToyModel() {};

            void set_config(const std::string& filename) { m_config_filename = filename; }
            void set_output(const std::string& filename) { m_output_filename = filename; }

            void run();

        private:

            std::string m_config_filename = "config.txt";
            std::string m_output_filename = "output.txt";

            static bool accept_pythia_event(std::vector<fastjet::PseudoJet> &jets, double pTmin, double pTmax);
    };

}

#endif
