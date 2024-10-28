#ifndef TOYMODELMERGER_H
#define TOYMODELMERGER_H

#include <string>
#include <vector>

namespace toymodel  
{
    class ToyModelMerger
    {
        public:
            ToyModelMerger() {};
            ~ToyModelMerger() {};

            void SetConfig(const std::string& filename) { m_config_filename = filename; }
            void SetInputSignal(const std::string& filename) { m_input_signal_filename = filename; }
            void SetInputBackground(const std::string& filename) { m_input_background_filename = filename; }
            void SetOutput(const std::string& filename) { m_output_filename = filename; }

            void Run();

        private:

            std::string m_config_filename{""};
            std::string m_output_filename{""};
            std::string m_input_signal_filename{""};
            std::string m_input_background_filename{""};

    };

}

#endif // TOYMODELMERGER_H
