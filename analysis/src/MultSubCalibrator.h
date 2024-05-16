#ifndef MULTSUBCALIBRATOR_H
#define MULTSUBCALIBRATOR_H

#include <string>
#include <vector>

class MultSubCalibrator
{
    public:

        MultSubCalibrator(const std::string &input_dir ) 
            : m_input_dir(input_dir)
        {
        }
        ~MultSubCalibrator() {}

        void OutputDir(const std::string &name) { m_output_dir = name; }
        std::string OutputDir() { return m_output_dir; }

        void SetPtBinning(const std::vector<double> &pt_bins) { m_pt_bins = pt_bins; m_user_pt_binning = true; }
        std::vector<double> PtBinning() { return m_pt_bins; }

        void SetCalibFileName(const std::string &name) { m_calib_file_name = name; }
        std::string CalibFile() { return m_calib_file_name; }

        int Run(bool overwrite=false) { return Calibrate(overwrite); }

    private:

        std::string m_input_dir{""};
        std::string m_input_list_name {""};
        std::string m_output_dir {""};
        std::string m_calib_file_name {""};

        bool m_user_pt_binning{false};
        std::vector<double> m_pt_bins;

        int Calibrate(bool overwrite);
};


#endif // MULTSUBCALIBRATOR_H