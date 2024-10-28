#ifndef SETTINGS_H
#define SETTINGS_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>

#include <string>

#include <utility>
#include <vector>

#include <cstdlib>
#include <algorithm>

namespace toymodel  
{


class ToyModelSettings 
{

public:

    ToyModelSettings() = default;
    ~ToyModelSettings() = default;

    struct MainSettings
    {
        int NumEvents = 1000;
        unsigned int RandomSeed = 0;
        int Verbosity = 1;
        
        struct ParticleSettings 
        {
            float MinPt = 0.0;
            float MaxEta = 1.1;
            void print(std::ostream& os) const
            {
                os << std::setw(iN) << std::left << "Main::Particle::MinPt" << "== " << std::setw(iM) << std::left << MinPt << std::setw(iM) << std::left << "# Minimum pT for all particles unless specified in Bkgd or Signal" << std::endl;
                os << std::setw(iN) << std::left << "Main::Particle::MaxEta" << "== " << std::setw(iM) << std::left << MaxEta << std::setw(iM) << std::left << "# Maximum eta for all particles unless specified in Bkgd or Signal" << std::endl;
            }

        } Particle;

        struct JetSettings
        {
            float MinPt = 0.0;
            float MinArea = 0.0;
            float R = 0.4;
            float MaxEta = -1.0; // not set -> defaults to Particle::MaxEta - Jet::R

            void print(std::ostream& os) const
            {
                os << std::setw(iN) << std::left << "Main::Jet::MinPt" << "== " << std::setw(iM) << std::left << MinPt << std::setw(iM) << std::left << "# Minimum pT for all jets unless specified in Signal" << std::endl;
                os << std::setw(iN) << std::left << "Main::Jet::MaxEta" << "== " << std::setw(iM) << std::left << MaxEta << std::setw(iM) << std::left << "# Maximum eta for all jets unless specified in Signal (-1.0 means default to Particle::MaxEta - Jet::R)" << std::endl;
                os << std::setw(iN) << std::left << "Main::Jet::MinArea" << "== " << std::setw(iM) << std::left << MinArea << std::setw(iM) << std::left << "# Minimum area for merged jet (0.0 means no cut)" << std::endl;
                os << std::setw(iN) << std::left << "Main::Jet::R" << "== " << std::setw(iM) << std::left << R << std::endl;
            }


        } Jet;


        void print(std::ostream& os) const
        {
            std::string title = " Main ToyModelSettings ";
            int titleLength = title.length();
            int diff = (formatLength - titleLength) / 2;
            os << std::string(diff, '#') << title << std::string(diff, '#') << std::endl;

            os << std::setw(iN) << std::left << "Main::NumEvents" << "== " << std::setw(iM) << std::left << NumEvents << std::endl;
            os << std::setw(iN) << std::left << "Main::RandomSeed" << "== " << std::setw(iM) << std::left << RandomSeed << std::setw(iM) << std::left << "# 0 means use time seed" << std::endl;
            os << std::setw(iN) << std::left << "Main::Verbosity" << "== " << std::setw(iM) << std::left << Verbosity << std::setw(iM) << std::left << "# 0 = quiet, 1 = normal, 2 = verbose" << std::endl;
            os << std::endl;
            
            Particle.print(os);
            os << std::endl;

            Jet.print(os);
            os << std::endl;
        }

    } Main;

    struct BkgdSettings 
    {
        int Multiplicity = 363;
        float ConstEventPlane1 = -1.0;
        float ConstEventPlane2 = -1.0;
        float ConstEventPlane3 = -1.0;
        float ConstEventPlane4 = -1.0;

        void print(std::ostream& os) const
        {
            std::string title = " Background ToyModelSettings ";
            int titleLength = title.length();
            int diff = (formatLength - titleLength) / 2;
            os << std::string(diff, '#') << title << std::string(diff, '#') << std::endl; 


            os << std::setw(iN) << std::left << "Bkgd::Multiplicity" << "== " << std::setw(iM) << std::left << Multiplicity << std::endl;
            
            // print bool as string
            os << std::endl;
            
            os << std::endl;
            os << std::setw(iN) << std::left << "Bkgd::ConstEventPlane1" << "== " << std::setw(iM) << std::left << ConstEventPlane1 << std::setw(iM) << std::left << "# Constant event plane angle for v1" << std::endl;
            os << std::setw(iN) << std::left << "Bkgd::ConstEventPlane2" << "== " << std::setw(iM) << std::left << ConstEventPlane2 << std::setw(iM) << std::left << "# Constant event plane angle for v2" << std::endl;
            os << std::setw(iN) << std::left << "Bkgd::ConstEventPlane3" << "== " << std::setw(iM) << std::left << ConstEventPlane3 << std::setw(iM) << std::left << "# Constant event plane angle for v3" << std::endl;
            os << std::setw(iN) << std::left << "Bkgd::ConstEventPlane4" << "== " << std::setw(iM) << std::left << ConstEventPlane4 << std::setw(iM) << std::left << "# Constant event plane angle for v4" << std::endl;

            os << std::endl;

        }

    } Bkgd;

    struct SignalSettings 
    {
        int NumJets = 1;

        struct JetSettings
        {
            float MinPt = 0.0;
            bool Rotate = true;
            std::string v2Function = "0.04";
            std::string v3Function = "0.01";
            std::string v4Function = "0.005";
           
            void print(std::ostream& os) const
            {
                
                os << std::setw(iN) << std::left << "Signal::Jet::MinPt" << "== " << std::setw(iM) << std::left << MinPt << std::setw(iM) << std::left << "# Minimum pT for Pythia jet (optional defaults to Main::Jet::MinPt)" << std::endl;
                os << std::endl;
                if (Rotate) os << std::setw(iN) << std::left << "Signal::Jet::Rotate" << "== " << std::setw(iM) << std::left << "true" << std::setw(iM) << std::left << "# Rotate jet to event plane" << std::endl;
                else os << std::setw(iN) << std::left << "Signal::Jet::Rotate" << "== " << std::setw(iM) << std::left << "false" << std::setw(iM) << std::left << "# Rotate jet to event plane" << std::endl;
                os << std::setw(iN) << std::left << "Signal::Jet::v2Function" << "== " << std::setw(iM) << std::left << v2Function << std::setw(iM) << std::left << "# v2 coefficient function form (pt is x)" << std::endl;
                os << std::setw(iN) << std::left << "Signal::Jet::v3Function" << "== " << std::setw(iM) << std::left << v3Function << std::setw(iM) << std::left << "# v3 coefficient function form (pt is x)" << std::endl;
                os << std::setw(iN) << std::left << "Signal::Jet::v4Function" << "== " << std::setw(iM) << std::left << v4Function << std::setw(iM) << std::left << "# v4 coefficient function form (pt is x)" << std::endl;
            }

        } Jet;

        struct PythiaSettings
        {
            std::string Config = "";
            std::vector<std::string> Commands = {};
            float PtHardMin = -1.0;
            float PtHardMax = -1.0;
            void print(std::ostream& os) const
            {
                os << std::setw(iN) << std::left << "Signal::Pythia::PythiaConfig" << "== " << std::setw(iM) << std::left << Config << std::setw(iM) << std::left << "# Pythia configuration file (optional)" << std::endl;
                
                os << std::string(iM, '#') << " Pythia commands are strings passed to Pythia8::Pythia::readString() " <<  std::endl;
                for (auto command : Commands)
                {
                    os << std::setw(iN) << std::left << "Signal::Pythia::PythiaCommand" << "== " << std::setw(iM) << std::left << command << std::endl;

                }
                os << std::setw(iN) << std::left << "Signal::Pythia::PtHardMin" << "== " << std::setw(iM) << std::left << PtHardMin << std::setw(iM) << std::left << "# Minimum pT hard bin (optional)" << std::endl;
                os << std::setw(iN) << std::left << "Signal::Pythia::PtHardMax" << "== " << std::setw(iM) << std::left << PtHardMax << std::setw(iM) << std::left << "# Maximum pT hard bin (optional)" << std::endl;
            }
        
        } Pythia;

        void print(std::ostream& os) const
        {
            std::string title = " Signal ToyModelSettings ";
            int titleLength = title.length();
            int diff = (formatLength - titleLength) / 2;
            os << std::string(diff, '#') << title << std::string(diff, '#') << std::endl;
        
            os << std::setw(iN) << std::left << "Signal::NumJets" << "== " << std::setw(iM) << std::left << NumJets << std::setw(iM) << std::left << "# Number of jets to passed back to merging" << std::endl;
            os << std::endl;
            Jet.print(os);
            os << std::endl;
            Pythia.print(os);
            os << std::endl;
        }

     
    } Signal;

    void SetPtHardRange(float min, float max)
    {
        this->Signal.Pythia.PtHardMin = min;
        this->Signal.Pythia.PtHardMax = max;
    }

    void print(std::ostream& os) const;
    void write(std::string filename) const;
    void read(const std::string& filename);


private:

    static const int formatLength = 60;
    static const int iN = 40;
    static const int iM = 40;

    std::string the_title = " Toy Model ToyModelSettings";

    void pound_line(std::ostream& os) const
    {
        os << std::string(formatLength, '#') << std::endl;
        return ;
    }

    void formatTitle(std::ostream& os, const std::string &title) const
    {
        int titleLength = title.length();
        int diff = (formatLength - titleLength) / 2;
        os << std::string(diff, '#') << title << std::string(diff, '#') << std::endl;
        return ;

    }

};
    
} // namespace toymodel

#endif // SETTINGS_H
