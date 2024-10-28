#include "Settings.h"


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


void ToyModelSettings::print(std::ostream& os) const
{   

    pound_line(os);
    formatTitle(os, the_title);
    pound_line(os);
    os << std::endl;

    this->Main.print(os);
    this->Bkgd.print(os);
    this->Signal.print(os);
}


void ToyModelSettings::write(std::string filename) const
{
    std::ofstream file;
    file.open(filename);
    
    this->Main.print(file);
    this->Bkgd.print(file);
    this->Signal.print(file);

    file.close();
}

void ToyModelSettings::read(const std::string& filename)
{

    std::ifstream configFile(filename);
    if (!configFile.is_open())
    {
        std::cerr << "Error: Could not open file " << filename << std::endl;
        exit(1);
    }

    std::string line;
    while (std::getline(configFile, line)) 
    {
        if (line.empty() || line[0] == '/' || line[0] == '#') 
        {
            continue; // Skip empty lines and comments
        }

        std::string Section;
        std::string Subsection;
        std::string Key;

        // remove comments at the end of the line
        size_t pos = line.find("#");
        if (pos != std::string::npos) 
        {
            line = line.substr(0, pos);
        }

        // Get section name
        // Section::Subsection::Key = Value
        
        
        pos = line.find("::"); // Find first occurrence of "::"


        if (pos != std::string::npos) 
        {
            Section = line.substr(0, pos); // Get section name
            // remove whitespace in front and back
            Section.erase(std::remove_if(Section.begin(), Section.end(), ::isspace), Section.end());


            line = line.substr(pos + 2); // Skip "::"
        }

        // See if there is a subsection (might be empty)
        pos = line.find("::");
        if (pos != std::string::npos) 
        {
            Subsection = line.substr(0, pos); // Get subsection name
            // remove whitespace if 
            Subsection.erase(std::remove_if(Subsection.begin(), Subsection.end(), ::isspace), Subsection.end());

            line = line.substr(pos + 2); // Skip "::"
        }

        pos = line.find("==");
        if (pos != std::string::npos) 
        {
            Key = line.substr(0, pos); // Get key name
            // remove whitespace if any
            Key.erase(std::remove_if(Key.begin(), Key.end(), ::isspace), Key.end());

            line = line.substr(pos + 2); // Skip "=="
        }

        std::string Value = line;
        // only remove whitespace in front and end
        pos = Value.find("#");
        if (pos != std::string::npos) 
        {
            Value = Value.substr(0, pos);
        }
        Value.erase(0, Value.find_first_not_of(" \t"));
        Value.erase(Value.find_last_not_of(" \t") + 1);


        // Convert value to the correct type
        if (Section == "Main") 
        {
            if (Subsection == "Particle") 
            {
                if (Key == "MinPt") 
                {
                    Main.Particle.MinPt = std::stof(Value);
                } 
                else if (Key == "MaxEta") 
                {
                    Main.Particle.MaxEta = std::stof(Value);
                }
                else 
                {
                    std::cerr << "Error: Unknown key " << Key << " in section " << Section << "::" << Subsection << std::endl;
                    exit(1);
                }
            } 
            else if (Subsection == "Jet") 
            {
                if (Key == "MinPt") 
                {
                    Main.Jet.MinPt = std::stof(Value);
                } 
                else if (Key == "MinArea") 
                {
                    Main.Jet.MinArea = std::stof(Value);
                } 
                else if (Key == "MaxEta") 
                {
                    Main.Jet.MaxEta = std::stof(Value);
                }
                else if (Key == "R") 
                {
                    Main.Jet.R = std::stof(Value);
                } 
                else 
                {
                    std::cerr << "Error: Unknown key " << Key << " in section " << Section << "::" << Subsection << std::endl;
                    exit(1);
                }
            } 
            else if (Key == "NumEvents") 
            {
                Main.NumEvents = std::stoi(Value);
            } 
            else if (Key == "RandomSeed") 
            {
                Main.RandomSeed = std::stoi(Value);
            } 
            else if (Key == "Verbosity") 
            {
                Main.Verbosity = std::stoi(Value);
            }
            else 
            {
                std::cerr << "Error: Unknown key " << Key << " in section " << Section << std::endl;
                exit(1);
            }
        } 
        else if (Section == "Bkgd") 
        {
            if (Key == "Multiplicity") 
            {
                Bkgd.Multiplicity = std::stoi(Value);
            } 
            else if (Key == "ConstEventPlane1") 
            {
                Bkgd.ConstEventPlane1 = std::stof(Value);
            }
            else if (Key == "ConstEventPlane2") 
            {
                Bkgd.ConstEventPlane2 = std::stof(Value);
            }
            else if (Key == "ConstEventPlane3") 
            {
                Bkgd.ConstEventPlane3 = std::stof(Value);
            }
            else if (Key == "ConstEventPlane4") 
            {
                Bkgd.ConstEventPlane4 = std::stof(Value);
            }
            else 
            {
                std::cerr << "Error: Unknown key " << Key << " in section " << Section << std::endl;
                exit(1);
            }
    

        }
        else if (Section == "Signal") 
        {
            if (Subsection == "Jet") 
            {
                if (Key == "MinPt") 
                {
                    Signal.Jet.MinPt = std::stof(Value);
                } 
                else if (Key == "Rotate") 
                {
                    Signal.Jet.Rotate = (Value == "true");
                } 
                else if (Key == "v2Function") 
                {
                    Signal.Jet.v2Function = Value;
                } 
                else if (Key == "v3Function") 
                {
                    Signal.Jet.v3Function = Value;
                } 
                else if (Key == "v4Function") 
                {
                    Signal.Jet.v4Function = Value;
                }
                else 
                {
                    std::cerr << "Error: Unknown key " << Key << " in section " << Section << "::" << Subsection << std::endl;
                    exit(1);
                }
            } 
            else if (Subsection == "Pythia") 
            {
                if (Key == "PythiaConfig") 
                {
                    Signal.Pythia.Config = Value;
                } 
                else if (Key == "PythiaCommand") 
                {
                    Signal.Pythia.Commands.push_back(Value);
                } 
                else if (Key == "PtHardMin") 
                {
                    Signal.Pythia.PtHardMin = std::stof(Value);
                } 
                else if (Key == "PtHardMax") 
                {
                    Signal.Pythia.PtHardMax = std::stof(Value);
                } 
                else 
                {
                    std::cerr << "Error: Unknown key " << Key << " in section " << Section << "::" << Subsection << std::endl;
                    exit(1);
                }
            } 
            else if (Key == "NumJets") 
            {
                Signal.NumJets = std::stoi(Value);
            }
            else 
            {
                std::cerr << "Error: Unknown key " << Key << " in section " << Section << std::endl;
                exit(1);
            }
        }


    }
    configFile.close();
}

} // namespace toymodel