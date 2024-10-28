import sys, os
import argparse

import ROOT
import numpy as np

def GoodRootFile(filename):
    
        # open the file
        f = ROOT.TFile(filename, "READ")
        if not f:
                print("Error opening file {0}".format(filename))
                return False
        t = f.Get("tree")
        if not t:
                print("Error getting tree from file {0}".format(filename))
                return False

        return True 

def ParseConfig(filename, config_vals):
        # open the file
        f = open(filename, "r")
        if not f:
                print("Error opening file {0}".format(filename))
                return {}
        
        lines = f.readlines()
        vals = {}
        for val in config_vals:
            for line in lines:
                var, value = line.split(":")
                if val == var.strip():
                    vals[val] = value.strip()
        return vals

def ParseQAFile(filename):
        QA_dict = {}
        with open(filename, "r") as f:
               # read the first 4 lines
                lines = f.readlines()
                config_file = lines[0].split(":")[1].strip()
                n_jobs = int(lines[1].split(":")[1].strip())
                n_events = int(lines[2].split(":")[1].strip())
                output_path = lines[3].split(":")[1].strip()

                QA_dict["config_file"] = config_file
                QA_dict["n_jobs"] = n_jobs
                QA_dict["n_events"] = n_events
                QA_dict["output_path"] = output_path

                for ijob in range(n_jobs):
                        # find the line with Job<ijob>:
                        found = False
                        for line in lines:
                                if "Job{0}".format(ijob) in line:
                                        found = True
                                        jobid = line.split(":")[0].strip()
                                        split_line = line.split(":")[1].split(",")
                                        jobname = split_line[0].strip()
                                        jobfile = split_line[1].strip()
                                        outfile_sig = split_line[2].strip()
                                        outfile_bkgd = split_line[3].strip()
                                        errfile = split_line[4].strip()
                                        logfile = split_line[5].strip()
                                        jobnumber = int(split_line[6].strip())
                                        QA_dict[jobid] = {"jobname": jobname, "jobfile": jobfile, "outfile_sig": outfile_sig, "outfile_bkgd": outfile_bkgd,
                                                           "errfile": errfile, "logfile": logfile, "jobnumber": jobnumber}
                        if not found:
                                print("Error: Job{0} not found".format(ijob))                
        return QA_dict

def ScanErrorFiles(errfiles):
        none_empty_errfiles = {}
        for key in errfiles:
                errfile = errfiles[key]
                
                delete_file = True
                if not os.path.exists(errfile):
                        continue
                with open(errfile, "r") as f:
                        lines = f.readlines()
                        if len(lines) > 0:
                                none_empty_errfiles[key] = errfile
                                delete_file = False
                if delete_file:
                        os.remove(errfile)

        return none_empty_errfiles

def ScanOutputFiles(output_files):
        bad_root_files = {}
        for key in output_files:
                outfile = output_files[key]
                if(os.path.exists(outfile)):
                    if not GoodRootFile(outfile):
                        bad_root_files[key] = outfile
                else:
                    bad_root_files[key] = outfile
        return bad_root_files

def ScanLogFiles(logfiles):
        # logfile_details = {}
        logfile_details = {}
        all_logfiles =[]
        keys_to_find = ["Main::NumEvents", "Main::RandomSeed", "Main::Particle::MinPt", "Main::Particle::MaxEta",
                                "Main::Jet::MinPt", "Main::Jet::MaxEta", "Main::Jet::MinArea", "Main::Jet::R",
                                "Bkgd::Multiplicity", "Bkgd::ConstEventPlane1", "Bkgd::ConstEventPlane2", "Bkgd::ConstEventPlane3", "Bkgd::ConstEventPlane4",
                                "Signal::Jet::MinPt", "Signal::Jet::Rotate", "Signal::Jet::v2Function", "Signal::Jet::v3Function", "Signal::Jet::v4Function",
                                "Signal::Pythia::PtHardMin", "Signal::Pythia::PtHardMax", "Elapsed time"]
        vars = ["n_events", "random_seed", "min_pt", "max_eta", "jet_min_pt", "jet_max_eta", "jet_min_area", "jet_R",
                "multiplicity", "event_plane1", "event_plane2", "event_plane3", "event_plane4",
                "jet_min_pt", "jet_rotate", "jet_v2", "jet_v3", "jet_v4",
                "pthard_min", "pthard_max", "elapsed_time"]
        for val in vars:
                logfile_details[val] = []

        for key in logfiles:
                all_logfiles.append(logfiles[key])
        
        for logfile in all_logfiles:
                with open(logfile, "r") as f:
                        lines = f.readlines()
                        for line in lines:
                                for key in keys_to_find:
                                        if key in line:
                                                if "Elapsed time" in key:
                                                        val = line.split(":")[1].strip()
                                                        val = val.replace("s", "")
                                                        logfile_details["elapsed_time"].append(float(val))
                                                else:
                                                        val = line.split("==")[1].strip()
                                                        if "#" in val:
                                                                val = val.split("#")[0].strip()
                                                        
                                                        #remove whitespace
                                                        logfile_details[vars[keys_to_find.index(key)]].append(val)

        keys_that_should_match = ["n_events", "min_pt", "max_eta", "jet_min_pt", "jet_max_eta",
                                  "jet_min_area", "jet_R", "multiplicity", "event_plane1", "event_plane2", "event_plane3", "event_plane4",
                                    "jet_min_pt", "jet_rotate", "jet_v2", "jet_v3", "jet_v4", "pthard_min", "pthard_max"]
        keys_that_should_diff = ["random_seed"]
        keys_that_dont_matter = ["elapsed_time"]

        dets = {}
        for key in logfile_details:
                dets[key] = np.unique(logfile_details[key])

        report={}
        has_error = False
        for key in dets:
                report_for_key = ""
                if key in keys_that_should_match:
                        if len(dets[key]) > 1:
                                report_for_key = "Mismatch in {0} for jobs: {1}".format(key, dets[key])
                                report_for_key += "\n"
                                report_for_key +="Files with mismatch: \n"
                                for val in dets[key]:
                                        for key2 in logfile_details:
                                                if logfile_details[key2] == val:
                                                        report_for_key += "{0}\n".format(key2)
                                report[key] = report_for_key
                                has_error = True
                if key in keys_that_should_diff:
                        if len(dets[key]) > len(all_logfiles):
                                has_error = True
                                matching_vals = []
                                for val in dets[key]:
                                        if len(logfile_details[key]) == logfile_details[key].count(val):
                                                matching_vals.append(val)
                                matching_files = []
                                for val in matching_vals:
                                        for key2 in logfile_details:
                                                if logfile_details[key2] == val:
                                                        matching_files.append(key2)

                                report_for_key = "match in {0} for jobs: {1}".format(key, matching_vals)
                                report_for_key += "\n"
                                report_for_key +="Files with match: \n"
                                for val in matching_files:
                                        report_for_key += "{0}\n".format(val)
                                
                                report[key] = report_for_key

                if key in keys_that_dont_matter:

                        report_for_key = "Average {0}: {1}\n".format(key, np.mean(dets[key]))
                        report_for_key += "Std Dev {0}: {1}\n".format(key, np.std(dets[key]))
                        report[key] = report_for_key
                



        return report, has_error

def CompareLogFiles(logfile_details):
        keys_that_should_match = ["n_events", "min_pt", "max_eta", "jet_min_pt", "jet_max_eta",
                                  "jet_min_area", "jet_R", "multiplicity", "event_plane1", "event_plane2", "event_plane3", "event_plane4",
                                    "jet_min_pt", "jet_rotate", "jet_v2", "jet_v3", "jet_v4", "pthard_min", "pthard_max"]
        keys_that_should_diff = ["random_seed"]
        keys_that_dont_matter = ["elapsed_time"]

        MASTER_KEYS = {}
        for match in keys_that_should_match:
                MASTER_KEYS[match] = [] 
        for diff in keys_that_should_diff:
                MASTER_KEYS[diff] = []
        for dont in keys_that_dont_matter:
                MASTER_KEYS[dont] = []

        for key in logfile_details:
                for match in keys_that_should_match:
                        if match in logfile_details[key]:
                                MASTER_KEYS[match].append(logfile_details[key][match])
                for diff in keys_that_should_diff:
                        if diff in logfile_details[key]:
                                MASTER_KEYS[diff].append(logfile_details[key][diff])
                for dont in keys_that_dont_matter:
                        if dont in logfile_details[key]:
                                MASTER_KEYS[dont].append(logfile_details[key][dont])
        
        NUNIQUE = {}
        for key in MASTER_KEYS:
                NUNIQUE[key] = len(np.unique(MASTER_KEYS[key]))

        AMALGAMATE = {}
        for match in keys_that_should_match:
                AMALGAMATE[match] = {}
                MATCH = {}
                MATCH["nunique"] = NUNIQUE[match]
                MATCH["values"] = {}
                for val in np.unique(MASTER_KEYS[match]):
                        MATCH["values"][val] = []
                        for key in logfile_details:
                                if logfile_details[key][match] == val:
                                        MATCH["values"][val].append(key)
                AMALGAMATE[match] = MATCH

        for diff in keys_that_should_diff:
                AMALGAMATE[diff] = {}
                DIFF = {}
                DIFF["nunique"] = NUNIQUE[diff]
                DIFF["values"] = {}
                for val in np.unique(MASTER_KEYS[diff]):
                        DIFF["values"][val] = []
                        for key in logfile_details:
                                if logfile_details[key][diff] == val:
                                        DIFF["values"][val].append(key)
                AMALGAMATE[diff] = DIFF

        for dont in keys_that_dont_matter:
                AMALGAMATE[dont] = {}
                DONT = {}
                vals = np.array(MASTER_KEYS[dont], dtype=np.float64)
                DONT["average"] = np.mean(vals)
                DONT["stddev"] = np.std(vals)
                print(MASTER_KEYS[dont])
                AMALGAMATE[dont] = DONT

        return AMALGAMATE

if __name__ == "__main__":
        
        # parse the command line arguments
        parser = argparse.ArgumentParser(description='Scan error files and report non-empty ones')
        parser.add_argument('-f', '--file', help='Job file to scan', required=True)
        args = parser.parse_args()

        # get all job files with the base name in the job directory
        QA_FILE = args.file
        QA_data = ParseQAFile(QA_FILE)
        FILEBASE= os.path.basename(QA_FILE).split(".")[0]
        pwd = os.getcwd()
        echo = os.system("echo {0}".format(pwd))
        REPORT_FILE= "{0}/reports/{1}_report.txt".format(pwd, FILEBASE)
        print("Scanning {0}".format(QA_FILE))
        print("Creating report {0}".format(REPORT_FILE))
        
        # if not os.path.exists("{0}/reports".format(pwd)):
        #         os.makedirs("{0}/reports".format(pwd), exist_ok=True)
        
        errfiles = {}
        logfiles = {}
        output_files_sig = {}
        output_files_bkgd = {}
        all_logfiles = []
        all_jobfiles = []
        for ijob in range(QA_data["n_jobs"]):
                jobid = "Job{0}".format(ijob)
                jobname = QA_data[jobid]["jobname"]
                jobfile = QA_data[jobid]["jobfile"]
                outfile_sig = QA_data[jobid]["outfile_sig"]
                outfile_bkgd = QA_data[jobid]["outfile_bkgd"]
                errfile = QA_data[jobid]["errfile"]
                logfile = QA_data[jobid]["logfile"]
                if '.e' in logfile:
                        logfile = logfile.replace(".e", ".o")
                jobnumber = QA_data[jobid]["jobnumber"]
                errfiles[jobid] = errfile
                logfiles[jobid] = logfile
                all_logfiles.append(logfile)
                all_jobfiles.append(jobfile)
                output_files_sig[jobid] = outfile_sig
                output_files_bkgd[jobid] = outfile_bkgd


        # scan the error files
        none_empty_errfiles = ScanErrorFiles(errfiles)

        # scan the output files
        bad_root_files = ScanOutputFiles(output_files_sig)
        bad_root_files.update(ScanOutputFiles(output_files_bkgd))

        logfile_details, dont_archive = ScanLogFiles(logfiles)
        with open(REPORT_FILE, "w") as f:
               
                f.write("=====================================================\n")
                f.write("Config File: {0}\n".format(QA_data["config_file"]))
                f.write("Number of jobs: {0}\n".format(QA_data["n_jobs"]))
                f.write("Number of events per job: {0}\n".format(QA_data["n_events"]))
                f.write("Output path: {0}\n".format(QA_data["output_path"]))
                f.write("=====================================================\n")
                f.write("None empty error files ({0}):\n".format(len(none_empty_errfiles)))
                for errfile in none_empty_errfiles:
                       f.write("{0}\n".format(errfile))
                f.write("=====================================================\n")
                f.write("Bad root files ({0}):\n".format(len(bad_root_files)))
                for rootfile in bad_root_files:
                        f.write("{0}\n".format(rootfile))
                f.write("=====================================================\n")

                # for key in DETS:
                #         subkeys = DETS[key].keys()
                #         if "nunique" in DETS[key]:
                #                 vals = DETS[key]["values"]
                #                 f.write("{0} ({1}): {2} \n".format(key, DETS[key]["nunique"],DETS[key]["values"]))
                #         if "average" in DETS[key]:
                #                 f.write("{0}\n".format(key))
                #                 f.write("\tAverage: {0}\n".format(DETS[key]["average"]))
                #                 f.write("\tStd Dev: {0}\n".format(DETS[key]["stddev"]))
                #         f.write("\n")
                # f.write("=====================================================\n")
                for key in logfile_details:
                        f.write("{0}\n".format(key))
                        f.write("{0}\n".format(logfile_details[key]))
                        f.write("\n")
                f.write("=====================================================\n")
        

        print("Report written to {0}".format(REPORT_FILE))
        # if dont_archive:
        #         sys.exit()
        # else:
        print("Archiving log files")
        log_archive = "{0}/arxiv/{1}_logs.tar".format(pwd, FILEBASE)
        if not os.path.exists("{0}/arxiv".format(pwd)):
                os.makedirs("{0}/arxiv".format(pwd), exist_ok=True)
        os.system("tar -cf {0} {1}".format(log_archive, " ".join(logfiles.values())))

        # rm log files
        for logfile in all_logfiles:
                os.remove(logfile)
        # rm job files
        for jobfile in all_jobfiles:
                os.remove(jobfile)

        # remove QA files
        os.remove(QA_FILE)

        print("Done")

        sys.exit()