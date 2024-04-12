#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include "omp.h"
#include <map>
#include <dirent.h>


int verbose = 1;

#define NUM_THREADS 12

struct AirQualityData {
    double latitude;
    double longitude;
    std::string timestamp;
    std::string pollutant;
    std::string unit;
    double concentration;
    double rawConcentration;
    int AQI;
    int category;
    std::string sName;
    std::string sAgency;
    std::string AQS;
    std::string fullAQS;
};

std::vector<std::string> find_files(const std::string& directory) {
    std::vector<std::string> filenames;
    DIR *dir;
    struct dirent *ent;
    if ((dir = opendir(directory.c_str())) != NULL) {
        while ((ent = readdir(dir)) != NULL) {
            if (ent->d_type == DT_DIR) {
                // Ignore "." and ".." directories
                if (std::strcmp(ent->d_name, ".") == 0 || std::strcmp(ent->d_name, "..") == 0) {
                    continue;
                }
                // Recursively search subdirectories
                std::string subdir = directory + "/" + ent->d_name;
                std::vector<std::string> subfiles = find_files(subdir);
                filenames.insert(filenames.end(), subfiles.begin(), subfiles.end());
            } else {
                std::string filename = ent->d_name;
                std::string foldername = directory.substr(directory.find_last_of('/') + 1);
                if (filename.find(foldername) == 0) {
                    filenames.push_back(directory + "/" + filename);
                }
            }
        }
        closedir(dir);
    } else {
        std::cerr << "Error opening directory: " << directory << std::endl;
    }
    
    
    
    return filenames;
}

std::vector<AirQualityData> read_csv_data(const std::string& filename) {
    std::vector<AirQualityData> data;

    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return data;
    }

    std::string line;

    while (std::getline(file, line)) {
        std::istringstream iss(line);
        std::string token;

        AirQualityData entry;

        // Read comma-separated values from the line


        std::getline(iss, token, ',');
        token = token.substr(1, token.size() - 2);
        try {
            entry.latitude = std::stod(token);
        } catch (const std::invalid_argument& e) {
            
            if (verbose > 2){
                std::cerr << token << " -> Invalid argument latitude: " << e.what() << std::endl;
            }
            entry.latitude =  0.0; 
        }
        

        std::getline(iss, token, ',');
        token = token.substr(1, token.size() - 2);
        try {
            entry.longitude = std::stod(token);
        } catch (const std::invalid_argument& e) {
            
            if (verbose > 2){
                std::cerr << token << " -> Invalid argument longitude: " << e.what() << std::endl;
            }
            entry.longitude =  0.0; 
        }

        std::getline(iss, entry.timestamp , ',');
        entry.timestamp = entry.timestamp.substr(1, entry.timestamp.size() - 2);

        std::getline(iss, entry.pollutant, ',');
        entry.pollutant = entry.pollutant.substr(1, entry.pollutant.size() - 2);
        

        std::getline(iss, token, ',');
        token = token.substr(1, token.size() - 2);
        try {
            entry.concentration = std::stod(token);
        } catch (const std::invalid_argument& e) {
            
            if (verbose > 2){
                std::cerr << token << " -> Invalid argument concentration: " << e.what() << std::endl;
            }
            entry.concentration =  0.0; 
        }
        if (entry.concentration == -999) {
            continue; // Skip the remaining part of the loop iteration
        }

        std::getline(iss, entry.unit, ',');
        entry.unit = entry.unit.substr(1, entry.unit.size() - 2);

        std::getline(iss, token, ',');
        token = token.substr(1, token.size() - 2);
         try {
            entry.rawConcentration = std::stod(token);
        } catch (const std::invalid_argument& e) {
            
            if (verbose > 2){
                std::cerr << token << " -> Invalid argument rawConcentration: " << e.what() << std::endl;
            }
            entry.rawConcentration =  0.0; 
        }
        if (entry.rawConcentration == -999) {
            continue; // Skip the remaining part of the loop iteration
        }

        std::getline(iss, token, ',');
        token = token.substr(1, token.size() - 2);
         try {
            entry.AQI = std::stod(token);
        } catch (const std::invalid_argument& e) {
            
            if (verbose > 2){
                std::cerr << token << " -> Invalid argument AQI: " << e.what() << std::endl;
            }
            entry.AQI =  0.0; 
        }
        if (entry.AQI ==999) {
            continue; // Skip the remaining part of the loop iteration
        }

        std::getline(iss, token, ',');
        token = token.substr(1, token.size() - 2);
         try {
            entry.category = std::stod(token);
        } catch (const std::invalid_argument& e) {
            
            if (verbose > 2){
                std::cerr << token << " -> Invalid argument category: " << e.what() << std::endl;
            }
            entry.category =  0.0; 
        }

        if (entry.category == -999) {
            continue; // Skip the remaining part of the loop iteration
        }

        std::getline(iss, entry.sName, ',');
        entry.sName = entry.sName.substr(1, entry.sName.size() - 2);


        std::getline(iss, entry.sAgency, ',');
        entry.sAgency = entry.sAgency.substr(1, entry.sAgency.size() - 2);


        std::getline(iss, entry.AQS, ',');
        entry.AQS = entry.AQS.substr(1, entry.AQS.size() - 2);


        std::getline(iss, entry.fullAQS);
        entry.fullAQS = entry.fullAQS.substr(1, entry.fullAQS.size() - 2);


        data.push_back(entry);
    }

    return data;
}

std::map<std::string, double> calculate_average(const std::vector<AirQualityData>& data) {
    std::map<std::string, double> sums;
    std::map<std::string, int> counts;
    std::map<std::string, double> averages;


    for (const auto& point : data) {
        sums[point.pollutant] += point.concentration;
        counts[point.pollutant]++;
    }

    for (const auto& pair : sums) {
        averages[pair.first] = pair.second / counts[pair.first];
    }

    return averages;
}

void write_to_file(std::string data){
    std::ofstream outfile("pollutant_avg.txt", std::ios::app);
    if (!outfile.is_open()) {
        std::cerr << "Error: Unable to open file for writing: pollutant_avg.txt" << std::endl;
        return;
    }
    //outfile<< "Average concentration for date "<< extract_date(filename)<< std::endl;
    
    outfile << data << std::endl;
    

    outfile.close();

}

void execution_time(int rank, double start, double end, int records, int totprocess) {
    std::string filename = "omp_execution_time_for_" + std::to_string(totprocess) + "_processes.txt";
    std::ofstream metrics(filename, std::ios::app);
    
    // Check if the file opened successfully
    if (metrics.is_open()) {
        // Calculate the duration of execution
        double duration = end - start;
        
        // Calculate records processed per second
        double recordsPerSecond = records / duration;
        
        // Write the metrics to the file
        metrics << "Thread ID: " << rank << ", Records processed: " << records <<", Execution time: " << duration << "s, Records per second: " << recordsPerSecond << std::endl;
        
        // Close the file
        metrics.close();
    } else {
        // Print an error message if file opening fails
        std::cerr << "Failed to open metrics file for process " << rank << std::endl;
    }
}

int main(int argc, char *argv[]) {

    std::map<std::string, double> partial_avgs;
   

    std::string directory = argv[1];
    omp_set_num_threads(NUM_THREADS);

    
    if (verbose > 2){
        std::cout<<"Directory : " <<directory<< std::endl;
    }
    std::vector<std::string> filenames;
    filenames = find_files(directory);

    if (verbose >4){
        for (int i = 0; i < filenames.size(); ++i){
            std::cout << "File " << i <<" "<< filenames[i]<<std::endl;
        }
    }
    
    std::ofstream outfile("pollutant_avg.txt", std::ios::out);
    outfile.close();
    

    #pragma omp parallel 
    {
        int id = omp_get_thread_num();

        int totRecords = 0;
        auto startTime = std::chrono::high_resolution_clock::now();
        for (int i = id; i < filenames.size(); i+=NUM_THREADS) {
            std::string filename = filenames[i];
            std::vector<AirQualityData> data = read_csv_data(std::string(filename));
            
            
            // Process the data (example: print the number of entries)
            if (verbose>1){
                std::cout << "ID " << id << " processed file " << filename << ". Number of entries: "<< data.size() << std::endl;
            }

            totRecords+=data.size();
            partial_avgs = calculate_average(data);
        }
        auto endTime = std::chrono::high_resolution_clock::now();
        execution_time(id,  
                std::chrono::duration<double>(startTime.time_since_epoch()).count(),  
                std::chrono::duration<double>(endTime.time_since_epoch()).count(), 
                totRecords, NUM_THREADS);   
    }

    double total_average = 0.0;
    for (std::map<std::string, double>::iterator it = partial_avgs.begin(); it != partial_avgs.end(); ++it) {
        total_average += it->second;
    }
    total_average /= filenames.size(); // Divide by the number of files

    std::cout << "Total average: " << total_average << std::endl;

    return 0;
    
}

