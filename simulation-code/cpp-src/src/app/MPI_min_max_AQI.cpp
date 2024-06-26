#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <mpi.h>
#include <map>
#include <dirent.h>


int verbose = 1;

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
    /*
    if ((dir = opendir(directory.c_str())) != NULL) {
        while ((ent = readdir(dir)) != NULL) {
            std::string filename = ent->d_name;
            if (filename.find("20200810-") == 0) {
                filenames.push_back(directory + "/" + filename);
            }
        }
        closedir(dir);
    } else {
        std::cerr << "Error opening directory: " << directory << std::endl;
    }
    */
    

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

        if (entry.AQI == 999) {
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
        if (entry.sName==""){
            std::cout<<"Empty";
            continue;
        }


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

std::string extract_date(const std::string& filepath) { 
    // Find the position of the last occurrence of '/'
    size_t last_slash_pos = filepath.rfind('/');

    size_t last_dot_pos = filepath.rfind('.');

    // Extract the substring representing the date-time
    std::string date_time = filepath.substr(last_slash_pos+1, last_dot_pos - last_slash_pos - 1);
    
    return date_time;
}


// int order: order =1  desc, order =0  asc
std::vector<AirQualityData> find_min_site_names(const std::vector<AirQualityData>& aqiData, int n, int order) {
    std::vector<AirQualityData> minSiteNames;
    std::vector<AirQualityData> sortedAQIData = aqiData;
    std::sort(sortedAQIData.begin(), sortedAQIData.end(), [order](const AirQualityData& a, const AirQualityData& b) {
        return (order  == 0)? a.AQI < b.AQI : a.AQI > b.AQI;
    });

    for (size_t i = 0; i < n && i < sortedAQIData.size(); ++i) {
        minSiteNames.push_back(sortedAQIData[i]);
    }

    return minSiteNames;
}

std::vector<char> serialize_air_quality_data(const std::vector<AirQualityData>& data) {
    std::vector<char> serialized_data;

    for (const auto& entry : data) {
        // Serialize AQI
        serialized_data.insert(serialized_data.end(), reinterpret_cast<const char*>(&entry.AQI), reinterpret_cast<const char*>(&entry.AQI) + sizeof(entry.AQI));

        // Serialize sName (including null terminator)
        serialized_data.insert(serialized_data.end(), entry.sName.begin(), entry.sName.end());
        serialized_data.push_back('\0');
    }

    return serialized_data;
}

// Function to deserialize a byte array into AQI and sName fields of AirQualityData objects
std::vector<AirQualityData> deserialize_air_quality_data(const std::vector<char>& serialized_data) {
    std::vector<AirQualityData> data;

    size_t offset = 0;
    while (offset < serialized_data.size()) {
        AirQualityData entry;

        // Deserialize AQI
        std::memcpy(&entry.AQI, &serialized_data[offset], sizeof(entry.AQI));
        offset += sizeof(entry.AQI);

        // Deserialize sName
        while (serialized_data[offset] != '\0') {
            entry.sName.push_back(serialized_data[offset]);
            ++offset;
        }
        // Move offset past the null terminator
        ++offset;

        // Add the deserialized entry to the data vector
        data.push_back(entry);
    }

    return data;
}

void execution_time(int rank, double start, double end, int records, int totprocess) {
    std::string filename = "min_values_exe_time_for_" + std::to_string(totprocess) + "_processes.txt";
    std::ofstream metrics(filename, std::ios::app);
    
    // Check if the file opened successfully
    if (metrics.is_open()) {
        // Calculate the duration of execution
        double duration = end - start;
        
        // Calculate records processed per second
        double recordsPerSecond = records / duration;
        
        // Write the metrics to the file
        metrics << "Process Rank: " << rank << ", Records processed: " << records <<", Execution time: " << duration << "s, Records per second: " << recordsPerSecond << std::endl;
        
        // Close the file
        metrics.close();
    } else {
        // Print an error message if file opening fails
        std::cerr << "Failed to open metrics file for process " << rank << std::endl;
    }
}


int main(int argc, char *argv[]) {
    MPI_Init(NULL, NULL);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    //std::cout<<size;

    if (size < 2) {
        std::cerr << "This program needs at least 2 MPI processes." << std::endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    std::vector<AirQualityData> all_minSiteNames; 

    int noofrecords = std::atoi(argv[2]);
    int order = std::atoi(argv[3]);


    if (rank == 0) {
        // Rank 0 reads the filenames and sends them to other ranks

        std::string directory = argv[1];
        
        if (verbose > 2){
            std::cout<<"Directory : " <<directory<< std::endl;
        }
        std::vector<std::string> filenames;
        filenames = find_files(directory);

        std::cout << "Rank " << rank << " processed files " << filenames.size()  << std::endl;

        if (verbose >4){
            for (int i = 0; i < filenames.size(); ++i){
                std::cout << "File " << i <<" "<< filenames[i]<<std::endl;
            }
        }

        // Send filenames to other ranks
        for (int i = 1; i < size ; ++i) {
            MPI_Send(filenames[i-1].c_str(), filenames[i-1].size() + 1, MPI_CHAR, i, 0, MPI_COMM_WORLD);
            if (verbose >3)
                std::cout << "Rank 0 sent file " << filenames[i-1] << "to  "<< i << std::endl;
        }
        
        for (int i = size-1; i < filenames.size(); ++i) {
            int ranktosend;
            MPI_Recv(&ranktosend, 1, MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            int vector_size;
            MPI_Recv(&vector_size, 1, MPI_INT, ranktosend, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            // Receive the entire vector
            std::vector<char> received_data(vector_size);
            MPI_Recv(received_data.data(), vector_size, MPI_CHAR, ranktosend, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            std::vector<AirQualityData> received_vector = deserialize_air_quality_data(received_data);
            all_minSiteNames.insert(all_minSiteNames.end(), received_vector.begin(), received_vector.end());

            
            if (verbose >3)
                std::cout << "Rank 0 remaining sending file " << filenames[i] << "to  "<< ranktosend << std::endl;
            MPI_Send(filenames[i].c_str(), filenames[i].size() + 1, MPI_CHAR, ranktosend, 0, MPI_COMM_WORLD);
        }
        // Signal termination to other processes
        for (int i = 1; i < size; ++i) {
            int ranktosend;
            MPI_Recv(&ranktosend, 1, MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);


            int vector_size;
            MPI_Recv(&vector_size, 1, MPI_INT, ranktosend, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            // Receive the entire vector
            std::vector<char> received_data(vector_size);
            MPI_Recv(received_data.data(), vector_size, MPI_CHAR, ranktosend, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            
            std::vector<AirQualityData> received_vector = deserialize_air_quality_data(received_data);
            all_minSiteNames.insert(all_minSiteNames.end(), received_vector.begin(), received_vector.end());

            if (verbose >3)
                std::cout << "Rank 0 sending empty file to  "<< ranktosend << std::endl;
            char empty_filename[1] = {'\0'};
            MPI_Send(empty_filename, 1, MPI_CHAR, ranktosend, 0, MPI_COMM_WORLD);
            if (verbose >3)
                std::cout << "Rank 0 sent empty file to  "<< ranktosend << std::endl;
        }

        
    } else{
        int i =0;
        int totRecords = 0;
        auto startTime = std::chrono::high_resolution_clock::now();
        // Receive filenames and process the corresponding files
        while(true){
            char filename[256];
            MPI_Recv(&filename, 256, MPI_CHAR, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            if (verbose> 3 && filename[0] == '\0'){
                std::cout << " recieved file " << filename << std::endl;
            }
            if (filename[0] == '\0') // Termination signal
            {
                auto endTime = std::chrono::high_resolution_clock::now();
                execution_time(rank,  
                                std::chrono::duration<double>(startTime.time_since_epoch()).count(),  
                                std::chrono::duration<double>(endTime.time_since_epoch()).count(), 
                                totRecords, size);
                break;
            }

            std::vector<AirQualityData> data = read_csv_data(std::string(filename));

            // Process the data (example: print the number of entries)
            if (verbose>1){
                std::cout << "Rank " << rank << " processed file " << filename << ". Number of entries: "<< data.size() << std::endl;
            }

            totRecords+=data.size();
            
            
            MPI_Send(&rank, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
            std::vector<AirQualityData> temp = find_min_site_names(data, noofrecords, order);

            std::vector<char> serialized_data = serialize_air_quality_data(temp);

            int vector_size = serialized_data.size();
            MPI_Send(&vector_size, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);

            // Send the entire vector
            MPI_Send(serialized_data.data(), vector_size, MPI_CHAR, 0, 0, MPI_COMM_WORLD);
            i++;
        }

    }

    if (rank == 0) {
        // Print the result
        all_minSiteNames = find_min_site_names(all_minSiteNames,noofrecords, order);
        std::cout << "5 Minimum Site Names for Maximum AQI Values:" << all_minSiteNames.size() << std::endl;
        std::ofstream outfile("min_max_AQI.txt", std::ios::app);
        if (outfile.is_open()) {
            for (const auto& siteName : all_minSiteNames) {
                std::cout << siteName.sName <<"  with AQI : " << siteName.AQI <<std::endl;
                outfile << siteName.sName<<"  with AQI : " << siteName.AQI << std::endl;
            }
        }
        else{
            std::cerr << "Error: Unable to open file for writing: min_max_AQI.txt" << std::endl;
        }

        outfile.close();

    }

    MPI_Finalize();
    return 0;
}
