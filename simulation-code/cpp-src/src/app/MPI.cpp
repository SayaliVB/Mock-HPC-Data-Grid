#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <mpi.h>
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
            
            if (verbose == 2){
                std::cerr << token << " -> Invalid argument latitude: " << e.what() << std::endl;
            }
            entry.latitude =  0.0; 
        }
        

        std::getline(iss, token, ',');
        token = token.substr(1, token.size() - 2);
        try {
            entry.longitude = std::stod(token);
        } catch (const std::invalid_argument& e) {
            
            if (verbose == 2){
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
            
            if (verbose == 2){
                std::cerr << token << " -> Invalid argument concentration: " << e.what() << std::endl;
            }
            entry.concentration =  0.0; 
        }


        std::getline(iss, entry.unit, ',');
        entry.unit = entry.unit.substr(1, entry.unit.size() - 2);

        std::getline(iss, token, ',');
        token = token.substr(1, token.size() - 2);
         try {
            entry.rawConcentration = std::stod(token);
        } catch (const std::invalid_argument& e) {
            
            if (verbose == 2){
                std::cerr << token << " -> Invalid argument rawConcentration: " << e.what() << std::endl;
            }
            entry.rawConcentration =  0.0; 
        }

        std::getline(iss, token, ',');
        token = token.substr(1, token.size() - 2);
         try {
            entry.AQI = std::stod(token);
        } catch (const std::invalid_argument& e) {
            
            if (verbose == 2){
                std::cerr << token << " -> Invalid argument AQI: " << e.what() << std::endl;
            }
            entry.AQI =  0.0; 
        }

        std::getline(iss, token, ',');
        token = token.substr(1, token.size() - 2);
         try {
            entry.category = std::stod(token);
        } catch (const std::invalid_argument& e) {
            
            if (verbose == 2){
                std::cerr << token << " -> Invalid argument category: " << e.what() << std::endl;
            }
            entry.category =  0.0; 
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


double calculate_average(const std::vector<AirQualityData>& data) {
    double total_concentration = 0.0;
    for (const auto& point : data) {
        total_concentration += point.concentration;
    }
    return total_concentration / data.size();
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

    std::vector<double> all_partial_averages(size);
    double partial_average;

    if (rank == 0) {
        // Rank 0 reads the filenames and sends them to other ranks
        std::vector<std::string> filenames;

        std::string directory = argv[1];
        
        if (verbose == 2){
            std::cout<<"Directory : " <<directory<< std::endl;
        }
        filenames = find_files(directory);

        std::cout << "Rank " << rank << " processed files " << filenames.size()  << std::endl;


        // Send filenames to other ranks
        for (int i = 1; i < size ; ++i) {
            //std::cout << "Rank 0 sending file " << filenames[i] << "to  "<< i << std::endl;
            MPI_Send(filenames[i].c_str(), filenames[i].size() + 1, MPI_CHAR, i, 0, MPI_COMM_WORLD);
            //std::cout << "Rank 0 sent file " << filenames[i] << "to  "<< i << std::endl;
        }
    } else{

        // Receive filenames and process the corresponding files
        char filename[256];
        MPI_Recv(filename, 256, MPI_CHAR, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        std::vector<AirQualityData> data = read_csv_data(std::string(filename));

        // Process the data (example: print the number of entries)
        std::cout << "Rank " << rank << " processed file " << filename << ". Number of entries: "<< data.size() << std::endl;

        // Calculate partial average
        partial_average = calculate_average(data);
    }

    MPI_Gather(&partial_average, 1, MPI_DOUBLE, all_partial_averages.data(), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        double total_concentration = 0.0;
        int total_data_points = 0;
        for (const auto& avg : all_partial_averages) {
            total_concentration += avg;
            total_data_points++;
        }
        double final_average = total_concentration / total_data_points;
        std::cout << "Final average concentration: " << final_average << std::endl;
    }

    MPI_Finalize();
    return 0;
}
