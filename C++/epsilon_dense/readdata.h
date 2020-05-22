#include <vector>
#include <array>
#include <string>
#include <fstream>
#ifndef READDATA_H
#define READDATA_H

std::vector<std::array<float,3>> readdata(std::string name_of_file);
void writeCSV2d(std::string name_of_file,std::vector<std::array<float,3>> vec);

#endif
