#include <vector>
#include <array>
#include <numeric> //for std::accumulate
#include <cmath>
#include <iostream>
#include <string>
#include <fstream>
#include <nanoflann.hpp>

#ifndef RANGESEARCH_H
#define RANGESEARCH_H
void RangeSearch(std::vector<std::vector<unsigned int>> &indices, std::vector<std::vector<float>> &dist, std::vector<std::vector<float>> &data, float radius);
void readdata(const int dim, std::vector<std::vector<float> >  &data, std::string name_of_file);
void writeCSV2d(std::string name_of_file,std::vector<std::vector<float>> vec);

#endif
