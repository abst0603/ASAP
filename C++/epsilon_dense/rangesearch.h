#include <vector>
#include <array>
#include <numeric> //for std::accumulate
#include <cmath>
#include <iostream>
#include <string>
#include <iostream>
#include <readdata.h>
#ifndef PREPROCESS_H
#define PREPROCESS_H

void pdist2(std::vector<std::array<float,3>> &data, std::array<float,3> &point, std::vector<float> &distmat);
void rangesearchfastV(std::vector<std::array<float,3>> &data, std::vector<std::vector<int>> &rmlist, float radius);
std::vector<std::array<float,3> > readdata(std::string name_of_file);

#endif
