#include <rangesearch.h>
#include "utils.h"
//#include <readdata.h>
using namespace nanoflann;

void pdist2(std::vector<std::array<float,3>> &data, std::array<float,3> &point, std::vector<float> &distmat){
    // instead of computing sqrt() of distance, compute the 2nd power of radius once and compare it again and again which is faster
    distmat.resize(data.size());
    for (unsigned int i=0;i<data.size();i++){
        distmat[i] = std::pow(data[i][0]-point[0],2) + std::pow(data[i][1]-point[1],2) + std::pow(data[i][2]-point[2],2);
    }
}

void convertTocloud(PointCloud<float> &point, std::vector<std::array<float,3>> &data)
{
	// Generating Random Point Cloud
	point.pts.resize(data.size());
	for (size_t i = 0; i < data.size(); i++)
	{
		point.pts[i].x = data[i][0];
		point.pts[i].y = data[i][1];
		point.pts[i].z = data[i][2];
	}
}

void RangeSearch(std::vector<std::vector<unsigned int>> &indices, std::vector<std::array<float,3>> &data, float radius)
{
    indices.resize(data.size());
    PointCloud<float> cloud;
    convertTocloud(cloud, data);
    // Nanoflann library gets r^2 as the radius
    radius = std::pow(radius,2);

    typedef KDTreeSingleIndexAdaptor<
		L2_Simple_Adaptor<float, PointCloud<float> > ,
		PointCloud<float>,
		3 /* dim */
		> my_kd_tree_t;

	my_kd_tree_t   index(3 /*dim*/, cloud, KDTreeSingleIndexAdaptorParams(10 /* max leaf */) );
	index.buildIndex();

    //const float radius = static_cast<float>(0.1);
    std::vector<std::pair<unsigned long int,float> >   ret_matches;
    //std::vector<std::vector<std::pair<unsigned long int,float> >> vec;

    nanoflann::SearchParams params;
    //params.sorted = false;
    int nMatches=0;
    for (unsigned int idx = 0; idx<data.size(); idx++)
    {
        nMatches = index.radiusSearch(&data[idx][0], radius, ret_matches, params);
        indices[idx].resize(nMatches);
        for (int j = 1; j < nMatches; j++)// I remove the first idx since it is the index of the selected point itself
            indices[idx][j] = ret_matches[j].first;

    }
}

std::vector<std::array<float,3> > readdata(std::string name_of_file){
    std::ifstream infile(name_of_file);
    if (!infile)
    {
       std::cerr << "Could not open file!"  <<  std::endl;
    }
    std::string line;
    int i = 0;
    std::array<float,3> v;
    int linelength, lines, linee;
    std::vector<std::array<float,3> > data;
    while(getline(infile,line)){
        linelength = line.length();
        lines = 0;
        i = 0;
        linee = line.find(',');
        while(linee>0){
            v[i] = stof(line.substr(lines,linee));
            i++;
            lines = linee+1;
            linee = line.find(',',linee+2);
        }
        v[2] = stof(line.substr(lines,linelength));
        data.push_back(v);
    }
    return data;
}

void writeCSV2d(std::string name_of_file,std::vector<std::array<float,3>> vec)
{
    std::ofstream outputfile;
    outputfile.open(name_of_file);
    for(unsigned int i=0;i<vec.size();i++){
        for(unsigned int j=0;j<2;j++){
            outputfile<<vec[i][j]<<',';
        }
        outputfile<<vec[i][2]<<'\n';
    }
    outputfile.close();
}
