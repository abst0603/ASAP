#include <rangesearch.h>
#include "utils.h"
#include "KDTreeVectorOfVectorsAdaptor.h"
using namespace nanoflann;

void pdist2(std::vector<std::vector<float>> &data, std::vector<float> &point, std::vector<float> &distmat){
    // instead of computing sqrt() of distance, compute the 2nd power of radius once and compare it again and again which is faster
    distmat.resize(data.size());
    for (unsigned int i=0;i<data.size();i++){
        distmat[i] = 0;
        for(unsigned int j=0;j<data[0].size();j++)
            distmat[i] += std::pow(data[i][j]-point[j],2);
    }
}

void RangeSearch(std::vector<std::vector<unsigned int>> &indices, std::vector<std::vector<float>> &data, float radius)
{
    indices.resize(data.size());
    PointCloud<float> cloud;
    // Nanoflann library gets r^2 as the radius
    radius = std::pow(radius,2);

    typedef KDTreeVectorOfVectorsAdaptor< std::vector<std::vector<float> >, float >  my_kd_tree_t;

	my_kd_tree_t   mat_index(data[0].size()  /*dim*/, data, 10 /* max leaf */ );
	mat_index.index->buildIndex();

    //const float radius = static_cast<float>(0.1);
    std::vector<std::pair<unsigned long int,float> >   ret_matches;
    //std::vector<std::vector<std::pair<unsigned long int,float> >> vec;

    nanoflann::SearchParams params;
    //params.sorted = false;
    int nMatches=0;
    for (unsigned int idx = 0; idx<data.size(); idx++)
    {
        nMatches = mat_index.index->radiusSearch(&data[idx][0], radius, ret_matches, params);
        indices[idx].resize(nMatches);
        for (int j = 1; j < nMatches; j++)// I remove the first idx since it is the index of the selected point itself
            indices[idx][j] = ret_matches[j].first;

    }
}

void readdata(const int dim, std::vector<std::vector<float> > &data,std::string name_of_file){
    std::ifstream infile(name_of_file);
    if (!infile)
    {
       std::cerr << "Could not open file!"  <<  std::endl;
    }
    std::string line;
    int i = 0;
    std::vector<float> v(dim);
    int linelength, lines, linee;
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
        v[dim-1] = stof(line.substr(lines,linelength));
        data.push_back(v);
    }

    if(data[0].size()!=dim)
        std::cerr<<"The number of dimensions of the input data is not equal to dim which is set in code.\n";
}

void writeCSV2d(std::string name_of_file,std::vector<std::vector<float>> vec)
{
    std::ofstream outputfile;
    unsigned int temp_size = vec[0].size();
    outputfile.open(name_of_file);
    for(unsigned int i=0;i<vec.size();i++){
        for(unsigned int j=0;j<temp_size-1;j++){
            outputfile<<vec[i][j]<<',';
        }
        outputfile<<vec[i][temp_size-1]<<'\n';
    }
    outputfile.close();
}
