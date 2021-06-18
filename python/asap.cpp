#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <Python.h>
namespace py = pybind11;

#include <iostream>
#include <array>
#include <vector>
#include <cmath>
#include <random>
#include <algorithm>
#include <ctime>
#include "KDTreeVectorOfVectorsAdaptor.h"
using namespace std;
typedef KDTreeVectorOfVectorsAdaptor< std::vector<std::vector<double> >, double >  my_kd_tree_t;
const int dim = 3; //You should set the value before compiling since nanoflann use this value for unrolling the loops.

using namespace nanoflann;

void RangeSearch(std::vector<std::vector<unsigned int>> &indices, std::vector<std::vector<double>> &dist, std::vector<std::vector<double>> &data, double radius)
{
    indices.resize(data.size());
    dist.resize(data.size());
    // Nanoflann library gets r^2 as the radius
    radius = std::pow(radius,2);
    typedef KDTreeVectorOfVectorsAdaptor< std::vector<std::vector<double> >, double >  my_kd_tree_t;
    my_kd_tree_t   mat_index(data[0].size()  /*dim*/, data, 10 /* max leaf */ );
    mat_index.index->buildIndex();

    std::vector<std::pair<unsigned long int,double> >   ret_matches;

    nanoflann::SearchParams params;
    //params.sorted = false;
    int nMatches=0;
    for (unsigned int idx = 0; idx<data.size(); idx++)
    {
        nMatches = mat_index.index->radiusSearch(&data[idx][0], radius, ret_matches, params);
        indices[idx].resize(nMatches-1);
        dist[idx].resize(nMatches-1);
        for (int j = 1; j < nMatches; j++)// I remove the first idx since it is the index of the selected point itself
            {indices[idx][j-1] = ret_matches[j].first;
            dist[idx][j-1] = ret_matches[j].second;
            }
    }
}

// On small scale datasets(a few ten thousands), the following version is faster than the next implementation, but on bigger one it may be different
void sample_by_covering_condition(std::vector<std::vector<double>> &data,std::vector<std::vector<double>> &samples, double radius2, int dim){
    std::vector<double> point(dim);
    while(data.size()>0){
        point = data[rand()%data.size()];
        samples.push_back(point);
        data.erase(std::remove_if(data.begin(),data.end(),
        //lambda expression for removing neighbors. This one works for Ndims data
        [&](std::vector<double> data_point){
        double temp_dist = 0;
        for(unsigned int cnt=0; cnt<dim;cnt++){temp_dist += std::pow(data_point[cnt]-point[cnt],2);}
        return temp_dist <= radius2;}
        ),data.end());
    }
    data.shrink_to_fit();
}

void find_nearest_neighbor(my_kd_tree_t &mat_index, std::vector<std::vector<double>> &data, std::vector<double> &newpos, int dim, int &ind){

    //    We only built the kdtree once and reuse it in all functions
    const size_t num_results = 1;
    std::vector<size_t>   ret_indexes(num_results);
    std::vector<double> out_dists_sqr(num_results);

    nanoflann::KNNResultSet<double> resultSet(num_results);

    resultSet.init(&ret_indexes[0], &out_dists_sqr[0] );
    mat_index.index->findNeighbors(resultSet, &newpos[0], nanoflann::SearchParams(10));
    ind = ret_indexes[0];

}

void check_covering_condition(my_kd_tree_t &mat_index, std::vector<std::vector<double>> &data,std::vector<std::vector<double>> &samples, double radius2, int dim){
    std::vector<unsigned int> indices;
    indices.resize(10*data.size());// just a big number so I now I don't need to resize again

//    We only built the kdtree once and reuse it in all functions
    std::vector<std::pair<unsigned long int,double> >   ret_matches;
    nanoflann::SearchParams params;
    params.sorted = false;
    int nMatches=0;
    unsigned int cnt = 0;
    for (unsigned int idx = 0; idx<samples.size(); idx++)
    {
        nMatches = mat_index.index->radiusSearch(&samples[idx][0], radius2, ret_matches, params);
        for (int j = 0; j < nMatches; j++)// I remove also the selected point itself --> sorted:false and starting from 0
            {indices[cnt] = ret_matches[j].first;
            cnt ++;}
    }
    indices.resize(cnt);
    std::sort(indices.begin(), indices.end());
    auto l_indices = std::unique(indices.begin(), indices.end());
    indices.erase(l_indices, indices.end());
    indices.shrink_to_fit();
    for (int i = indices.size()-1;i >= 0; --i)
    {
        unsigned int it = indices[i];
        data.erase(data.begin()+it);
    }
    data.shrink_to_fit();
}

void coulomb_law_movements(my_kd_tree_t &mat_index, std::vector<std::vector<double>> &data,std::vector<std::vector<double>> &samples, double radius, int dim, double lr){
    std::vector<std::vector<unsigned int>> idx;
    std::vector<std::vector<double>> dist;
    std::vector<std::vector<double>> nsamples(dim);
    std::vector<double> distmat;
    std::vector<double> sumtmp(dim,0);
    double tmpdist2 = 0;
    int tmpidx = 0;
    bool flag = false;
    double radius2 = std::pow(radius,2);
    RangeSearch(idx, dist, samples, 2*radius);
    for (unsigned int i=0; i<samples.size(); i++){
        if(idx[i].size() == 0)
            continue;
        nsamples.resize(idx[i].size(),sumtmp);
        for(unsigned int j = 0; j<idx[i].size();j++)
            nsamples[j] = samples[idx[i][j]];
        // Distance between neighboring samples
        distmat = dist[i];
        // Compute the displacement
        for(unsigned int j = 0; j<idx[i].size();j++)
            for(unsigned int cnt = 0; cnt<dim;cnt++)
                sumtmp[cnt] = sumtmp[cnt] + ((samples[i][cnt] - nsamples[j][cnt])/sqrt(distmat[j]))*(lr/distmat[j]);
        // Compute the displaced point
        for(unsigned int cnt = 0; cnt<dim;cnt++)
            sumtmp[cnt] = sumtmp[cnt] + samples[i][cnt];

        //finding the nearest sample to the new position which is computed based on coulomb law
        tmpidx = 0;
        find_nearest_neighbor(mat_index, data, sumtmp, dim, tmpidx);

        sumtmp.assign(dim,0);
        //condition check
        flag = false;
        for(unsigned int j=0;j<nsamples.size();j++){
            tmpdist2 = 0;
            for(unsigned int cnt = 0; cnt<dim;cnt++)
                tmpdist2 += std::pow(nsamples[j][cnt]-data[tmpidx][cnt],2);
            //Check the packing condition after moving the sample based on coulomb law. If it is true, it means that new sample is in distance less than radius to another sample which is not desirable
            if(tmpdist2 <= radius2){
                flag = true;
                break;
            }
        }
        if(!flag)
            samples[i] = data[tmpidx];
        std::vector<std::vector<double>>().swap(nsamples);
    }
}

py::array sampling(py::array arr, double radius, double lr, double tau){

    auto arr_obj_prop = arr.request();
    //initialize values
    double *vals = (double*) arr_obj_prop.ptr;

    unsigned int shape_1 = arr_obj_prop.shape[0];
    unsigned int shape_2 = arr_obj_prop.shape[1];


    std::vector<std::vector <double>> data( shape_1, std::vector<double> (shape_2));

    for(unsigned int i = 0; i < shape_1; i++){
      for(unsigned int j = 0; j < shape_2; j++){
        data[i][j] = vals[i*shape_2 + j];
      }
    }   

///////////////////////////////////////////////////////////////////////////
    std::srand ( unsigned ( std::time(0) ) );
    std::vector<std::vector<double>> datatmp;
    std::vector<std::vector<double>> samples,samples2;

    datatmp = data;
    // double radius = 0.15;

    struct parameters{
    double lr;
    double tau;
    double radius;
    double radius2;
    } params;

    params.lr = lr;
    params.tau = tau;
    params.radius = radius;
    params.radius2 = std::pow(radius,2);

    // Building the kdtree on the data once
    typedef KDTreeVectorOfVectorsAdaptor< std::vector<std::vector<double> >, double >  my_kd_tree_t;
    my_kd_tree_t   mat_index(dim  /*dim*/, data, 10 /* max leaf */ );
    mat_index.index->buildIndex();

    sample_by_covering_condition(data, samples, params.radius2, dim);


    // in this parts points repel each other so they will fill the whole space
    // and more space will be provided so we can sample again in this space. The
    // final thing will be more dense and the unwanted wholes will be removed

    int counter = 1;
    while(params.lr > 0.1*pow(params.radius,3)){
        params.lr = (pow(params.radius,3))*exp(-1*counter/params.tau);
        data = datatmp;
        std::random_shuffle ( samples.begin(), samples.end() );

        counter ++;
        coulomb_law_movements(mat_index, data,samples,params.radius,dim, params.lr);

        // check if all part of the data is covered by at least one sample
        check_covering_condition(mat_index, data, samples, params.radius2, dim);

        // If a part of "data" is not sampled, we sample it here.
        sample_by_covering_condition(data, samples, params.radius2, dim);
        std::cout<<counter<<"\n";
    }
    sample_by_covering_condition(samples, samples2, params.radius2, dim);

///////////////////////////////////////////////////////////////////////////////
    py::array ret =  py::cast(samples2);
    return ret;

}

PYBIND11_MODULE(asap, m) {

    m.doc() = "pybind11 module for iterating over generations";

    m.def("sampling", &sampling,
      "the function which loops over a numpy array", py::arg("data"), py::arg("radius") = 0.15, py::arg("lr") = 1, py::arg("tau") = 3.5);
}