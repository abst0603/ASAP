#include <iostream>
#include <array>
#include <vector>
#include <rangesearch.h>
#include <cmath>
#include <random>
#include <algorithm>
#include <ctime>
using namespace std;

void sample_by_covering_condition(std::vector<std::vector<float>> &data,std::vector<std::vector<float>> &samples, float radius2, int dim){
    std::vector<float> point(dim);
    while(data.size()>0){
        point = data[rand()%data.size()];
        samples.push_back(point);
        data.erase(std::remove_if(data.begin(),data.end(),
        //lambda expression for removing neighbors. This one works for Ndims data
        [&](std::vector<float> data_point){
        float temp_dist = 0;
        for(unsigned int cnt=0; cnt<data_point.size();cnt++){temp_dist += std::pow(data_point[cnt]-point[cnt],2);}
        return temp_dist <= radius2;}
        ),data.end());
    }
    data.shrink_to_fit();
}

void check_covering_condition(std::vector<std::vector<float>> &data,std::vector<std::vector<float>> &samples, float radius2, int dim){
    std::vector<float> point(dim);
    for(unsigned int i=0; i<samples.size();i++){
        point = samples[i];
        data.erase(std::remove_if(data.begin(),data.end(),
        //lambda expression for removing neighbors. This one works for Ndims data
        [&](std::vector<float> data_point){
        float temp_dist = 0;
        for(unsigned int cnt=0; cnt<data_point.size();cnt++){temp_dist += std::pow(data_point[cnt]-point[cnt],2);}
        return temp_dist <= radius2;}
        ),data.end());
    }
    data.shrink_to_fit();
}

void coulomb_law_movements(std::vector<std::vector<float>> &data,std::vector<std::vector<float>> &samples, float radius, int dim, float lr){
    std::vector<std::vector<unsigned int>> idx;
    std::vector<std::vector<float>> nsamples(dim);
    std::vector<float> distmat;
    std::vector<float> sumtmp(dim,0);
    float tmpdist1 = 0;
    float tmpdist2 = 0;
    int tmpidx = 0;
    bool flag = false;
    float radius2 = std::pow(radius,2);
    RangeSearch(idx, samples, 2*radius);
    for (unsigned int i=0; i<samples.size(); i++){
        if(idx[i].size() == 0)
            continue;
        nsamples.resize(idx[i].size(),sumtmp);
        for(unsigned int j = 0; j<idx[i].size();j++)
            nsamples[j] = samples[idx[i][j]];
        pdist2(nsamples, samples[i], distmat);
        for(unsigned int j = 0; j<idx[i].size();j++)
            for(unsigned int cnt = 0; cnt<dim;cnt++)
                sumtmp[cnt] = sumtmp[cnt] + ((samples[i][cnt] - nsamples[j][cnt])/sqrt(distmat[j]))*(lr/distmat[j]);

        for(unsigned int cnt = 0; cnt<dim;cnt++)
            sumtmp[cnt] = sumtmp[cnt] + samples[i][cnt];

        //finding the nearest sample to the new position which is computed based on coulomb law
        tmpdist1 = 1000; //a big number
        tmpidx = 0;
        for(unsigned int j = 0; j<data.size();j++){
            tmpdist2 = 0;
            for(unsigned int cnt = 0; cnt<dim;cnt++)
                tmpdist2 += std::pow(data[j][cnt]-sumtmp[cnt],2);
            if(tmpdist1 > tmpdist2){
                tmpdist1 = tmpdist2;
                tmpidx = j;
            }
        }
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
        std::vector<std::vector<float>>().swap(nsamples);
    }
}

int main()
{
	const int dim = 5; //You should set the value before compiling since nanoflann use this value for unrolling the loops.
    std::srand ( unsigned ( std::time(0) ) );
    std::vector<std::vector<float>> data,datatmp;
    std::vector<std::vector<float>> samples;
    std::string nameofinp;
    std::cout<<"Please insert the name of data file. Note that it should be a 'csv' file and a n*D table.\n You should insert the D as dim in the code.\n";
    std::cin>>nameofinp;
    readdata(dim, data, nameofinp);
    datatmp = data;
    std::cout<<"Please insert the radius value.\n";
    float radius = 0.2;
    std::cin>>radius;

    struct parameters{
    float lr = 1;
    float tau = 3.5;
    float radius;
    float radius2;
    } params;

    params.radius = radius;
    params.radius2 = std::pow(radius,2);
    sample_by_covering_condition(data, samples, params.radius2, dim);
    writeCSV2d("dataout.csv",samples);

    //in this parts points repel each other so they will fill the whole space
    //and more space will be provided so we can sample again in this space. The
    //final thing will be more dense and the unwanted wholes will be removed

    int counter = 1;
    while(params.lr > 0.1*pow(params.radius,3)){
        params.lr = (pow(params.radius,3))*exp(-1*counter/params.tau);
        data = datatmp;
        std::random_shuffle ( samples.begin(), samples.end() );

        counter ++;
        coulomb_law_movements(data,samples,params.radius,dim, params.lr);

    // check if all part of the data is covered by at least one sample
    check_covering_condition(data, samples, params.radius2, dim);

    // If a part of "data" is not sampled, we sample it here.
    sample_by_covering_condition(data, samples, params.radius2, dim);
    }
    std::vector<std::vector<float>> samples2;
    sample_by_covering_condition(samples, samples2, params.radius2, dim);

    writeCSV2d("dataout_columb.csv",samples2);
    return 0;
}

