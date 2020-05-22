#include <iostream>
#include <array>
#include <vector>
#include <readdata.h>
#include <rangesearch.h>
#include <cmath>
#include <random>
#include <algorithm>
#include <ctime>

using namespace std;

int main()
{

    std::srand ( unsigned ( std::time(0) ) );
    std::vector<std::array<float,3>> data,datatmp;
    std::vector<std::array<float,3>> samples;
    std::string nameofinp;
    std::cout<<"Please insert the name of data file. Note that it should be a 'csv' file and a n*3 table.\n";
    std::cin>>nameofinp;
    data = readdata(nameofinp);
    datatmp = data;
    std::cout<<"Please insert the radius value.\n";
    float radius = 0.2;
    std::cin>>radius;
    float radius2 = std::pow(radius,2);
    std::array<float,3> point;
    std::vector<std::vector<int>> idx;
    int i = 0;
    while(data.size()>0){
        point = data[rand()%data.size()];
        samples.push_back(point);
        data.erase(std::remove_if(data.begin(),data.end(),[&](std::array<float,3> j){return std::pow(j[0]-point[0],2) + std::pow(j[1]-point[1],2) + std::pow(j[2]-point[2],2) <= radius2;}),data.end());
        if(i%100 == 0){
        data.shrink_to_fit();
        std::cout<<i/100<<" "<<data.size()<<"\n";
        }
        i++;
    }
    writeCSV2d("dataout.csv",samples);

    //in this parts points repel each other so they will fill the whole space
    //and more space will be provided so we can sample again in this space. The
    //final thing will be more dense and the unwanted wholes will be removed
    float lr = 1;
    float tau = 3.5;
    int counter = 1;
    std::vector<std::array<float,3>> nsamples;
    std::vector<float> distmat;
    std::array<float,3> sumtmp = {0,0,0};
    float tmpdist1 = 0;
    float tmpdist2 = 0;
    int tmpidx = 0;
    bool flag = false;
    while(lr > 0.1*pow(radius,3)){
        lr = (pow(radius,3))*exp(-1*counter/tau);
        data = datatmp;
        std::random_shuffle ( samples.begin(), samples.end() );

        counter ++;
        rangesearchfastV(samples, idx, 2*radius);
        for (unsigned int i=0; i<samples.size(); i++){
            if(idx[i].size() == 0)
                continue;
            nsamples.resize(idx[i].size(),sumtmp);
            for(unsigned int j = 0; j<idx[i].size();j++)
                nsamples[j] = samples[idx[i][j]];
            pdist2(nsamples, samples[i], distmat);
            for(unsigned int j = 0; j<idx[i].size();j++){
                sumtmp[0] = sumtmp[0] + ((samples[i][0] - nsamples[j][0])/sqrt(distmat[j]))*(lr/distmat[j]);
                sumtmp[1] = sumtmp[1] + ((samples[i][1] - nsamples[j][1])/sqrt(distmat[j]))*(lr/distmat[j]);
                sumtmp[2] = sumtmp[2] + ((samples[i][2] - nsamples[j][2])/sqrt(distmat[j]))*(lr/distmat[j]);
            }
            sumtmp[0] = sumtmp[0] + samples[i][0];
            sumtmp[1] = sumtmp[1] + samples[i][1];
            sumtmp[2] = sumtmp[2] + samples[i][2];

            tmpdist1 = 1000;
            tmpidx = 0;
            for(unsigned int j = 0; j<data.size();j++){
            tmpdist2 = std::pow(data[j][0]-sumtmp[0],2) + std::pow(data[j][1]-sumtmp[1],2) + std::pow(data[j][2]-sumtmp[2],2);
                if(tmpdist1 > tmpdist2){
                    tmpdist1 = tmpdist2;
                    tmpidx = j;
                }
            }
            sumtmp = {0,0,0};
            //condition check
            flag = false;
            for(unsigned int j=0;j<nsamples.size();j++){
                if(std::pow(nsamples[j][0]-data[tmpidx][0],2) + std::pow(nsamples[j][1]-data[tmpidx][1],2) + std::pow(nsamples[j][2]-data[tmpidx][2],2) <= radius2){
                    flag = true;
                    break;
                }
            }
            if(!flag)
                samples[i] = data[tmpidx];
            std::vector<std::array<float,3>>().swap(nsamples);
        }
    for(unsigned int i=0; i<samples.size();i++){
        point = samples[i];
        data.erase(std::remove_if(data.begin(),data.end(),[&](std::array<float,3> j){return std::pow(j[0]-point[0],2) + std::pow(j[1]-point[1],2) + std::pow(j[2]-point[2],2) <= radius2;}),data.end());
        if(i%100 == 0){
        data.shrink_to_fit();
        }
    }
    while(data.size()>0){
        point = data[rand()%data.size()];
        samples.push_back(point);
        data.erase(std::remove_if(data.begin(),data.end(),[&](std::array<float,3> j){return std::pow(j[0]-point[0],2) + std::pow(j[1]-point[1],2) + std::pow(j[2]-point[2],2) <= radius2;}),data.end());
        if(i%100 == 0){
        data.shrink_to_fit();
        }
        i++;
    }

    }
    std::vector<std::array<float,3>> samples2;
    while(samples.size()>0){
        point = samples[rand()%samples.size()];
        samples2.push_back(point);
        samples.erase(std::remove_if(samples.begin(),samples.end(),[&](std::array<float,3> j){return std::pow(j[0]-point[0],2) + std::pow(j[1]-point[1],2) + std::pow(j[2]-point[2],2) <= radius2;}),samples.end());
        if(i%100 == 0){
        samples.shrink_to_fit();
        }
        i++;
    }
    writeCSV2d("dataout_columb.csv",samples2);
    return 0;
}
