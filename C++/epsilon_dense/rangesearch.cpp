#include <rangesearch.h>
#include <readdata.h>

void pdist2(std::vector<std::array<float,3>> &data, std::array<float,3> &point, std::vector<float> &distmat){
    // instead of computing sqrt() of distance, compute the 2nd power of radius once and compare it again and again which is faster
    distmat.resize(data.size());
    for (unsigned int i=0;i<data.size();i++){
        distmat[i] = std::pow(data[i][0]-point[0],2) + std::pow(data[i][1]-point[1],2) + std::pow(data[i][2]-point[2],2);
    }
}

void rangesearchfastV(std::vector<std::array<float,3>> &data, std::vector<std::vector<int>> &rmlist, float radius){
    radius = std::pow(radius,2);
    rmlist.resize(data.size());
    // instead of computing sqrt() of distance, compute the 2nd power of radius once and compare it again and again which is faster
    // If i is a neighbour of j then j is a neighbour of i too. Thus the first loop goes till size-1 and second loop start from i+1
    // Because of this setting of loops i==j never is true.
    for (unsigned int i=0;i<data.size()-1;i++){
        for (unsigned int j=i+1;j<data.size();j++){//I already remove point 'i' be neighbour with point 'i'
            if(std::pow(data[i][0]-data[j][0],2) + std::pow(data[i][1]-data[j][1],2) + std::pow(data[i][2]-data[j][2],2)<=radius){
                rmlist[i].push_back(j);
                rmlist[j].push_back(i);}
        }
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
