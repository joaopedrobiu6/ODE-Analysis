#include "DataReader.h"

DataReader::DataReader(std::string filename)
{
    std::fstream rFile(filename); // read mode
    if (rFile.fail())
    {
        std::cout << "Error: Unable to load data file" << std::endl;
        exit(1);
    }

    std::string line;
    while (getline(rFile, line))
    {                               // loop on file lines
        std::stringstream ss(line); // build object stringstream
        float d;
        std::vector<float> temp;
        while (ss >> d)
        { // parse line words to numbers (empty space separated)
            temp.push_back(d);
        }
        ss.clear(); // erase stringstream contents
        data.push_back(temp);
    }
    rFile.close();
};
DataReader::DataReader(const DataReader &D)
{
    data = D.data;
}

void DataReader::dump()
{
    for (int i = 0; i < data.size(); i++)
    {
        for (int j = 0; j < data[0].size(); j++)
        {
            std::cout << data[i][j] << ", ";
        }
        std::cout << "\n";
    }
};

const std::vector<std::vector<float>> &DataReader::GetData()
{
    return data;
}

std::ostream &operator<<(std::ostream &s, const DataReader &DR)
{
    for (int i = 0; i < DR.data.size(); i++)
    {
        for (int j = 0; j < DR.data[0].size(); j++)
        {
            s << DR.data[i][j] << "\t";
        }
        std::cout << "\n";
    }
    return s;
}