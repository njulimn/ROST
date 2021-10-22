#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <sstream>
#include <string>
#include "./skiplist2.hpp"
using namespace std;

unsigned int readFromCSV(vector<snode> &data)
{
    ifstream inFile("../static/exp_link50.csv", ios::in);
    if (!inFile)
    {
        cout << "打开文件失败！" << endl;
        exit(1);
    }
    vector<snode> data_x;
    // int i = 0;
    string line;
    string field;
    unsigned int res=0;
    while (getline(inFile, line))//getline(inFile, line)表示按行读取CSV文件中的数据
    {
        string field1,field2;
        istringstream sin(line); //将整行字符串line读入到字符串流sin中

        getline(sin, field1, ','); //将字符串流sin中的字符读入到field字符串中，以逗号为分隔符 
        // cout<<atoi(field.c_str())<<" ";//将刚刚读取的字符串转换成int

        getline(sin, field2, ','); //将字符串流sin中的字符读入到field字符串中，以逗号为分隔符 
        // cout << atoi(field.c_str()) << " ";//将刚刚读取的字符串转换成int
        // cout<<endl;
        // i++;
        snode x;
        x.key = atoi(field1.c_str());
        x.level = atoi(field2.c_str());
        data_x.push_back(x);
        res = max(res,x.level);
    }
    inFile.close();

    data.assign(data_x.begin(),data_x.end());
    return res;
    // cout << "read " << i << "lines" << endl;
    // cout << "读取数据完成" << endl;
}