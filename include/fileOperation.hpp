#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <sstream>
#include <string>
#include "./skiplist.hpp"
using namespace std;


bool check_file_test(char const *fileName)
{
    // 用 ifstream 来判断文件是否存在
    ifstream testFile(fileName);
    if(!testFile)
    {
        cerr << "file not exit" << endl;
        return false;
    }
    else
    {
        cerr << "file exits" << endl;
        return true;
    }
    return false;
}

bool write_into_file(char const *fileName, char const *content)
{
    ofstream out(fileName,ios::app);
    // ofstream out;
    // out.open(fileName);
    if(!out.is_open())
    {
        cerr << "file not exit" << endl;
        return false;
    }
    else
    {
        out << content;
        // cerr << "write succeed" << endl;
        out.close();
        // usleep(100);
        return true;
    }
    
    return false;
}

bool generate_file_test(char const *fileName)
{
    ofstream out;
    out.open(fileName);
    // 判断文件是否已经打开
    if(out.is_open())
    {
        cerr << "file created succeed" << endl;
        out.close();
        return true;
    }
    else
    {
        cerr << "file created failed" << endl;
        out.close();
        return false;
    }
    
    out.close();
    return false;
}

unsigned int readFromCSV(vector<snode> &data)
{
    ifstream inFile("../static/exp_log2.csv", ios::in);
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