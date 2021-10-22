#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <sstream>
#include <string>
using namespace std;

int main()
{
    ifstream inFile("../static/exp_link.csv", ios::in);
    if (!inFile)
    {
        cout << "打开文件失败！" << endl;
        exit(1);
    }
    vector<unsigned int> data_x;
    vector<unsigned int> data_l;
    int i = 0;
    string line;
    string field;
    while (getline(inFile, line))//getline(inFile, line)表示按行读取CSV文件中的数据
    {
        string field;
        istringstream sin(line); //将整行字符串line读入到字符串流sin中

        getline(sin, field, ','); //将字符串流sin中的字符读入到field字符串中，以逗号为分隔符 
        data_x.push_back(atoi(field.c_str()));
        cout<<atoi(field.c_str())<<" ";//将刚刚读取的字符串转换成int

        getline(sin, field, ','); //将字符串流sin中的字符读入到field字符串中，以逗号为分隔符 
        data_l.push_back(atoi(field.c_str()));
        cout << atoi(field.c_str()) << " ";//将刚刚读取的字符串转换成int
        cout<<endl;
        i++;
    }
    inFile.close();
    cout << "read " << i << "lines" << endl;
    // cout << "读取数据完成" << endl;
    return 0;
}