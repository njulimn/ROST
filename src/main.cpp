#pragma once
#include "../include/skiplist.hpp"
#include <fstream>
#include <string>
#include <algorithm>
#include <time.h>

#define NUMBERDATA 1000000
#define MaxL (int)(log(NUMBERDATA)/log(2))
#define Gma 128
#define Verify 0

using namespace std;

int length;
vector<snode> dataInput;
int MaxLevel = 0;
vector<snode> nofind;

unsigned int readFromCSV(vector<snode> &data)
{
    ifstream inFile("../static/exp_log2.csv", ios::in);
    if (!inFile)
    {
        cout << "打开文件失败！" << endl;
        exit(1);
    }
    srand((int)time(0));
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

        // getline(sin, field2, ','); //将字符串流sin中的字符读入到field字符串中，以逗号为分隔符 
        // cout << atoi(field.c_str()) << " ";//将刚刚读取的字符串转换成int
        // cout<<endl;
        // i++;
        snode x;
        x.key = atoi(field1.c_str());
        int NewLevel = 1;
        int P = 50;
        while(1){
            int t = rand() % 101;
            if(t<P)
                NewLevel++;
            else
                break;
        }
        NewLevel =  min(NewLevel,MaxL);
        x.level = NewLevel;
        data_x.push_back(x);
        res = max(res,x.level);
    }
    inFile.close();

    data.assign(data_x.begin(),data_x.end());
    return res;
    // cout << "read " << i << "lines" << endl;
    // cout << "读取数据完成" << endl;
}

void GetData()
{
    int length = 0;
    srand((int)time(0));
    unsigned int *dataIoT;
    std::ifstream inF("E:\\skiplistExp\\iotdeleteZ.bin", std::ios::binary);
    inF.seekg(0, inF.end);
    length = inF.tellg();
    inF.seekg(0, inF.beg);

    std::cerr << length << std::endl;
    dataIoT = (unsigned int*)malloc(sizeof(unsigned int) * length);

    inF.read(reinterpret_cast<char *>(dataIoT), sizeof(unsigned int) * length);
    inF.close();
    vector<unsigned int> dataq0(length);
    int j=0;
    for(int i = 0;i<length;i++){
        if(dataIoT[i]!=0){
            dataq0[j] = dataIoT[i];
            j++;
        }
    }
    dataq0.resize(j);
    length = j;

    int bound = min(NUMBERDATA,length);
    for(int i=0;i<bound;i++){
        snode x;
        x.key = dataq0[i];
        int NewLevel = 1;
        int P = 50;
        while(1){
            int t = rand() % 101;
            if(t<P)
                NewLevel++;
            else
                break;
        }
        NewLevel =  min(NewLevel,MaxL);
        x.level = NewLevel;
        dataInput.push_back(x);
        MaxLevel = max(MaxLevel,NewLevel);
    }
}

void ExpSearch(skiplist* list){
    clock_t start,end;
    double sumTime = 0;
    int bound = NUMBERDATA;
    cerr<<"bound"<<bound<<endl;
    int nofindcnt = 0;
    start = clock();
    for (int i = 0; i < bound; i++) {
        node* res = list->Search(dataInput[i].key);
    }
    end = clock();
    sumTime =(double(end-start)/CLOCKS_PER_SEC);
    // cerr<<"no find cnt:"<<nofindcnt<<endl;
    cerr<<"Search time = "<<sumTime<<"s"<<endl;  //输出时间（单位：ｓ）
}

int main(){
    // MaxLevel = readFromCSV(dataInput);
    GetData();
    length = dataInput.size();
    cerr<<"length:"<<length<<endl;
    // cerr<<"Max Level:"<<MaxLevel<<endl;
    skiplist* list = new skiplist(MaxL,Gma);
    list->setup(dataInput);
    list->ShowNodeDis();

    //show space size
    list->ComputeSpace();
    //search test
    ExpSearch(list);
#if Verify
    cerr<<"verify:"<<endl;
    int error = 0;
    Segment_pt* seg1 = list->header;
    for(int i = 0;i<list->size;i++){
        seg1 = seg1->forward[1];
        // cerr<<"start index:"<<seg1->fisrt_index<<endl;
    int bound = seg1->nodes.size();//min(4,(int)nofind.size());
    for (int j = 0; j < bound; j++) {
        unsigned int key = seg1->nodes[j].key;
        // cerr<<"data: "<<key;
        int pred = seg1->slope*key+seg1->intercept;
        // cerr<<"\tpred:"<<pred<<"\treal:"<<j<<endl;
        if(abs(pred-j)> 2*Gma){
            error++;
        }
        // Segment_pt* x = list->header;
        // bool locate = false;
        // unsigned int pred;
        // for(int i = list->level;i>=1;i--){
        //     while(nofind[j].key > x->forward[i]->stop){
        //         cerr<<"->\t";
        //         x = x->forward[i];
        //     }
        //     if(nofind[j].key >= x->forward[i]->start){
        //         x = x->forward[i];
        //         cerr<<"pred :"<<pred<<"\treal:"<<nofind[j].level<<endl;
        //         break;
        //     }
        //     cerr<<"down\t";
        // }
        // cerr<<endl;
    }
    cerr<<"error:"<<error<<endl;
    }
#endif
    return 0;
}