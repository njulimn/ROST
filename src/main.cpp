#pragma once
#include "../include/skiplist.hpp"
#include <fstream>
#include <string>
#include <algorithm>
#include <time.h>

#define NUMBERDATA 5000000
#define MaxL (int)(log(NUMBERDATA)/log(2))
#define Gma 128
#define Verify 0
#define NOFINDDEBUG 0

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
    std::ifstream inF("../../../datasets/iotdeleteZ.bin", std::ios::binary);
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
    cerr<<"search bound"<<bound<<endl;
    int nofindcnt = 0;
    start = clock();
    for (int i = 0; i < bound; i++) {
        node* res = list->Search(dataInput[i].key);
#if NOFINDDEBUG
        if(!res || res->key!=dataInput[i].key){
            nofindcnt++;
            snode x;
            x.key = dataInput[i].key;
            x.level = i;
            nofind.push_back(x);
        }
#endif
    }
    end = clock();
    sumTime =(double(end-start)/CLOCKS_PER_SEC);
#if NOFINDDEBUG
    cerr<<"no find cnt:"<<nofindcnt<<endl;
#endif
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
    for(int i = 0;i<1;i++){
        seg1 = seg1->forward[1];
        seg1->show(0);
        cerr<<"start index:"<<seg1->nodes[0].key<<endl;
        int bound = seg1->node_size;//min(4,(int)nofind.size());
        for (int j = 0; j < min(10,bound); j++) {
            unsigned int key = seg1->nodes[j].key;
            cerr<<"data: "<<key;
            double pred = seg1->slope*key;
            cerr<<"\tslope"<<seg1->slope<<"*key = "<<pred<<"\tintercept:"<<seg1->intercept;
            pred+=seg1->intercept;
            cerr<<"\tpred:"<<pred<<"\treal:"<<j<<endl;
            if(abs(pred-j)> 2*Gma){
                error++;
                cerr<<abs(pred-j)<<"\t";
            }
        }
    }
    cerr<<"error:"<<error<<endl;

#endif

#if NOFINDDEBUG
    cerr<<"NOFINDDEBUG:"<<endl;
    int error = 0;
    Segment_pt* seg1 = list->header;

    for(int i = 0;i<min(4,(int)nofind.size());i++){
        Segment_pt* x = list->header;
        // unsigned int pred;
        unsigned key = nofind[i].key;
        for(int i = list->level;i>=1;i--){
            int f = 0;
            while(key > x->forward[i]->stop){
                x = x->forward[i];
            }
            if(key >= x->forward[i]->start){
                cerr<<"key "<<key<<"locate in"<<endl;
                x = x->forward[i];
                x->show(0);
                int pred = x->slope*key+x->intercept;
                int l=max((pred-list->gamma),0);
                // node* data = Seg2Data[x];
                int r=min(x->node_size-1,pred+list->gamma),mid;
                
                while(l<=r){
                    mid = l+(r-l)/2;
                    if(x->nodes[mid].key == key){
                        cerr<<"find \n";
                        f = 1;
                        break;
                    }
                    else if(x->nodes[mid].key > key){
                        r = mid-1;
                    }
                    else{
                        l = mid+1;
                    }
                }

                if(f){
                    break;
                }
            }
        }    
    }
#endif
    return 0;
}