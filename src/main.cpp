#pragma once
#include "../include/fileOperation.hpp"
#include <fstream>
#include <string>
#include <algorithm>
#include <time.h>

using namespace std;

int length;
vector<snode> dataInput;

void ExpSearch(skiplist* list){
    clock_t start,end;
    double sumTime = 0;
    int bound = min(50000,length);
    cerr<<"bound"<<bound<<endl;
    for (int i = 0; i < bound; i++) {
        start = clock();
        list->Search(dataInput[i].key);
        end = clock();
        sumTime +=(double(end-start)/CLOCKS_PER_SEC);
    }
    cerr<<"Search time = "<<sumTime<<"s"<<endl;  //输出时间（单位：ｓ）
}

int main(){
    srand((int)time(0));
    
    unsigned int MaxLevel = 0;
    MaxLevel = readFromCSV(dataInput);
    length = dataInput.size();
    cerr<<"length:"<<length<<endl;
    // cerr<<"Max Level:"<<MaxLevel<<endl;
    skiplist* list = new skiplist(MaxLevel);
    list->setup(dataInput);
    list->ShowNodeDis();

    //show space size
    list->ComputeSpace();
    //search test
    ExpSearch(list);
    return 0;
}