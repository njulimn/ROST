#pragma once
#include "../include/fileOperation.hpp"
#include <fstream>
#include <string>
#include <algorithm>
#include <time.h>

using namespace std;



int main(){
    srand((int)time(0));
    vector<snode> data;
    unsigned int MaxLevel = 0;
    MaxLevel = readFromCSV(data);
    // cerr<<"Max Level:"<<MaxLevel<<endl;
    skiplist* list = new skiplist(MaxLevel);
    list->setup(data);
    list->ShowNodeDis();

    //search test
    node* res = nullptr;
    res = list->Search(1487269949);
    if(res){
        cerr<<"find key:"<<res->value<<endl;
    }
    else{
        cerr<<"not find"<<endl;
    }
    return 0;
}