#include "../include/Conskiplist.hpp"
#include <iostream>
#include <thread>

#define MM 1000000
#define NUMBERDATA (64*MM)
#define PREINSERT 4*MM
#define SkiplistMaxLevel 9//(int)(log(NUMBERDATA)/log(2))
#define THREAD_NUMBER 16
#define NOFINDDEBUG 0
#define QUERY_ 0

int key_dis = (NUMBERDATA-PREINSERT)/THREAD_NUMBER;
unsigned int *dataq0 = new unsigned int[NUMBERDATA];
skiplist<unsigned int,int> *list = new skiplist<unsigned int,int>(SkiplistMaxLevel,Gm);

using namespace chrono;

void GetData(){
    //ycsb64M
    // char unique_dadta_file[] = "/home/yhzhou/datasets/ycsb64M.csv";//unique_iot_web_unique_shuffle
    // // //unique_iot_web_unique
    // ifstream fp(unique_dadta_file);
    // string line;
    // for(int i =0;i<NUMBERDATA;i++){
    //     getline(fp,line);
    //     dataq0[i] = atoi(line.c_str());
    //     // cerr<<"key["<<i<<"]"<<dataq0[i]<<endl;    
    // }
    // fp.close();
    for(int i = 0;i<NUMBERDATA;i++){
        dataq0[i] = i+1;
    }
    std::random_device rd;
    std::mt19937 g(rd());
    std::shuffle(dataq0, dataq0+NUMBERDATA, g);
}

void test(const int id,const int bound_l,const int bound_r ){
    for(int i = bound_l;i<bound_r;i++){
        if(dataq0[i] == 0 || dataq0[i] == UNINT_MAX)
            continue;
        // cout<<i<<endl;
        list->Add(dataq0[i],i);
    }
}

int main(){
    GetData();
    srand((int)time(0));
    std::cout<<"NUMBERDATA:"<<NUMBERDATA<<std::endl;
    std::vector<thread> threads;

    const auto start_time = std::chrono::steady_clock::now();
    for(int i = 0;i<PREINSERT;i++){
        if(dataq0[i] == 0 || dataq0[i] == UNINT_MAX)
            continue;
        list->Add(dataq0[i],i);
    }
    std::cout<<"PREINSERT finish"<<std::endl;
    int kk = PREINSERT;
    for(int idx=0; idx < THREAD_NUMBER; idx++){
        threads.push_back(thread(test, idx,kk,kk+key_dis));
        kk+=key_dis;
    } 
    for(int idx=0; idx <THREAD_NUMBER; idx++){
        threads[idx].join();
    }
    const auto end_time = std::chrono::steady_clock::now();
    const auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);
    std::cout << "insert time: " << duration.count() << "us"<<std::endl;
    int no_find = 0;
    for(int i = 0;i<NUMBERDATA;i++){
        std::pair<int,int> res = list->Lookup(dataq0[i]);
        if(! res.first || res.second != i){
            no_find++;
        }
    }
    std::cout<<"no_find:"<<no_find<<std::endl;
    // std::pair<int,int> res = list->Lookup(1);
    // std::cout<<res.first<<","<<res.second<<std::endl;
}