#include "../include/Conskiplist.hpp"
#include <iostream>
#include <thread>

#define PRFO 0

#if PRFO
#include<gperftools/profiler.h>
#endif

#define MM 1000000
#define NUMBERDATA (64*MM)
#define SkiplistMaxLevel (int)(log(1500)/log(2))
#define THREAD_NUMBER 32
#define NOFINDDEBUG 0
#define QUERY_TEST 1

int key_dis = (NUMBERDATA)/THREAD_NUMBER;
unsigned int *dataq0 = new unsigned int[NUMBERDATA];
skiplist<unsigned int,int> *list = new skiplist<unsigned int,int>(SkiplistMaxLevel,Gm);

// int max_depths[THREAD_NUMBER];
// long long path_depths[THREAD_NUMBER];
// long long scan_cnts[THREAD_NUMBER];
// long long partial_r[THREAD_NUMBER];

using namespace chrono;

void GetData(){
    //ycsb64M
    char unique_dadta_file[] = "/home/yhzhou/datasets/ycsb64M.csv";//unique_iot_web_unique_shuffle
    // //unique_iot_web_unique
    ifstream fp(unique_dadta_file);
    string line;
    for(int i =0;i<NUMBERDATA;i++){
        getline(fp,line);
        dataq0[i] = atoi(line.c_str());   
    }
    fp.close();
    // for(int i = 0;i<NUMBERDATA;i++){
    //     dataq0[i] = i+1;
    // }
    // std::random_device rd;
    // std::mt19937 g(rd());
    // std::shuffle(dataq0, dataq0+NUMBERDATA, g);
}

void test(const int id,const int bound_l,const int bound_r ){
    // long long caslock_acquire = 0,bottom_move = 0,split_blk = 0; 
    // std::pair<int,int> res;
    // long long scan_cnt = 0;
    // long long depth = 0;
    // int max_depth_ = 0;
    int cnt = 0;
    for(int i = bound_l;i<bound_r;i++){
        if(dataq0[i] == 0 || dataq0[i] == UNINT_MAX)
            continue;
        list->Add(dataq0[i],i,cnt);
        // int cnt = 0;
        // if(list->Add(dataq0[i],i,cnt)){
        //     scan_cnt++;
        // }
        // depth+=cnt;
        // max_depth_ = max(max_depth_,cnt);
    }
    // path_depths[id] = depth;
    // scan_cnts[id] = scan_cnt;
    // max_depths[id] = max_depth_;
}

void test_query(const int id,const int bound_l,const int bound_r ){
    std::pair<int,int> res;
    int nofind = 0;
    for(int i = bound_l;i<bound_r;i++){
        if(dataq0[i] == 0 || dataq0[i] == UNINT_MAX)
            continue;
        // list->Lookup(dataq0[i]);
        res = list->Lookup(dataq0[i]);
        if(!(res.first) || res.second != i){
            nofind++;
        }
    }
    if(nofind)
        std::cout<<"nofind:\t"<<nofind<<std::endl;
}

void Insert_Half(int &kk){
    std::vector<thread> threads;
    #if PRFO
    ProfilerStart("test.prof");
    #endif
    const auto start_time = std::chrono::steady_clock::now();
    for(int idx=0; idx < THREAD_NUMBER; idx++){
        threads.push_back(thread(test, idx,kk,kk+key_dis));
        kk+=key_dis;
    } 
    for(int idx=0; idx <THREAD_NUMBER; idx++){
        threads[idx].join();
    }
    #if PRFO
    ProfilerStop();
    #endif
    const auto end_time = std::chrono::steady_clock::now();
    const auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);
    std::cout << "insert time: " << duration.count() << "us"<<std::endl;
}

void Query_Half(int st){
    std::vector<thread> threads_2;
    const auto start_time = std::chrono::steady_clock::now();
    for(int idx=0; idx < THREAD_NUMBER; idx++){
        threads_2.push_back(thread(test_query, idx,st,st+key_dis));
        st+=key_dis;
    } 
    for(int idx=0; idx <THREAD_NUMBER; idx++){
        threads_2[idx].join();
    }
    const auto end_time = std::chrono::steady_clock::now();
    const auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);
    std::cout << "query time: " << duration.count() << "us"<<std::endl;
}

int main(){
    GetData();
    srand((int)time(0));
    std::cout<<"NUMBERDATA:"<<NUMBERDATA<<std::endl;
    std::cout<<"THREAD_NUMBER:"<<THREAD_NUMBER<<std::endl;
    std::cout<<"SkiplistMaxLevel:"<<SkiplistMaxLevel<<std::endl;
    std::cout<<"list max segment size:"<<list->segment_max_size<<std::endl;
    std::cout<<"DETA_INSERT:"<<DETA_INSERT<<std::endl;
    int kk = 0;
    Insert_Half(kk);
    Query_Half(0);   

    // long long scan = 0;
    // long long depth_sum = 0;
    // long long partial_ = 0;
    // int max_depth_ = 0;
    // for(int i = 0;i<THREAD_NUMBER;i++){
    //     scan+=scan_cnts[i];
    //     depth_sum+=path_depths[i];
    //     max_depth_ = max(max_depth_,max_depths[i]);
    //     partial_+=partial_r[i];
    // }
    // std::cout<<"scan rebuild/split cnt: "<<scan<<std::endl;
    // std::cout<<"depth of insert: "<<depth_sum<<"\t average: "<<depth_sum*1.0/NUMBERDATA<<"\tmax:"<<max_depth_<<std::endl;
    // std::cout<<"partial_rebuild cnt: "<<partial_<<std::endl;
    // list->ShowSegmentNumber();
    #if 0
    std::vector<thread> threads_2;
    kk = 0;
    for(int idx=0; idx < THREAD_NUMBER; idx++){
        threads_2.push_back(thread(test_query, idx,kk,kk+key_dis));
        kk+=key_dis;
    } 
    for(int idx=0; idx <THREAD_NUMBER; idx++){
        threads_2[idx].join();
    }
    #endif
}