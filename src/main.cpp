#include "../include/Conskiplist.hpp"
#include <iostream>
#include <thread>

#define PRFO 0
#define PRFOINSERT 0
#define PRFOQUERY 0

#if PRFO
#include<gperftools/profiler.h>
#endif

#define MM 1000000
#define NUMBERDATA (64*MM)
#define SkiplistMaxLevel (int)(log(1500)/log(2))
#define THREAD_NUMBER 32
#define NOFINDDEBUG 0
#define QUERY_TEST 1
#define TEST_4 0

#if TEST_4
int key_dis = (4*MM)/THREAD_NUMBER;
#else
int key_dis = (NUMBERDATA)/THREAD_NUMBER;
#endif
unsigned int *input_data = new unsigned int[NUMBERDATA];
skiplist<unsigned int,int> *list = new skiplist<unsigned int,int>(SkiplistMaxLevel,Gm);

long long partial_r[THREAD_NUMBER];
long long MAX_DEPTH_array[THREAD_NUMBER];
long long scan_time_on_check[THREAD_NUMBER];
int split_total[THREAD_NUMBER];
long long collision_total[THREAD_NUMBER];
long long insert_time_total[THREAD_NUMBER];
char search_time[] = "./search_time.csv";
char insert_time[] = "./insert_time.csv";

using namespace chrono;

void GetData(){
    //ycsb64M
    char unique_dadta_file[] = "/home/yhzhou/datasets/ycsb64M.csv";
    ifstream fp(unique_dadta_file);
    string line;
    for(int i =0;i<NUMBERDATA;i++){
        getline(fp,line);
        input_data[i] = atoi(line.c_str());  
    }
    fp.close();
    // std::random_device rd;
    // std::mt19937 g(rd());
    // std::shuffle(input_data, input_data+NUMBERDATA, g);
}

void test_insert(const int id,const int bound_l,const int bound_r ,int dijiduan){
    long long rebuild_cnt = 0;
    int MaxDepth = 1;
    int SplitCnt = 0;
    long long scan_on_rebuild_topk = 0;
    long long collision_ = 0;
    for(int i = bound_l;i<bound_r;i++){
        if(input_data[i] == 0 || input_data[i] == UNINT_MAX)
            continue;
        // list->Add(input_data[i],i,cnt);
        int cnt = 0;
        long long scan = 0;
        long long collision_cnt = 0;
        if(list->Add(input_data[i],i,cnt,scan,SplitCnt,collision_cnt)){
            rebuild_cnt++;
        }
        scan_on_rebuild_topk+=scan;
        collision_+=collision_cnt;
    }
    partial_r[id] = rebuild_cnt;
    scan_time_on_check[id] = scan_on_rebuild_topk;
    split_total[id] = SplitCnt;
    collision_total[id] = collision_;
}

void test_query(const int id,const int bound_l,const int bound_r ){
    std::pair<int,int> res;
    int nofind = 0;
    for(int i = bound_l;i<bound_r;i++){
        if(input_data[i] == 0 || input_data[i] == UNINT_MAX)
            continue;
        // list->Lookup(input_data[i]);
        res = list->Lookup(input_data[i]);
        if(!(res.first) || res.second != i){
            nofind++;
        }
    }
    if(nofind)
        std::cout<<"nofind:\t"<<nofind<<std::endl;
}

void Insert_Part(int &kk,int id){
    std::vector<thread> threads;
    #if PRFOINSERT
    string name = "test"+std::to_string(id)+".prof";
    ProfilerStart(name.c_str());
    #endif
    for(int i = 0;i<THREAD_NUMBER;i++){
        partial_r[i] = 0;
    }
    const auto start_time = std::chrono::steady_clock::now();
    for(int idx=0; idx < THREAD_NUMBER; idx++){
        threads.push_back(thread(test_insert, idx,kk,kk+key_dis,id));
        kk+=key_dis;
    } 
    for(int idx=0; idx <THREAD_NUMBER; idx++){
        threads[idx].join();
    }
    #if PRFOINSERT
    ProfilerStop();
    #endif
    const auto end_time = std::chrono::steady_clock::now();
    const auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);
    long long cnt = 0;
    long long split_cnt = 0;
    long long time_scan = 0;
    long long collision_cnt = 0;
    for(int i = 0;i<THREAD_NUMBER;i++){
        cnt+=partial_r[i];
        time_scan+=scan_time_on_check[i];
        // std::cout<<i<<" time on split/scan:"<<time_scan<<std::endl;
        split_cnt+=split_total[i];
        collision_cnt+=collision_total[i];
    }
    insert_time_total[id] = duration.count();
    string outs = "insert time: "+ std::to_string( duration.count()) + "us\trebuild cnt: " + std::to_string(cnt)+ 
        "\ttime on split/rebuild: " + std::to_string(time_scan) +"\tsplit cnt:"+std::to_string(split_cnt)+
        "\tcollision cnt:"+std::to_string(collision_cnt)+"\t";
    std::cout<<outs;
    list->ShowSegmentNumber(false);
    string k = std::to_string(duration.count())+"\n";
    write_into_file(insert_time,k.c_str());
}

void Query_Part(int &st){
    std::vector<thread> threads_2;
    const auto start_time = std::chrono::steady_clock::now();
    int query_key_dis = NUMBERDATA/THREAD_NUMBER;
    #if PRFOQUERY
    string name = "query.prof";
    ProfilerStart(name.c_str());
    #endif
    for(int idx=0; idx < THREAD_NUMBER; idx++){
        threads_2.push_back(thread(test_query, idx,st,st+query_key_dis));
        st+=query_key_dis;
    } 
    for(int idx=0; idx <THREAD_NUMBER; idx++){
        threads_2[idx].join();
    }
    #if PRFOQUERY
    ProfilerStop();
    #endif
    const auto end_time = std::chrono::steady_clock::now();
    const auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);
    std::cout << "query time: " << duration.count() << "us"<<std::endl;
    string k = std::to_string(duration.count())+"\n";
    write_into_file(search_time,k.c_str());
}

void RangeQuery_Test(){
    vector<std::pair<unsigned int,int>> rg_res;
    unsigned int st = 500,ed = 1480000;
    list->Lookup(st,ed,rg_res);
    for(auto it:rg_res){
        std::cout<<it.first<<","<<it.second<<std::endl;
    }
    
}

int main(){
    const auto start_time = std::chrono::steady_clock::now();
    GetData();
    srand((int)time(0));
    std::cout<<"NUMBERDATA:"<<NUMBERDATA<<"\tTHREAD_NUMBER:"<<THREAD_NUMBER<<std::endl;
    std::cout<<"SkiplistMaxLevel:"<<SkiplistMaxLevel<<"\tmax segment size:"<<list->segment_max_size<<"\tDELTA_INSERT:"<<DELTA_INSERT<<std::endl;
    const auto end_time = std::chrono::steady_clock::now();
    const auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);
    std::cout << "prepare time: " << duration.count() << "us"<<std::endl;
    
    int n = 0;
    
    #if TEST_4
    for(int i = 0;i<16;i++){
        Insert_Part(n,i);
    }
    long long time_sum = 0;
    for(int i = 0;i<THREAD_NUMBER;i++){
        time_sum+=insert_time_total[i];
    }
    std::cout<<"sum of time:"<<time_sum<<std::endl;
    #else
    Insert_Part(n,0);
    #endif
    #if QUERY_TEST
    n = 0;
    Query_Part(n);
    #endif   

    std::cout<<"LearnedSkiplist space size:"<<list->SpaceSize() - NUMBERDATA * (sizeof(unsigned int)+ sizeof(int)) <<"\tByte"<<std::endl;
}