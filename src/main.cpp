#include "../include/Conskiplist.hpp"
#include <iostream>
#include <thread>
#include <set>

#define PRFO 0
#define PRFOINSERT 0
#define PRFOQUERY 0

#if PRFO
#include<gperftools/profiler.h>
#endif

#define MM 1000000
#define NUMBERDATA (64*MM)
#define THREAD_NUMBER 32
#define NOFINDDEBUG 0
#define QUERY_TEST 1
#define WriteRatio (1.0)
#define FOLDKEY 0
#define SPACEPRINT 0
#define SegmentGamma 16

int read_offset = NUMBERDATA * WriteRatio;
int write_dis = (read_offset)/THREAD_NUMBER;
int read_dis = (NUMBERDATA * (1 - WriteRatio)) / THREAD_NUMBER;
char data_file[] = "/root/LSDataset/twitter/twitter64M01.txt";

KeyType *input_data = new KeyType[NUMBERDATA];
skiplist<KeyType,VaueType,ModelType> *list = new skiplist<KeyType,VaueType,ModelType>(SkiplistMaxLevel,SegmentGamma);

long long partial_r[THREAD_NUMBER];
long long MAX_DEPTH_array[THREAD_NUMBER];
long long scan_time_on_check[THREAD_NUMBER];
int split_total[THREAD_NUMBER];
long long collision_total[THREAD_NUMBER];
long long insert_time_total[THREAD_NUMBER];
char search_time[] = "./search_time.csv";
char insert_time[] = "./insert_time.csv";
char workload_time[] = "./workload01_time.csv";

using namespace chrono;

void GetDataFromText(char path[]){
    // set<KeyType> unique;
    ifstream myfile(path);
    if (!myfile.is_open())
	{
		cout << "can not open this file" << endl;
		exit(-1);
	}
    for(int i = 0;i<NUMBERDATA;i++){
        myfile >> input_data[i];
        int k;
        myfile >> k;
    }
    myfile.close();
    // std::random_device rd;
    // std::mt19937 g(rd());
    // std::shuffle(input_data, input_data+NUMBERDATA, g);
}

void test_query(const int id,const int bound_l,const int bound_r ){
    std::pair<int,VaueType> res;
    int nofind = 0;
    for(int i = bound_l;i<bound_r;i++){
        if(input_data[i] == 0 || input_data[i] == KeyMax)
            continue;
        // list->Lookup(input_data[i]);
        res = list->Lookup(id,input_data[i]);
        if(!(res.first)){
            nofind++;
        }
    }
    if(nofind)
        std::cout<<"nofind:\t"<<nofind<<std::endl;
}

void WorkloadTest(int id){
    int bound_l = id * (write_dis);
    int bound_r = (id+1)*(write_dis);
    // printf("range[%d,%d)\n",bound_l,bound_r);
    int SplitCnt = 0;
    skiplist<KeyType,VaueType,ModelType>::subtree** route = (skiplist<KeyType,VaueType,ModelType>::subtree**)malloc(sizeof(skiplist<KeyType,VaueType,ModelType>::subtree*)*MAX_DEPTH*2);
    for(int i = bound_l;i<bound_r;i++){
        if(input_data[i] == 0 || input_data[i] == KeyMax)
            continue;
        int cnt = 0;
        long long collision_cnt = 0;
        // if((i - bound_l) % MM == 0)
            // cout<<i<<std::endl;
        list->Add(id,input_data[i],i,route,cnt,SplitCnt,collision_cnt);
    }
    // std::cout<<"insert finish"<<std::endl;
    bound_l = read_offset + id*(read_dis);
    bound_r = read_offset + (id+1)*(read_dis);
    int no_find = 0;
    for(int i = bound_l;i<bound_r;i++){
        if(input_data[i] == 0 || input_data[i] == KeyMax)
            continue;
        auto res = list->Lookup(id,input_data[i]);
        if(!res.first){
            no_find++;
        }
    }
    free(route);
    // std::cout<<"no find:"<<no_find<<std::endl;
}

void WorkloadInserPart(int id){
    int bound_l = id * (write_dis);
    int bound_r = (id+1)*(write_dis);
    // printf("range[%d,%d)\n",bound_l,bound_r);
    int SplitCnt = 0;
    skiplist<KeyType,VaueType,ModelType>::subtree** route = (skiplist<KeyType,VaueType,ModelType>::subtree**)malloc(sizeof(skiplist<KeyType,VaueType,ModelType>::subtree*)*MAX_DEPTH*2);
    for(int i = bound_l;i<bound_r;i++){
        if(input_data[i] == 0 || input_data[i] == KeyMax)
            continue;
        int cnt = 0;
        long long collision_cnt = 0;
        list->Add(id,input_data[i],i,route,cnt,SplitCnt,collision_cnt);
    }
    free(route);
}

void WorkloadQueryPart(int id){
    int bound_l = read_offset + id*(read_dis);
    int bound_r = read_offset + (id+1)*(read_dis);
    for(int i = bound_l;i<bound_r;i++){
        if(input_data[i] == 0 || input_data[i] == KeyMax)
            continue;
        auto res = list->query(id,input_data[i]);
    }
}

void Insert_Part(int &kk,int id){
    #if PRFOINSERT
    string name = "workload73.prof";
    ProfilerStart(name.c_str());
    #endif
    std::vector<thread> threads;
    for(int i = 0;i<THREAD_NUMBER;i++){
        partial_r[i] = 0;
    }
    const auto start_time = std::chrono::steady_clock::now();
    // std::cout<<"start time:\t"<<start_time.time_since_epoch().count()<<std::endl;
    for(int idx=0; idx < THREAD_NUMBER; idx++){
        threads.push_back(thread(WorkloadTest,idx));
    } 
    for(int idx=0; idx <THREAD_NUMBER; idx++){
        threads[idx].join();
    }
    #if PRFOINSERT
    ProfilerStop();
    #endif
    const auto end_time = std::chrono::steady_clock::now();
    const auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);
    string outs = "operation time: "+ std::to_string( duration.count()) +"\t";
    std::cout<<outs<<std::endl;
    list->ShowSegmentNumber(false);
    string k = std::to_string(duration.count())+"\n";
    write_into_file(workload_time,k.c_str());
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
    vector<std::pair<KeyType,VaueType>> rg_res;
    KeyType st = 500,ed = 1480000;
    list->Lookup(0,st,ed,rg_res);
    for(auto it:rg_res){
        std::cout<<it.first<<","<<it.second<<std::endl;
    }
    
}

int main(){
    const auto start_time = std::chrono::steady_clock::now();
    GetDataFromText(data_file);
    srand((int)time(0));
    std::cout<<"NUMBERDATA:"<<NUMBERDATA<<"\tTHREAD_NUMBER:"<<THREAD_NUMBER<<std::endl;
    std::cout<<"SkiplistMaxLevel:"<<SkiplistMaxLevel<<"\tmax segment size:"<<SEGMENT_MAX_SIZE<<std::endl
        <<"DELTA_INSERT:"<<DELTA_INSERT<<"\tLinearity:"<<LINAERITY<<std::endl;
    const auto end_time = std::chrono::steady_clock::now();
    const auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);
    std::cout << "prepare time: " << duration.count() << "us"<<std::endl;
    rcu_init(THREAD_NUMBER);
    int n = 0;
    Insert_Part(n,0);  
    #if QUERY_TEST
    #if FOLDKEY
    list->FoldSegment();
    #endif
    n = 0;
    Query_Part(n);
    #endif
    // list->ShowIndexLayer();
    #if SPACEPRINT
    string outs = std::to_string(list->SpaceSize())+"\n";
    write_into_file("./space_size.txt",outs.c_str());
    #endif
}