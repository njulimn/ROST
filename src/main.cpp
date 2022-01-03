#include "skiplist.hpp"
#include <iostream>
#include <thread>
// #include<gperftools/profiler.h>

#define MM 1000000
#define NUMBERDATA (32*MM)
#define PREINSERT 4*MM
#define SkiplistMaxLevel 9//(int)(log(NUMBERDATA)/log(2))
#define THREAD_NUMBER 16
#define NOFINDDEBUG 0

int key_dis = (NUMBERDATA-PREINSERT)/THREAD_NUMBER;
// unsigned int dataq0[] = {1,2,3,4,5,6,7,8,9,10,11,12};
unsigned int *dataq0 = new unsigned int[NUMBERDATA];
skiplist *list = new skiplist(SkiplistMaxLevel,Gm);
int conflicts_thread[THREAD_NUMBER];
 
using namespace chrono;

void GetData2(){
    //ycsb64M
    char unique_dadta_file[] = "/home/yhzhou/datasets/ycsb64M.csv";//unique_iot_web_unique_shuffle
    //unique_iot_web_unique
    // vector<unsigned int> dataq0(26900528);
    // long long length = 0;
    /*数据集*/
    ifstream fp(unique_dadta_file);
    string line;
    for(int i =0;i<NUMBERDATA;i++){
        getline(fp,line);
        dataq0[i] = atoi(line.c_str());
        // cerr<<"key["<<i<<"]"<<dataq0[i]<<endl;    
    }
    fp.close();
    // for(int i = 0;i<NUMBERDATA;i++){
    //     if(i%2 == 0)
    //         dataq0[i] = i+1;
    //     else
    //         dataq0[i] = NUMBERDATA+i;
    // }
    // std::random_device rd;
    // std::mt19937 g(rd());
    // std::shuffle(dataq0, dataq0+NUMBERDATA, g);
}

void test(const int id,const int bound_l,const int bound_r ){
    // cerr<<bound_r<<endl;
    skiplist::State mystate;
    // mystate.ppreds = (skiplist::Segment_pt**)malloc(sizeof(skiplist::Segment_pt*)*SkiplistMaxLevel);
    mystate.preds = (skiplist::Segment_pt**)malloc(sizeof(skiplist::Segment_pt*)*SkiplistMaxLevel);
    // mystate.currs = (skiplist::Segment_pt**)malloc(sizeof(skiplist::Segment_pt*)*SkiplistMaxLevel);
    mystate.succs = (skiplist::Segment_pt**)malloc(sizeof(skiplist::Segment_pt*)*SkiplistMaxLevel);
    // memset(mystate.ppreds,0,sizeof(skiplist::Segment_pt*)*SkiplistMaxLevel);
    memset(mystate.preds,0,sizeof(skiplist::Segment_pt*)*SkiplistMaxLevel);
    // memset(mystate.currs,0,sizeof(skiplist::Segment_pt*)*SkiplistMaxLevel);
    memset(mystate.succs,0,sizeof(skiplist::Segment_pt*)*SkiplistMaxLevel);
    int nofind = 0;
    int conflicts = 0;
    int max_conflict = 0;
    for(int i = bound_l;i<bound_r;i++){
        if(dataq0[i] == 0 || dataq0[i] == UNINT_MAX)
            continue;
        mystate.key = dataq0[i];
        mystate.value = i;
        int t = 0;
        // cout<<i<<endl;
        list->Insert(&mystate,t);
        max_conflict = max(max_conflict,t);
        conflicts+=t;
        // cout<<i<<endl;
    }
    conflicts_thread[id] = conflicts;
    cerr<<"max conflict:"<<max_conflict<<endl;
    // free(mystate.ppreds);
    free(mystate.preds);
    // free(mystate.currs);
    free(mystate.succs);
}

void test_query(const int id,const int bound_l,const int bound_r ){
    // cerr<<bound_r<<endl;
    skiplist::State mystate;
    // mystate.ppreds = (skiplist::Segment_pt**)malloc(sizeof(skiplist::Segment_pt*)*SkiplistMaxLevel);
    mystate.preds = (skiplist::Segment_pt**)malloc(sizeof(skiplist::Segment_pt*)*SkiplistMaxLevel);
    // mystate.currs = (skiplist::Segment_pt**)malloc(sizeof(skiplist::Segment_pt*)*SkiplistMaxLevel);
    mystate.succs = (skiplist::Segment_pt**)malloc(sizeof(skiplist::Segment_pt*)*SkiplistMaxLevel);
    // memset(mystate.ppreds,0,sizeof(skiplist::Segment_pt*)*SkiplistMaxLevel);
    memset(mystate.preds,0,sizeof(skiplist::Segment_pt*)*SkiplistMaxLevel);
    // memset(mystate.currs,0,sizeof(skiplist::Segment_pt*)*SkiplistMaxLevel);
    memset(mystate.succs,0,sizeof(skiplist::Segment_pt*)*SkiplistMaxLevel);
    // int nofind = 0;
    for(int i = bound_l;i<bound_r;i++){
        if(dataq0[i] == 0 || dataq0[i] == UNINT_MAX)
            continue;
        mystate.key = dataq0[i];
        mystate.value = i;
        list->Lookup(&mystate);
        // cout<<i<<endl;
        std::pair<int,int> res = list->Lookup(&mystate);
        // if(!res.first || res.second != i){
        //     // cerr<<dataq0[i];
        //     // string kk = to_string(dataq0[i])+"\n";
        //     // write_into_file("./nofind.txt",kk.c_str());
        //     nofind++;
        //     // list->Lookup(&mystate);
        // }
    }
    // cerr<<"nofind:"<<nofind<<endl;
    // free(mystate.ppreds);
    free(mystate.preds);
    // free(mystate.currs);
    free(mystate.succs);
}

int main(){
    GetData2();
    skiplist::State mystate;
    // mystate.ppreds = (skiplist::Segment_pt**)malloc(sizeof(skiplist::Segment_pt*)*SkiplistMaxLevel);
    mystate.preds = (skiplist::Segment_pt**)malloc(sizeof(skiplist::Segment_pt*)*SkiplistMaxLevel);
    // mystate.currs = (skiplist::Segment_pt**)malloc(sizeof(skiplist::Segment_pt*)*SkiplistMaxLevel);
    mystate.succs = (skiplist::Segment_pt**)malloc(sizeof(skiplist::Segment_pt*)*SkiplistMaxLevel);
    // memset(mystate.ppreds,0,sizeof(skiplist::Segment_pt*)*SkiplistMaxLevel);
    memset(mystate.preds,0,sizeof(skiplist::Segment_pt*)*SkiplistMaxLevel);
    // memset(mystate.currs,0,sizeof(skiplist::Segment_pt*)*SkiplistMaxLevel);
    memset(mystate.succs,0,sizeof(skiplist::Segment_pt*)*SkiplistMaxLevel);
    srand((int)time(0));
    cerr<<"NUMBERDATA:"<<NUMBERDATA<<endl;
    cerr<<"THREAD_NUMBER:"<<THREAD_NUMBER<<endl;
    cerr<<"SkiplistMaxLevel:"<<SkiplistMaxLevel<<endl;
    cerr<<"PREINSERT:"<<PREINSERT<<endl;
    cerr<<"segment max size:"<<list->segment_max_size<<endl;
    std::vector<thread> threads;//(THREAD_NUMBER);
    if(THREAD_NUMBER == 16){
        list->segment_max_size = 1e4;
    }
    const auto start_time = std::chrono::steady_clock::now();
    // ProfilerStart("test.prof");
    int t = 0;
    for(int i = 0;i<PREINSERT;i++){
        if(dataq0[i] == 0 || dataq0[i] == UNINT_MAX)
            continue;
        mystate.key = dataq0[i];
        mystate.value = i;
        list->Insert(&mystate,t);
    }
    cerr<<"skip list segment:"<<list->GetSegCnt()<<endl;
    cerr<<"preinsert finish"<<endl;
    const auto end_time1 = std::chrono::steady_clock::now();
    const auto duration1 = std::chrono::duration_cast<std::chrono::microseconds>(end_time1 - start_time);
    std::cout << "insert time: " << duration1.count() << "us"<<std::endl;
    if(THREAD_NUMBER == 16){
        list->segment_max_size = 1e5;
    }
    cerr<<"segment max size:"<<list->segment_max_size<<endl;
    int kk = PREINSERT;
    for(int idx=0; idx < THREAD_NUMBER; idx++){
        // cerr<<kk<<","<<kk+key_dis<<endl;
        threads.push_back(thread(test, idx,kk,kk+key_dis));
        kk+=key_dis;
    } 
    for(int idx=0; idx <THREAD_NUMBER; idx++){
        threads[idx].join();
    }
	// ProfilerStop();
    const auto end_time = std::chrono::steady_clock::now();
    // ProfilerStop();
    const auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);
    std::cout << "insert time: " << duration.count() << "us"<<std::endl;
    cerr<<"skip list segment total:"<<list->GetSegCnt()<<endl;
    long long conflicts = 0;
    for(int i = 0;i<THREAD_NUMBER;i++){
        conflicts+=conflicts_thread[i];
    }
    cerr<<"total conflicts:"<<conflicts<<"\taverage conflicts:"<<conflicts*1.0/(NUMBERDATA-PREINSERT)<<endl;
    // std::vector<thread> threads2;//(THREAD_NUMBER);
    // const auto start_time2 = std::chrono::steady_clock::now();
    // // int nofind = 0;
    // // for(int i =0;i<NUMBERDATA;i++){
    // //     if(dataq0[i] == 0 || dataq0[i] == UNINT_MAX)
    // //         continue;
    // //     mystate.key = dataq0[i];
    // //     std::pair<int,int> res = list->Lookup(&mystate);
    // //     if(!res.first || res.second != i){
    // //         // cerr<<dataq0[i];
    // //         // string kk = to_string(dataq0[i])+"\n";
    // //         // write_into_file("./nofind.txt",kk.c_str());
    // //         nofind++;
    // //     }
    // // }
    // kk = 0;
    // for(int idx=0; idx < THREAD_NUMBER; idx++){
    //     // cerr<<kk<<","<<kk+key_dis<<endl;
    //     threads2.push_back(thread(test_query, idx,kk,kk+key_dis));
    //     kk+=key_dis;
    // } 
    // for(int idx=0; idx <THREAD_NUMBER; idx++){
    //     threads2[idx].join();
    // }
    // const auto end_time2 = std::chrono::steady_clock::now();
    // const auto duration2 = std::chrono::duration_cast<std::chrono::microseconds>(end_time2 - start_time2);
    // std::cout << "query time: " << duration2.count() << "us"<<std::endl;
    // free(mystate.ppreds);
    free(mystate.preds);
    // free(mystate.currs);
    free(mystate.succs);
    // list->showList();
}