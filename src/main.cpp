#include "skiplist.hpp"
#include <iostream>
#include <thread>

#define MM 1000000
#define NUMBERDATA (1000000)
#define PREINSERT 0//1*MM//(4*MM)
#define SkiplistMaxLevel 8//(int)(log(NUMBERDATA)/log(2))
#define THREAD_NUMBER 2
#define NOFINDDEBUG 0

int key_dis = NUMBERDATA/THREAD_NUMBER;
// unsigned int dataq0[] = {1,2,3,4,5,6,7,8,9,10,11,12};
unsigned int *dataq0 = new unsigned int[NUMBERDATA];
skiplist *list = new skiplist(SkiplistMaxLevel,Gm);

using namespace chrono;

void GetData2(){
    //ycsb64M
    // char unique_dadta_file[] = "/home/yhzhou/datasets/unique_iot_web_unique_shuffle.csv";//unique_iot_web_unique_shuffle
    //unique_iot_web_unique
    // vector<unsigned int> dataq0(26900528);
    // long long length = 0;
    /*数据集*/
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

// void ExpSearch(){
//     clock_t start,end;
//     double sumTime = 0;
// #if NOFINDDEBUG
//     int nofindcnt = 0;
// #endif
//     start = clock();
//     for (int i = 0; i < NUMBERDATA; i++) {
//         if(dataq0[i] == UNINT_MAX || dataq0[i] == 0){
//             continue;
//         }
//         // std::pair<int,int> res = list->Lookup(dataq0[i]);  
// #if NOFINDDEBUG
//         if(!res.first || (res.second)->key != dataq0[i]){
//             cerr<<"not find "<<dataq0[i]<<" "<<i<<endl;
//             nofindcnt++;
//         }
//
// #endif
//     }
//     end = clock();
//     sumTime =(double(end-start)/CLOCKS_PER_SEC);
// #if NOFINDDEBUG
//     cerr<<"no find cnt:"<<nofindcnt<<endl;
// #endif
//     cerr<<"Search time = "<<sumTime<<"s"<<endl;  //输出时间（单位：ｓ）
// }

// void ExpInsert(){
//     clock_t start,end;
//     double sumTime = 0;
//     // int indx = 0;
//     start = clock();
//     for(int i = 0;i<NUMBERDATA;i++){
//         if(dataq0[i] == UNINT_MAX || dataq0[i] == 0){
//             continue;
//         }
//         list->Insert(dataq0[i],i);
//     }
//     end = clock();
//     sumTime =(double(end-start)/CLOCKS_PER_SEC);
//     cerr<<"Insert time = "<<sumTime<<"s"<<endl;  //输出时间（单位：ｓ）
//     // list->show();
// }

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
    for(int i = bound_l;i<bound_r;i++){
        // cerr<<i<<endl;
        mystate.key = dataq0[i];
        mystate.value = i;
        // cout<<i<<endl;
        list->Insert(&mystate);
    }
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
    cerr<<NUMBERDATA<<endl;
    cerr<<"THREAD_NUMBER:"<<THREAD_NUMBER<<endl;
    cerr<<"SkiplistMaxLevel:"<<SkiplistMaxLevel<<endl;
    std::vector<thread> threads;//(THREAD_NUMBER);
    const auto start_time = std::chrono::steady_clock::now();
    //time_point<system_clock> start = system_clock::now();
    // ProfilerStart("test.prof");
    // start = clock();
    for(int i = 0;i<PREINSERT;i++){
        mystate.key = dataq0[i];
        mystate.value = i;
        cerr<<i<<endl;
        list->Insert(&mystate);
    }
    cerr<<"preinsert finish"<<endl;
    int kk = PREINSERT;
    for(int idx=0; idx < THREAD_NUMBER; idx++){
        threads.push_back(thread(test, idx,kk,kk+key_dis));
        kk+=key_dis;
    } 
    for(int idx=0; idx <THREAD_NUMBER; idx++){
        threads[idx].join();
    }
	// ProfilerStop();
    const auto end_time = std::chrono::steady_clock::now();
    const auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);
    std::cout << "insert time: " << duration.count() << "us"<<std::endl;
    for(int i =0;i<NUMBERDATA;i++){
        mystate.key = dataq0[i];
        std::pair<int,int> res = list->Lookup(&mystate);
        if(!res.first || res.second != i){
            cerr<<dataq0[i]<<" "<<i<<endl;
        }
    }
    // free(mystate.ppreds);
    free(mystate.preds);
    // free(mystate.currs);
    free(mystate.succs);
    // list->showList();
}