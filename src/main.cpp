/*
this version is based on subtree + pccs split + split with (greedy plr / random / greedy pccs)
*/

#include "../include/skiplist.hpp"
#include <time.h>
#include <set>


#define MM 1000000
#define NUMBERDATA (24*MM)
#define MaxL 8//(int)(log(NUMBERDATA)/log(2))
#define Verify 0
#define NOFINDDEBUG 0
#define PPROF 0

#if PPROF
    #include<gperftools/profiler.h>
#endif

using namespace std;

int length;
// vector<snode> dataInput;
// vector<unsigned int> dataq0(NUMBERDATA);
int MaxLevel = 0;
vector<snode> nofind;
unsigned int *dataq0 = new unsigned int[NUMBERDATA];
char search_time[] = "./search_time.csv";
char insert_time[] = "./insert_time.csv";
// set<unsigned> keys_unique;

long long real = 0;


void GetData2(){
    char unique_dadta_file[] = "/home/yhzhou/datasets/unique_iot_web_unique_shuffle.csv";//unique_iot_web_unique_shuffle
    //unique_iot_web_unique
    ifstream fp(unique_dadta_file);
    string line;
    for(int i =0;i<NUMBERDATA;i++){
        getline(fp,line);
        dataq0[i] = atoi(line.c_str());
        // cerr<<"key["<<i<<"]"<<dataq0[i]<<endl;    
    }
    fp.close();
}

// void GetData()
// {
//     int length = 0;
//     srand((int)time(0));
//     unsigned int *dataIoT;
//     std::ifstream inF("/home/yhzhou/datasets/iotShuffle.bin", std::ios::binary);
//     inF.seekg(0, inF.end);
//     length = inF.tellg();
//     inF.seekg(0, inF.beg);
//
//     // std::cerr << length << std::endl;
//     dataIoT = (unsigned int*)malloc(sizeof(unsigned int) * length);
//
//     inF.read(reinterpret_cast<char *>(dataIoT), sizeof(unsigned int) * length);
//     inF.close();
//     int j=0;
//     for(int i = 0;i<length && j<NUMBERDATA;i++){
//         if(dataIoT[i]!=0 && !keys_unique.count(dataIoT[i])){
//             dataq0[j] = dataIoT[i];
//             keys_unique.insert(dataIoT[i]);
//             j++;
//         }
//     }
//     cerr<<"key unique:"<<keys_unique.size()<<endl;
//     dataq0.resize(j);
//     length = j;
//     free(dataIoT);
// }

void ExpSearch(skiplist* list){
    clock_t start,end;
    double sumTime = 0;
#if NOFINDDEBUG
    int nofindcnt = 0;
#endif
    start = clock();
    for (int i = 0; i < NUMBERDATA; i++) {
        if(!dataq0[i]) continue;
        std::pair<bool,node*> res = list->Search(dataq0[i]);  
#if NOFINDDEBUG
        if(!res.first || (res.second)->key != dataq0[i]){
            cerr<<"not find "<<dataq0[i]<<" "<<i<<endl;
            nofindcnt++;
        }
#endif
    }
    end = clock();
    sumTime =(double(end-start)/CLOCKS_PER_SEC);
#if NOFINDDEBUG
    cerr<<"no find cnt:"<<nofindcnt<<endl;
#endif
    cerr<<"Search time = "<<sumTime<<"s"<<endl;  //输出时间（单位：ｓ）
    string k = to_string(sumTime)+"\n";
    write_into_file(search_time,k.c_str());
    cerr<<"Average jmp seg:"<<jmp_seg*1.0/real<<"\tjmp subtree:"<<jmp_subtree*1.0/real<<endl;
}

void ExpInsert(skiplist* list){
    clock_t start,end;
    double sumTime = 0;
    // int indx = 0;
#if PPROF
    ProfilerStart("lipp2.prof");
#endif
    start = clock();
    for(int i = 0;i<NUMBERDATA;i++){
        if(!dataq0[i]) continue;
        list->Insert(dataq0[i],i);
        real++;
#if NOFINDDEBUG  
        std::pair<int,node*> res = list->Search(dataq0[i]);

        if(!res.first || (res.second)->value != i){
            cerr<<"not find "<<dataq0[i]<<" "<<i<<endl;
        }  

  
#endif
    }
    end = clock();
#if PPROF
    ProfilerStop();
#endif
    sumTime =(double(end-start)/CLOCKS_PER_SEC);
    cerr<<"Insert time = "<<sumTime<<"s"<<endl;  //输出时间（单位：ｓ）
    string k = to_string(sumTime)+"\n";
    write_into_file(insert_time,k.c_str());
    cerr<<"scan ele cnt:"<<scan_<<"\tcollision:"<<collsion<<"\trebuild:"<<rebuild_cnt<<"\tsplit"<<split_cnt<<endl;
    // list->show();
}

int main(){
    cerr<<"split according to pccs"<<endl;
    cerr<<"USEPLR:"<<USEPLR<<endl;
    cerr<<"USEGREEDYPCCS:"<<USEGREEDYPCCS<<endl;
    cerr<<"LINAERITY:"<<LINAERITY<<endl;
    cerr<<"MAX_DEPTH:"<<MAX_DEPTH<<endl;
    cerr<<"MaxL:"<<MaxL<<endl;
    GetData2();
    cerr<<"total count:"<<NUMBERDATA<<endl;
    // cerr<<"Max Level:"<<MaxLevel<<endl;
    skiplist* list = new skiplist(MaxL,Gm);

    ExpInsert(list);

    ExpSearch(list);

    // list->setup(dataInput);
    // list->ShowNodeDis();

    //show space size
    list->ComputeSpace();
    
    list->show();

    cerr<<"list segments:"<<list->segCnt<<endl;

    return 0;
}