#include "../include/skiplist.hpp"
#include <fstream>
#include <string>
#include <algorithm>
#include <time.h>

#define NUMBERDATA 50000
#define MaxL 9//(int)(log(NUMBERDATA)/log(2))
#define Gma 128
#define Verify 0
#define NOFINDDEBUG 0

using namespace std;

int length;
vector<snode> dataInput;
vector<unsigned int> dataq0(NUMBERDATA);
int MaxLevel = 0;
vector<snode> nofind;

extern int collsion;
extern int rebuild_cnt;
extern int split_cnt;

void GetData()
{
    int length = 0;
    srand((int)time(0));
    unsigned int *dataIoT;
    std::ifstream inF("../../iotShuffle.bin", std::ios::binary);
    inF.seekg(0, inF.end);
    length = inF.tellg();
    inF.seekg(0, inF.beg);

    // std::cerr << length << std::endl;
    dataIoT = (unsigned int*)malloc(sizeof(unsigned int) * length);

    inF.read(reinterpret_cast<char *>(dataIoT), sizeof(unsigned int) * NUMBERDATA);
    inF.close();
    int j=0;
    for(int i = 0;i<length;i++){
        if(dataIoT[i]!=0){
            dataq0[j] = dataIoT[i];
            j++;
        }
    }
    dataq0.resize(j);
    length = j;
    free(dataIoT);

    // int bound = min(NUMBERDATA,length);

    // for(int i=0;i<bound;i++){
    //     snode x;
    //     x.key = dataq0[i];
    //     int NewLevel = 1;
    //     int P = 50;
    //     while(1){
    //         int t = rand() % 101;
    //         if(t<P)
    //             NewLevel++;
    //         else
    //             break;
    //     }
    //     NewLevel =  min(NewLevel,MaxL);
    //     x.level = NewLevel;
    //     dataInput.push_back(x);
    //     MaxLevel = max(MaxLevel,NewLevel);
    // }
}

// void ExpSearch(skiplist* list){
//     clock_t start,end;
//     double sumTime = 0;
//     int bound = NUMBERDATA;
//     cerr<<"search bound"<<bound<<endl;
//     int nofindcnt = 0;
//     start = clock();
//     for (int i = 0; i < bound; i++) {
//         node* res = list->Search(dataInput[i].key);
// #if NOFINDDEBUG
//         if(!res || res->key!=dataInput[i].key){
//             nofindcnt++;
//             snode x;
//             x.key = dataInput[i].key;
//             x.level = i;
//             nofind.push_back(x);
//         }
// #endif
//     }
//     end = clock();
//     sumTime =(double(end-start)/CLOCKS_PER_SEC);
// #if NOFINDDEBUG
//     cerr<<"no find cnt:"<<nofindcnt<<endl;
// #endif
//     cerr<<"Search time = "<<sumTime<<"s"<<endl;  //输出时间（单位：ｓ）
// }

void ExpInsert(skiplist* list){
    clock_t start,end;
    double sumTime = 0;
    start = clock();
    for(int i = 0;i<NUMBERDATA;i++){
        if(!dataq0[i]) continue;
        list->Insert(dataq0[i],i);
        std::pair<int,Segment_pt::Item*> res = list->Search(dataq0[i]);
        if(!res.first || (res.second)->comp.data.value != i){
            cerr<<"not find"<<endl;
        }
    }
    end = clock();
    sumTime =(double(end-start)/CLOCKS_PER_SEC);
    cerr<<"Insert time = "<<sumTime<<"s"<<endl;  //输出时间（单位：ｓ）
    cerr<<"collision:"<<collsion<<"\trebuild:"<<rebuild_cnt<<"\tsplit"<<split_cnt<<endl;
    // list->show();
}

int main(){
    // MaxLevel = readFromCSV(dataInput);
    GetData();
    length = dataInput.size();
    cerr<<"length:"<<length<<endl;
    // cerr<<"Max Level:"<<MaxLevel<<endl;
    skiplist* list = new skiplist(MaxL,Gma);

    ExpInsert(list);

    // list->setup(dataInput);
    // list->ShowNodeDis();

    //show space size
    // list->ComputeSpace();
    //search test
    // ExpSearch(list);

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