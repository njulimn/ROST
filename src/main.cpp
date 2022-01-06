#include "../include/Conskiplist.hpp"
#include <iostream>
#include <thread>

#define MM 1000000
#define NUMBERDATA (64*MM)
#define PREINSERT 4*MM
#define SkiplistMaxLevel 9//(int)(log(NUMBERDATA)/log(2))
#define THREAD_NUMBER 32
#define NOFINDDEBUG 0
#define QUERY_ 0

skiplist<unsigned int,int> *list = new skiplist<unsigned int,int>(SkiplistMaxLevel,Gm);

int main(){
    list->Add(1,2);
    std::pair<int,int> res = list->Lookup(1);
    std::cout<<res.first<<","<<res.second<<std::endl;
}