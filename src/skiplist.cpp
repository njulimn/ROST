#include "../include/skiplist.hpp"
#include <fstream>
#include <string>
#include <algorithm>
#include <time.h>
#include <time.h>
#include <stack>
#include <queue>
// #include <unistd.h>

#define DEBUG 0
#define FMCD 0
// #define Gm 128

int collsion = 0;
int rebuild_cnt = 0;
int split_cnt = 0;


bool check_file_test(char const *fileName)
{
    // 用 ifstream 来判断文件是否存在
    ifstream testFile(fileName);
    if(!testFile)
    {
        cerr << "file not exit" << endl;
        return false;
    }
    else
    {
        cerr << "file exits" << endl;
        return true;
    }
    return false;
}

bool write_into_file(char const *fileName, char const *content)
{
    ofstream out(fileName,ios::app);
    // ofstream out;
    // out.open(fileName);
    if(!out.is_open())
    {
        cerr << "file not exit" << endl;
        return false;
    }
    else
    {
        out << content;
        // cerr << "write succeed" << endl;
        out.close();
        // usleep(100);
        return true;
    }
    
    return false;
}

bool generate_file_test(char const *fileName)
{
    ofstream out;
    out.open(fileName);
    // 判断文件是否已经打开
    if(out.is_open())
    {
        cerr << "file created succeed" << endl;
        out.close();
        return true;
    }
    else
    {
        cerr << "file created failed" << endl;
        out.close();
        return false;
    }
    
    out.close();
    return false;
}

/////////////////////////////////////////////

bool Segment_pt::SplitSegment(Segment_pt** Update,skiplist* list){
    // int factor = 2;
    split_cnt++;
    subtree* x = DataArray;
    int HalfCumItem = 0;
    int posStone[2] = {x->num_items>>1,x->num_items};
    int pt = 0;
    //compute the number of ele(including subtree) in [0,x->num_items>>1)
    for(int i = 0;i<x->num_items;i++){
        if(i>=posStone[pt]){
            break;
        }
        if(BITMAP_GET(x->none_bitmap, i) == 0){
            if(BITMAP_GET(x->child_bitmap, i) == 0){
                HalfCumItem++;
            }
            else{
                HalfCumItem+=(x->items[i].comp.child->size);
            }
        }
    }
    int SmallHalf = min(HalfCumItem,x->size-HalfCumItem);
    int LargerHalf = max(HalfCumItem,x->size-HalfCumItem);
    if(SmallHalf*2 >= LargerHalf)
        return false;
#if DEBUG
    cerr<<"split subtree"<<x<<" into two part:"<<HalfCumItem<<","<<x->size-HalfCumItem<<endl;
    cerr<<"split"<<x;
#endif
    
    const int ESIZE = x->size;
    unsigned int* keys = new unsigned int[ESIZE];
    int* values = new int[ESIZE];
    memset(keys,0,ESIZE);
    memset(values,0,ESIZE);
    unsigned int start_x = x->start,stop_x = x->stop;

    scan_and_destroy_subtree(x,keys,values);
#if DEBUG
    cerr<<"range "<<keys[0]<<"to"<<keys[ESIZE-1]<<endl;
    
    for(int i=0;i<ESIZE;i++){
        if(keys[i]==0){
            cerr<<"scan get key0 from "<<x<<endl;
        }
    }
#endif
    /*分布case 
    左半 all 右半 0--->需要进一步限制stop
    左半 all-1 右半 1--->右半不需要rebuild，直接build_none
    左半 <= all-2 右半 >=2--->右半正常rebuild
    另外三种情况与上述左右互换
    */
   
    subtree* new_subtree1 = nullptr;
    subtree* new_subtree2 = nullptr;
   if(HalfCumItem == LargerHalf){
        if(HalfCumItem <= ESIZE-2){
            new_subtree1 = rebuild_tree(keys,values,HalfCumItem,start_x,keys[HalfCumItem]);
            new_subtree2  = rebuild_tree(keys+HalfCumItem,values+HalfCumItem,ESIZE-HalfCumItem,keys[HalfCumItem],stop_x);
        }
        else if(HalfCumItem == ESIZE -1){
            new_subtree1 = rebuild_tree(keys,values,HalfCumItem,start_x,keys[HalfCumItem]);
            new_subtree2  = build_tree_none();
            BITMAP_CLEAR(new_subtree2->none_bitmap,0);
            new_subtree2->size++;
            new_subtree2->num_inserts++;
            new_subtree2->items[0].comp.data.key = keys[HalfCumItem];
            new_subtree2->items[0].comp.data.value = values[HalfCumItem];
            new_subtree2->start = keys[HalfCumItem];
            new_subtree2->stop = stop_x;
        }
        else{
            new_subtree1 = rebuild_tree(keys,values,ESIZE,start_x,keys[ESIZE-1]);
            this->DataArray = new_subtree1;
            cerr<<"split end"<<endl;
            return true;
        }
#if DEBUG
        cerr<<"replace "<<this->DataArray<<"with "<<new_subtree1<<",insert a new top subtree "<<new_subtree2<<endl;
#endif
        this->DataArray = new_subtree1;
        Segment_pt* newSeg = new Segment_pt(list->RandLevel(),new_subtree2);
        if(newSeg->level > list->skip_level){
            for(int i = newSeg->level;i>list->skip_level;i--){
                Update[i] = list->header;
            }
        }
        if(newSeg->level > this->level){
            for(int i = newSeg->level;i>this->level;i--){
                newSeg->forward[i] = Update[i]->forward[i];
                Update[i]->forward[i] = newSeg;
            }
        }
        for(int i =min(this->level,newSeg->level);i>=1;i--){
            newSeg->forward[i] = Update[i]->forward[i]->forward[i];
            Update[i]->forward[i]->forward[i] = newSeg;
        }
   }
   else{
        if(HalfCumItem >= 2){
            new_subtree1 = rebuild_tree(keys,values,HalfCumItem,start_x,keys[HalfCumItem]);

            new_subtree2  = rebuild_tree(keys+HalfCumItem,values+HalfCumItem,ESIZE-HalfCumItem,keys[HalfCumItem],stop_x);
        }
        else if(HalfCumItem == 1){
            new_subtree1 = build_tree_none();
            new_subtree1->start = start_x;
            new_subtree1->stop = keys[1];
            BITMAP_CLEAR(new_subtree1->none_bitmap,0);
            new_subtree1->size++;
            new_subtree1->num_inserts++;
            new_subtree1->items[0].comp.data.key = keys[0];
            new_subtree1->items[0].comp.data.value = values[0];

            new_subtree2  = rebuild_tree(keys+1,values+1,ESIZE-1,keys[1],stop_x);        
        }
        else{
            new_subtree1 = rebuild_tree(keys,values,ESIZE,keys[0],stop_x);
            this->DataArray = new_subtree1;
#if DEBUG
            cerr<<"split end"<<endl;
#endif
            return true;
        }
#if DEBUG
        cerr<<"replace "<<this->DataArray<<"with "<<new_subtree1<<",insert a new top subtree "<<new_subtree2<<endl;this->DataArray = new_subtree1;
#endif
        this->DataArray = new_subtree1;
        Segment_pt* newSeg = new Segment_pt(list->RandLevel(),new_subtree2);
        if(newSeg->level > list->skip_level){
            for(int i = newSeg->level;i>list->skip_level;i--){
                Update[i] = list->header;
            }
            list->skip_level = newSeg->level;
        }
        if(newSeg->level > this->level){
            for(int i = newSeg->level;i>this->level;i--){
                newSeg->forward[i] = Update[i]->forward[i];
                Update[i]->forward[i] = newSeg;
            }
        }
        for(int i =min(this->level,newSeg->level);i>=1;i--){
            newSeg->forward[i] = Update[i]->forward[i]->forward[i];
            Update[i]->forward[i]->forward[i] = newSeg;
        }
   }
#if DEBUG
    if(new_subtree1){
        // int shouldbe = new_subtree1->size;
        // int real = 0;
        // for(int i = 0;i<new_subtree1->num_items;i++){
        //     if(!BITMAP_GET(new_subtree1->none_bitmap,i)){
        //         if(BITMAP_GET(new_subtree1->child_bitmap,i)){
        //             real+=(new_subtree1->items[i].comp.child->size);
        //         }
        //         else{
        //             real++;
        //         }
        //     }
        // }
        // if(shouldbe!=real){
        //     cerr<<"rebuild error"<<endl;
        // }
        if(!checkTree(new_subtree1)){
            cerr<<"checkTree false"<<endl;
            RT_ASSERT(0);
        }
    }
    if(new_subtree2){
        // int shouldbe = new_subtree2->size;
        // int real = 0;
        // for(int i = 0;i<new_subtree2->num_items;i++){
        //     if(!BITMAP_GET(new_subtree2->none_bitmap,i)){
        //         if(BITMAP_GET(new_subtree2->child_bitmap,i)){
        //             real+=(new_subtree2->items[i].comp.child->size);
        //         }
        //         else{
        //             real++;
        //         }
        //     }
        // }
        // if(shouldbe!=real){
        //     cerr<<"rebuild error"<<endl;
        // }
        if(!checkTree(new_subtree2)){
            cerr<<"checkTree false"<<endl;
            RT_ASSERT(0);
        }
    }
#endif 
    list->segCnt++;
#if DEBUG
    cerr<<"split end"<<endl;
#endif
    return true;   
}

Segment_pt::subtree* Segment_pt::build_tree_none(){
    subtree* n = new_subtree(1);
    n->is_two = 0;
    n->build_size = 0;
    n->size = 0;
    n->fixed = 0;
    n->child_ptr = 0;
    n->num_inserts = n->num_insert_to_data = 0;
    n->num_items = 1;
    n->slope = n->intercept = 0;
    n->start = UNINT_MAX;
    n->stop = UNINT_MAX;
    n->items = new_items(1);
    memset(n->items,0,1);
    n->none_bitmap = new_bitmap(1);
    n->none_bitmap[0] = 0;
    BITMAP_SET(n->none_bitmap, 0);
    n->child_bitmap = new_bitmap(1);
    n->child_bitmap[0] = 0;
    return n;
}

Segment_pt::subtree* Segment_pt::build_tree_two(unsigned int key1, int value1, unsigned int key2, int value2,unsigned int start,unsigned int stop,subtree* x)
{
    if (key1 > key2) {
        std::swap(key1, key2);
        std::swap(value1, value2);
    }
    // RT_ASSERT(key1 < key2);
    // static_assert(BITMAP_WIDTH == 8);

    subtree* n = NULL;
    if (pending_two.empty()) {
        n = new_subtree(1);
        n->is_two = 1;
        n->build_size = 2;
        n->size = 2;
        n->fixed = 0;
        n->child_ptr = 0;
        n->num_inserts = n->num_insert_to_data = 0;
        n->num_items = 8;
        n->items = new_items(n->num_items);
        memset(n->items,0,n->num_items);
        n->none_bitmap = new_bitmap(1);
        n->child_bitmap = new_bitmap(1);
        n->none_bitmap[0] = 0xff;
        n->child_bitmap[0] = 0;
    } else {
        n = pending_two.top();
        n->child_ptr = 0;
        memset(n->child_bitmap,0,n->num_items);
        memset(n->none_bitmap,0xff,n->num_items);
        pending_two.pop();
    }
    if(!x){
        n->start = start;
        n->stop = stop;  
    }else{
    //从size=1增加元素重新生成subtree,[start,stop)
        n->start = max(x->start,2*key1-key2);
        n->stop = min(x->stop,2*key2-key1);
    }
    
    const long double mid1_key = static_cast<long double>(key1)-static_cast<long double>(n->start);
    const long double mid2_key = static_cast<long double>(key2)-static_cast<long double>(n->start);

    const double mid1_target = n->num_items / 3;
    const double mid2_target = n->num_items * 2 / 3;

    n->slope = (mid2_target - mid1_target) / (mid2_key - mid1_key);
    if(isinf(n->slope)){
        cerr<<"197slope chushu:(mid2_key - mid1_key)"<<(mid2_key - mid1_key)<<endl;
        cerr<<"inf!"<<endl;
        RT_ASSERT(0);
    }
    n->intercept = mid1_target - n->slope * mid1_key;
    

    // RT_ASSERT(isfinite(n->slope));
    // RT_ASSERT(isfinite(n->intercept));

    { // insert key1&value1
        int pos = predict(n, key1);
        // RT_ASSERT(BITMAP_GET(n->none_bitmap, pos) == 1);
        BITMAP_CLEAR(n->none_bitmap, pos);
        n->items[pos].comp.data.key = key1;
        n->items[pos].comp.data.value = value1;
    }
    { // insert key2&value2
        int pos = predict(n, key2);
        // RT_ASSERT(BITMAP_GET(n->none_bitmap, pos) == 1);
        BITMAP_CLEAR(n->none_bitmap, pos);
        n->items[pos].comp.data.key = key2;
        n->items[pos].comp.data.value = value2;
    }
    RT_ASSERT(n!=NULL);

    return n;
}

Segment_pt::subtree* Segment_pt::rebuild_tree(unsigned int *_keys,int*  _values,int _size,unsigned int start,unsigned int stop){
    RT_ASSERT(_size > 1);
    
    typedef struct {
            int begin;
            int end;
            int lvl; // top lvl = 1
            subtree* n;
    } Seg;
    std::queue<Seg> s;
    subtree* ret = new_subtree(1);
    ret->child_ptr = 0;
    ret->start = start;
    ret->stop = stop;
#if DEBUG
    cerr<<"rebuild tree new subtree=ret"<<ret<<endl;
#endif
    s.push({0, _size, 1, ret});
    while(!s.empty()){
        rebuild_cnt++;
        const int begin = s.front().begin;
        const int end = s.front().end;
        const int lvl = s.front().lvl;
#if DEBUG
        cerr<<"rebuild level"<<lvl<<endl;
#endif
        subtree* n = s.front().n;
        s.pop();
        RT_ASSERT(end - begin >= 2);
        if(end-begin == 2){
            //n要不是root要不是子subtree，start和stop都事先确定了
            subtree* _ = build_tree_two(_keys[begin], _values[begin], _keys[begin+1], _values[begin+1],n->start,n->stop);
            memcpy(n, _, sizeof(subtree));
            delete_subtree(_, 1);
            if(lvl == 1){
                ret = n;
            }
#if DEBUG
            cerr<<"real put ele :2 with btw"<<endl;
#endif
        }
        else{
            unsigned int* keys = _keys + begin;
            int* values = _values + begin;
            const int size = end - begin;
            const int BUILD_GAP_CNT = compute_gap_count(size);

            n->is_two = 0;
            n->build_size = size;
            n->size = size;
            n->fixed = 0;
            n->num_inserts = n->num_insert_to_data = 0;
#if FMCD
        {
            const int L = size * static_cast<int>(BUILD_GAP_CNT + 1);
            int i = 0;
            int D = 1;
            RT_ASSERT(D <= size-1-D);
            double Ut = (static_cast<long double>(keys[size - 1 - D]-static_cast<long double>(n->start)) - static_cast<long double>(keys[D]-static_cast<long double>(n->start))) /
                        (static_cast<double>(L - 2)) + 1e-6;
            while (i < size - 1 - D) {
                while (i + D < size && keys[i + D] - keys[i] >= Ut) {
                    i ++;
                }
                if (i + D >= size) {
                    break;
                }
                D = D + 1;
                if (D * 3 > size) break;
                RT_ASSERT(D <= size-1-D);
                Ut = (static_cast<long double>(keys[size - 1 - D]-static_cast<long double>(n->start)) - static_cast<long double>(keys[D]-static_cast<long double>(n->start))) /
                        (static_cast<double>(L - 2)) + 1e-6;
            }
            if (D * 3 <= size) {
                // stats.fmcd_success_times ++;
                if(abs(Ut)<1e-8){
                    cerr<<"slope may not a number"<<endl;
                }
                n->slope = 1.0 / Ut;
                if(isinf(n->slope)){
                    cerr<<"slope chushu Ut:"<<Ut<<endl;
                    cerr<<"inf!"<<endl;
                    RT_ASSERT(0);
                }
                n->intercept = (L - n->slope * (static_cast<long double>(keys[size - 1 - D]-static_cast<long double>(n->start)) +
                                                        static_cast<long double>(keys[D]-static_cast<long double>(n->start)))) / 2;
                // RT_ASSERT(isfinite(node->slope));
                // RT_ASSERT(isfinite(node->intercept));
                n->num_items = L;
                

            } else {
                // stats.fmcd_broken_times ++;
                
                int mid1_pos = (size - 1) / 3;
                int mid2_pos = (size - 1) * 2 / 3;

                // RT_ASSERT(0 <= mid1_pos);
                // RT_ASSERT(mid1_pos < mid2_pos);
                // RT_ASSERT(mid2_pos < size - 1);

                const long double mid1_key = (static_cast<long double>(keys[mid1_pos]-static_cast<long double>(n->start)) +
                                                static_cast<long double>(keys[mid1_pos + 1]-static_cast<long double>(n->start))) / 2;
                const long double mid2_key = (static_cast<long double>(keys[mid2_pos]-static_cast<long double>(n->start)) +
                                                static_cast<long double>(keys[mid2_pos + 1]-static_cast<long double>(n->start))) / 2;

                n->num_items = size * static_cast<int>(BUILD_GAP_CNT + 1);
                const double mid1_target = mid1_pos * static_cast<int>(BUILD_GAP_CNT + 1) + static_cast<int>(BUILD_GAP_CNT + 1) / 2;
                const double mid2_target = mid2_pos * static_cast<int>(BUILD_GAP_CNT + 1) + static_cast<int>(BUILD_GAP_CNT + 1) / 2;
                
                if(abs(mid2_key - mid1_key) < 1e-6){
                    cerr<<"slope may not a number"<<endl;
                }
                n->slope = (mid2_target - mid1_target) / (mid2_key - mid1_key);
                if(isinf(n->slope)){
                    cerr<<"slope chushu:(mid2_key - mid1_key)"<<(mid2_key - mid1_key)<<endl;
                    cerr<<"inf!"<<endl;
                    RT_ASSERT(0);
                }
                n->intercept = mid1_target - n->slope * mid1_key;
                // RT_ASSERT(isfinite(node->slope));
                // RT_ASSERT(isfinite(node->intercept));
                
            }
        }
#else
        {
            int mid1_pos = (size - 1) / 3;
            int mid2_pos = (size - 1) * 2 / 3;

            RT_ASSERT(0 <= mid1_pos);
            RT_ASSERT(mid1_pos < mid2_pos);
            RT_ASSERT(mid2_pos < size - 1);

            const long double mid1_key =
                    (static_cast<long double>(keys[mid1_pos])-static_cast<long double>(n->start) + static_cast<long double>(keys[mid1_pos + 1])-static_cast<long double>(n->start) ) / 2;
            const long double mid2_key =
                    (static_cast<long double>(keys[mid2_pos]) -static_cast<long double>(n->start)+ static_cast<long double>(keys[mid2_pos + 1]) -static_cast<long double>(n->start)) / 2;

            n->num_items = size * static_cast<int>(BUILD_GAP_CNT + 1);
            //有点不理解
            const double mid1_target = mid1_pos * static_cast<int>(BUILD_GAP_CNT + 1) + static_cast<int>(BUILD_GAP_CNT + 1) / 2;
            const double mid2_target = mid2_pos * static_cast<int>(BUILD_GAP_CNT + 1) + static_cast<int>(BUILD_GAP_CNT + 1) / 2;
                        RT_ASSERT(0 <= mid1_pos);
            RT_ASSERT(abs(mid2_key - mid1_key)>1e-8);
            n->slope = (mid2_target - mid1_target) / (mid2_key - mid1_key);
            n->intercept = mid1_target - n->slope * mid1_key;
            RT_ASSERT(isfinite(n->slope));
            RT_ASSERT(isfinite(n->intercept));
        }        
#endif
            
#if DEBUG           
            if(n->slope < 0){
                cerr<<"slope < 0"<<endl;
            }
#else
            RT_ASSERT(n->slope >= 0);
#endif
            const int lr_remains = static_cast<int>(size * n->BUILD_LR_REMAIN);
            n->intercept += lr_remains;
            n->num_items += lr_remains * 2;
            //fix size?
            if (size > 2e4) {
                n->fixed = 1;
            }

            n->items = new_items(n->num_items);
            memset(n->items,0,n->num_items);
            const int bitmap_size = BITMAP_SIZE(n->num_items);
            n->none_bitmap = new_bitmap(bitmap_size);
            n->child_bitmap = new_bitmap(bitmap_size);
            memset(n->none_bitmap, 0xff, sizeof(bitmap_t) * bitmap_size);
            memset(n->child_bitmap, 0, sizeof(bitmap_t) * bitmap_size);
            int real_put_ele = 0;
            for (int item_i = predict(n, keys[0]), offset = 0; offset < size; ) {
                //item_i 表示当前元素应插入的位置
                //next_i 表示下一个元素应插入的位置
                int next = offset + 1, next_i = -1;
                while (next < size) {     
                    next_i = predict(n, keys[next]);
                 
#if DEBUG
                    cerr<<"keys["<<next<<"]"<<keys[next]<<"predict:"<<next_i<<endl;
#endif
                    if (next_i == item_i) {
                        next ++;
                    } else {
                        break;
                    }
                }
                if (next == offset + 1) {
                    
#if DEBUG
                    cerr<<"set key in pos:"<<item_i<<endl;
#endif
                    real_put_ele++;
#if DEBUG
                    if(!(BITMAP_GET(n->none_bitmap,item_i))){
                        cerr<<"fugai"<<endl;
                    }
#else
                    RT_ASSERT(BITMAP_GET(n->none_bitmap,item_i) == 1);
                    RT_ASSERT(BITMAP_GET(n->child_bitmap,item_i) == 0);
#endif
                    BITMAP_CLEAR(n->none_bitmap, item_i);
                    n->items[item_i].comp.data.key = keys[offset];
                    n->items[item_i].comp.data.value = values[offset];
                } else {
                    // ASSERT(next - offset <= (size+2) / 3);
                    RT_ASSERT(BITMAP_GET(n->none_bitmap,item_i) == 1);
                    RT_ASSERT(BITMAP_GET(n->child_bitmap,item_i) == 0);
                    BITMAP_CLEAR(n->none_bitmap, item_i);
                    BITMAP_SET(n->child_bitmap, item_i);
                    n->child_ptr++;
                    n->items[item_i].comp.child = new_subtree(1);
                    
                    n->items[item_i].comp.child->child_ptr = 0;
                    std::pair<unsigned int,unsigned int> line = computeRange(item_i,n);
                    n->items[item_i].comp.child->start = min(line.first,keys[offset]);
                    n->items[item_i].comp.child->stop = line.second;
                    real_put_ele+=(next-offset);
#if DEBUG
                    cerr<<"set child "<<item_i<<" with:"<<n->items[item_i].comp.child<<endl;
#endif
                    s.push({begin + offset, begin + next, lvl + 1, n->items[item_i].comp.child});
                }
                if (next >= size) {
                    break;
                } else {
                    item_i = next_i;
#if DEBUG
                    cerr<<"item i:"<<item_i<<",next:"<<next<<endl;
#endif
                    offset = next;
                }
            }

#if DEBUG
            cerr<<"real put ele:"<<real_put_ele<<"in "<<n<<"n:size:"<<n->size<<"diff:"<<n->size-real_put_ele<<endl;
#endif
        }   
    }
#if DEBUG
    cerr<<"return ret:"<<ret<<endl;


    if(!checkTree(ret)){
            cerr<<"checkTree false"<<endl;
            RT_ASSERT(0);
    }
#endif
    return ret;

}

void Segment_pt::destroy_subtree(subtree* root){
    std::stack<subtree*> s;
    s.push(root);
    while (!s.empty()) {
        subtree* n = s.top(); s.pop();

        for (int i = 0; i < n->num_items; i ++) {
            if (BITMAP_GET(n->child_bitmap, i) == 1) {
                s.push(n->items[i].comp.child);
            }
        }

        if (n->is_two) {
            // RT_ASSERT(n->build_size == 2);
            // RT_ASSERT(n->num_items == 8);
            n->size = 2;
            n->num_inserts = n->num_insert_to_data = 0;
            n->none_bitmap[0] = 0xff;
            n->child_bitmap[0] = 0;
            pending_two.push(n);
        } else {
            delete_items(n->items, n->num_items);
            const int bitmap_size = BITMAP_SIZE(n->num_items);
            delete_bitmap(n->none_bitmap, bitmap_size);
            delete_bitmap(n->child_bitmap, bitmap_size);
            delete_subtree(n, 1);
        }
    }
}

void Segment_pt::scan_and_destroy_subtree(subtree* _root,unsigned int*keys,int* values, bool destory){
    typedef std::pair<int, subtree*> Seg; // <begin, subtree*>
    std::stack<Seg> s;

    s.push(Seg(0, _root));
    int keylen = _root->size;
    subtree* prec = _root;
    while (!s.empty()) {
        int begin = s.top().first;
        subtree* n = s.top().second;
        const int SHOULD_END_POS = begin + n->size;
        s.pop();
#if DEBUG
        cerr<<"\n"<<n<<"\'s items fill keys from"<<begin<<" to "<<SHOULD_END_POS<<endl;
#endif
        for (int i = 0; i < n->num_items; i ++) {
            if (BITMAP_GET(n->none_bitmap, i) == 0) {
                if (BITMAP_GET(n->child_bitmap, i) == 0) {
                    keys[begin] = n->items[i].comp.data.key;
#if DEBUG
                    cerr<<i<<" "<<keys[begin]<<",";
#endif
                    values[begin] = n->items[i].comp.data.value;
                    begin ++;
                } else {
                    Seg next_seg;
                    next_seg.first = begin;
                    next_seg.second = n->items[i].comp.child;
                    s.push(next_seg);
                    prec = n;
#if DEBUG
                    cerr<<"\n"<<i<<" s.push(next_seg);"<<n->items[i].comp.child<<" size:"<<n->items[i].comp.child->size<<endl;
#endif
                    begin += n->items[i].comp.child->size;
                }
            }
        }
#if DEBUG
        if(SHOULD_END_POS != begin){
            int a = 5;
            cerr<<"SHOULD_END_POS != begin"<<endl;
        }
#else
        RT_ASSERT(SHOULD_END_POS == begin);
#endif
        if (destory) {
            if (n->is_two) {
                RT_ASSERT(n->build_size == 2);
                RT_ASSERT(n->num_items == 8);
                n->size = 2;
                n->num_inserts = n->num_insert_to_data = 0;
                n->none_bitmap[0] = 0xff;
                n->child_bitmap[0] = 0;
                pending_two.push(n);
            } 
            else {
                delete_items(n->items, n->num_items);
                const int bitmap_size = BITMAP_SIZE(n->num_items);
                delete_bitmap(n->none_bitmap, bitmap_size);
                delete_bitmap(n->child_bitmap, bitmap_size);
                delete_subtree(n, 1);
            }
        }
    } 
}

void Segment_pt::insert_subtree(unsigned int key,int value,Segment_pt** Update,skiplist* list){
    constexpr int MAX_DEPTH = 10;
    subtree* path[MAX_DEPTH];
    int path_size = 0;
    int insert_to_data = 0;
    subtree* n = this->DataArray;
    while(1){
        path[path_size++] = n;
        int pos = predict(n,key);
        if (BITMAP_GET(n->none_bitmap, pos) == 1) {
#if DEBUG
                if(!checkTree(n)){
                    cerr<<"checkTree false"<<endl;
                    cerr<<"key:"<<key<<endl;
                    RT_ASSERT(0);
                }
#endif
                RT_ASSERT(BITMAP_GET(n->child_bitmap, pos) == 0);
                BITMAP_CLEAR(n->none_bitmap, pos);
                n->items[pos].comp.data.key = key;
                n->items[pos].comp.data.value = value;
                n->size++;
                n->num_inserts++;                
                break;
        } else if (BITMAP_GET(n->child_bitmap, pos) == 0) {
            //key冲突
            if(n->items[pos].comp.data.key == key){
                n->items[pos].comp.data.value = value;
                n->size++;
                n->num_inserts++; 
                break;
            }
            collsion++;
            
            if(path_size == 1 && n->size == 1){
                //从tree_one到tree_two，将tree_one作为参数传入后限制start、stop范围
                subtree* newtree =  build_tree_two(key, value, n->items[pos].comp.data.key, n->items[pos].comp.data.value,0,0,n);
                memcpy(this->DataArray, newtree, sizeof(subtree));
                path[0] = this->DataArray;
                delete_subtree(newtree, 1);
#if DEBUG
                if(!checkTree(this->DataArray)){
                    cerr<<"checkTree false"<<endl;
                    RT_ASSERT(0);
                }
#endif
                break;
            }
            std::pair<unsigned int,unsigned int>line = computeRange(pos,n);
            BITMAP_SET(n->child_bitmap, pos);
            n->child_ptr++;
            //产生非顶层subtree，此情况下star stop考虑其父tree
            if(n->is_two){
                n->is_two = 0;
            }
            n->size++;
            n->num_inserts++;
            n->items[pos].comp.child  = build_tree_two(key, value, n->items[pos].comp.data.key, n->items[pos].comp.data.value,min(line.first,min(key,n->items[pos].comp.data.key)),line.second);
            insert_to_data = 1;
            
            break;
        } else {
            n->size++;
            n->num_inserts++;
#if !DEBUG
            RT_ASSERT(n->items[pos].comp.child != nullptr);
#else
            if(! BITMAP_GET(n->none_bitmap, pos) && !n->items[pos].comp.child)
            {
                cerr<<n->none_bitmap[(pos) / BITMAP_WIDTH]<<endl;
                cerr<<n->child_bitmap[(pos) / BITMAP_WIDTH]<<endl;
                cerr<<"shenqi chuxianle!!"<<endl;
                RT_ASSERT(0);
            }
#endif
            n = n->items[pos].comp.child;
        }
    }
    for(int i = 0;i<path_size;i++){
        path[i]->num_insert_to_data += insert_to_data;
    }

    for(int i = 0;i<path_size;i++){
        n = path[i];
        const int num_inserts = n->num_inserts;
        const int num_insert_to_data = n->num_insert_to_data;
        //待调整
        const bool need_rebuild = n->fixed == 0 && n->size >= n->build_size * 4 && n->size >= 64 && num_insert_to_data * 10 >= num_inserts;
        // int T;
        if(need_rebuild || (!i && path_size >= MAX_DEPTH)){
            if(!i){
                bool split = SplitSegment(Update,list);
                if(split)
                    return;
            }
#if DEBUG
            cerr<<"rebuild"<<n<<"it has "<<n->child_ptr<<" child ptr"<<endl;
#endif
            const int ESIZE = n->size;
            unsigned int* keys = new unsigned int[ESIZE];
            int* values = new int[ESIZE];
            unsigned int n_start = n->start,n_stop = n->stop;
            memset(keys,0,ESIZE);
            memset(values,0,ESIZE);
#if DEBUG
            if(!checkTree(n)){
                cerr<<"checkTree false"<<endl;
                RT_ASSERT(0);
            }
#endif
    
            scan_and_destroy_subtree(n,keys,values);
            subtree* newsubtree = rebuild_tree(keys,values,ESIZE,n_start,n_stop);
#if DEBUG
            cerr<<"(level:"<<i<<")replace it with"<<newsubtree<<" it's child ptr:"<<newsubtree->child_ptr<<endl;
#endif
            delete[] keys;
            delete[] values;
#if DEBUG
            if(!checkTree(newsubtree)){
                cerr<<"checkTree false"<<endl;
                RT_ASSERT(0);
            }
#endif
            if(i>0){
                int pos = predict(path[i-1],key);
                RT_ASSERT(BITMAP_GET(path[i-1]->none_bitmap, pos) == 0);
                RT_ASSERT(BITMAP_GET(path[i-1]->child_bitmap, pos) == 1);
                path[i-1]->items[pos].comp.child = newsubtree;
            }
            else{
                this->DataArray = newsubtree;
            }
            return;

        }
    }

}

std::pair<int,Segment_pt::Item*> Segment_pt::find_subtree(unsigned int key){
    subtree* n = this->DataArray;
    while(1){
        int pos = predict(n,key);
        if (BITMAP_GET(n->none_bitmap, pos) == 1) {
            //checkTypeEqualNull
            break;
            // return {false,0};
        }
        else if(BITMAP_GET(n->child_bitmap, pos) == 0){
            //checkTypeEqualData
            return {n->items[pos].comp.data.key == key,&(n->items[pos])};
        }
        else{
            n = n->items[pos].comp.child;
        }
    }
    return {false,0};
}

//////////////////////////

// void verify_gamma(double gamma, vector<unsigned int> data, vector<Segment*>seg) {
// 	int j = 0;
// 	Segment* cur = seg[j];
// 	for (int i = 0; i < data.size(); i++)
// 	{
// 		double pred = 0;
// 		while (j < seg.size()) {
// 			if (data[i] >= cur->start) {
// 				if (data[i] < cur->stop) {
// 					break;
// 				}
// 				j++;
// 				cur = seg[j];
// 			}
// 			else {
// 				cerr << "data out of bound" << endl;
// 				break;
// 			}
// 		}
// 		if (j < seg.size()) {
// 			pred = cur->slope * data[i] + cur->intercept;
// 			if (fabs(pred - i) <= gamma) {
// 				cerr << "within error tolerence ";
// 			}
// 			else {
// 				cerr << "out of error tolerence ";
// 			}
// 			cerr << "data [" << i << "]" << data[i] << " pred:" << pred << endl;
// 		}
// 	}
// }

// void skiplist::insert_static(vector<Segment_pt*> &Update,Segment* seg,unsigned int st,unsigned int ed,node* input,int size,int level){
//     Segment_pt* newSeg = new Segment_pt(size,input,level,st,ed,seg->slope,seg->intercept);
//     //insert
//     if(level > this->level){
//         for(int l = level;l>this->level;l--){
//             Update[l] = this->header;
//         }
//         this->level = level;
//     }
//     for(int l = 1;l<=level;l++){
//         newSeg->forward[l] = Update[l]->forward[l];
//         Update[l]->forward[l] = newSeg;
//         Update[l] = newSeg;
//     }
//     // Seg2Data[newSeg] = input;
//     this->size++;
// }

// void skiplist::setup(vector<snode> input){
//     //init
//     vector<Segment_pt*> Update(this->MaxLevel+1);
//     for(int i = 0;i<=this->MaxLevel;i++){
//         Update[i] = this->header;
//     }
//     int st = 0,ed =0;
// 
//     srand((int)time(0));
// 
//     GreedyPLR* plr = new GreedyPLR(this->gamma);
//     vector<Segment*> seg;
//     vector<int> segment_stIndex;
//     for (int i = 0; i < input.size(); i++) {
// 		Segment* seg_res = nullptr;
// 		seg_res = plr->Process(input[i].key, i-st);
// 		if(seg_res) {
//             segment_stIndex.push_back(st);
//             seg.push_back(seg_res);
//             st = ed = i;
// 		}
//         else{
//             ed = i;
//         }
// 	}
//     Segment* seg_res = nullptr;
// 	seg_res = plr->finish();
// 	if (seg_res) {
//         segment_stIndex.push_back(st);
//         seg.push_back(seg_res);
// 		// vector<node> inputPart(ed-st+1);
//         // for(int l=0;l<ed-st+1;l++){
//         //     inputPart[l].key = input[l+st].key;
//         //     inputPart[l].value = inputPart[l].key;
//         // }
//         // this->insert_static(Update,seg_res,input[st].key,UNINT_MAX,inputPart,input[st].level);
// 	}
//     cerr<<"seg size:"<<seg.size()<<endl;
//     this->MaxLevel = (log(seg.size())/log(2));
// 
//     //segment_data = (int**)malloc(sizeof(int*)*seg.size());//申请seg.size个空间存储每段数据的数组首地址
// 
//     for(int i=0;i<seg.size();i++){
//         st = segment_stIndex[i];
//         if(i == seg.size()-1)
//             ed = input.size()-1;
//         else
//             ed = segment_stIndex[i+1]-1;
//         node* data_segement = nullptr;
//         data_segement = (node*)malloc(sizeof(node)*(ed-st+1));
//         // vector<node> inputPart(ed-st+1);
//         for(int l=0;l<ed-st+1;l++){
//             data_segement[l].key = input[l+st].key;
//             data_segement[l].value = data_segement[l].key;
//         }
//         int level = this->RandLevel();
//         this->insert_static(Update,seg[i],input[st].key,input[ed].key,data_segement,ed-st+1,level);
//     }
// 
//     // unsigned int nodeLevel = 1;
//
// }

// node* skiplist::binarySearch(Segment_pt* x,unsigned int key){
//     int pred = x->slope*key+x->intercept;
//     int l=max((pred-this->gamma),0);
//     // node* data = Seg2Data[x];
//     int r=min(x->node_size-1,pred+this->gamma),mid;
//     while(l<=r){
//         mid = l+(r-l)/2;
//         if(x->nodes[mid].key == key){
//             return &(x->nodes[mid]);
//         }
//         else if(x->nodes[mid].key > key){
//             r = mid-1;
//         }
//         else{
//             l = mid+1;
//         }
//     }
//     return nullptr;
// }

std::pair<int,Segment_pt::Item*> skiplist::Search(unsigned int key){
    Segment_pt* x = this->header;
    // unsigned int pred;
    for(int i = this->skip_level;i>=1;i--){
        while(key >= x->forward[i]->DataArray->stop){
            x = x->forward[i];
        }
        if( (key > x->forward[i]->DataArray->start && x->forward[i]->DataArray->size == 1)||
        (key >= x->forward[i]->DataArray->start)
        ){
            x = x->forward[i];
            std::pair<int,Segment_pt::Item*> res = x->find_subtree(key);
            //this->binarySearch(x,key);
           return res;
        }
    }    
    return {0,0};

}

void skiplist::Insert(unsigned int key,int value){
    Segment_pt* x = this->header;
    Segment_pt** Update = nullptr;//store the predecessor in different level of locating segment
    Update = (Segment_pt**)malloc(sizeof(Segment_pt*)*(this->MaxLevel+1));
    // for(int i = 1;i<=this->MaxLevel;i++){
    //     Update[i] = this->header;
    // }
    Segment_pt* locate = nullptr;
    for(int i = this->skip_level;i>=1;i--){
        while(key >= x->forward[i]->DataArray->stop){
        //if segment's number of element is 1,x turn forward only when key >= stop
        //if segment's number of element > 1,x turn forward only when key > stop
            x = x->forward[i];
        }
        
        //we find the rightmost of segment in level i which is smaller than key
        Update[i] = x;

        if(!locate && 
        ((key > x->forward[i]->DataArray->start && x->forward[i]->DataArray->size == 1)||
        (key >= x->forward[i]->DataArray->start)
        )){//if segment's number of element is 1,locate segment only when key > start
        //if segment's number of element >1,locate segment only when key>=start
            locate = x->forward[i];
        }
    }
    if(locate){
        std::pair<int,Segment_pt::Item*> find_ = locate->find_subtree(key);
        if(find_.first){
            (find_.second)->comp.data.value = value;
            return;
        }
        locate->insert_subtree(key,value,Update,this);
    }
    else{
        //create a new segment
        Segment_pt* newSeg = new Segment_pt(this->RandLevel());
        newSeg->DataArray->size++;
        newSeg->DataArray->num_inserts++;
        newSeg->DataArray->items[0].comp.data.key = key;
        newSeg->DataArray->items[0].comp.data.value = value;
        BITMAP_CLEAR(newSeg->DataArray->none_bitmap, 0);
        if(newSeg->level > this->skip_level){
            for(int i = newSeg->level;i>this->skip_level;i--){
                Update[i] = this->header;
            }
            this->skip_level = newSeg->level;
        }
        newSeg->DataArray->stop = Update[1]->forward[1]->DataArray->start;
        for(int i = 1;i<=newSeg->level;i++){
            newSeg->forward[i] = Update[i]->forward[i];
            Update[i]->forward[i] = newSeg;
        }
        if(Update[1] == this->header){
            newSeg->DataArray->start = 0;
        }else{
            newSeg->DataArray->start = Update[1]->DataArray->stop;
        }
        this->segCnt++;
    }
}

void skiplist::show(){
    Segment_pt* x = header;
    for(int i = 0;i<segCnt;i++){
        x = x->forward[1];
        cerr<<"Segment_pt level:"<<x->level<<" start:"<<x->DataArray->start<<" stop:"<<x->DataArray->stop<<" slope:"<<x->DataArray->slope<<" intertectpt:"<<x->DataArray->intercept<<"\n";
        x->showArray(x->DataArray);
    }
}

// void skiplist::ComputeSpace(){
//     unsigned int space = 0;
//     Segment_pt *x = this->header;
//     space+=sizeof(this);
//     while (x && x->forward[1] != this->header) {
//         // cerr<<"node:"<<sizeof(*(x->forward[1]))<<endl;
//         space+=(sizeof(*(x->forward[1])));
//         space+=((this->level+1) * sizeof(Segment_pt*));
//         space+=((x->forward[1]->node_size) * sizeof(node));
//         x = x->forward[1];
//     }
//     // space+=( (sizeof(Segment_pt*)+sizeof(node*) )*Seg2Data.size() );
//
//     cerr<<"space size:"<<space<<endl;
// }

// void skiplist::ShowNodeDis(){
//     Segment_pt *x = this->header;
//     char filePath[] = "./log/exp_log_dis.txt";
//     char filePath2[] = "./log/exp_log_plr.csv";
//     // char filePath2[] = "./log/exp_log_dis2.txt";
//     vector<int> Count(this->level+1,0);
//     map<int,int> nodemap;
//     vector<vector<int>> jump(this->level+1,vector<int>(1));
//     int i = 1;
//     while (x && x->forward[1] != this->header) {
//         x->forward[1]->show(1);
//         if(!Count[x->forward[1]->level]){
//             jump[x->forward[1]->level][0] = i;
//         }
//         else{
//             jump[x->forward[1]->level].push_back(i);
//         }
//         Count[x->forward[1]->level]++;
//         if(!nodemap.count(x->forward[1]->node_size)){
//             nodemap[x->forward[1]->node_size] = 1;
//         }
//         else{
//             nodemap[x->forward[1]->node_size]++;
//         }
//         i++;
//         x = x->forward[1];
//     }
//    
//     string s1 = "NodeDistribution---------------------\nthe number of segment:"+to_string (this->size)+"\n";
//     write_into_file(filePath,s1.c_str());
//
//     map<int,int>::iterator it;
//     for(it = nodemap.begin();it!=nodemap.end();it++){
//         string sit = to_string((*it).first)+","+to_string((*it).second)+"\n";
//         write_into_file(filePath2,sit.c_str());
//     }
//
//     for(int i=1;i<=this->level;i++){
//         string s2 = "level "+to_string (i)+",segment count "+to_string(Count[i])+"\n";
//         write_into_file(filePath,s2.c_str());
//         for(int j = 0;j<jump[i].size();j++){
//             string s3 = to_string(jump[i][j])+"\t";
//             write_into_file(filePath,s3.c_str());
//         }
//         write_into_file(filePath,"\n");
//     }
// }

// void Segment_pt::show(int w){
//     char filePath[] = "./log/exp_log.txt";
//     string s1 ="Segment_pt level:"+to_string(this->level)+" start:"+to_string(this->start)+" stop:"+to_string(this->stop)+" slope:"+to_string(this->slope)+" intertectpt:"+to_string(this->intercept)+"\n";
//     string s2 = "node count:"+to_string(this->node_size)+" node range:"+to_string(this->nodes[0].key)+" "+to_string(this->nodes[this->node_size-1].key)+"\n";
//     if(w==1){
//         write_into_file(filePath,s1.c_str());
//         write_into_file(filePath,s2.c_str());
//     }
//     else if(w==0){
//         cerr<<s1;
//         cerr<<s2<<endl;
//     }
//     else if(w == 2){
//         string s3 ="Segment_pt level:"+to_string(this->level)+" start:"+to_string(this->start)+" stop:"+to_string(this->stop)+" node count:"+to_string(this->node_size)+"\n";
//         cerr<<s3;
//     }
// }
