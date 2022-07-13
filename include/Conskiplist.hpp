#include <iostream>
#include <algorithm>
#include <vector>
#include <limits.h>
#include <math.h>
#include <cfloat>
#include <fstream>
#include <iomanip>
#include <vector>
#include <sstream>
#include <map>
#include <stack>
#include <cstring>
#include <time.h>
#include <queue>
#include <atomic>
#include <assert.h>
#include <random>
#include <mutex>
#include <condition_variable>
#include <folly/MicroLock.h>
#include <shared_mutex>

#define PREFETCH(addr, rw, locality) __builtin_prefetch
typedef uint8_t bitmap_t;

typedef unsigned long long KeyType;
typedef int VaueType;
typedef long double ModelType; //long double
#define KeyMax (ULLONG_MAX)

#define BITMAP_WIDTH (sizeof(bitmap_t) * 8)
#define BITMAP_SIZE(num_items) (((num_items) + BITMAP_WIDTH - 1) / BITMAP_WIDTH)
#define BITMAP_GET(bitmap, pos) (((bitmap)[(pos) / BITMAP_WIDTH] >> ((pos) % BITMAP_WIDTH)) & 1)
#define BITMAP_SET(bitmap, pos) ((bitmap)[(pos) / BITMAP_WIDTH] |= 1 << ((pos) % BITMAP_WIDTH))
#define BITMAP_CLEAR(bitmap, pos) ((bitmap)[(pos) / BITMAP_WIDTH] &= ~bitmap_t(1 << ((pos) % BITMAP_WIDTH)))
#define BITMAP_NEXT_1(bitmap_item) __builtin_ctz((bitmap_item))

#define p 0.5
#define UNINT_MAX 0xffffffff
#define SkiplistMaxLevel (int)(log(5882)/log(2))//20//(int)(log(1500)/log(2))
#define Epslion 1e-8
#define LINAERITY 0.8
#define USEPLR 0//now means not use plr's model to rebuild segment 
#define Gm 16
#define MAX_DEPTH 6
#define INIT_DEPTH 6//3
#define INSERT_ROUTE 0
#define SEGMENT_MAX_SIZE 1e5
#define WRITESEG 0
#define DELTA_INSERT 2e5
#define IsREAD_THREAD 1
#define StateType int
#define DEBUG_ASSERT 1
#define PLR_DATA_PREPROCESS 0 //control whether data of plr is preprocessed
#define USINGLinearity 0
#define MEMORY_FETCH 0

// runtime assert
#define RT_ASSERT(expr) \
{ \
    if (!(expr)) { \
        fprintf(stderr, "RT_ASSERT Error at %s:%d, `%s`\n", __FILE__, __LINE__, #expr); \
        exit(0); \
    } \
}
using namespace folly;
using namespace std;

int file_num = 0;
std::mt19937 rd(time(0));

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

enum GreedyState {
	Need2 = 0, Need1, Ready
};

enum ItemState {
    Empty = 0, Element, Subtree
};

typedef struct RCUStatus {
    std::atomic<uint64_t> status;
    std::atomic<bool> waiting;
}RCUStatus;

//--------------------------rcu operation--------------------------//
std::unique_ptr<RCUStatus[]> rcu_status;
uint32_t worker_n;
inline void rcu_init(uint32_t worker_n_){
    worker_n = worker_n_;
    rcu_status = std::make_unique<RCUStatus[]>(worker_n);
    for (uint32_t worker_i = 0; worker_i < worker_n; worker_i++) {
        rcu_status[worker_i].status = 0;
        rcu_status[worker_i].waiting = false;
    }
}

inline void rcu_progress(const uint32_t worker_id) {
    // rcu_status[worker_id].status++;
    rcu_status[worker_id].status++;
}

void rcu_barrier(const uint32_t worker_id){
    uint64_t prev_status[worker_n];
    for (uint32_t w_i = 0; w_i < worker_n; w_i++) {
        prev_status[w_i] = rcu_status[w_i].status;
    }
    for (uint32_t w_i = 0; w_i < worker_n; w_i++) {
        if(worker_id != w_i){
            while (rcu_status[w_i].status <= prev_status[w_i]){;}
        }
    }
}

template<class K, class V,class M>
class skiplist {
    public:
        class Point {
        public:
            M x;
            M y;

            Point();

            Point(M x0, M y0) {
                x = x0;
                y = y0;
            }

            Point* upper_bound(M gamma) {
                Point* res = new Point(this->x, this->y + gamma);
                return res;
            }

            Point* lower_bound(M gamma) {
                Point* res = new Point(this->x, this->y - gamma);
                return res;
            }
        };

        class Line {
        public:
            M slope;
            M intercept;

            Line();

            Line(Point* a, Point* b) {
                this->slope = (b->y - a->y) / (b->x - a->x);
                this->intercept = b->y - b->x * this->slope;
            }

            Point* Intersection(Line* b) {
                M x,y;
                M deta = this->slope - b->slope;
                x = (b->intercept - this->intercept)/deta;
                y = (this->slope * b->intercept - this->intercept*b->slope)/deta;
                Point* res = new Point(x, y);
                return res;
            }

            bool AboveAssert(Point* k) {
                return k->y > k->x * this->slope + this->intercept;
            }

            bool BeblowAssert(Point* k) {
                return k->y < k->x * this->slope + this->intercept;
            }

        };

        class Segment {
            public:
                K start;
                K stop;
                M slope;
                M intercept;
                Segment();
                Segment(K start, K stop, M slope, M intercept) {
                    this->start = start;
                    this->stop = stop;
                    this->slope = slope;
                    this->intercept = intercept;
                }
        };

        class GreedyPLR {
        public:
            GreedyState state = GreedyState::Need2;
            int gamma;
            Point* s0 = nullptr;
            Point* s1 = nullptr;
            Point* sint = nullptr;//point of intersection
            Point* s_last = nullptr;
            Line* rho_lower = nullptr;
            Line* rho_upper = nullptr;

            GreedyPLR();

            GreedyPLR(int ga):gamma(ga){}

            void setup(Point* s0, Point* s1)
            {
                this->s0 = s0;
                this->s1 = s1;

                this->rho_lower = new Line(s0->upper_bound(this->gamma), s1->lower_bound(this->gamma));

                this->rho_upper = new  Line(s0->lower_bound(this->gamma), s1->upper_bound(this->gamma));
                this->sint = nullptr;
                this->sint = this->rho_upper->Intersection(this->rho_lower);
                if (this->sint == nullptr) {
                    cerr << "there is no intersection between upper line and lower line " << endl;
                    RT_ASSERT(this->sint != nullptr);
                }
            };

            Segment* CurrentSegment(M end) {
                if (this->state != GreedyState::Ready) {
                    return nullptr;
                }
                Segment* res = nullptr;
                M s_start = this->s0->x;
                M s_stop = end;
                M s_slope;
                s_slope = (this->rho_lower->slope + this->rho_upper->slope) / 2.0;
                M s_intercept = this->sint->y - this->sint->x * s_slope;
                res = new Segment(s_start, s_stop, s_slope, s_intercept);
                return res;
            }

            Segment* Process_pt(Point* k) {
                if (this->state != GreedyState::Ready) {
                    return nullptr;
                }
                if (!(this->rho_lower->AboveAssert(k) && this->rho_upper->BeblowAssert(k))) {
                    //重新开一个段
                    Segment* current = this->CurrentSegment(this->s_last->x);//(k->x);
                    // delete this->s0;
                    k->y = 0;
                    k->x = 0;
                    // this->s0 = nullptr;
                    this->s0 = k;
                    // this->s_last = nullptr;
                    this->s_last = k;
                    this->state = GreedyState::Need1;

                    return current;
                }
                Point* s_u = k->upper_bound(this->gamma);
                Point* s_l = k->lower_bound(this->gamma);
                if (this->rho_upper->BeblowAssert(s_u)) {
                    delete this->rho_upper;
                    // this->rho_upper = nullptr;
                    this->rho_upper = new Line(this->sint, s_u);
                }

                if (this->rho_lower->AboveAssert(s_l)) {
                    delete this->rho_lower;
                    // this->rho_lower = nullptr;
                    this->rho_lower = new Line(this->sint, s_l);
                }
                return nullptr;

            }

            Segment* Process(M x, M y) {
                //delete this->s_last;
                Segment* res = nullptr;
                Point* newp = new Point(x, y);
                // GreedyState newState = this->state;
                if (this->state == GreedyState::Need2) {
                    this->s0 = nullptr;
                    this->s_last = nullptr;
                    newp->y = 0;
                    this->s_last = newp;
                    this->s0 = newp;
                    // newState = GreedyState::Need1;
                    this->state = GreedyState::Need1;
                }
                else if (this->state == GreedyState::Need1) {
                    //delete this->s1;
                    this->s1 = nullptr;
                    this->s_last = nullptr;
                    this->s_last = newp;
                    this->s1 = newp;
                    this->setup(this->s0, this->s1);
                    // newState = GreedyState::Ready;
                    this->state = GreedyState::Ready;
                }
                else if (this->state == GreedyState::Ready) {
                    res = this->Process_pt(newp);
                    if (res != nullptr) {
                        // newState = GreedyState::Need1;
                        this->state = GreedyState::Need1;
                    }
                    else {
                        // newState = GreedyState::Ready;
                        this->state = GreedyState::Ready;
                    }
                }
                // this->state = newState;
                return res;
            }

            Segment* finish() {
                if (this->state == GreedyState::Need2) {
                    return nullptr;
                }
                else if (this->state == GreedyState::Need1) {
                    Segment* curnt = new Segment(this->s0->x, KeyMax, 0.0, this->s0->y);
                    // cerr <<"finish slopt:"<< curnt->slope << " "<<endl;
                    return curnt;
                }
                else if (this->state == GreedyState::Ready) {
                    return this->CurrentSegment(KeyMax);
                }
                return nullptr;
            }
        };

        class Semaphore;
        typedef struct spinlock {
            std::atomic<bool> lock_ = {0};

            void lock() noexcept {
                for (;;) {
                // Optimistically assume the lock is free on the first try
                if (!lock_.exchange(true, std::memory_order_acquire)) {
                    return;
                }
                // Wait for lock to be released without generating cache misses
                while (lock_.load(std::memory_order_acquire)) {
                    // Issue X86 PAUSE or ARM YIELD instruction to reduce contention between
                    // hyper-threads
                    __builtin_ia32_pause();
                }
                }
            }

            bool try_lock() noexcept {
                // First do a relaxed load to check if lock is free in order to prevent
                // unnecessary cache misses if someone does while(!try_lock())
                return !lock_.load(std::memory_order_acquire) &&
                    !lock_.exchange(true, std::memory_order_acquire);
            }

            void unlock() noexcept {
                lock_.store(false, std::memory_order_release);
            }

        }slock;
        struct node{
            K key;
            V value;
        };
        struct SNode{            
            SNode* Next() {
                // Use an 'acquire load' so that we observe a fully initialized
                // version of the returned Node.
                return (next_.load(std::memory_order_acquire));
            }

            void SetNext(SNode* x) {
                next_.store(x, std::memory_order_release);
            }

            bool CASNext(SNode* expected, SNode* x) {
                return next_.compare_exchange_strong(expected, x);
            }

            bool isIndexNode(){
                return isIndex;
            }

            K Key;//the lower bound
            bool isIndex; //read only, true-> index node, false -> sgement node
            std::atomic<SNode*> next_; // next 指针用于构建链表，直接继承 
            SNode(K base,bool type):Key(base),isIndex(type),next_(nullptr){}  
        };
        struct inode{
            K key;
            SNode *address;
        };      
        static bool cmp(inode a,inode b){
            return a.key <= b.key;
        }
        struct inodeArray{
            std::allocator<inode> inodes_allocator;
            inode* new_inodes(int n){
                inode* x = inodes_allocator.allocate(n);
                if(x!=nullptr){
                    return x;
                }
                std::cout<<"allocate failed"<<std::endl;
                return nullptr;
            }
            void delete_inodes(inode* x, int n){
                inodes_allocator.deallocate(x, n);
            }

            inline int predict(K key){
                M v = (slope) * (static_cast<long double>(key)) + (intercept);
                if(v > std::numeric_limits<int>::max() / 2) {
                    return size - 1;
                }
                if(v < 0 || static_cast<int>(v)<0){
                    return 0;
                }
                return std::min(size - 1, static_cast<int>(v));
            }

            M slope = 0,intercept = 0;
            int size;
            inode* iArray;
            inodeArray(int size_, M slope_,M intercept_,inode *iArray_) :slope(slope_), intercept(intercept_),size(size_), iArray(iArray_){}
            inodeArray(int size_, M slope_,M intercept_):slope(slope_), intercept(intercept_),size(size_){
                iArray = nullptr;
                iArray = new_inodes(size_);
            }
            inodeArray(M slope_,M intercept_,K key,SNode *value):slope(slope_), intercept(intercept_),size(1){
                iArray = nullptr;
                iArray = new_inodes(1);
                iArray[0] = {key,value};
            }
            ~inodeArray(){
                if(iArray){
                    delete_inodes(iArray,size);
                    iArray = nullptr;
                }
            }
        };
        struct subtree;
        struct Item{
            union {
                node data;
                subtree* child;
            }comp;
        };
        struct subtree{
            //--------------------------subtree operation--------------------------//
            inline int predict(K key){
                M v = (static_cast<long double>(this->slope)) * (static_cast<long double>(key)) + (static_cast<long double>(this->intercept));
                if(v > std::numeric_limits<int>::max() / 2) {
                    return this->num_items - 1;
                }
                if(v < 0 || static_cast<int>(v)<0){
                    return 0;
                }
                return std::min(this->num_items - 1, static_cast<int>(v));
            }

            std::pair<bool,V> find_key_in_subtree(K key){
                subtree* n = this;
                int pos = n->predict(key);
                while(!n->TryLockItem(pos)){
                    // printf("try lock %lld 's %d item failed\n",n,pos);
                }
                // n->LockItem(pos);
                while(1){
                    if (BITMAP_GET(n->none_bitmap, pos) == 1) {
                        //checkTypeEqualNull
                        n->ReleaseItem(pos);
                        break;
                    }
                    else if(BITMAP_GET(n->child_bitmap, pos) == 0){
                        //checkTypeEqualData
                        bool res_1 = n->items[pos].comp.data.key == key;
                        V res_2 = n->items[pos].comp.data.value;
                        n->ReleaseItem(pos);
                        return {res_1,res_2};
                    }
                    else{
                        subtree* next_subtree = n->items[pos].comp.child;
                        n->ReleaseItem(pos);
                        int pos_next = next_subtree->predict(key);
                        while(!next_subtree->TryLockItem(pos_next)){
                            // printf("try lock %lld 's %d item failed\n",next_subtree,pos_next);
                        }
                        n = next_subtree;
                        pos = pos_next;
                    }
                }
                return {false,0};
            }

            std::pair<bool,V> GetKeyValueNoMutex(K key){
                subtree* n = this;
                int pos = n->predict(key);
                while(1){
                    if (BITMAP_GET(n->none_bitmap, pos) == 1) {
                        break;
                    }
                    else if(BITMAP_GET(n->child_bitmap, pos) == 0){
                        //checkTypeEqualData
                        bool res_1 = n->items[pos].comp.data.key == key;
                        V res_2 = n->items[pos].comp.data.value;
                        return {res_1,res_2};
                    }
                    else{
                        subtree* next_subtree = n->items[pos].comp.child;
                        int pos_next = next_subtree->predict(key);
                        n = next_subtree;
                        pos = pos_next;
                    }
                }
                return {false,0};
            }

            //TODO:层锁获取后面再看一下
            void find_key_in_subtree(K key1,K key2,std::vector<std::pair<K,V>> &result){
                //目前的做法是扫描某一层时逐渐将其下一层锁住，扫描本层结束释放本层的锁
                struct seg{
                    int pos_start,pos_end;
                    subtree *array; 
                };
                std::queue<seg> check_list;
                subtree *n = this;
                int pos_start_ = n->predict(key1);
                int pos_end_ = n->predict(key2);
                n->LockItem(pos_start_);
                check_list.push({pos_start_,pos_end_,n});
                while (!check_list.empty()){
                    pos_start_ = check_list.front().pos_start;
                    pos_end_ = check_list.front().pos_end;
                    n = check_list.front().array;
                    check_list.pop();
                    for(int pos = pos_start_ ; pos <= pos_end_ ; pos++){
                        if (BITMAP_GET(n->none_bitmap, pos) == 0) {
                            if(BITMAP_GET(n->child_bitmap, pos) == 0){
                                K key_p = n->items[pos].comp.data.key;
                                if(key_p >= key1 && key_p < key2){
                                    result.push_back({key_p,n->items[pos].comp.data.value});
                                }
                            }else{
                                subtree *next = n->items[pos].comp.child;
                                int p_st = next->predict(key1);
                                int p_ed = next->predict(key2);
                                next->LockItem(p_st);
                                check_list.push({p_st,p_ed,next});
                            }
                        }
                    }
                    n->ReleaseItem(pos_start_);
                }
                
            }

            inline bool TryLockItem(int n){
                // printf("try to lock %lld 's %d item\n",this,n);
                return this->mlock[n/BITMAP_WIDTH].try_lock();
            }

            //get the lock of n->items[pos]
            inline void LockItem(int n){
                this->mlock[n/BITMAP_WIDTH].lock();
                // printf("lock %lld 's %d item\n",this,n);
            }

            //release the lock of n->items[pos]
            inline void ReleaseItem(int n){
                this->mlock[n/BITMAP_WIDTH].unlock();
                // printf("unlock %lld 's %d item\n",this,n);
            }

            long long SpaceSize(){
                long long res = sizeof(is_two) + sizeof(size) + sizeof(num_items) + 
                    sizeof(M)*2 + sizeof(bitmap_t*)*2;
                long long child_cnt = 0;
                const int bitmap_size_ = BITMAP_SIZE(num_items);
                res += (sizeof(bitmap_t) * bitmap_size_ * 2);
                // res += (sizeof(Item) * num_items);
                for(int i = 0;i<num_items;i++){
                    if(BITMAP_GET(none_bitmap,i) == 0){
                        if(BITMAP_GET(child_bitmap,i) == 1){
                            res += ((items[i].comp.child)->SpaceSize());
                            child_cnt++;
                        }
                    }
                }
                return res;
            }

            bool is_two = 0; // is special node for only two keys
            std::atomic<int> size;
            int num_items = 0; // size of items
            M slope = 0,intercept = 0;
            // K start = KeyMax,stop=KeyMax;
            Item* items = nullptr;
            slock *mlock = nullptr;
            bitmap_t* none_bitmap = nullptr; // 1 means None, 0 means Data or Child
            bitmap_t* child_bitmap = nullptr; // 1 means Child. will always be 0 when none_bitmap is 1
        };
        struct Segment_pt: SNode{
            //--------------------------allocator---------------------------//
            std::allocator<subtree> subtree_allocator;
            subtree* new_subtree(int n){
                subtree* x = subtree_allocator.allocate(n);
                if(x!=nullptr){
                    return x;
                }
                std::cout<<"allocate failed"<<std::endl;
                return nullptr;
            }
            void delete_subtree(subtree* x,int n){
                subtree_allocator.deallocate(x,n);
            }
            void destroy_subtree(subtree* root){
                std::stack<subtree*> s;
                s.push(root);
                while (!s.empty()) {
                    subtree* n = s.top(); s.pop();
                    for (int i = 0; i < n->num_items; i ++) {
                        if (BITMAP_GET(n->child_bitmap, i) == 1) {
                            s.push(n->items[i].comp.child);
                        }
                    }
                    delete_items(n->items, n->num_items);
                    const int bitmap_size = BITMAP_SIZE(n->num_items);
                    delete_bitmap(n->none_bitmap, bitmap_size);
                    delete_bitmap(n->child_bitmap, bitmap_size);
                    delete_slock(n->mlock,bitmap_size);
                    delete_subtree(n, 1);
                }
            }

            std::allocator<Item> item_allocator;
            Item* new_items(int n)
            {
                Item* x = item_allocator.allocate(n);
                if(x!=nullptr){
                    return x;
                }
                std::cout<<"allocate failed"<<std::endl;
                return nullptr;
            }
            void delete_items(Item* x, int n)
            {
                item_allocator.deallocate(x, n);
            }

            std::allocator<bitmap_t> bitmap_allocator;
            bitmap_t* new_bitmap(int n){
                bitmap_t* x = bitmap_allocator.allocate(n);
                if(x!=nullptr){
                    return x;
                }
                std::cout<<"allocate failed"<<std::endl;
                return nullptr;
            }
            void delete_bitmap(bitmap_t* x, int n){
                bitmap_allocator.deallocate(x, n);
            }
                       
            std::allocator<slock> slock_allocator;
            slock* new_slock(int n){
                slock *x = slock_allocator.allocate(n);
                if(x!=nullptr){
                    return x;
                }
                std::cout<<"allocate failed"<<std::endl;
                return nullptr;
            }
            void delete_slock(slock* x,int n){
                slock_allocator.deallocate(x,n);
            }

            subtree* build_tree_none(K st,K ed){
                subtree* n = new_subtree(1);
                n->is_two = 0;
                n->size.store(0, std::memory_order_release);
                n->num_items = 1;
                n->slope = n->intercept = 0;
                n->items = new_items(1);
                memset(n->items,0,sizeof(Item));
                n->mlock = new_slock(1);
                n->TryLockItem(0);
                n->ReleaseItem(0);
                n->none_bitmap = new_bitmap(1);
                n->none_bitmap[0] = 0;
                BITMAP_SET(n->none_bitmap, 0);
                n->child_bitmap = new_bitmap(1);
                n->child_bitmap[0] = 0;
                return n;
            }

            subtree* build_tree_two_nokey(K start,K stop){
                subtree* n = NULL;
                // while(!CASPool()){;}
                // if(tree_pool.empty()){
                    n = new_subtree(1);
                    n->is_two = 0;
                    n->size.store(0, std::memory_order_release);
                    n->num_items = 512;
                    n->items = new_items(n->num_items);
                    memset(n->items,0,(n->num_items)*sizeof(Item));
                    
                    const int bitmap_size = BITMAP_SIZE(n->num_items);
                    n->mlock = new_slock(bitmap_size);
                    for(int i = 0;i<bitmap_size;i++){
                        n->TryLockItem(i*8);
                        n->ReleaseItem(i*8);
                    }
                    n->none_bitmap = new_bitmap(bitmap_size);
                    n->child_bitmap = new_bitmap(bitmap_size);
                    for(int i = 0;i<bitmap_size;i++){
                        n->none_bitmap[i] = 0xff;
                        n->child_bitmap[i] = 0;
                    }
                    
                // }else{
                //     n = tree_pool.top();
                //     n->size.store(0, std::memory_order_release);
                //     tree_pool.pop();
                // }
                // ReleasePool();
                
                const long double mid1_key = (static_cast<long double>(stop) - static_cast<long double>(start))/3.0;
                const long double mid2_key = mid1_key*2;

                const long double mid1_target = n->num_items / 3;
                const long double mid2_target = n->num_items * 2 / 3;

                n->slope = (mid2_target - mid1_target) / (mid2_key - mid1_key);
                #if DEBUG_ASSERT
                RT_ASSERT(isinf(n->slope)==0);
                #endif
                n->intercept = mid1_target - n->slope * mid1_key;
                #if DEBUG_ASSERT
                RT_ASSERT(isfinite(n->slope));
                RT_ASSERT(isfinite(n->intercept));
                RT_ASSERT(n!=NULL);
                #endif
                return n;
            }

            subtree* build_tree_two(K key1, V value1, K key2, V value2){
                if (key1 > key2) {
                    std::swap(key1, key2);
                    std::swap(value1, value2);
                }
                subtree* n = NULL;
                while(!CASPool()){;}
                if(tree_pool.empty()){
                    n = new_subtree(1);
                    n->is_two = 1;
                    n->size.store(2, std::memory_order_release);
                    n->num_items = 8;
                    n->items = new_items(n->num_items);
                    memset(n->items,0,(n->num_items)*sizeof(Item));
                    n->mlock = new_slock(1);
                    n->TryLockItem(0);
                    n->ReleaseItem(0);
                    n->none_bitmap = new_bitmap(1);
                    n->child_bitmap = new_bitmap(1);
                    n->none_bitmap[0] = 0xff;
                    n->child_bitmap[0] = 0;
                }else{
                    n = tree_pool.top();
                    // n->TryLockItem(0);
                    // n->ReleaseItem(0);
                    tree_pool.pop();
                }
                ReleasePool();
                
                const long double mid1_key = static_cast<long double>(key1);
                const long double mid2_key = static_cast<long double>(key2);

                const long double mid1_target = n->num_items / 3;
                const long double mid2_target = n->num_items * 2 / 3;

                n->slope = (mid2_target - mid1_target) / (mid2_key - mid1_key);
                #if DEBUG_ASSERT
                if(isinf(n->slope) != 0){
                    std::cout<<"???"<<std::endl;
                }
                // RT_ASSERT(isinf(n->slope)==0);
                #endif
                n->intercept = mid1_target - n->slope * mid1_key;
                #if DEBUG_ASSERT
                RT_ASSERT(isfinite(n->slope));
                RT_ASSERT(isfinite(n->intercept));
                #endif
                { // insert key1&value1
                    int pos = n->predict(key1);
                    #if DEBUG_ASSERT
                    RT_ASSERT(BITMAP_GET(n->none_bitmap, pos) == 1);
                    #endif
                    BITMAP_CLEAR(n->none_bitmap, pos);
                    n->items[pos].comp.data.key = key1;
                    n->items[pos].comp.data.value = value1;
                }
                { // insert key2&value2
                    int pos = n->predict(key2);
                    #if DEBUG_ASSERT
                    RT_ASSERT(BITMAP_GET(n->none_bitmap, pos) == 1);
                    #endif
                    BITMAP_CLEAR(n->none_bitmap, pos);
                    n->items[pos].comp.data.key = key2;
                    n->items[pos].comp.data.value = value2;
                }
                #if DEBUG_ASSERT
                RT_ASSERT(n!=NULL);
                #endif
                return n;
            }

            inline int compute_gap_count(int size) {
                if (size >= 1000000) return 1;
                if (size >= 100000) return 2;
                return 5;
            }

            inline std::pair<bool,V> Lookup(K key){
                return DataArray->find_key_in_subtree(key);
            }

            inline void Lookup(K key1,K key2,std::vector<std::pair<K,V>> &result){
                this->DataArray->find_key_in_subtree(key1,key2,result);
            }
        
            void insert_subtree(K key,V value,int &path_size,subtree** route_,int &element_cnt,long long &collision){
                subtree* n = DataArray;
                int pos = n->predict(key);
                StateType state_raw;
                // n->mlock->lock();
                n->LockItem(pos);
                element_cnt = n->size.load(std::memory_order_acquire) + 1;
                state_raw = BITMAP_GET(n->none_bitmap, pos) == 1?0:BITMAP_GET(n->child_bitmap, pos)+1;
                // bool duplicate = false;
                while(1){
                    path_size++;
                    if(state_raw == ItemState::Empty){
                        BITMAP_CLEAR(n->none_bitmap, pos);
                        n->items[pos].comp.data.key = key;
                        n->items[pos].comp.data.value = value;
                        n->size.fetch_add(1, std::memory_order_acquire);
                        // n->mlock->unlock();
                        n->ReleaseItem(pos);
                        break;
                    }else if(state_raw == ItemState::Element){
                        if(n->items[pos].comp.data.key == key){
                            n->items[pos].comp.data.value = value;
                            // n->mlock->unlock();
                            n->ReleaseItem(pos);
                            for(int i = 0;i<path_size-1;i++){
                                route_[i]->size.fetch_sub(1,std::memory_order_acquire);
                            }
                            break;
                        }
                        subtree* next_subtree = nullptr;
                        next_subtree = build_tree_two(key, value,n->items[pos].comp.data.key, 
                        n->items[pos].comp.data.value);   
                        //replace
                        BITMAP_SET(n->child_bitmap, pos);
                        n->items[pos].comp.child = next_subtree;
                        n->size.fetch_add(1, std::memory_order_acquire);
                        path_size++;
                        // n->mlock->unlock();
                        n->ReleaseItem(pos);
                        collision++;
                        break;
                    }else{//ItemState::Subtree
                        n->size.fetch_add(1, std::memory_order_acquire);
                        RT_ASSERT(path_size < MAX_DEPTH*2);
                        route_[path_size-1] = n;
                        subtree* next_n = n->items[pos].comp.child;//n->items[pos].comp.child;
                        // n->mlock->unlock();
                        n->ReleaseItem(pos);
                        int next_pos = next_n->predict(key);
                        // next_n->mlock->lock();
                        next_n->LockItem(next_pos);
                        state_raw = BITMAP_GET(next_n->none_bitmap, next_pos) == 1?0:BITMAP_GET(next_n->child_bitmap, next_pos)+1;
                        n = next_n;
                        pos = next_pos;
                        collision++;
                    }
                }

            }
            
            void scan_and_destroy_subtree(subtree* root,K *keys,V *values,int ESIZE, bool destory = true){
                typedef std::pair<int, subtree*> Seg; // <begin, subtree*>
                std::stack<Seg> s;
                s.push(Seg(0, root));
                while (!s.empty()) {
                    int begin = s.top().first;
                    subtree* n = s.top().second;
                    s.pop();
                    for (int i = 0; i < n->num_items;i++) {
                        int none_byte = n->none_bitmap[i/BITMAP_WIDTH];
                        if(none_byte == 0xff){
                            i+=7;
                            continue;
                        }else if (((none_byte >> ((i) % BITMAP_WIDTH)) & 1) == 0) {
                            if (BITMAP_GET(n->child_bitmap, i) == 0) {
                                keys[begin] = n->items[i].comp.data.key;
                                values[begin] = n->items[i].comp.data.value;
                                begin ++;
                            } else {
                                Seg next_seg;
                                next_seg.first = begin;
                                next_seg.second = n->items[i].comp.child;
                                s.push(next_seg);
                                begin += n->items[i].comp.child->size.load(std::memory_order_acquire);
                            }
                        }
                    }
                    if (destory) {
                        if(n->is_two){
                            n->size.store(2, std::memory_order_release);
                            n->num_items = 8;
                            memset(n->items,0,(n->num_items)*sizeof(Item));
                            n->none_bitmap[0] = 0xff;
                            n->child_bitmap[0] = 0;
                            n->TryLockItem(0);
                            n->ReleaseItem(0);
                            while(!CASPool()){;}
                            tree_pool.push(n);
                            ReleasePool();
                        }else{
                            delete_items(n->items, n->num_items);
                            const int bitmap_size = BITMAP_SIZE(n->num_items);
                            delete_bitmap(n->none_bitmap, bitmap_size);
                            delete_bitmap(n->child_bitmap, bitmap_size);
                            delete_slock(n->mlock,bitmap_size);
                            delete_subtree(n, 1);
                        }
                    }
                }
            }

            long double ScanAndDestroySubtreeWithPCCs(subtree* root,K *keys,V *values,int ESIZE, bool destory = true){
                typedef std::pair<int, subtree*> Seg; // <begin, subtree*>
                std::stack<Seg> s;
                s.push(Seg(0, root));
                long double keylen = static_cast<long double>(ESIZE);
                long double sum_of_key = 0;//sigma(x-x0)
                long double sum_of_y = (keylen-1)* keylen/2.0;//simga(y)
                long double E_Y = ( keylen-1)/2.0;//E(Y)
                long double sum_of_keyMuly = 0;//sigma((x-x0)y)
                long double sum_of_key2 = 0;//sigma((x-x0)^2)
                long double sum_of_y2 = (keylen-1) * keylen * (2*keylen - 1)/6.0;//sigma(y^2)
                long double  key_0 = -1;//x0
                while (!s.empty()) {
                    int begin = s.top().first;
                    subtree* n = s.top().second;
                    s.pop();
                    for (int i = 0; i < n->num_items;i++) {
                        int none_byte = n->none_bitmap[i/BITMAP_WIDTH];
                        if(none_byte == 0xff){
                            i+=7;
                            continue;
                        }else if (((none_byte >> ((i) % BITMAP_WIDTH)) & 1) == 0) {
                            if (BITMAP_GET(n->child_bitmap, i) == 0) {
                                keys[begin] = n->items[i].comp.data.key;
                                values[begin] = n->items[i].comp.data.value;
                                if(key_0 == -1)
                                    key_0 = static_cast<long double>(keys[begin]);
                                long double xi = static_cast<long double>(keys[begin])-key_0;
                                sum_of_key+=(xi);
                                sum_of_keyMuly+=(xi*static_cast<long double>(begin));
                                sum_of_key2+=(xi*xi);
                                begin ++;
                            } else {
                                Seg next_seg;
                                next_seg.first = begin;
                                next_seg.second = n->items[i].comp.child;
                                s.push(next_seg);
                                begin += n->items[i].comp.child->size.load(std::memory_order_acquire);
                            }
                        }
                    }
                    if (destory) {
                        if(n->is_two){
                            n->size.store(2, std::memory_order_release);
                            n->num_items = 8;
                            memset(n->items,0,(n->num_items)*sizeof(Item));
                            n->none_bitmap[0] = 0xff;
                            n->child_bitmap[0] = 0;
                            n->TryLockItem(0);
                            n->ReleaseItem(0);
                            while(!CASPool()){;}
                            tree_pool.push(n);
                            ReleasePool();
                        }else{
                            delete_items(n->items, n->num_items);
                            const int bitmap_size = BITMAP_SIZE(n->num_items);
                            delete_bitmap(n->none_bitmap, bitmap_size);
                            delete_bitmap(n->child_bitmap, bitmap_size);
                            delete_slock(n->mlock,bitmap_size);
                            delete_subtree(n, 1);
                        }
                    }
                }

                long double E_key = sum_of_key/keylen;
                long double st_x = sum_of_key2 - sum_of_key*E_key;
                long double st_y = sum_of_y2 - sum_of_y*E_Y;
                if(st_x > 1e-6 && st_y > 1e-6 ){
                    long double PCCs = sum_of_keyMuly - sum_of_key*E_Y;
                    PCCs/=sqrt(st_x);
                    PCCs/=sqrt(st_y);
                    return PCCs;
                }
                return 1;
            }

            std::pair<subtree*,int> rebuild_tree(K *_keys,V * _values,int _size,bool plr = false,M top_slope = 0,M top_intercept = 0){
                //plr means whether to use the model from plr
                //top_slope , top_intercept are the model from the plr
                RT_ASSERT(_size > 1);
                typedef struct {
                        int begin;
                        int end;
                        int lvl; // top lvl = 1
                        subtree* n;
                } Seg;
                std::queue<Seg> s;
                subtree* ret = new_subtree(1);
                s.push({0, _size, 1, ret});
                int max_depth = 0;
                while(!s.empty()){
                    const int begin = s.front().begin;
                    const int end = s.front().end;
                    const int lvl = s.front().lvl;
                    max_depth = max(max_depth,lvl);
                    subtree* n = s.front().n;
                    s.pop();
                    #if DEBUG_ASSERT
                    RT_ASSERT(end - begin >= 2);
                    #endif
                    if(end-begin == 2){
                        //start and stop are set up in advance
                        subtree* _ = build_tree_two(_keys[begin], _values[begin], _keys[begin+1], _values[begin+1]);
                        memcpy(n, _, sizeof(subtree));//shalow copy
                        delete_subtree(_, 1);
                        if(lvl == 1){
                            ret = n;
                        }
                    }
                    else{
                        K * keys = _keys + begin;
                        V * values = _values + begin;
                        const int size = end - begin;
                        const int BUILD_GAP_CNT = compute_gap_count(size);
                        n->is_two = 0;
                        n->size.store(size, std::memory_order_release);
                        if(plr && lvl == 1){
                            n->num_items = size;
                            n->slope = top_slope;
                            n->intercept = top_intercept;
                            #if DEBUG_ASSERT
                            RT_ASSERT(isfinite(n->slope));
                            RT_ASSERT(isfinite(n->intercept));
                            #endif
                        }
                        else{
                            int mid1_pos = (size - 1) / 3;
                            int mid2_pos = (size - 1) * 2 / 3;
                            #if DEBUG_ASSERT
                            RT_ASSERT(0 <= mid1_pos);
                            RT_ASSERT(mid1_pos < mid2_pos);
                            RT_ASSERT(mid2_pos < size - 1);
                            #endif
                            const long double mid1_key =
                                    (static_cast<long double>(keys[mid1_pos]) + static_cast<long double>(keys[mid1_pos + 1]) ) / 2;
                            const long double mid2_key =
                                    (static_cast<long double>(keys[mid2_pos]) + static_cast<long double>(keys[mid2_pos + 1]) ) / 2;

                            n->num_items = size * static_cast<int>(BUILD_GAP_CNT + 1);
                            
                            const long double mid1_target = mid1_pos * static_cast<int>(BUILD_GAP_CNT + 1) + static_cast<int>(BUILD_GAP_CNT + 1) / 2;
                            const long double mid2_target = mid2_pos * static_cast<int>(BUILD_GAP_CNT + 1) + static_cast<int>(BUILD_GAP_CNT + 1) / 2;
                            #if DEBUG_ASSERT
                            RT_ASSERT(abs(mid2_key - mid1_key)>1e-8);
                            #endif
                            n->slope = (mid2_target - mid1_target) / (mid2_key - mid1_key);
                            n->intercept = mid1_target - n->slope * mid1_key;
                            #if DEBUG_ASSERT
                            RT_ASSERT(isfinite(n->slope));
                            RT_ASSERT(isfinite(n->intercept));
                            #endif
                        }        
                        #if DEBUG_ASSERT
                        RT_ASSERT(n->slope >= 0);
                        #endif
                        n->items = new_items(n->num_items);
                        memset(n->items,0,n->num_items*sizeof(Item));
                        
                        const int bitmap_size = BITMAP_SIZE(n->num_items);
                        n->mlock = new_slock(bitmap_size);
                        for(int i = 0;i<bitmap_size;i++){
                            n->TryLockItem(i*8);
                            n->ReleaseItem(i*8);
                        }
                        n->none_bitmap = new_bitmap(bitmap_size);
                        n->child_bitmap = new_bitmap(bitmap_size);
                        memset(n->none_bitmap, 0xff, sizeof(bitmap_t) * bitmap_size);
                        memset(n->child_bitmap, 0, sizeof(bitmap_t) * bitmap_size);
                        int item_i = n->predict(keys[0]);
                        for (int offset = 0; offset < size; ) {
                            int next = offset + 1, next_i = -1;
                            while (next < size) {     
                                next_i = n->predict( keys[next]);
                                if (next_i == item_i) {
                                    next ++;
                                } else {
                                    break;
                                }
                            }
                            #if DEBUG_ASSERT
                            RT_ASSERT(item_i < n->num_items);
                            #endif
                            if (next == offset + 1) {
                                #if DEBUG_ASSERT
                                RT_ASSERT(BITMAP_GET(n->none_bitmap,item_i) == 1);
                                RT_ASSERT(BITMAP_GET(n->child_bitmap,item_i) == 0);
                                #endif
                                BITMAP_CLEAR(n->none_bitmap, item_i);
                                n->items[item_i].comp.data.key = keys[offset];
                                n->items[item_i].comp.data.value = values[offset];
                            } else {
                                #if DEBUG_ASSERT
                                RT_ASSERT(BITMAP_GET(n->none_bitmap,item_i) == 1);
                                RT_ASSERT(BITMAP_GET(n->child_bitmap,item_i) == 0);
                                #endif
                                BITMAP_CLEAR(n->none_bitmap, item_i);
                                BITMAP_SET(n->child_bitmap, item_i);
                                n->items[item_i].comp.child = new_subtree(1);
                                s.push({begin + offset, begin + next, lvl + 1, n->items[item_i].comp.child});
                            }
                            if (next >= size) {
                                break;
                            } else {
                                item_i = next_i;
                                offset = next;
                            }
                        }
                    }   
                }
                return {ret,max_depth};
            }

            void CreateMemPool(int size_ = 200){
                for(int i = 0 ; i < size_ ; i++){
                    subtree* n = new_subtree(1);
                    n->is_two = 1;
                    n->size.store(2, std::memory_order_release);
                    n->num_items = 8;
                    n->items = new_items(n->num_items);
                    memset(n->items,0,(n->num_items)*sizeof(Item));
                    n->mlock = new_slock(1);
                    n->TryLockItem(0);
                    n->ReleaseItem(0);
                    n->none_bitmap = new_bitmap(1);
                    n->child_bitmap = new_bitmap(1);
                    n->none_bitmap[0] = 0xff;
                    n->child_bitmap[0] = 0;
                    n->intercept = 0;
                    n->slope = 0;
                    tree_pool.push(n);
                }
            }

            inline bool ReadPool(){
                return lock_pool.load(std::memory_order_acquire);
            }

            inline bool CASPool(bool free = false,bool block=true){
                return lock_pool.compare_exchange_strong(free, block);
            }

            inline void ReleasePool(){
                lock_pool.store(false, std::memory_order_release);
            }
            
            //reset segment max size when build model or split
            inline void SetSegMax(int slot_,int esize,K key_space){
                SegmentMaxSize = min((K)(SEGMENT_MAX_SIZE),key_space);
                // SegmentMaxSize = ComputeSegmentMaxSize(slot_,esize,key_space);
            }

            long long SpaceSize(){
                long long space_sum = sizeof(Segment_pt);
                //DataArray
                space_sum += DataArray->SpaceSize();
                //tree_pool
                return space_sum;
            }

            K bound;//right guard,the upper bound
            // int level;
            subtree *DataArray;
            Semaphore delta_inserts;//concurrency control
            int SegmentMaxSize;
            std::atomic<bool> lock_pool;
            std::stack<subtree*> tree_pool;
            
            Segment_pt(K base,K bnd,int deta_insert_init = INIT_DEPTH):SNode(base,false),bound(bnd),
            delta_inserts(Semaphore(deta_insert_init)),lock_pool(false){}
        };
        struct Index: SNode{
            inline void InitIArray(Index *root,K *keys_array,SNode **values_array,const int n_size,M slope_,M intercept_){
                root->SetDown(values_array[0]);
                
                inodeArray *iarray_ = new inodeArray(n_size,slope_,intercept_);
                inode *ary = iarray_->iArray;
                for(int i = 0;i<n_size;i++){
                    ary[i] = {keys_array[i],values_array[i]};
                }
                root->IArray.store((void*)iarray_,std::memory_order_release);
            }
            
            //get downward_
            SNode* Down() {
                return (downward_.load(std::memory_order_acquire));
                // RT_ASSERT(IArray);
                // return IArray[0].address;
            }

            //set downward_
            void SetDown(SNode* x) {
                downward_.store(x, std::memory_order_release);
            }

            inline void AppendNext(SNode* next_left,SNode* next_right){
                while (true){
                    SNode *pr_next = this->Next();//1.暂存当前Index的next
                    next_right->SetNext(pr_next);//2.将next的后继修改为pr_next
                    if(this->CASNext(pr_next,next_left)){//3.compare and set
                        return;
                    }
                    //4.false then loop
                }
            }

            // 这里node 有可能是index，也有可能是segment，统一这样处理了
            // inline void AppendBelow(SNode* below) {
            //     this->SetDown(below);
            // }

            inline bool CASBufferLock(bool free = false,bool block=true){
                return lock_buffer.compare_exchange_strong(free, block);
            }
            
            inline void ReleaseBufferLock(){
                lock_buffer.store(false, std::memory_order_release);
            }

            bool Find(const uint32_t worker_id,K key){
                rcu_progress(worker_id);
                inodeArray *array_address  = (inodeArray *)(IArray.load(std::memory_order_acquire));
                inode *ary = array_address->iArray;
                int size = array_address->size;
                if(!ary)
                    return false;
                if(buffer){
                    while(!CASBufferLock()){;}
                    for(auto it:*buffer){
                        if(it.key == key){
                            ReleaseBufferLock();
                            return true;
                        }
                    }
                    ReleaseBufferLock();
                } 
                int pos = array_address->predict(key);
                if(ary[pos].key == key){
                    return true;
                }
                int l = (pos-Gm + size)%size, r = (pos+Gm)%size;
                while (l <= r)
                {
                    int mid = l + ((r-l)>>1);
                    if(ary[mid].key == key){
                        return true;
                    }else if(ary[mid].key > key){
                        r = mid -1;
                    }else{
                        l = mid+1;
                    }
                }
                return false;

            }

            SNode* FindPrecursor(const uint32_t worker_id, K key){
                inodeArray *array_address  = (inodeArray *)(IArray.load(std::memory_order_acquire));
                inode *ary = array_address->iArray;
                int size = array_address->size;
                // RT_ASSERT(array_address);
                if(!ary)
                    return nullptr;
                PREFETCH(ary[0],0,1);
                int pos = array_address->predict(key);
                if(ary[pos].key == key){
                    RT_ASSERT(ary[pos].address);
                    return ary[pos].address;
                }
                //buffer中的不去读
                int l = max(0,pos - 2*Gm), r = min(size-1,pos + 2*Gm);
                SNode *res = nullptr;
                while (l <= r){
                    int mid = l + (r-l)/2;
                    K mid_key  = ary[mid].key;
                    if(mid_key == key){
                        return ary[mid].address;
                    }
                    else if(key < mid_key){
                        r = mid-1;
                    }else{
                        res = ary[mid].address;
                        l = mid + 1;
                    }
                }
                return res;
            }

            inline Index* FindParent(K key,SNode* pre,int level){
                SNode *next_ = pre->Next();
                while(next_ && key >= next_->Key){
                    pre = next_;
                    next_ = pre->Next();
                }
                return reinterpret_cast<Index*>(pre);
            }

            //假定new_index_node 的父节点就是本index node，那么Keys和new_index_node中的第一个元素不需要参与插入
            void Insert(const uint32_t worker_id,std::vector<K> &Keys,std::vector<SNode*> &new_index_node,SNode* preds[],int level,skiplist *list){
                while(!CASBufferLock()){;}
                int b_size = buffer->size(),n_size = new_index_node.size()-1;
                if(b_size + n_size < 10){
                    //lock_buffer
                    for(int i = 1;i<n_size;i++){
                        buffer->push_back({Keys[i],new_index_node[i]});
                    }
                    ReleaseBufferLock();
                    return;
                }else{
                    inodeArray *array_address  = (inodeArray *)(IArray.load(std::memory_order_acquire));
                    inode *ary = array_address->iArray;
                    int size = array_address->size;
                    const int new_size = b_size + size + n_size - 1;
                    K *keys_array = new K[new_size];
                    SNode* values_array[new_size];
                    inodeArray *old_array = array_address;
                    const int old_size = size;

                    MergeArray(ary,Keys,new_index_node,keys_array,values_array,size,b_size,n_size);
                    //using greedy plr split
                    std::vector<Index*> group;
                    int split_cnt = SplitArray(keys_array,values_array,new_size,group);//change this->IArray,this->buffer
                    //wait reader thread

                    rcu_barrier(worker_id);
                    delete old_array;

                    //insert into the upper layer index node
                    if(split_cnt > 1 && level < SkiplistMaxLevel){
                        K key_inode = reinterpret_cast<SNode*>(this)->Key;
                        Index *pre_inode = FindParent(key_inode,preds[level+1],level);//得到该层的最大前驱(snode.Key<=key_inode)
                        Index *parent_inode = pre_inode;
                        std::vector<K> split_inode_keys;
                        std::vector<SNode*> split_inodes_;
                        for(int i = 0;i<split_cnt;i++){
                            split_inode_keys.push_back(reinterpret_cast<SNode*>(group[i])->Key);
                            split_inodes_.push_back(reinterpret_cast<SNode*>(group[i]));
                        }
                        if(!pre_inode->Find(worker_id,key_inode)){//若最大前驱inode中无key_inode节点，则认为无父节点，需要手动创建
                            //create a new index
                            parent_inode =  new Index(key_inode,this);
                            pre_inode->AppendNext(parent_inode,parent_inode);
                            list->UpdateHeight(level+1);
                        }
                        ReleaseBufferLock();
                        parent_inode->Insert(worker_id,split_inode_keys,split_inodes_,preds,level+1,list);
                    }else{
                        ReleaseBufferLock();
                    }
                    delete []keys_array;
                    return;
                }
                // ReleaseBufferLock();
            }

            void Merge(inode *ary,int size,std::vector<inode> &middle_array,int b_size){
                int p1 = 0,p2 = 0,p_ = 0;
                while (p1 < size && p2 < b_size){
                    if(ary[p1].key <= buffer->at(p2).key){
                        middle_array[p_] = {ary[p1].key,ary[p1].address};
                        p1++;
                    }else{
                        middle_array[p_] = {buffer->at(p2).key,buffer->at(p2).address};
                        p2++;
                    }
                    p_++;
                }

                while(p1 < size){
                    middle_array[p_] = {ary[p1].key,ary[p1].address};
                    p1++;
                    p_++;
                }

                while(p2 < b_size){
                    middle_array[p_] = {buffer->at(p2).key,buffer->at(p2).address};
                    p2++;
                    p_++;
                }
                return;
            }

            int SplitArray(K *keys_array,SNode **values_array,const int n_size,std::vector<Index*> &group){
                GreedyPLR* plr = new GreedyPLR(Gm);
                vector<Segment*> seg;
                vector<int> segment_stIndex;
                int st = 0,ed = 0;
                for (int i = 0; i < n_size; i++) {
                    Segment* seg_res = nullptr;
                    seg_res = plr->Process(static_cast<M>(keys_array[i]), static_cast<M>(i-st));       
                    if(seg_res) {
                        segment_stIndex.push_back(st);
                        seg.push_back(seg_res);
                        st = ed = i;
                    }
                    else{
                        ed = i;
                    }
                }
                Segment* seg_res = nullptr;
                seg_res = plr->finish();
                if (seg_res) {
                    segment_stIndex.push_back(st);
                    seg.push_back(seg_res);
                }

                //create new index node
                int new_index_cnt = seg.size();
                RT_ASSERT(segment_stIndex.size() == new_index_cnt);
                if(new_index_cnt == 1){
                    delete plr;
                    InitIArray(this,keys_array,values_array,n_size,seg[0]->slope,seg[0]->intercept);//atomic change IArray
                    buffer->clear();
                    return 1;
                }
                group.resize(new_index_cnt);
                group[0] = this;
                for(int i = 1;i<new_index_cnt;i++){
                    K *keys_array_ = keys_array + segment_stIndex[i];
                    SNode **values_array_ = values_array + segment_stIndex[i];
                    group[i] = new Index(keys_array_[0]);
                    int i_size;
                    if(i == new_index_cnt -1 ){
                        i_size = n_size - segment_stIndex[i];
                    }else{
                        i_size = segment_stIndex[i+1] - segment_stIndex[i];
                    }
                    InitIArray(group[i],keys_array_,values_array_,i_size,seg[i]->slope,seg[i]->intercept);
                    group[i]->buffer = new std::vector<inode>();
                    //link 1...new_index_cnt-1 at current layer
                    if(i>1)
                        reinterpret_cast<SNode*>(group[i-1])->SetNext(reinterpret_cast<SNode*>(group[i]));
                }
                this->AppendNext(group[1],group[new_index_cnt-1]);
                //atomic change address
                InitIArray(this,keys_array,values_array,segment_stIndex[1],seg[0]->slope,seg[0]->intercept);
                buffer->clear(); 
                delete plr;      
                return new_index_cnt;     
            }

            void MergeArray(inode *ary,std::vector<K> &Keys,std::vector<SNode*> &new_index_node,K *keys_array,SNode **values_array,
                int size,int b_size,int n_size,bool including_first = false){
                sort(buffer->begin(),buffer->end(),cmp);
                std::vector<inode> middle_array(b_size + size);
                Merge(ary,size,middle_array,b_size);
                int p1 = 0,p2 = 0,p_ = 0;
                if(!including_first){
                    p2 = 1;
                }
                int merge_size = b_size + size;
                while (p1 < merge_size && p2 < n_size){
                    if(middle_array[p1].key <= Keys[p2]){
                        keys_array[p_] = middle_array[p1].key;
                        values_array[p_] = middle_array[p1].address;
                        p1++;
                    }else{
                        keys_array[p_] = Keys[p2];
                        values_array[p_] = new_index_node[p2];
                        p2++;
                    }
                    p_++;
                }

                while(p1 < merge_size){
                    keys_array[p_] = middle_array[p1].key;
                    values_array[p_] = middle_array[p1].address;
                    p1++;
                    p_++;
                }

                while(p2 < n_size){
                    keys_array[p_] = Keys[p2];
                    values_array[p_] = new_index_node[p2];
                    p2++;
                    p_++;
                }

            }

            void Show(){
                //start key
                inodeArray *iary = (inodeArray *)(IArray.load(std::memory_order_acquire));
                int size = iary->size;
                inode *ary = iary->iArray;
                cout<<"size:\t"<<size<<std::endl;
                cout<<ary[0].key;
                //键对数量，键对打出来
                for(int i = 1;i<size;i++){
                    cout<<","<<ary[i].key;
                }
                cout<<endl;
            }

            long long SpaceSize(){
                long long space_size = sizeof(Index) - sizeof(lock_buffer);
                //IArray
                void *t = IArray;
                inodeArray *array_s = (inodeArray*)t;
                space_size += array_s->SpaceSize();
                space_size += ((buffer->capacity()) * sizeof(inode));
                return space_size;
            }

            std::atomic<SNode*> downward_; // 下一层SNode,可能是Index, Segment_pt,
            std::atomic<void*> IArray;//inodeArray*,// inode *IArray;//compacted no gap
            std::atomic<bool> lock_buffer;
            std::vector<inode> *buffer;
            Index(K base,SNode *init_add = nullptr):SNode(base,true),downward_(nullptr),IArray(nullptr),lock_buffer(false),buffer(nullptr){
                if(init_add){
                    SetDown(init_add);
                    inodeArray *ary = new inodeArray(0,0,base,init_add);
                    IArray.store((void*)ary,std::memory_order_release);
                    buffer = new std::vector<inode>();
                }
            }
        };      
        class Semaphore {
        public:
            //wait_cnt(0),
            explicit Semaphore(int h_):state(false),max_height(h_),work_cnt(0){
            }
            ~Semaphore(){}
            //true means thread need to rebuild
            bool Signal(int h_,int size = 10,int id=0){
                std::unique_lock<std::mutex> lock(mutex_);

                if(h_ > max_height){
                    max_height = h_;
                }
                if(max_height > MAX_DEPTH || size > SEGMENT_MAX_SIZE){
                    state = true;//busy for rebuild
                }
                work_cnt--;
                if(work_cnt == 0 && state){
                    // std::cout<<max_height<<std::endl;
                    return true;
                }
                else{
                    return false;
                }
            }

            //用于rebuild/split结束 和 新建节点完成Index layer连接后
            void SignalForState(int h_=0,int esize=10,int id=0){
                std::unique_lock<std::mutex> lock(mutex_);
                state = false;
                if(h_)
                    max_height = h_;
                // if(esize > 5e4)
                    // depth_cap = MAX_DEPTH;
                state_cv.notify_all();
            }

            void Block(){
                std::unique_lock<std::mutex> lock(mutex_);
                state = true;
            }

            //false means skip,true means work_cnt++
            bool Wait(K key,Segment_pt *root,int id=0){
                std::unique_lock<std::mutex> lock(mutex_);
                while(state){
                    state_cv.wait(lock);
                }
                if(root->bound > key){
                    work_cnt++;
                    return true;
                }else{
                    return false;
                }                  
            }
        
        private:
            // std::mutex wait_mutex;
            // std::condition_variable wait_cv;
            // int wait_cnt;//记录当前正在等待的after insert线程数

            std::mutex mutex_;
            std::condition_variable state_cv;
            bool state;
            int max_height;
            int work_cnt;//Records the current number of threads that are working
            // int depth_cap;
        };
                
        //--------------------------new segment_pt--------------------------//

        SNode* NewSegment_pt(K base,K bound,subtree* n){
            //data array = n,if no data,n = nullptr
            if(!n)
                return NewSegment_pt(base,bound,false);
            SNode *newseg = nullptr;
            // int slots = (n->num_items);
            int ele_size = n->size.load(std::memory_order_acquire);
            if(ele_size < 2){
                newseg = new Segment_pt(base,bound);
            }
            else{
                newseg = new Segment_pt(base,bound,MAX_DEPTH);
            }
            reinterpret_cast<Segment_pt*>(newseg)->DataArray = n;
            reinterpret_cast<Segment_pt*>(newseg)->CreateMemPool();
            reinterpret_cast<Segment_pt*>(newseg)->SetNext(nullptr);
            // reinterpret_cast<Segment_pt*>(newseg)->SetSegMax(slots,ele_size,bound - base);
            return newseg;
        }
        
        SNode* NewSegment_pt(K base,K bound,bool ht){
            //ht:is head node or tail node
            SNode *newseg = new Segment_pt(base,bound,INIT_DEPTH);
            reinterpret_cast<Segment_pt*>(newseg)->CreateMemPool();
            //data array = build_tree_none
            if(ht)
                reinterpret_cast<Segment_pt*>(newseg)->DataArray =  reinterpret_cast<Segment_pt*>(newseg)->build_tree_none(base,bound);
            else
                reinterpret_cast<Segment_pt*>(newseg)->DataArray =  reinterpret_cast<Segment_pt*>(newseg)->build_tree_two_nokey(base,bound);
            reinterpret_cast<Segment_pt*>(newseg)->SetNext(nullptr);
            // reinterpret_cast<Segment_pt*>(newseg)->SetSegMax(8,0,bound - base);
            return newseg;
        }
        

        //--------------------------new index--------------------------//
        //create the multi layer of head/tail 
        SNode* NewIndexBoundArray(int level,K base,K bound,vector<SNode*> &SnodeArray){
            SNode* pr_ = NewSegment_pt(base,bound,true);//bottom segment
            SnodeArray[0] = pr_;
            for(int l = 1;l<=level;l++){
                Index *inode_n = new Index(base,pr_);
                SNode* curr = reinterpret_cast<SNode*>(inode_n);
                curr->SetNext(nullptr);//TODO to delete
                pr_ = curr;
                SnodeArray[l] = curr;
            }
            return pr_;//top index node
        }
        
        //--------------------------skiplist operation--------------------------//
        
        inline int GetMaxHeight() const {
            return max_height_.load(std::memory_order_acquire);
        }
        
        inline int RandLevel(){
            int lvl = 1;
            default_random_engine e(rd());
            uniform_int_distribution<int>u(1,100);
            for(int i = 0;i<MaxLevel-1;i++){
                if(u(e) > 50){
                    break;
                }
                lvl++;
            }
            return lvl;
        }

        inline void UpdateHeight(int new_height){
            int old_max_height = GetMaxHeight();
            while (new_height > old_max_height) {
                if (max_height_.compare_exchange_weak(old_max_height, new_height)) {
                    // successfully updated it
                    break;
                }
            }
        }

        //for Lookup
        //TODO
        SNode* Scan(const uint32_t worker_id,K key,SNode* preds[]){
            SNode* pred = head_;
            SNode* curr = nullptr;
            SNode* succ = nullptr;
            SNode* locate = nullptr;
            int top_level = GetMaxHeight();//Get the highest INDEX height at this time
            //down
            //get predecessor at [MaxLevel,top_level) layer 
            for(int l = MaxLevel;l>top_level;l--){
                preds[l] = pred;
                pred = reinterpret_cast<Index*>(pred)->Down();
            }

            //get the preds from [top_level,0)
            for(int l = top_level;l>=0;l--){
                //not sure key > curr's key
                curr = reinterpret_cast<Index*>(pred)->Next();
                #if DEBUG_ASSERT
                RT_ASSERT(curr!=nullptr);
                #endif
                succ = reinterpret_cast<Index*>(curr)->Next();
                while(succ && succ->Key <= key){//pred proceeds on, only when the left bound of succ is not greater than key
                    pred = curr;
                    curr = succ;
                    succ = reinterpret_cast<Index*>(succ)->Next();
                }
                //case1 pred <= key < curr ;
                //case2 pred < curr <= key < succ
                if(key >= curr->Key){
                    //key is bwtween curr and succ,case 2
                    pred = curr;//curr != tail_
                }
                preds[l] = pred;
                RT_ASSERT(preds[l]->Key >= preds[l+1]->Key);
                //TODO TEST
                //最大前驱查找需要修改
                if(l){
                    pred = reinterpret_cast<Index*>(pred)->FindPrecursor(worker_id,key);
                }else if(pred->Key == 0){
                    pred = reinterpret_cast<Index*>(pred)->Down();
                }
                // if(l)
                //     pred = reinterpret_cast<Index*>(pred)->Down();
                RT_ASSERT(pred);
            }
            //return segment_pt
            locate = pred;
            RT_ASSERT(reinterpret_cast<Segment_pt*>(locate)->bound > key);
            return locate;
        }

        SNode* Scan(const uint32_t worker_id,K key){
            SNode* pred = head_;
            SNode* curr = nullptr;
            SNode* succ = nullptr;
            SNode* locate = nullptr;
            int top_level = GetMaxHeight();
            //down
            for(int l = MaxLevel;l>top_level;l--){
                pred = reinterpret_cast<Index*>(pred)->Down();
            }
            //get the preds from [top_level,0)
            for(int l = top_level;l>=0;l--){
                //not sure key > curr's key
                curr = reinterpret_cast<Index*>(pred)->Next();
                #if DEBUG_ASSERT
                RT_ASSERT(curr!=nullptr);
                #endif
                succ = reinterpret_cast<Index*>(curr)->Next();   
                while(succ && succ->Key <= key){
                    pred = curr;
                    curr = succ;
                    succ = reinterpret_cast<Index*>(succ)->Next();
                }
                if(key >= curr->Key){
                    //key is bwtween curr and succ
                    pred = curr;
                }
                if(l){
                    pred = reinterpret_cast<Index*>(pred)->FindPrecursor(worker_id,key);
                }else if(pred->Key == 0){
                    pred = reinterpret_cast<Index*>(pred)->Down();
                }
            }
            //return segment_pt
            locate = pred;//reinterpret_cast<Seg*>(pred)->Down();
            return locate;
        }

        std::pair<bool,V> query(const uint32_t worker_id,K key){
            rcu_progress(worker_id);
            Segment_pt* locate = reinterpret_cast<Segment_pt*>(Scan(worker_id,key));
            // return locate->Lookup(key);
            // return locate->DataArray->GetKeyValueNoMutex(key);
            while(true){
                bool  noskip = locate->delta_inserts.Wait(key,locate);
                if(noskip){
                    std::pair<bool,V> res = locate->DataArray->GetKeyValueNoMutex(key);
                    bool rebuild = locate->delta_inserts.Signal(1);
                    if(rebuild){
                        locate->delta_inserts.SignalForState(1);
                    }
                    return res;
                }
                locate = reinterpret_cast<Segment_pt*>(locate->Next());
            }
            return {false,0};
        }

        //skiplist lookup a key
        std::pair<bool,V> Lookup(const uint32_t worker_id,K key){
            rcu_progress(worker_id);
            SNode* preds[MaxLevel+1];
            for(int i =0;i<=MaxLevel;i++){
                preds[i] = nullptr;
            }
            Segment_pt* locate = reinterpret_cast<Segment_pt*>(Scan(worker_id,key,preds));
            while(true){
                bool  noskip = locate->delta_inserts.Wait(key,locate);
                if(noskip){
                    std::pair<bool,V> res = locate->Lookup(key);
                    bool rebuild = locate->delta_inserts.Signal(1);
                    if(rebuild){
                        locate->delta_inserts.SignalForState(1);
                    }
                    return res;
                }
                locate = reinterpret_cast<Segment_pt*>(locate->Next());
            }
            return {false,0};
        }
        
        //return the corresponding values whose keys satisfy [key1,key2)
        void Lookup(const uint32_t worker_id,K key1,K key2,std::vector<std::pair<K,V>> &result){
            if(key2 < key1){
                return;
            }
            rcu_progress(worker_id);
            SNode* preds[MaxLevel+1];
            for(int i =0;i<=MaxLevel;i++){
                preds[i] = nullptr;
            }
            Segment_pt *locate = reinterpret_cast<Segment_pt*>(Scan(worker_id,key1,preds));
            while(true){
                bool  noskip = locate->delta_inserts.Wait(key1,locate);//范围查询不参与rebuild
                if(noskip){
                    locate->Lookup(key1,key2,result);
                    if(key2 < locate->bound){
                        //read compelete
                        bool rebuild = locate->delta_inserts.Signal(1);
                        if(rebuild)
                            locate->delta_inserts.SignalForState();
                        return;
                    }else{//traverse the next node
                        Segment_pt *next_seg = reinterpret_cast<Segment_pt*>(locate->Next());//void repeatedly reading
                        bool rebuild = locate->delta_inserts.Signal(1);
                        if(rebuild)
                            locate->delta_inserts.SignalForState();
                        locate = next_seg;
                        continue;
                    }
                }
                locate = reinterpret_cast<Segment_pt*>(locate->Next());
            }
            return;
        }

        //skiplist add a key 
        bool Add(const uint32_t worker_id,K key,V value,subtree** route_,int &path_depth,int &split_cnt,long long &collision){
            rcu_progress(worker_id);
            SNode* preds[MaxLevel+1];
            for(int i =0;i<=MaxLevel;i++){
                preds[i] = nullptr;
            }
            Segment_pt* locate = reinterpret_cast<Segment_pt*>(Scan(worker_id,key,preds));
            while(true){    
                if(key < locate->bound){//skip when the bound of locate segment is smaller than key because of splitting
                    bool noskip = locate->delta_inserts.Wait(key,locate);
                    if(noskip){
                        //work_cnt++ so won't rebuild/split
                        int after_insert_size = 0;
                        locate->insert_subtree(key,value,path_depth,route_,after_insert_size,collision);
                        bool rebuild = locate->delta_inserts.Signal(path_depth,after_insert_size);
                        if(rebuild){
                            //rebuild or split
                            #if USINGLinearity
                            int esize = RebuildSegmentWithLinearity(preds,locate,split_cnt);
                            #else
                            bool not_signal = RebuildSegment(worker_id,preds,locate,split_cnt);
                            #endif
                            if(not_signal)
                                locate->delta_inserts.SignalForState(INIT_DEPTH);
                            return true;
                        }
                        return false;
                    }
                }
                //skip 
                SNode* next = locate->Next();
                locate = reinterpret_cast<Segment_pt*>(next);
            }
            return false;	
        }

        //rebuild or split segment
        bool RebuildSegment(const uint32_t worker_id,SNode** preds,Segment_pt* locate,int &split_cnt){
            subtree *n = locate->DataArray;
            const int ESIZE = n->size.load(std::memory_order_acquire);
            K n_start = reinterpret_cast<SNode*>(locate)->Key,n_stop = locate->bound;
            //new_segment_max_size
            if(ESIZE>SEGMENT_MAX_SIZE){//split
                K* keys = new K[ESIZE];
                V* values = new V[ESIZE];
                locate->scan_and_destroy_subtree(n,keys,values,ESIZE);
                locate->DataArray = nullptr;
                std::vector<SNode*> split_nodes;
                std::vector<K> inode_keys;
                // std::cout<<"split "<<ESIZE<<std::endl;
                SplitSegment(locate,keys,values,ESIZE,n_start,n_stop,split_nodes,inode_keys);
                int new_size = locate->DataArray->size.load(std::memory_order_acquire);
                delete[] keys;
                delete[] values;
                split_cnt++;
                Index *pre_inode = reinterpret_cast<Index*>(preds[1])->FindParent(n_start,preds[1],1);
                Index *parent_inode = pre_inode;
                if(!pre_inode->Find(worker_id,n_start)){
                    //create a new index
                    parent_inode =  new Index(n_start,split_nodes[0]);
                    pre_inode->AppendNext(parent_inode,parent_inode);
                }
                locate->delta_inserts.SignalForState(INIT_DEPTH);
                parent_inode->Insert(worker_id,inode_keys,split_nodes,preds,1,this);
                return false;
            }else{
                K *keys = new K[ESIZE];
                V *values = new V[ESIZE];
                locate->scan_and_destroy_subtree(n,keys,values,ESIZE);
                locate->DataArray = nullptr;

                std::pair<subtree*,int> res_ = locate->rebuild_tree(keys,values,ESIZE);
                locate->DataArray = res_.first;
                // int new_slot = locate->DataArray->num_items;
                // locate->SetSegMax(new_slot,ESIZE,n_stop - n_start);
                delete[] keys;
                delete[] values;
                // new_segment_max_size = locate->SegmentMaxSize;
                return true;
            }
            return true;
        }

        int SplitSegment(Segment_pt *root,K *_keys,V *_values,int _size,K _start,K _stop,
            std::vector<SNode*> &split_nodes,std::vector<K> &inode_keys){
            //1. partition array
            int st = 0,ed =0;//the rank of the first/last+1 element of one segment
            GreedyPLR* plr = new GreedyPLR(gamma);
            vector<Segment*> seg;
            vector<int> segment_stIndex;
            for (int i = 0; i < _size; i++) {
                Segment* seg_res = nullptr;
                #if PLR_DATA_PREPROCESS
                seg_res = plr->Process(static_cast<M>(_keys[i]) - static_cast<M>(_keys[st]), static_cast<M>(i-st));
                #else
                // if(_stop - _start > 0x7fffffff)//TODO:需要随着key类型及数据范围调整
                    seg_res = plr->Process(static_cast<M>(_keys[i]), static_cast<M>(i-st));
                // else
                    // seg_res = plr->Process(static_cast<M>(_keys[i]) - static_cast<M>(_keys[st]), static_cast<M>(i-st));
                #endif
                                
                if(seg_res) {
                    segment_stIndex.push_back(st);
                    seg.push_back(seg_res);
                    st = ed = i;
                }
                else{
                    ed = i;
                }
            }
        
            Segment* seg_res = nullptr;
            seg_res = plr->finish();
            if (seg_res) {
                segment_stIndex.push_back(st);
                seg.push_back(seg_res);
            }
            
            //2. generate segments group
            int seg_size = static_cast<int>(seg.size());
            split_nodes.resize(seg_size);
            inode_keys.resize(seg_size);

            //create new segment
            split_nodes[0] = root;
            subtree *root_dataarray = nullptr;
            K split_lower_bound = 0;//the right bound of the first segment node after split
            for(int i = 0;i<seg_size;i++){
                int ed,st=segment_stIndex[i];
                K base,bound;
                //2.1 set base
                base = i==0?_start:_keys[st];
                inode_keys[i] = base;
                //2.2 set bound,ed
                if(i == seg_size-1){
                    bound = _stop;
                    ed = _size;
                }else{
                    ed = segment_stIndex[i+1];
                    bound = _keys[ed];
                }
                #if DEBUG_ASSERT
                RT_ASSERT(ed-st>0);
                #endif
                //2.3 rebuild subtree
                subtree* new_subtree = nullptr;
                if(ed-st == 1){
                    new_subtree = root->build_tree_none(base,bound);
                    BITMAP_CLEAR(new_subtree->none_bitmap,0);
                    new_subtree->size.store(1, std::memory_order_release);
                    new_subtree->items[0].comp.data.key = _keys[st];
                    new_subtree->items[0].comp.data.value = _values[st];
                }
                else{
                    std::pair<subtree*,int> res_;
                    #if USEPLR
                        if(seg[i]->slope > 0)
                            res_ = root->rebuild_tree(_keys+st,_values+st,ed-st,true,seg[i]->slope,seg[i]->intercept);
                        else
                            res_ = root->rebuild_tree(_keys+st,_values+st,ed-st);
                    #else
                        res_ = root->rebuild_tree(_keys+st,_values+st,ed-st);
                    #endif
                    new_subtree = res_.first;
                }

                //2.4 create segment 
                if(i){
                    //2.4.1 create index array(stored in SNodeArray) with subtree
                    SNode* next_seg = NewSegment_pt(base,bound,new_subtree);
                    split_nodes[i] = next_seg;
                    if(i>1){
                        split_nodes[i-1]->SetNext(next_seg);
                    }
                }else{
                    root_dataarray = new_subtree;
                    split_lower_bound = bound;
                }
                  
            }

            if(seg_size == 1){
                root->DataArray = root_dataarray;
                delete plr;
                return 1;
            }else{
                //2.6 link the bottom(segment_pt)
                split_nodes[seg_size-1]->SetNext(root->Next());
                root->SetNext(split_nodes[1]);
                //2.7 update root(segment_pt)
                root->bound = split_lower_bound;
                root->DataArray = root_dataarray;
                //TODO not sure whether the new node can split 
                //release the split locks of segment group 
                delete plr;
                return seg_size;
            }
        }

        /*bool SplitSegmentRaw(Segment_pt *root,SNode** preds,K *_keys,V *_values,int _size,K _start,K _stop){
            //1. partition array
            int st = 0,ed =0;//the rank of the first/last+1 element of one segment
            GreedyPLR* plr = new GreedyPLR(gamma);
            vector<Segment*> seg;
            vector<int> segment_stIndex;
            vector<int> height_seq;
            int segment_max_height = 1;//the max height of split segment
            for (int i = 0; i < _size; i++) {
                Segment* seg_res = nullptr;
                #if PLR_DATA_PREPROCESS
                seg_res = plr->Process(static_cast<M>(_keys[i]) - static_cast<M>(_keys[st]), static_cast<M>(i-st));
                #else
                // if(_stop - _start > 0x7fffffff)//TODO:需要随着key类型及数据范围调整
                    seg_res = plr->Process(static_cast<M>(_keys[i]), static_cast<M>(i-st));
                // else
                    // seg_res = plr->Process(static_cast<M>(_keys[i]) - static_cast<M>(_keys[st]), static_cast<M>(i-st));
                #endif
                                
                if(seg_res) {
                    segment_stIndex.push_back(st);
                    seg.push_back(seg_res);
                    if(st != 0){
                        int h = RandLevel();
                        segment_max_height = max(segment_max_height,h);
                        height_seq.push_back(h);
                    }
                    st = ed = i;
                }
                else{
                    ed = i;
                }
            }
        
            Segment* seg_res = nullptr;
            seg_res = plr->finish();
            if (seg_res) {
                segment_stIndex.push_back(st);
                seg.push_back(seg_res);
                if(st != 0){
                    int h = RandLevel();
                    height_seq.push_back(h);
                    segment_max_height = max(segment_max_height,h);
                }
            }
            
            //2. generate segments group
            int seg_size = static_cast<int>(seg.size());
            // std::cout<<"split segement number:"<<seg_size<<std::endl;
            //create new segment and index
            vector<Segment_pt*> split_segments(seg_size);//store the all segment_pt
            vector<SNode*> SNode_left(MaxLevel+1,nullptr);//segment group leftest snode
            vector<SNode*> SNode_right(MaxLevel+1,nullptr);//segment group rightest snode

            split_segments[0] = root;
            subtree *root_dataarray = nullptr;
            K split_lower_bound = 0;//the right bound of the first segment node after split
            for(int i = 0;i<seg_size;i++){
                int height_tmp=1,ed,st=segment_stIndex[i];
                K base,bound;
                //2.1 set height_tmp,base
                if(i == 0){
                    base = _start;
                }else{
                    height_tmp = height_seq[i-1];
                    base  = _keys[st];
                }
                //2.2 set bound,ed
                if(i == seg_size-1){
                    bound = _stop;
                    ed = _size;
                }else{
                    ed = segment_stIndex[i+1];
                    bound = _keys[ed];
                }
                #if DEBUG_ASSERT
                RT_ASSERT(ed-st>0);
                #endif
                //2.3 rebuild subtree
                subtree* new_subtree = nullptr;
                if(ed-st == 1){
                    new_subtree = root->build_tree_none(base,bound);
                    BITMAP_CLEAR(new_subtree->none_bitmap,0);
                    new_subtree->size.store(1, std::memory_order_release);
                    new_subtree->items[0].comp.data.key = _keys[st];
                    new_subtree->items[0].comp.data.value = _values[st];
                }
                else{
                    std::pair<subtree*,int> res_;
                    #if USEPLR
                        if(seg[i]->slope > 0)
                            res_ = root->rebuild_tree(_keys+st,_values+st,ed-st,true,seg[i]->slope,seg[i]->intercept);
                        else
                            res_ = root->rebuild_tree(_keys+st,_values+st,ed-st);
                    #else
                        res_ = root->rebuild_tree(_keys+st,_values+st,ed-st);
                    #endif
                    new_subtree = res_.first;
                }

                //2.4 create segment ,index without data and link inner part
                SNode* next_seg = nullptr;//the bottom segment_pt
                vector<SNode*> SNodeArray(height_tmp+1,nullptr);
                if(i){
                    //2.4.1 create index array(stored in SNodeArray) with subtree
                    NewIndexArrayWithData(height_tmp,base,bound,new_subtree,SNodeArray);
                    next_seg = SNodeArray[0];
                    split_segments[i] = reinterpret_cast<Segment_pt*>(next_seg);//store segment in split_segments
                    reinterpret_cast<Segment_pt*>(next_seg)->delta_inserts.Block();//restrain split
                    //2.5.1 link segment_pt if pre_segment exsits
                    if(SNode_right[0]){
                        reinterpret_cast<Segment_pt*>(SNode_right[0])->SetNext(next_seg);   
                    }
                    if(i == 1)
                        SNode_left[0] = next_seg;
                    SNode_right[0] = next_seg;
                    //2.5.2 link index in [height_tmp,0) layer if pre_segment exsits
                    for(int l = height_tmp;l>0;l--){
                        SNode* pr_ = SNodeArray[l];//get the l level index 
                        if(SNode_right[l]){
                            reinterpret_cast<Index*>(SNode_right[l])->SetNext(pr_);
                        }else{
                            SNode_left[l] = pr_;
                        }
                        SNode_right[l] = pr_;
                    }

                }else{
                    root_dataarray = new_subtree;
                    split_lower_bound = bound;
                    // root->SetSegMax(new_subtree->num_items,new_subtree->size.load(std::memory_order_acquire),bound - base);
                }
                  
            }

            // std::cout<<"split number:"<<seg_size<<std::endl;
            if(seg_size == 1){
                root->DataArray = root_dataarray;
                delete plr;
                return true;
            }else{
                //2.6 link the bottom(segment_pt)
                reinterpret_cast<Segment_pt*>(SNode_right[0])->SetNext(root->Next());//the right link
                root->SetNext(SNode_left[0]);//the left link
                //2.7 update root(segment_pt)
                root->bound = split_lower_bound;
                root->DataArray = root_dataarray;


                //2.8 update maxheight
                int max_height = GetMaxHeight();
                while (segment_max_height > max_height) {
                    if (max_height_.compare_exchange_weak(max_height, segment_max_height)) {
                        // successfully updated it
                        break;
                    }
                }

                //2.9 insert index from 1 to segment_max_height, reusing the information of preds array 
                SNode *pred = nullptr;
                SNode *curr = nullptr;
                SNode *succ = nullptr;
                for(int i = 1;i<=segment_max_height;i++){
                    pred = preds[i];
                    //move forward
                    curr = pred->Next();
                    succ = curr->Next();
                    while(succ && succ->Key <= split_lower_bound){
                        pred = curr;
                        curr = succ;
                        succ = succ->Next();
                    }
                    if(split_lower_bound >= curr->Key){
                        pred = curr;
                    }
                    if(SNode_left[i]){
                        reinterpret_cast<Index*>(pred)->AppendNext(SNode_left[i],SNode_right[i]);
                    }
                }
                //2.11 release the split locks of segment group 
                for(int i = 1;i<seg_size;i++){
                    split_segments[i]->delta_inserts.SignalForState(MAX_DEPTH);
                }
                delete plr;
                return true;
            }
            
            return true;
        }
        */

        void ShowIndexLayer(){
            int cur_height = GetMaxHeight();
            SNode *pr = head_;
            for(int i = MaxLevel;i>cur_height;i--){
                pr = reinterpret_cast<Index*>(pr)->Down();
            }
            SNode *cur = pr;
            for(int i = cur_height;i>0;i--){
                cout<<"layer "<<i<<":"<<std::endl;
                SNode *n = cur;
                while(n->Key != KeyMax){
                    reinterpret_cast<Index*>(n)->Show();
                    n = n->Next();
                }
                cur = reinterpret_cast<Index*>(cur)->Down();
            }
        }

        void ShowSegmentNumber(bool write_ = true){
            Segment_pt* head_seg = nullptr;
            SNode* pr_ = head_;
            for(int i = MaxLevel;i>0;i--){
                pr_ = reinterpret_cast<Index*>(pr_)->Down();
            }
            head_seg = reinterpret_cast<Segment_pt*>(pr_);
            int cnt = 0;
            SNode* cur = head_seg->Next();
            while(cur->Key != KeyMax){
                if(write_){
                    auto ESIZE = reinterpret_cast<Segment_pt*>(cur)->DataArray->size.load(std::memory_order_acquire);
                    string kk = to_string(cur->Key)+"\t"+to_string(reinterpret_cast<Segment_pt*>(cur)->bound)+"\t"+to_string(ESIZE)+"\t"+
                    to_string(ESIZE*1.0/(reinterpret_cast<Segment_pt*>(cur)->DataArray->num_items))+"\n";
                    write_into_file("./segment.csv",kk.c_str());
                }
                cnt++;
                cur = cur->Next();
            }
            std::cout<<"SegmentNumber:"<<cnt<<std::endl;
        }

        void FoldSegment(){
            SNode *x = head_;
            while(x->isIndexNode()){
                x = reinterpret_cast<Index*>(x)->Down();
            }
            x = x->Next();
            Segment_pt *n=nullptr;
            while(x->Key != KeyMax){
                n = reinterpret_cast<Segment_pt*>(x);
                subtree *seg = n->DataArray;
                const int ESIZE = seg->size.load(std::memory_order_acquire);
                if(ESIZE>1){
                    K *keys = new K[ESIZE];
                    V *values = new V[ESIZE];
                    n->scan_and_destroy_subtree(seg,keys,values,ESIZE,false);
                    n->DataArray = nullptr;
                    std::pair<subtree*,int> res_ = n->rebuild_tree(keys,values,ESIZE);
                    n->DataArray = res_.first;
                    delete[] keys;
                    delete[] values;
                }
                x = x->Next();
            }

        }

        long long SpaceSize(){
            long long space_sum  = 0;
            //skiplist class
            space_sum += sizeof(skiplist);
            SNode *x = head_;
            SNode *pr = x;
            for(int i = MaxLevel;i>0;i--){
                pr = x;
                while(pr){
                    space_sum += reinterpret_cast<Index*>(pr)->SpaceSize();
                    pr = pr->Next();
                }
                x = reinterpret_cast<Index*>(x)->Down();
            }
            cout<<"index layer size:"<<space_sum<<endl;
            pr = x;
            while (pr){
                space_sum += reinterpret_cast<Segment_pt*>(pr)->SpaceSize();
                pr = pr->Next();
            }
            
            return space_sum;
        }

        //--------------------------skiplist member variable--------------------------//
        std::atomic<int> max_height_; 
        int MaxLevel;
		int gamma;
        double linearity;
        SNode* head_;
        SNode* tail_;
        skiplist();
        skiplist(int MaxLevel,int gamma):max_height_(1),MaxLevel(MaxLevel),gamma(gamma),linearity(LINAERITY){
            /*
            头尾节点 各自有1-Maxlevel的index node，作指示作用
            中间节点，一开始不设置index node，仅有data layer,插入只会发生在split后
            */
            vector<SNode*> left_index(MaxLevel+1,nullptr);//store the all level SNode of head_
            vector<SNode*> right_index(MaxLevel+1,nullptr);//store the all level SNode of tail_
            head_ = NewIndexBoundArray(MaxLevel,0,0,left_index);
            tail_ = NewIndexBoundArray(MaxLevel,KeyMax,KeyMax,right_index);

            SNode *middle_node = NewSegment_pt(1,KeyMax,false);
            {//segment link
                reinterpret_cast<Segment_pt*>(left_index[0])->SetNext(middle_node);
                reinterpret_cast<Segment_pt*>(middle_node)->SetNext(right_index[0]);
            }
            {//index link
                for(int i = 1;i<=MaxLevel;i++){
                    reinterpret_cast<Index*>(left_index[i])->SetNext(right_index[i]);
                }
            }
        }
};