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

#define BITMAP_WIDTH (sizeof(bitmap_t) * 8)
#define BITMAP_SIZE(num_items) (((num_items) + BITMAP_WIDTH - 1) / BITMAP_WIDTH)
#define BITMAP_GET(bitmap, pos) (((bitmap)[(pos) / BITMAP_WIDTH] >> ((pos) % BITMAP_WIDTH)) & 1)
#define BITMAP_SET(bitmap, pos) ((bitmap)[(pos) / BITMAP_WIDTH] |= 1 << ((pos) % BITMAP_WIDTH))
#define BITMAP_CLEAR(bitmap, pos) ((bitmap)[(pos) / BITMAP_WIDTH] &= ~bitmap_t(1 << ((pos) % BITMAP_WIDTH)))
#define BITMAP_NEXT_1(bitmap_item) __builtin_ctz((bitmap_item))

#define p 0.5
#define UNINT_MAX 0xffffffff
#define Epslion 1e-8
#define LINAERITY 0.98
#define USEPLR 0//now means not use plr's model to rebuild segment 
#define USEGREEDYPCCS 0
#define Gm 64
#define MAX_DEPTH 6
#define INSERT_ROUTE 0
#define SEGMENT_MAX_SIZE 2e5
#define WRITESEG 0
#define DETA_INSERT 2e5
#define IsREAD_THREAD 1
#define StateType int
#define INCREASE_MODEL 3
#define W 4
#define Q 0.5
#define SLOT_MAX_SIZE 4096
#define NOREBUILD_th 0.02
#define PARTIALREBUILD_th 0.2

#define SHOWLINEAR 0
#define DEBUG 0
#define SEMA_DEBUG 0
#define MICROLOCK_DEBUG 0
#define TOPK_BUILD 0
#define OUTPUT_TOPK 0
#define SHOWBLOCK 0
#define SHOWREBUILD 0

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

class Point {
public:
	double x;
	double y;

	Point();

	Point(double x0, double y0) {
		x = x0;
		y = y0;
	}

	Point* upper_bound(double gamma) {
		Point* res = new Point(this->x, this->y + gamma);
		return res;
	}

	Point* lower_bound(double gamma) {
		Point* res = new Point(this->x, this->y - gamma);
		return res;
	}
};

class Line {
public:
	double slope;
	double intercept;

	Line();

	Line(Point* a, Point* b) {
		this->slope = (b->y - a->y) / (b->x - a->x);
		this->intercept = b->y - b->x * this->slope;
	}

	Point* Intersection(Line* b) {
		double x,y;
		double deta = this->slope - b->slope;
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
		unsigned int start;
		unsigned int stop;
		double slope;
		double intercept;
		Segment();
		Segment(unsigned int start, unsigned int stop, double slope, double intercept) {
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

	GreedyPLR(int ga) { this->gamma = ga; }

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
		}
	};

	Segment* CurrentSegment(double end) {
		if (this->state != GreedyState::Ready) {
			return nullptr;
		}
		Segment* res = nullptr;
		double s_start = this->s0->x;
		double s_stop = end;
		double s_slope;
		s_slope = (this->rho_lower->slope + this->rho_upper->slope) / 2.0;
		// if(s_slope < 0 && this->rho_lower->slope >= 0){
		// 	s_slope = this->rho_lower->slope;
		// }
		// RT_ASSERT(s_slope>=0);
		// cerr << "current slope:" << s_slope << endl;
		double s_intercept = this->sint->y - this->sint->x * s_slope;
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

	Segment* Process(double x, double y) {
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
			Segment* curnt = new Segment(this->s0->x, UNINT_MAX, 0.0, this->s0->y);
			// cerr <<"finish slopt:"<< curnt->slope << " "<<endl;
			return curnt;
		}
		else if (this->state == GreedyState::Ready) {
			return this->CurrentSegment(UNINT_MAX);
		}
		return nullptr;
	}
};

enum ItemState {
    Empty = 0, Element, Subtree
};

int ComputeNextCount(const int seg_max_size,int esize,int slot_){
    int remain_size = max(1,seg_max_size - esize);
    int count = 0;
    if(slot_ < 4096){
        count  = int(slot_*2);
    }
    else if(slot_ < 4e5){
        count = min(int(slot_*4),remain_size);
    }else{
        count = min(int(slot_*4),remain_size);
    }
    count = min(int(DETA_INSERT),count);
    return count;
}

int ComputeSegmentMaxSize(int slot_,int esize,unsigned int key_space){
    int res = 0; 
    res = min((unsigned int)(2e5),key_space);
    // if(esize < 5e4 && key_space > 0x3FFFFFF){// 
    //     res = min((unsigned int)(2e5),key_space);
    // }
    // else if(esize < 1e5){
    //     res = min((unsigned int)(2e5),key_space);
    // }
    // else{
    //     res = min((unsigned int)(1e6),key_space);
    // }
    return res;
}

template<class K, class V>
class skiplist {
    public:
        class Semaphore {
        public:
            explicit Semaphore(int count ,int upbound ) : count_(count),threshold(upbound),work_(0) {
            }
            ~Semaphore(){}
            //for delta_inserts V,esize=the current slot number
            void SignalAll(const int seg_max_size,int esize,int slot_=DETA_INSERT,int id=0){
                std::unique_lock<std::mutex> lock(mutex_);
                // std::cout<<"blocked:"<<count_<<std::endl;
                #if SHOWBLOCK
                std::cout<<"blocked:"<<count_<<std::endl;
                #endif 
                count_ = ComputeNextCount(seg_max_size,esize,slot_);
                #if SEMA_DEBUG
                std::cout<<"thread "<<id<<" notify all"<<std::endl;
                #endif
                // std::cout<<count_<<std::endl;
                cv_.notify_all();
            }
            //for delta_inserts V
            void SignalOne(int id=0){
                std::unique_lock<std::mutex> lock(mutex_);
                count_ = 1;
                #if SEMA_DEBUG
                std::cout<<"thread "<<id<<" notify one blocked thread"<<std::endl;
                #endif
                cv_.notify_one();
            }

            //for delta_inserts,当count_<0时，应该循环被挂起，返回值一定为非负数
            int Wait(int read_thread = 0,int id=0) {
                std::unique_lock<std::mutex> lock(mutex_);
                --count_;
                while(count_ < 0){
                     #if SEMA_DEBUG
                    std::cout<<"thread "<<id<<" ->wait"<<" state:"<<count_<<std::endl;
                    #endif
                    cv_.wait(lock);//只有收到notify信号才会结束wait
                     #if SEMA_DEBUG
                    std::cout<<"thread "<<id<<" ->wake up"<<std::endl;
                    #endif
                    --count_;
                }
                //从wait到--count，一定是notify all，所有相应的线程被唤醒，轮到本线程执行
                //--count后其他线程应该仍阻塞，因此当读线程执行++count，读线程只会被写线程阻塞，但不会使得写线程的count<0，令写线程在该信号量下被阻塞
                if(read_thread){// && count_ == 0
                    ++count_;
                    if(count_ == 1) 
                        cv_.notify_one();//读线程执行到此，一种是没进入while【当前没有blocked的线程，没必要notify_one】
                                        //一种可能是notify all的唤醒的【当前没有blocked的线程，没必要notify one】
                                        //一种是notify_one唤醒的【可能被state=0的误入写线程，可能被state=0的误入写线程唤醒的state=1的读线程唤醒】
                                        //由于无法区别以上情况，统一notify_one,因此可能会notify 空
                }
                #if SEMA_DEBUG
                std::cout<<"thread "<<id<<" state:"<<count_<<std::endl; 
                #endif
                std::unique_lock<std::mutex> lock_(mutex_w);
                ++work_;
                return count_; 
            }

            //for work threads
            void AddWork(int id=0){
                std::unique_lock<std::mutex> lock(mutex_w);
                ++work_;
                 #if SEMA_DEBUG
                std::cout<<"thread "<<id<<" work"<<std::endl;
                #endif
            }
            //for work threads V
            void SignalWork(int id=0){
                std::unique_lock<std::mutex> lock(mutex_w);
                --work_;
                 #if SEMA_DEBUG
                std::cout<<"thread "<<id<<" finish"<<std::endl;
                #endif
                if (work_ == 0){
                     #if SEMA_DEBUG
                    std::cout<<"thread "<<id<<" notify_one to work"<<std::endl;
                    #endif
                    cv_w.notify_one();
                }
            }
            //for work threads P
            void WaitWork(int id=0){
                std::unique_lock<std::mutex> lock(mutex_w);
                --work_;
                if(work_>0){
                     #if SEMA_DEBUG
                    std::cout<<"rebuild thread "<<id<<" wait working threads"<<std::endl;
                    #endif
                    // std::cout<<"wait to rebuild\n";
                    cv_w.wait(lock);
                }
                if(work_ == 0){
                    work_ = 1;
                }
                #if SEMA_DEBUG
                std::cout<<"rebuild thread "<<id<<" start rebuild"<<std::endl;
                #endif
            }
            
            //only test 
            void UpdateCount(int cnt){
                count_ = cnt;
            }

        private:
            std::mutex mutex_;
            std::condition_variable cv_;
            int count_;
            int threshold;
            std::mutex mutex_w;
            std::condition_variable cv_w;
            int work_;
        };
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

            K Key;//the lower bound
            bool isIndex; //只读变量，指示是否是index node, true-> index node, false sgement node
            std::atomic<SNode*> next_; // next 指针用于构建链表，直接继承 
            SNode(K base,bool type):Key(base),isIndex(type){}  
        };
        struct subtree;
        struct HashTable{
            std::pair<int,int> HashFunction(int pos){
                return {pos*293%capacity,pos*673%capacity};
            }

            void Update(int pos,int cnt){
                std::pair<int,int> predict_ = HashFunction(pos);
                if(index_[predict_.first] == pos){
                    count_[predict_.first] = cnt;
                }else if(index_[predict_.second] == pos){
                    count_[predict_.second] = cnt;
                }else{
                    int replace = count_[predict_.first] <= count_[predict_.second]?predict_.first:predict_.second;
                    if(count_[replace] < cnt){
                        index_[replace] = pos;
                        count_[replace] = cnt;
                    }
                }
            }

            std::pair<int,int> SumUP(){
                int sum = 0;
                int max_cnt = 0;
                for(int i = 0;i<capacity;i++){
                    if(index_[i] != -1){
                        sum+=(count_[i]);
                        max_cnt = max(max_cnt,count_[i]);
                    }
                }
                return {sum,max_cnt};
            }

            void Init(){
                for(int i = 0;i<capacity;i++){
                    index_[i] = -1;
                    count_[i] = 0;
                }
            }

            void Create(int capa){
                capacity = capa;
                index_ = new int[capacity];
                count_ = new int[capacity];
                for(int i = 0;i<capacity;i++){
                    index_[i] = -1;
                    count_[i] = 0;
                }
            }

            int capacity;//read only
            int *index_ = nullptr;//存放index 缺省值-1
            int *count_ = nullptr;//存放对应的元素个数，缺省值0
        };
        struct Item{
            union {
                node data;
                subtree* child;
            }comp;
        };
        struct subtree{
            //--------------------------subtree operation--------------------------//
            inline int predict(K key){
                double v = this->slope * (static_cast<long double>(key)) + this->intercept;//-static_cast<long double>(this->start)
                // if (std::isnan(v)){
                // 	return this->num_items - 1;
                // }
                if(v > std::numeric_limits<int>::max() / 2) {
                    return this->num_items - 1;
                }
                if(v < 0 || static_cast<int>(v)<0){
                    return 0;
                }
                return std::min(this->num_items - 1, static_cast<int>(v));
            }

            std::pair<unsigned int,unsigned int> computeRange(int pos){
                unsigned int start_ = (pos - this->intercept)/(this->slope * 1000.0)*1000.0;
                unsigned int stop_ = (pos+1 - this->intercept)/(this->slope * 1000.0)*1000.0;
                if(pos == 0)
                    start_ = this->start;
                if(pos == this->num_items-1)
                    stop_ = this->stop;
                return {start_,stop_};//+this->start
            }

            std::pair<bool,V> find_key_in_subtree(K key){
                subtree* n = this;
                vector<subtree*> paths;
                vector<int> poss;
                int pos = n->predict(key);
                n->LockItem(pos);
                while(1){
                    paths.push_back(n);
                    // int pos = n->predict(key);
                    poss.push_back(pos);
                    if (BITMAP_GET(n->none_bitmap, pos) == 1) {
                        //checkTypeEqualNull
                        std::cout<<"the pos is null "<<n->items[pos].comp.data.key<<std::endl;
                        n->ReleaseItem(pos);
                        break;
                        // return {false,0};
                    }
                    else if(BITMAP_GET(n->child_bitmap, pos) == 0){
                        //checkTypeEqualData
                        bool res_1 = n->items[pos].comp.data.key == key;
                        V res_2 = n->items[pos].comp.data.value;
                        n->ReleaseItem(pos);
                        return {res_1,res_2};
                        // return {n->items[pos].comp.data.key == key,n->items[pos].comp.data.value};
                    }
                    else{
                        subtree* next_subtree = n->items[pos].comp.child;
                        int pos_next = next_subtree->predict(key);
                        next_subtree->LockItem(pos_next);
                        n->ReleaseItem(pos);
                        // n = n->items[pos].comp.child;
                        n = next_subtree;
                        pos = pos_next;
                    }
                }
                return {false,0};
            }

            //for reader thread to get the lock of n->items[pos]
            inline void LockItem(int n){
                // while(this->mlock->try_lock() == false){;}
                this->mlock->lock();
                // this->mlock->lock();
            }

            //for reader thread to release the lock of n->items[pos]
            inline void ReleaseItem(int n){
                this->mlock->unlock();
            }
           
            //for writer thread to atomically get the state and value of n->items[pos], atomic through the microlock
            void GetStateItem(int n,int &state,Item &raw_item){
            //state 0:empty 1:element 2:subtree
                // while(this->mlock->try_lock() == false){;}
                this->mlock->lock();
                #if MICROLOCK_DEBUG
                std::cout<<this<<"\t"<<n<<"\tlock"<<std::endl;
                #endif
                state = BITMAP_GET(this->none_bitmap, n) == 1?0:BITMAP_GET(this->child_bitmap, n)+1;
                memcpy(&raw_item,this->items+n,sizeof(Item));
            }

            //for writer thread to atomically change the state and value of n->items[pos], atomic through the microlock
            bool CASStateItem(int n,int &raw_state,Item &raw_item,int new_state,Item &new_item){
                // if(this->mlock->try_lock()){
                    #if MICROLOCK_DEBUG
                    std::cout<<this<<"\t"<<n<<"\tlock"<<std::endl;
                    #endif
                    int real_state = BITMAP_GET(this->none_bitmap, n) == 1?0:BITMAP_GET(this->child_bitmap, n)+1;
                    if(raw_state != real_state){
                        #if MICROLOCK_DEBUG
                        std::cout<<this<<"\t"<<n<<"\tunlock"<<std::endl;
                        #endif
                        return false;
                    }
                    if(this->items[n].comp.child != raw_item.comp.child){
                        #if MICROLOCK_DEBUG
                        std::cout<<this<<"\t"<<n<<"\tunlock"<<std::endl;
                        #endif
                        return false;
                    }
                    memcpy(this->items+n,&new_item,sizeof(Item));
                    //new state can not be 0
                    if(new_state == ItemState::Subtree){
                        BITMAP_SET(this->child_bitmap, n);
                    }else{
                        BITMAP_CLEAR(this->none_bitmap, n);
                        // BITMAP_CLEAR(this->child_bitmap, n);
                    }
                    int check = BITMAP_GET(this->none_bitmap, n) == 1?0:BITMAP_GET(this->child_bitmap, n)+1;
                    this->mlock->unlock();
                    #if MICROLOCK_DEBUG
                    std::cout<<this<<"\t"<<n<<"\tunlock"<<std::endl;
                    #endif
                    return true;
                // }else{
                //     return false;
                // }
            }


            int is_two = 0; // is special node for only two keys
            int build_size = 0; // tree size (include sub nodes) when node created
            std::atomic<int> size;
            int fixed = 0;//fixed node will not trigger rebuild
            int confliction = 0;//the total insert confliction in this subtree
            int num_items = 0; // size of items
            double slope = 0,intercept = 0;
            K start = UNINT_MAX,stop=UNINT_MAX;
            Item* items = nullptr;
            std::mutex *mlock = nullptr;
            bitmap_t* none_bitmap = nullptr; // 1 means None, 0 means Data or Child
            bitmap_t* child_bitmap = nullptr; // 1 means Child. will always be 0 when none_bitmap is 1
            //HashTable TopK;//a hash table stored topk
        };
        struct Segment_pt: SNode{
            //--------------------------allocator---------------------------//
            std::allocator<subtree> subtree_allocator;
            subtree* new_subtree(int n){
                subtree* x = subtree_allocator.allocate(n);
                if(x!=nullptr){
                    return x;
                }
                cerr<<"allocate failed"<<endl;
                return nullptr;
            }
            void delete_subtree(subtree* x,int n){
                // delete[] x->TopK.count_;
                // delete[] x->TopK.index_;
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
                        delete_subtree(n, 1);
                }
            }

            std::allocator<Item> item_allocator;
            Item* new_items(int n)
            {
                Item* x = item_allocator.allocate(n);
                // RT_ASSERT(p != NULL && p != (Item*)(-1));
                return x;
            }
            void delete_items(Item* x, int n)
            {
                item_allocator.deallocate(x, n);
            }

            std::allocator<bitmap_t> bitmap_allocator;
            bitmap_t* new_bitmap(int n){
                bitmap_t* x = bitmap_allocator.allocate(n);
                // RT_ASSERT(x != NULL && x != (bitmap_t*)(-1));
                return x;
            }
            void delete_bitmap(bitmap_t* x, int n){
                bitmap_allocator.deallocate(x, n);
            }
            
            #if 0
            std::allocator<MicroLock> microlock_allocator;
            MicroLock* new_MicroLock(int n){
                MicroLock* x = microlock_allocator.allocate(n);
                return x;
            }
            void delete_MicroLock(MicroLock* x, int n){
                microlock_allocator.deallocate(x, n);
            }
            #endif  
           
            subtree* build_tree_none(){
                subtree* n = new_subtree(1);
                n->is_two = 0;
                n->build_size = 0;
                n->size.store(0, std::memory_order_release);
                n->fixed = 0;
                n->num_items = 1;
                n->confliction = 0;
                n->slope = n->intercept = 0;
                n->start = UNINT_MAX;
                n->stop = UNINT_MAX;
                n->items = new_items(1);
                memset(n->items,0,sizeof(Item));
                n->mlock = new mutex();//new_MicroLock(1);
                // n->mlock->try_lock();
                // n->mlock->unlock();
                n->none_bitmap = new_bitmap(1);
                n->none_bitmap[0] = 0;
                BITMAP_SET(n->none_bitmap, 0);
                n->child_bitmap = new_bitmap(1);
                n->child_bitmap[0] = 0;
                // n->TopK.Create(10);
                return n;
            }

            subtree* build_tree_none(K st,K ed){
                subtree* n = new_subtree(1);
                n->is_two = 0;
                n->build_size = 0;
                n->size.store(0, std::memory_order_release);
                n->fixed = 0;
                n->num_items = 1;
                n->confliction = 0;
                n->slope = n->intercept = 0;
                n->start = st;
                n->stop = ed;
                n->items = new_items(1);
                memset(n->items,0,sizeof(Item));
                n->mlock = new mutex();//new_MicroLock(1);
                // n->mlock->try_lock();
                // n->mlock->unlock();
                n->none_bitmap = new_bitmap(1);
                n->none_bitmap[0] = 0;
                BITMAP_SET(n->none_bitmap, 0);
                n->child_bitmap = new_bitmap(1);
                n->child_bitmap[0] = 0;
                // n->TopK.Create(10);
                return n;
            }

            subtree* build_tree_two(K key1, int value1, K key2, int value2,K start=0,K stop=0,subtree* x = nullptr){
                if (key1 > key2) {
                    std::swap(key1, key2);
                    std::swap(value1, value2);
                }
                // RT_ASSERT(key1 < key2);
                // static_assert(BITMAP_WIDTH == 8);
                subtree* n = NULL;
                while(!CASPool()){;}
                if(tree_pool.empty()){
                    n = new_subtree(1);
                    n->is_two = 1;
                    n->build_size = 2;
                    n->size.store(2, std::memory_order_release);
                    n->fixed = 0;
                    n->confliction = 0;
                    n->num_items = 8;
                    n->items = new_items(n->num_items);
                    memset(n->items,0,(n->num_items)*sizeof(Item));
                    n->mlock = new mutex();
                    n->none_bitmap = new_bitmap(1);
                    n->child_bitmap = new_bitmap(1);
                    n->none_bitmap[0] = 0xff;
                    n->child_bitmap[0] = 0;
                    // n->TopK.Create(10);
                }else{
                    n = tree_pool.top();
                    tree_pool.pop();
                }
                ReleasePool();

                if(!x){
                    n->start = start;
                    n->stop = stop;  
                }else{
                //从size=1增加元素重新生成subtree,[start,stop)
                    n->start = x->start;//max(x->start,2*key1-key2);
                    n->stop = x->stop;//min(x->stop,2*key2-key1);
                }
                
                const long double mid1_key = static_cast<long double>(key1);//-static_cast<long double>(n->start)
                const long double mid2_key = static_cast<long double>(key2);//-static_cast<long double>(n->start)

                const double mid1_target = n->num_items / 3;
                const double mid2_target = n->num_items * 2 / 3;

                n->slope = (mid2_target - mid1_target) / (mid2_key - mid1_key);
                // RT_ASSERT(isinf(n->slope)==0);
                n->intercept = mid1_target - n->slope * mid1_key;
                

                // RT_ASSERT(isfinite(n->slope));
                // RT_ASSERT(isfinite(n->intercept));

                { // insert key1&value1
                    int pos = n->predict(key1);
                    // RT_ASSERT(BITMAP_GET(n->none_bitmap, pos) == 1);
                    BITMAP_CLEAR(n->none_bitmap, pos);
                    n->items[pos].comp.data.key = key1;
                    n->items[pos].comp.data.value = value1;
                    // n->TopK.Update(pos,1);
                }
                { // insert key2&value2
                    int pos = n->predict(key2);
                    // RT_ASSERT(BITMAP_GET(n->none_bitmap, pos) == 1);
                    BITMAP_CLEAR(n->none_bitmap, pos);
                    n->items[pos].comp.data.key = key2;
                    n->items[pos].comp.data.value = value2;
                    // n->TopK.Update(pos,1);
                }
                RT_ASSERT(n!=NULL);

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
        
            void insert_subtree(K key,V value,subtree* path[],int &path_size,long long &collision){
                subtree* n = DataArray;//此时不会变动dataarray
                int pos = n->predict(key);
                StateType state_raw;
                n->mlock->lock();
                state_raw = BITMAP_GET(n->none_bitmap, pos) == 1?0:BITMAP_GET(n->child_bitmap, pos)+1;
                #if INSERT_ROUTE
                int path_capacity = MAX_DEPTH*2;
                #endif
                while(1){
                    #if INSERT_ROUTE
                    if(path_size == path_capacity){
                        path.reserve(2*path_capacity);
                        path.resize(2*path_capacity);
                        path_capacity*=2;
                    }
                    path[path_size] = n;
                    #endif
                    path_size++;
                    if(state_raw == ItemState::Empty){
                        BITMAP_CLEAR(n->none_bitmap, pos);
                        n->items[pos].comp.data.key = key;
                        n->items[pos].comp.data.value = value;
                        n->size.fetch_add(1, std::memory_order_acquire);
                        // n->TopK.Update(pos,1);
                        n->mlock->unlock();
                        break;
                    }else if(state_raw == ItemState::Element){
                        //none_tree由delta_insert控制写数量，因此当另一key被写入后，需要更新DataArray
                        if(path_size == 1 &&  n->size.load(std::memory_order_acquire) == 1){
                            subtree* newtree =  build_tree_two(key, value, n->items[pos].comp.data.key, n->items[pos].comp.data.value,0,0,n);
                            memcpy(DataArray, newtree, sizeof(subtree));//shallow copy,so just free newtree
                            delete_subtree(newtree, 1);
                            #if INSERT_ROUTE
                            path[0] = DataArray;
                            #endif
                            if(n->mlock->try_lock())
                                n->mlock->unlock();
                            break;
                        }else{
                            subtree* next_subtree = nullptr;
                            //compute
                            std::pair<unsigned int,unsigned int>line = n->computeRange(pos);
                            next_subtree = build_tree_two(key, value,n->items[pos].comp.data.key, 
                            n->items[pos].comp.data.value,min(line.first,min(key,n->items[pos].comp.data.key)),
                            line.second);   
                            //replace
                            BITMAP_SET(n->child_bitmap, pos);
                            n->items[pos].comp.child = next_subtree;
                            n->size.fetch_add(1, std::memory_order_acquire);
                            // n->TopK.Update(pos,2);
                            n->confliction++;
                            #if INSERT_ROUTE
                            if(path_size == path_capacity){
                                path.reserve(path_capacity+1);
                                path.resize(path_capacity+1);
                            }
                            path[path_size] = n->items[pos].comp.child;
                            #endif
                            path_size++;
                            n->mlock->unlock();
                            collision++;
                            break;
                        }
                    }else{//ItemState::Subtree
                        n->size.fetch_add(1, std::memory_order_acquire);
                        n->confliction++;
                        subtree* next_n = n->items[pos].comp.child;//n->items[pos].comp.child;
                        int next_pos = next_n->predict(key);
                        next_n->mlock->lock();
                        state_raw = BITMAP_GET(next_n->none_bitmap, next_pos) == 1?0:BITMAP_GET(next_n->child_bitmap, next_pos)+1;
                        // n->TopK.Update(pos,1+next_n->size.load(std::memory_order_acquire));
                        n->mlock->unlock();
                        n = next_n;
                        pos = next_pos;
                        collision++;
                    }
                }
            }
            
            std::pair<long double,int> scan_and_destroy_subtree(subtree* root,K *keys,V *values,int ESIZE, bool destory = true){
                typedef std::pair<int, subtree*> Seg; // <begin, subtree*>
                std::stack<Seg> s;
                s.push(Seg(0, root));
                int keylen = ESIZE;
                long double sum_of_key = 0;//sigma(x-x0)
                long double sum_of_y = (static_cast<long double>( keylen)-1)* keylen/2;//simga(y)
                long double E_Y = ( static_cast<long double>(keylen)-1)/2.0;//E(Y)
                long double sum_of_keyMuly = 0;//sigma((x-x0)y)
                long double sum_of_key2 = 0;//sigma((x-x0)^2)
                long double sum_of_y2 = (static_cast<long double>(keylen)-1) * keylen * (2*keylen - 1)/6;//sigma(y^2)
                long double  key_0 = -1;//x0
                int half_stone = root->num_items/2;
                int mark = -1;
                while (!s.empty()) {
                    int begin = s.top().first;
                    subtree* n = s.top().second;
                    // const int SHOULD_END_POS = begin + n->size.load(std::memory_order_acquire);
                    s.pop();
                    for (int i = 0; i < n->num_items;i++) {
                        if(mark == -1 && i>=half_stone){
                            mark = begin;
                        }
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
            #if SHOWLINEAR
                                string tmp = to_string(xi)+","+to_string(xi*static_cast<long double>(begin))+"\n";
                                string name = "./test"+to_string(file_num)+".csv";
                                write_into_file(name.c_str(),tmp.c_str());
            #endif
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
                    // RT_ASSERT(SHOULD_END_POS == begin);
                    if (destory) {
                        if(n->is_two){
                            n->build_size = 2;
                            n->size.store(2, std::memory_order_release);
                            n->fixed = 0;
                            n->confliction = 0;
                            n->num_items = 8;
                            memset(n->items,0,(n->num_items)*sizeof(Item));
                            n->none_bitmap[0] = 0xff;
                            n->child_bitmap[0] = 0;
                            // n->TopK.Init();
                            while(!CASPool()){;}
                            tree_pool.push(n);
                            ReleasePool();
                        }else{
                            delete_items(n->items, n->num_items);
                            const int bitmap_size = BITMAP_SIZE(n->num_items);
                            delete_bitmap(n->none_bitmap, bitmap_size);
                            delete_bitmap(n->child_bitmap, bitmap_size);
                            delete_subtree(n, 1);
                        }
                    }
                } 
                file_num++;
                long double E_key = sum_of_key/(static_cast<long double>(keylen));
                long double st_x = sum_of_key2 - sum_of_key*E_key;
                long double st_y = sum_of_y2 - sum_of_y*E_Y;
                if(st_x > 1e-6 && st_y > 1e-6 ){
                    long double PCCs = sum_of_keyMuly - sum_of_key*E_Y;
                    PCCs/=sqrt(st_x);
                    PCCs/=sqrt(st_y);
                    return {PCCs,mark};
                }
                return {1,mark};
            }

            std::pair<subtree*,int> rebuild_tree(K *_keys,int* _values,int _size,K start,K stop,bool plr = false,double top_slope = 0,double top_intercept = 0){
                RT_ASSERT(_size > 1);
                typedef struct {
                        int begin;
                        int end;
                        int lvl; // top lvl = 1
                        subtree* n;
                } Seg;
                std::queue<Seg> s;
                subtree* ret = new_subtree(1);
                // ret->child_ptr = 0;
                ret->start = start;//_keys[0];
                ret->stop = stop;//_keys[_size-1]+1;
                s.push({0, _size, 1, ret});
                int max_depth = 0;
                while(!s.empty()){
                    const int begin = s.front().begin;
                    const int end = s.front().end;
                    const int lvl = s.front().lvl;
                    //TODO:test
                    max_depth = max(max_depth,lvl);
                    subtree* n = s.front().n;
                    s.pop();
                    RT_ASSERT(end - begin >= 2);
                    if(end-begin == 2){
                        //n要不是root要不是子subtree，start和stop都事先确定了
                        subtree* _ = build_tree_two(_keys[begin], _values[begin], _keys[begin+1], _values[begin+1],n->start,n->stop);
                        memcpy(n, _, sizeof(subtree));//shalow copy
                        delete_subtree(_, 1);
                        if(lvl == 1){
                            ret = n;
                        }
                    }
                    else{
                        unsigned int* keys = _keys + begin;
                        int* values = _values + begin;
                        const int size = end - begin;
                        const int BUILD_GAP_CNT = compute_gap_count(size);

                        n->is_two = 0;
                        n->build_size = size;
                        n->size.store(size, std::memory_order_release);
                        n->fixed = 0;
                        n->confliction = 0;
                        if(plr && lvl == 1){
                            n->num_items = size;
                            n->slope = top_slope;
                            n->intercept = top_intercept;
                            RT_ASSERT(isfinite(n->slope));
                            RT_ASSERT(isfinite(n->intercept));
                        }
                        else{
                            int mid1_pos = (size - 1) / 3;
                            int mid2_pos = (size - 1) * 2 / 3;
                            if(mid2_pos >= size){
                                std::cout<<"out of bound mid2_pos\n";
                            }
                            RT_ASSERT(0 <= mid1_pos);
                            RT_ASSERT(mid1_pos < mid2_pos);
                            RT_ASSERT(mid2_pos < size - 1);
                            const long double mid1_key =
                                    (static_cast<long double>(keys[mid1_pos]) + static_cast<long double>(keys[mid1_pos + 1]) ) / 2;
                            const long double mid2_key =
                                    (static_cast<long double>(keys[mid2_pos]) + static_cast<long double>(keys[mid2_pos + 1]) ) / 2;
                            // const long double mid1_key =
                            //         (static_cast<long double>(keys[mid1_pos])-static_cast<long double>(n->start) + static_cast<long double>(keys[mid1_pos + 1])-static_cast<long double>(n->start) ) / 2;
                            // const long double mid2_key =
                            //         (static_cast<long double>(keys[mid2_pos]) -static_cast<long double>(n->start)+ static_cast<long double>(keys[mid2_pos + 1]) -static_cast<long double>(n->start)) / 2;

                            n->num_items = size * static_cast<int>(BUILD_GAP_CNT + 1);
                            
                            const double mid1_target = mid1_pos * static_cast<int>(BUILD_GAP_CNT + 1) + static_cast<int>(BUILD_GAP_CNT + 1) / 2;
                            const double mid2_target = mid2_pos * static_cast<int>(BUILD_GAP_CNT + 1) + static_cast<int>(BUILD_GAP_CNT + 1) / 2;
                                        RT_ASSERT(0 <= mid1_pos);
                            RT_ASSERT(abs(mid2_key - mid1_key)>1e-8);
                            n->slope = (mid2_target - mid1_target) / (mid2_key - mid1_key);
                            n->intercept = mid1_target - n->slope * mid1_key;
                            RT_ASSERT(isfinite(n->slope));
                            RT_ASSERT(isfinite(n->intercept));
                        }        
                        RT_ASSERT(n->slope >= 0);
                        // if(n->num_items >= SLOT_MAX_SIZE){
                        //     n->fixed = 1;
                        // }

                        n->items = new_items(n->num_items);
                        memset(n->items,0,n->num_items*sizeof(Item));
                        n->mlock = new mutex();
                        // n->TopK.Create(10);
                        const int bitmap_size = BITMAP_SIZE(n->num_items);
                        n->none_bitmap = new_bitmap(bitmap_size);
                        n->child_bitmap = new_bitmap(bitmap_size);
                        memset(n->none_bitmap, 0xff, sizeof(bitmap_t) * bitmap_size);
                        memset(n->child_bitmap, 0, sizeof(bitmap_t) * bitmap_size);
                        int real_put_ele = 0;
                        int item_i = n->predict(keys[0]);
                        for (int offset = 0; offset < size; ) {
                            //item_i 表示当前元素应插入的位置,next_i 表示下一个元素应插入的位置
                            int next = offset + 1, next_i = -1;
                            while (next < size) {     
                                next_i = n->predict( keys[next]);
                                if (next_i == item_i) {
                                    next ++;
                                } else {
                                    break;
                                }
                            }
                            RT_ASSERT(item_i < n->num_items);//"out of bound item_i\n";
                            if (next == offset + 1) {
                                real_put_ele++;
                                RT_ASSERT(BITMAP_GET(n->none_bitmap,item_i) == 1);
                                RT_ASSERT(BITMAP_GET(n->child_bitmap,item_i) == 0);
                                BITMAP_CLEAR(n->none_bitmap, item_i);
                                n->items[item_i].comp.data.key = keys[offset];
                                n->items[item_i].comp.data.value = values[offset];
                                // n->TopK.Update(item_i,1);
                            } else {
                                // ASSERT(next - offset <= (size+2) / 3);
                                RT_ASSERT(BITMAP_GET(n->none_bitmap,item_i) == 1);
                                RT_ASSERT(BITMAP_GET(n->child_bitmap,item_i) == 0);
                                BITMAP_CLEAR(n->none_bitmap, item_i);
                                BITMAP_SET(n->child_bitmap, item_i);
                                // n->child_ptr++;
                                n->items[item_i].comp.child = new_subtree(1);
                                std::pair<unsigned int,unsigned int> line = n->computeRange(item_i);
                                n->items[item_i].comp.child->start = min(line.first,keys[offset]);
                                n->items[item_i].comp.child->stop = line.second;
                                real_put_ele+=(next-offset);
                                n->confliction+=(next-offset-1);
                                s.push({begin + offset, begin + next, lvl + 1, n->items[item_i].comp.child});
                                // n->TopK.Update(item_i,next-offset);
                            }
                            if (next >= size) {
                                break;
                            } else {
                                item_i = next_i;
                                offset = next;
                            }
                        }
                        RT_ASSERT(real_put_ele == size);
                    }   
                }
                RT_ASSERT(ret->start == start);
                RT_ASSERT(ret->stop == stop);
                return {ret,max_depth};
            }

            void CreateMemPool(int size_ = 200){
                for(int i = 0 ; i < size_ ; i++){
                    subtree* n = new_subtree(1);
                    n->is_two = 1;
                    n->build_size = 2;
                    n->size.store(2, std::memory_order_release);
                    n->fixed = 0;
                    n->confliction = 0;
                    n->num_items = 8;
                    n->items = new_items(n->num_items);
                    memset(n->items,0,(n->num_items)*sizeof(Item));
                    n->mlock = new mutex();
                    // n->TopK.Create(10);
                    n->none_bitmap = new_bitmap(1);
                    n->child_bitmap = new_bitmap(1);
                    n->none_bitmap[0] = 0xff;
                    n->child_bitmap[0] = 0;
                    n->start = 0;
                    n->stop = 0;
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
            
            //第一次建模/每一次rebuild/split都重新设置segment max size
            void SetSegMax(int slot_,int esize,unsigned int key_space){
                SegmentMaxSize = ComputeSegmentMaxSize(slot_,esize,key_space);
                #if SHOWBLOCK
                std::cout<<"slot:"<<slot_<<"\tSegmentMaxSize:"<<SegmentMaxSize<<std::endl;
                #endif
                // std::cout<<"SegmentMaxSize:"<<SegmentMaxSize<<"\tslot:"<<slot_<<"\tkey_space"<<key_space<<std::endl;
            }

            K bound;//right guard,the upper bound
            int level;
            subtree *DataArray;
            Semaphore delta_inserts;
            int SegmentMaxSize;
            std::atomic<bool> lock_pool;
            std::stack<subtree*> tree_pool;
            
            Segment_pt(K base,K bnd,int lvl,int deta_insert_init = DETA_INSERT,int upbound_ = DETA_INSERT):SNode(base,false),bound(bnd),level(lvl),
            delta_inserts(Semaphore(deta_insert_init,upbound_)),lock_pool(false){}
        };
        struct Index: SNode{
            //get downward_
            SNode* Down() {
                // Use an 'acquire load' so that we observe a fully initialized
                // version of the returned Node.
                return (downward_.load(std::memory_order_acquire));
            }

            //set downward_
            void SetDown(SNode* x) {
                downward_.store(x, std::memory_order_release);
            }

            //comapre downward_ with expected, if equal,then replace downward_ with x
            //return true,otherwise return false
            bool CASDown(SNode* expected, SNode* x) {
                return downward_.compare_exchange_strong(expected, x);
            }

            // //链表添加后置元素
            // inline void AppendNext(SNode* next) {
            //    while (true){
            //         SNode *pr_next = this->Next();//1.暂存当前Index的next
            //         next->SetNext(pr_next);//2.将next的后继修改为pr_next
            //         if(this->CASNext(pr_next,next)){//3.compare and set
            //             return;
            //         }
            //         //4.false then loop
            //     }
            // }

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
            inline void AppendBelow(SNode* below) {
                this->SetDown(below);
            }

            std::atomic<SNode*> downward_; // 下一层SNode,可能是Index, Segment_pt
            // Index(K base){isIndex = true;Key = base;};  // 设置指示变量
            Index(K base):SNode(base,true){}
        };      
        typedef struct {
            int begin;
            int end;
            int lvl; // top lvl = 1
            subtree* n;
        } Seg;

        //--------------------------new segment_pt--------------------------//

        //data array = n,if no data,n = nullptr
        SNode* NewSegment_pt(K base,K bound,int level,subtree* n){
            if(!n)
                return NewSegment_pt(base,bound,level);
            SNode *newseg;
            //TODO:size or slot decide delta_insert
            //now choose slot
            int slots = (n->num_items);
            int ele_size = n->size.load(std::memory_order_relaxed);
            if(ele_size < 2){//because of the 
                newseg = new Segment_pt(base,bound,level,1,16);//if the number of ele in subtree is smaller than 2,
                                                            //then set deta_insert = 2,else set deta_insert = DETA_INSERT
            }
            else{
                int seg_max_size = ComputeSegmentMaxSize(slots,ele_size,bound - base);
                int count_set = ComputeNextCount(seg_max_size,ele_size,slots);
                newseg = new Segment_pt(base,bound,level,count_set,DETA_INSERT);
            }
            reinterpret_cast<Segment_pt*>(newseg)->DataArray = n;
            reinterpret_cast<Segment_pt*>(newseg)->CreateMemPool();
            reinterpret_cast<Segment_pt*>(newseg)->SetNext(nullptr);
            reinterpret_cast<Segment_pt*>(newseg)->SetSegMax(slots,ele_size,bound - base);
            return newseg;
        }
        
        //data array = build_tree_none
        SNode* NewSegment_pt(K base,K bound,int level){
            SNode *newseg = new Segment_pt(base,bound,level,2,16);
            reinterpret_cast<Segment_pt*>(newseg)->CreateMemPool();
            reinterpret_cast<Segment_pt*>(newseg)->DataArray =  reinterpret_cast<Segment_pt*>(newseg)->build_tree_none(base,bound);
            // reinterpret_cast<Segment_pt*>(newseg)->InitSplit();
            reinterpret_cast<Segment_pt*>(newseg)->SetNext(nullptr);
            reinterpret_cast<Segment_pt*>(newseg)->SetSegMax(8,0,bound - base);//build_tree_none待定
            return newseg;
        }
        

        //--------------------------new index--------------------------//
        
        //return the top index,SnodeArray[0] is segment layer,1-level are index
        SNode* NewIndexArray(int level,K base,K bound,vector<SNode*> &SnodeArray){
            // SNode* bottom_seg = new NewSegment_pt(base,bound);
            SNode* pr_ = NewSegment_pt(base,bound,level);
            SnodeArray[0] = pr_;
            for(int l = 1;l<=level;l++){
                SNode* curr = new Index(base);
                reinterpret_cast<Index*>(curr)->SetNext(nullptr);//no need to loop
                reinterpret_cast<Index*>(curr)->AppendBelow(pr_);
                SnodeArray[l] = curr;
                pr_ = curr;
            }
            return pr_;
        }

        SNode* NewIndexArrayWithData(int level,K base,K bound,subtree* n,vector<SNode*> &SnodeArray){
            // SNode* bottom_seg = new NewSegment_pt(base,bound);
            SNode* pr_ = NewSegment_pt(base,bound,level,n);
            SnodeArray[0] = pr_;
            for(int l = 1;l<=level;l++){
                SNode* curr = new Index(base);
                reinterpret_cast<Index*>(curr)->SetNext(nullptr);
                reinterpret_cast<Index*>(curr)->AppendBelow(pr_);
                SnodeArray[l] = curr;
                pr_ = curr;
            }
            return pr_;
        }

        SNode* NewIndex(int level,K base,K bound){
            // SNode* bottom_seg = new NewSegment_pt(base,bound);
            SNode* pr_ = NewSegment_pt(base,bound,level);
            for(int l = 0;l<level;l++){
                SNode* curr = new Index(base);
                reinterpret_cast<Index*>(curr)->SetNext(nullptr);
                reinterpret_cast<Index*>(curr)->AppendBelow(pr_);
                pr_ = curr;
            }
            return pr_;
        }
        
        //--------------------------skiplist operation--------------------------//
        
        inline int GetMaxHeight() const {
            return max_height_.load(std::memory_order_relaxed);
        }
        
        inline int RandLevel(){
            int lvl = 1;
            default_random_engine e(rd());
            uniform_int_distribution<int>u(0,9);
            for(int i = 0;i<MaxLevel-1;i++){
                if(u(e) > 4){
                    break;
                }
                lvl++;
            }
            return lvl;
        }

        //for Lookup
        SNode* Scan(K key){
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
            for(int l = top_level;l>0;l--){
                //not sure key > curr's key
                curr = reinterpret_cast<Index*>(pred)->Next();
                if(curr == nullptr){
                    curr = nullptr;
                }
                RT_ASSERT(curr!=nullptr);
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
                pred = reinterpret_cast<Index*>(pred)->Down();
            }
            //return segment_pt
            locate = pred;//reinterpret_cast<Seg*>(pred)->Down();
            return locate;
        }

        //for Add,preds[0,MaxLevel]
        SNode* Scan(K key, SNode* preds[]){
            // SNode* ppred = nullptr;
            SNode* pred = head_;
            SNode* curr = nullptr;
            SNode* succ = nullptr;
            SNode* locate = nullptr;
            int top_level = GetMaxHeight();//获取此时最高的index高度
            //down
            //获取[MaxLevel,top_level)区间高度的前驱，即为head_的不同高度的index
            for(int l = MaxLevel;l>top_level;l--){
                preds[l] = pred;
                // ppred = pred;
                pred = reinterpret_cast<Index*>(pred)->Down();
            }

            //get the preds from [top_level,0)
            for(int l = top_level;l>0;l--){//0 是因为 每个底层节点保证一定有一层index
                //not sure key > curr's key
                RT_ASSERT(pred->isIndex);
                curr = reinterpret_cast<Index*>(pred)->Next();//pred右侧index
                RT_ASSERT(curr!=nullptr);
                succ = reinterpret_cast<Index*>(curr)->Next();//curr右侧index
                while(succ && succ->Key <= key){//只有succ的左边界(key) 不超过 当前key时，pred右移
                    // ppred = pred;
                    pred = curr;
                    curr = succ;
                    succ = reinterpret_cast<Index*>(succ)->Next();
                }
                
                //case1 pred <= key < curr ;
                //case2 pred < curr <= key < succ
                if(key >= curr->Key){
                    //key is bwtween curr and succ,case 2
                    // ppred = pred;
                    pred = curr;//curr != tail_
                }
                // else{
                //     //key is before or located in curr
                //     if(!locate && key >= curr->Key){
                //         //key is located in curr
                //         locate = curr;
                //         preds[l] = curr;
                //         break;
                //     }
                // }
                preds[l] = pred;
                // ppred = pred;
                pred = reinterpret_cast<Index*>(pred)->Down();//当l = 1时，pred更新为底层节点
            }
            //return segment_pt
            locate = pred;
            preds[0] = locate;
            return locate;
        }

        //skiplist lookup a key
        std::pair<int,V> Lookup(K key){
            Segment_pt* locate = reinterpret_cast<Segment_pt*>(Scan(key));
            while(true){
                // int state = locate->delta_inserts.Wait(IsREAD_THREAD);
                locate->delta_inserts.Wait(IsREAD_THREAD);
                if(key < locate->bound){
                    // locate->work_threads.AddWork();
                    // //TODO:to check read thread中该部分应该在lookup结束后，还是addwork结束与lookup之间
                    // if(state == 0)
                    //     locate->delta_inserts.SignalOne();//避免threads一直无法被唤醒
                    std::pair<int,V> res = locate->Lookup(key);
                    locate->delta_inserts.SignalWork();
                    return res;
                }
                // if(state == 0)//同Add中对应code
                //     locate->delta_inserts.SignalOne();
                //skip
                locate->delta_inserts.SignalWork();
                locate = reinterpret_cast<Segment_pt*>(locate->Next());
            }
            return {-1,0};
        }

        //skiplist add a key 
        bool Add(K key,V value,int &path_depth,long long &scan_time,int &split_cnt,long long &collision){
            SNode* preds[MaxLevel+1];
            for(int i =0;i<=MaxLevel;i++){
                preds[i] = nullptr;
            }
            Segment_pt* locate = reinterpret_cast<Segment_pt*>(Scan(key,preds));
            while(true){    
                if(key < locate->bound){//locate可能由于split，bound缩小了，需要跳过
#if DEBUG
                    cerr<<"try to acquire"<<locate<<endl;
#endif
                    /*
                    从if 到Wait之间可能异常的情况：
					1.locate刚好split结束，边界值变化了，需要skip，只能先delta_inserts.Wait，然后判断locate边界， 
					2.前面有线程要对locate进行split操作，但还未结束，这时，Wait会使得本线程挂起 
					*/
                    int state = locate->delta_inserts.Wait();
                    // locate->ExclusiveLock();
                    if(key < locate->bound){
                        return Add_key_in_Segment(key,value,preds,locate,state,path_depth,scan_time,split_cnt,collision);//will ReleaseExclusive and finish operations on delta_inserts and work_threads
                        // return false;
                    }
                    // locate->ReleaseExclusive();
                    if(state == 0){//由于刚好前面的线程对locate进行了rebuild/split，本该当前线程notify all，现在要执行SignalOne再skip,
                        locate->delta_inserts.SignalWork();
                        locate->delta_inserts.SignalOne();
                    }//deta_insert=1 and wake up one blocked thread to do split/rebuild 
                    else                                    //which will notify the other threads 
                        locate->delta_inserts.SignalWork();//其他情况state > 0，可以乐观的认为插了一个重复的键，这里不做其他处理
                    #if DEBUG
                    cerr<<"Release "<<locate<<endl;
                    #endif
                }
                //skip 
                // while(locate->ReadSplit());
                SNode* next = locate->Next();
                //split时，底层是先连接后面的段，然后才改变第一个段的bound
                //可能出现locate的next在这两步被改的情况，如[2,10) [5,8) [8,10)
                //但key 一定大于第一个段分裂前的bound，因此直接取其next即可
                locate = reinterpret_cast<Segment_pt*>(next);
            }
            return false;	
        }

        //add a key-value into segment_pt
        bool Add_key_in_Segment(K key,V value,SNode** preds,Segment_pt* locate,int state,int &insert_depth,long long &scan_time,int &split_cnt,long long &collision){
            #if INSERT_ROUTE
            subtree* path[MAX_DEPTH*2];
            #else
            subtree* path[0];
            #endif
            int path_size = 0;
            int flag = 0;
            if(state == 0){
                int ele_size = locate->DataArray->size.load(std::memory_order_acquire);
                if(ele_size <= 1){
                    locate->delta_inserts.WaitWork();//wait working threads before this,work_thread --count
                    flag = 1;
                }
            }
            locate->insert_subtree(key,value,path,path_size,collision);
            insert_depth = path_size;
            #if INSERT_ROUTE
            delete path;
            #endif
            subtree *n = locate->DataArray;
            if(flag){
                //flag = 1,means state = 0 and the size of n <=2 after insert_subtree, DataArray was updated
                int new_size = n->size.load(std::memory_order_acquire);
                int new_slot_ = n->num_items;
                int key_space = n->stop - n->start;
                locate->SetSegMax(new_slot_,new_size, key_space);
                int new_segment_max_size = locate->SegmentMaxSize;
                locate->delta_inserts.SignalWork();
                #if SHOWBLOCK
                std::cout<<"rebuild model\t";
                #endif
                locate->delta_inserts.SignalAll(new_segment_max_size,new_size,new_slot_);//wake up blocked threads
                return false;
            }
            const bool need_rebuild = ( state == 0);
            if(!need_rebuild){
                //结束插入不需要修改delta_inserts
                locate->delta_inserts.SignalWork();//working threads number --
                return false;
            }else{
                locate->delta_inserts.WaitWork();//wait working threads before this
                const int ESIZE = n->size.load(std::memory_order_acquire);
                int new_segment_max_size = locate->SegmentMaxSize;
                int n_slots = n->num_items;
                const int REBUILD_SIZE = 64;
                const int SLOTS_ = 256;
                if(ESIZE < REBUILD_SIZE){
                    locate->delta_inserts.SignalWork();
                    //TODO:强制增大count
                    locate->delta_inserts.SignalAll(new_segment_max_size,ESIZE,SLOTS_);//wake up blocked threads
                    return false;
                }
                const auto start_time = std::chrono::steady_clock::now();
                if(ESIZE >= new_segment_max_size ){//split
                    K* keys = new K[ESIZE];
                    V* values = new V[ESIZE];
                    unsigned int n_start = n->start,n_stop = n->stop;
                    // SplitSegmentNoScan(locate,preds,ESIZE,n_start,n_stop);
                    locate->scan_and_destroy_subtree(n,keys,values,ESIZE);//std::pair<long double,int> res = 
                    locate->DataArray = nullptr;
                    // std::cout<<ESIZE<<",\t"<<n_slots<<std::endl;
                    SplitSegment(locate,preds,keys,values,ESIZE,n_start,n_stop);
                    //SplitSegment already released the write lock of locate
                    int new_slot = locate->DataArray->num_items;
                    int new_size = locate->DataArray->size.load(std::memory_order_acquire);
                    new_segment_max_size = locate->SegmentMaxSize;
                    delete[] keys;
                    delete[] values;
                    split_cnt++;
                    const auto end_time = std::chrono::steady_clock::now();
                    const auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);
                    locate->delta_inserts.SignalWork();
                    #if SHOWREBUILD
                    std::cout<<std::to_string(duration.count())<<std::endl;
                    #endif
                    locate->delta_inserts.SignalAll(new_segment_max_size,new_size,new_slot);//wake up blocked threads
                    scan_time+=duration.count();
                    // std::cout<<n_start<<"\t"<<n_stop<<"\t"<<locate->DataArray->stop<<endl;
                    return false;
                }else{
                    //rebuild每次将深度的树展开成扁平的数组，是否会影响并发性能
                    // if(ESIZE > 2e5){
                    //     locate->delta_inserts.SignalWork();
                    //     #if SHOWBLOCK
                    //     const auto end_time = std::chrono::steady_clock::now();
                    //     const auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);
                    //     std::cout<<"no build\t"<<std::to_string(duration.count())<<"\t";
                    //     #endif
                    //     locate->delta_inserts.SignalAll(new_segment_max_size,ESIZE,n_slots);//wake up blocked threads
                    //     return false;
                    // }
                    //TODO sum_of_top_k 若低于NOREBUILD_th 则不rebuild
                    #if TOPK_BUILD
                    const int n_Top = 10;//max(10,int(n->num_items * 0.01));
                    int *Top_k_pos = new int[n_Top];
                    int *Top_k_size = new int[n_Top];
                    int cnt_ = 0;
                    // std::pair<int,int> TopK_res = n->TopK.SumUP();
                    // double sum_topk = TopK_res.first*1.0 / ESIZE;
                    double sum_topk = GetTopk(n_Top,ESIZE,n_slots,n,NOREBUILD_th,Top_k_pos,Top_k_size,cnt_,scan_time);
                    // std::cout<<sum_topk<<std::endl;
                    if(sum_topk < NOREBUILD_th){
                        delete[] Top_k_pos;
                        delete[] Top_k_size;
                        locate->delta_inserts.SignalWork();
                        //激进的插入
                        const int next_count = n_slots; 
                        locate->delta_inserts.SignalAll(new_segment_max_size,ESIZE,next_count);//wake up blocked threads
                        return false;
                    }
                    delete[] Top_k_pos;
                    delete[] Top_k_size;
                    #endif
                    K* keys = new K[ESIZE];
                    V* values = new V[ESIZE];
                    unsigned int n_start = n->start,n_stop = n->stop;
                    locate->scan_and_destroy_subtree(n,keys,values,ESIZE);//std::pair<long double,int> res = 
                    locate->DataArray = nullptr;
                    std::pair<subtree*,int> res_ = locate->rebuild_tree(keys,values,ESIZE,n_start,n_stop);
                    locate->DataArray = res_.first;
                    int new_slot = locate->DataArray->num_items;
                    // std::cout<<ESIZE<<"\t"<<ESIZE*1.0/n_slots<<"\tafter rebuild:"<<ESIZE*1.0/new_slot<<std::endl;
                    locate->SetSegMax(new_slot,ESIZE,n_stop - n_start);
                    delete[] keys;
                    delete[] values;
                    const auto end_time = std::chrono::steady_clock::now();
                    const auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);
                    new_segment_max_size = locate->SegmentMaxSize;
                    locate->delta_inserts.SignalWork();
                    #if SHOWREBUILD
                    std::cout<<"rebuild new_segment_max_size:\t"<<new_segment_max_size<<"\t"<<duration.count()<<std::endl;
                    #endif
                    locate->delta_inserts.SignalAll(new_segment_max_size,ESIZE,new_slot);//wake up blocked threads
                    scan_time+=duration.count();
                    return true;
                }
            }
            return false;
        }

        double GetTopk(const int TOPK,const int ESIZE,const int n_slots,subtree *n,const double sum_k,int *Top_k_pos,int *Top_k_size,int &cnt_,long long &scan_time){
            struct pos{
                int x,y;//x,index;y,esize
                pos(){}
                pos(int x,int y):x(x),y(y){}
                bool operator< (const pos &a)const {
                    return y < a.y;
                }
            };
            struct cmp
            {
                bool operator()(const pos &a,const pos &b)
                {
                    return a.y < b.y;
                }
            };
            priority_queue<pos,vector<pos>,cmp > top_k_set;
            const auto start_time = std::chrono::steady_clock::now();
            // int TOPK = 10;
            // vector<double> distribution(n_slots);
            for(int i = 0;i<n_slots;i++){
                if(BITMAP_GET(n->none_bitmap,i) == 0){
                    if(BITMAP_GET(n->child_bitmap,i) == 1){
                        int num_ = n->items[i].comp.child->size.load(std::memory_order_acquire);
                        if(top_k_set.size() == TOPK){
                            if(top_k_set.top().y < num_){
                                top_k_set.pop();
                                top_k_set.push({i,num_});
                            }
                        }else{
                            top_k_set.push({i,num_});
                        }
                        // distribution[i] = num_*1.0/ESIZE;
                    }
                }
            }
            // sort(distribution.begin(),distribution.end(),greater<double>());
            double sum_topk = 0;
            // for(int i = 0;i<min(n_slots,TOPK);i++){
            //     sum_topk+=distribution[i];
            // }
            cnt_ = 0;
            int sum_esize = 0;
            while(!top_k_set.empty()){
                Top_k_pos[cnt_] = top_k_set.top().x;
                Top_k_size[cnt_] = top_k_set.top().y;
                sum_esize+=top_k_set.top().y;
                cnt_++;
                top_k_set.pop();
            }
            sum_topk = sum_esize*1.0/ESIZE;
            const auto end_time = std::chrono::steady_clock::now();
            const auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);
            scan_time+=duration.count();
            #if OUTPUT_TOPK
            string out = std::to_string(n->start)+"\ttotal size:"+std::to_string(ESIZE)+
                "\ttotal slots:"+std::to_string(n->num_items) + "\ttop10 sum:"+std::to_string(sum_topk)+"\n";
            std::cout<<out;
            // for(int i = 0;i<min(TOPK,n_slots);i++){
            //     out+="\t";
            //     out+=(std::to_string(distribution[i]));
            // }
            #endif
            // std::cout << "overhead of scan top array:" << duration.count() << "us"<<std::endl;
            return sum_topk;
        }

        bool SplitSegmentNoScan(Segment_pt *root,SNode** preds,const int ESIZE,K _start,K _stop){
            struct Seg{
                int begin;
                int height;
                subtree *array;
            };
            std::stack<Seg> s;
            subtree *head_prepare = nullptr;
            subtree *tail_prepare = nullptr;
            int keylen = ESIZE;
            const int n_slot = root->DataArray->num_items;
            int seg_size = 0;
            //1. partition array
            //get head tail
            if(BITMAP_GET(root->DataArray->none_bitmap,0) == 0){
                if(BITMAP_GET(root->DataArray->child_bitmap,0) == 1){
                    int num_ = root->DataArray->items[0].comp.child->size.load(std::memory_order_acquire);
                    if(num_ > 1e2){
                        head_prepare = root->DataArray->items[0].comp.child;
                        keylen-=num_;
                        seg_size++;
                    }
                }
            }
            if(BITMAP_GET(root->DataArray->none_bitmap,n_slot-1) == 0){
                if(BITMAP_GET(root->DataArray->child_bitmap,n_slot-1) == 1){
                    int num_ = root->DataArray->items[n_slot-1].comp.child->size.load(std::memory_order_acquire);
                    if(num_ > 1e2){
                        tail_prepare = root->DataArray->items[n_slot-1].comp.child;
                        keylen-=num_;
                        seg_size++;
                    }
                }
            }
            RT_ASSERT(keylen > 0);
            K *keys = new K[keylen];
            V * values = new V[keylen];
            s.push({0,0, root->DataArray});
            while (!s.empty()) {
                int begin = s.top().begin;
                int height = s.top().height;
                subtree* n = s.top().array;
                s.pop();
                for (int i = 0; i < n->num_items;i++) {
                    if(i == 0 && height == 0 && head_prepare != nullptr){
                        continue;
                    }
                    if(i == n->num_items-1 && height == 0 && tail_prepare != nullptr){
                        continue;
                    }
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
                            s.push({begin,height+1,n->items[i].comp.child});
                            begin += n->items[i].comp.child->size.load(std::memory_order_acquire);
                        }
                    }
                }
                if(n->is_two){
                    n->build_size = 2;
                    n->size.store(2, std::memory_order_release);
                    n->fixed = 0;
                    n->confliction = 0;
                    n->num_items = 8;
                    memset(n->items,0,(n->num_items)*sizeof(Item));
                    n->none_bitmap[0] = 0xff;
                    n->child_bitmap[0] = 0;
                    // n->TopK.Init();
                    while(!root->CASPool()){;}
                    root->tree_pool.push(n);
                    root->ReleasePool();
                }else{
                    root->delete_items(n->items, n->num_items);
                    const int bitmap_size = BITMAP_SIZE(n->num_items);
                    root->delete_bitmap(n->none_bitmap, bitmap_size);
                    root->delete_bitmap(n->child_bitmap, bitmap_size);
                    root->delete_subtree(n, 1);
                }
            }
        
            int st = 0,ed =0;//the rank of the first/last+1 ele of one segment
            GreedyPLR* plr = new GreedyPLR(gamma);
            vector<Segment*> seg;
            vector<int> segment_stIndex;
            for (int i = 0; i < keylen; i++) {
                Segment* seg_res = nullptr;
                seg_res = plr->Process(static_cast<double>(keys[i])-static_cast<double>(keys[st]), static_cast<double>(i-st));
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
            int segment_max_height = 1;//the max height of split segment[]
            seg_size += static_cast<int>(seg.size());
            
            //create new segment and index
            vector<Segment_pt*> split_segments(seg_size);//store the all segment_pt
            vector<SNode*> SNode_inner_pred(MaxLevel+1,nullptr);//段group当前右边界
            vector<SNode*> SNode_left(MaxLevel+1,nullptr);//段group左边界
            vector<SNode*> SNode_right(MaxLevel+1,nullptr);//段group右边界
            split_segments[0] = root;
            subtree *root_dataarray = nullptr;
            unsigned int split_lower_bound = 0;
            // std::cout<<"split range:"<<_start<<"\t"<<_stop<<std::endl;
            for(int i = 0,j = 0 ;i < seg_size ; i++) {
                int height_tmp=1;
                int st=segment_stIndex[j],ed;//st ed代表keys下标
                unsigned int base,bound;
                subtree* new_subtree = nullptr;
                if(i == 0 && head_prepare){
                    new_subtree = head_prepare;
                    base = new_subtree->start,bound = new_subtree->stop;
                }else if(i == seg_size-1 && tail_prepare){
                    new_subtree = tail_prepare;
                    base = new_subtree->start;//keys[keylen-1]
                    bound = new_subtree->stop;
                    height_tmp = RandLevel();
                    segment_max_height = max(segment_max_height,height_tmp);
                }else{
                    if(i == 0){
                        base = _start;//as head,no head_prepare
                    }else{//middle segment
                        height_tmp = RandLevel();
                        segment_max_height = max(segment_max_height,height_tmp);
                        if(i == 1 &&head_prepare){
                            base = head_prepare->stop;
                        }else{
                            base  = keys[st];
                        }
                    }

                    if(i == seg_size-1){//as tail,no tail_prepare
                        ed = keylen;
                        bound = _stop;
                    }else{//middle segment,maybe exsit tail_prepare
                        if(j+1 == segment_stIndex.size()){
                            ed = keylen-1;
                            bound = tail_prepare->start;
                        }else{
                            ed = segment_stIndex[j+1];
                            bound = keys[ed];
                        }
                    }
                    j++;


                    if(ed-st == 1){
                        new_subtree = root->build_tree_none(base,bound);
                        BITMAP_CLEAR(new_subtree->none_bitmap,0);
                        new_subtree->size.store(1, std::memory_order_release);
                        new_subtree->items[0].comp.data.key = keys[st];
                        new_subtree->items[0].comp.data.value = values[st];
                    }else{
                        std::pair<subtree*,int> res_;
                        #if USEPLR
                        if(seg[i]->slope > 0)
                            res_ = root->rebuild_tree(keys+st,values+st,ed-st,base,bound,true,seg[i]->slope,seg[i]->intercept);
                        else
                            res_ = root->rebuild_tree(keys+st,values+st,ed-st,base,bound);
                        #else
                            res_ = root->rebuild_tree(keys+st,values+st,ed-st,base,bound);
                        #endif
                        new_subtree = res_.first;
                    }
                }
                
                //2.4 create segment ,index without data and link inner part
                SNode* next_seg = nullptr;//the bottom segment_pt
                // SNode* top_index = nullptr;//the top level of index
                vector<SNode*> SNodeArray(height_tmp+1,nullptr);
                if(i){
                    //create index array(stored in SNodeArray) with subtree
                    NewIndexArrayWithData(height_tmp,base,bound,new_subtree,SNodeArray);
                    next_seg = SNodeArray[0];
                    split_segments[i] = reinterpret_cast<Segment_pt*>(next_seg);//store segment in split_segments
                    reinterpret_cast<Segment_pt*>(next_seg)->delta_inserts.AddWork();//限制split
                    //2.5.1 link segment_pt if pre_segment exsits
                    if(SNode_inner_pred[0]){
                        reinterpret_cast<Segment_pt*>(SNode_inner_pred[0])->SetNext(next_seg);//no need to loop   
                    }
                    SNode_inner_pred[0] = next_seg;
                    if(i == 1)
                        SNode_left[0] = next_seg;
                    SNode_right[0] = next_seg;
                    //link [height_tmp,0) index if pre_segment exsits
                    for(int l = height_tmp;l>0;l--){
                        SNode* pr_ = SNodeArray[l];//get the l level index 
                        // reinterpret_cast<Index*>(pr_)->SetNext(nullptr);//this step finished in NewIndexArrayWithData  function
                        if(SNode_inner_pred[l]){
                            reinterpret_cast<Index*>(SNode_inner_pred[l])->SetNext(pr_);//no need to loop
                            SNode_inner_pred[l] = pr_;
                        }else{
                            SNode_inner_pred[l] = pr_;
                            SNode_left[l] = pr_;
                        }
                        SNode_right[l] = pr_;
                    }
                }else{
                    root_dataarray = new_subtree;
                    split_lower_bound = bound;
                    root->SetSegMax(new_subtree->num_items,new_subtree->size.load(std::memory_order_acquire),bound - base);
                }
            }

            if(seg_size == 1){
                root->DataArray = root_dataarray;
                delete plr;
                #if SHOWREBUILD
                std::cout<<_size<<" split new segment number:"<<seg_size-1<<"\t";
                #endif
                return true;
            }else{
                //link the bottom(segment_pt)
                reinterpret_cast<Segment_pt*>(SNode_right[0])->SetNext(root->Next());//the right link
                // split_segments[seg_size-1]->SetNext(root->Next());
                root->SetNext(SNode_left[0]);//the left link
                //2.7 update root(segment_pt)
                root->bound = split_lower_bound;
                root->DataArray = root_dataarray;
                //2.8 release the write lock of the new segment
                // for(int i = 1;i<seg_size;i++){
                //     split_segments[i]->ReleaseExclusive();
                //     #if DEBUG
                //     cerr<<"Release "<<split_segments[i]<<endl;
                //     #endif
                // }
                //2.9 update maxheight
                int max_height = GetMaxHeight();
                while (segment_max_height > max_height) {
                    if (max_height_.compare_exchange_weak(max_height, segment_max_height)) {
                        // successfully updated it
                        max_height = segment_max_height;
                        break;
                    }
                }

                //2.10 insert index from 1 to segment_max_height, reusing the information of preds array 
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
                    split_segments[i]->delta_inserts.SignalWork();
                }
                delete plr;
                // std::cout<<"raw size:"<<_size<<"key_space:"<<_start<<","<<_stop<<" split new segment number:"<<seg_size-1<<std::endl;
                #if SHOWREBUILD
                std::cout<<_size<<" split new segment number:"<<seg_size-1<<"\t";
                #endif
                return true;
            }       
        }

        bool SplitSegment(Segment_pt *root,SNode** preds,K *_keys,int*  _values,int _size,K _start,K _stop){
            //1. partition array
            int st = 0,ed =0;//the rank of the first/last+1 ele of one segment
            GreedyPLR* plr = new GreedyPLR(gamma);
            vector<Segment*> seg;
            vector<int> segment_stIndex;
            vector<int> height_seq;
            int segment_max_height = 1;//the max height of split segment
            for (int i = 0; i < _size; i++) {
                Segment* seg_res = nullptr;
                seg_res = plr->Process(static_cast<double>(_keys[i])-static_cast<double>(_keys[st]), static_cast<double>(i-st));
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
            //create new segment and index
            vector<Segment_pt*> split_segments(seg_size);//store the all segment_pt
            vector<SNode*> SNode_inner_pred(MaxLevel+1,nullptr);//段group当前右边界
            vector<SNode*> SNode_left(MaxLevel+1,nullptr);//段group左边界
            vector<SNode*> SNode_right(MaxLevel+1,nullptr);//段group右边界

            split_segments[0] = root;
            subtree *root_dataarray = nullptr;
            unsigned int split_lower_bound = 0;
            for(int i = 0;i<seg_size;i++){
                int height_tmp=1,ed,st=segment_stIndex[i];
                unsigned int base,bound;
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
                
                RT_ASSERT(ed-st>0);
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
                    //TODO:test plr influence
                    std::pair<subtree*,int> res_;
                    #if USEPLR
                        if(seg[i]->slope > 0)
                            res_ = root->rebuild_tree(_keys+st,_values+st,ed-st,base,bound,true,seg[i]->slope,seg[i]->intercept);
                        else
                            res_ = root->rebuild_tree(_keys+st,_values+st,ed-st,base,bound);
                    #else
                        res_ = root->rebuild_tree(_keys+st,_values+st,ed-st,base,bound);
                    #endif
                    new_subtree = res_.first;
                    #if 0
                    string strs = to_string(res_.second)+"\t";
                    std::cout<<strs;
                    #endif
                }

                //2.4 create segment ,index without data and link inner part
                SNode* next_seg = nullptr;//the bottom segment_pt
                // SNode* top_index = nullptr;//the top level of index
                vector<SNode*> SNodeArray(height_tmp+1,nullptr);
                if(i){
                    //2.4.1 create index array(stored in SNodeArray) with subtree
                    NewIndexArrayWithData(height_tmp,base,bound,new_subtree,SNodeArray);
                    next_seg = SNodeArray[0];
                    split_segments[i] = reinterpret_cast<Segment_pt*>(next_seg);//store segment in split_segments
                    reinterpret_cast<Segment_pt*>(next_seg)->delta_inserts.AddWork();//限制split
                    //2.5.1 link segment_pt if pre_segment exsits
                    if(SNode_inner_pred[0]){
                        reinterpret_cast<Segment_pt*>(SNode_inner_pred[0])->SetNext(next_seg);//no need to loop   
                    }
                    SNode_inner_pred[0] = next_seg;
                    if(i == 1)
                        SNode_left[0] = next_seg;
                    SNode_right[0] = next_seg;
                    //2.5.2 link [height_tmp,0) index if pre_segment exsits
                    for(int l = height_tmp;l>0;l--){
                        SNode* pr_ = SNodeArray[l];//get the l level index 
                        // reinterpret_cast<Index*>(pr_)->SetNext(nullptr);//this step finished in NewIndexArrayWithData  function
                        if(SNode_inner_pred[l]){
                            reinterpret_cast<Index*>(SNode_inner_pred[l])->SetNext(pr_);//no need to loop
                            SNode_inner_pred[l] = pr_;
                        }else{
                            SNode_inner_pred[l] = pr_;
                            SNode_left[l] = pr_;
                        }
                        SNode_right[l] = pr_;
                    }

                }else{
                    root_dataarray = new_subtree;
                    split_lower_bound = bound;
                    root->SetSegMax(new_subtree->num_items,new_subtree->size.load(std::memory_order_acquire),bound - base);
                }
                  
            }

            if(seg_size == 1){
                root->DataArray = root_dataarray;
                delete plr;
                #if SHOWREBUILD
                std::cout<<_size<<" split new segment number:"<<seg_size-1<<"\t";
                #endif
                return true;
            }else{
                //2.6 link the bottom(segment_pt)
                reinterpret_cast<Segment_pt*>(SNode_right[0])->SetNext(root->Next());//the right link
                // split_segments[seg_size-1]->SetNext(root->Next());
                root->SetNext(SNode_left[0]);//the left link
                //2.7 update root(segment_pt)
                root->bound = split_lower_bound;
                root->DataArray = root_dataarray;
                // root->ReleaseSplit();
                // root->ReleaseExclusive();
                // #if DEBUG
                // cerr<<"Release "<<root<<endl;
                // #endif
                //2.8 release the write lock of the new segment
                // for(int i = 1;i<seg_size;i++){
                //     split_segments[i]->ReleaseExclusive();
                //     #if DEBUG
                //     cerr<<"Release "<<split_segments[i]<<endl;
                //     #endif
                // }
                //2.9 update maxheight
                int max_height = GetMaxHeight();
                while (segment_max_height > max_height) {
                    if (max_height_.compare_exchange_weak(max_height, segment_max_height)) {
                        // successfully updated it
                        max_height = segment_max_height;
                        break;
                    }
                }

                //2.10 insert index from 1 to segment_max_height, reusing the information of preds array 
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
                    split_segments[i]->delta_inserts.SignalWork();
                }
                delete plr;
                // std::cout<<"raw size:"<<_size<<"key_space:"<<_start<<","<<_stop<<" split new segment number:"<<seg_size-1<<std::endl;
                #if SHOWREBUILD
                std::cout<<_size<<" split new segment number:"<<seg_size-1<<"\t";
                #endif
                return true;
            }
            return true;
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
            while(cur->Key != UNINT_MAX){
                if(write_){
                    auto ESIZE = reinterpret_cast<Segment_pt*>(cur)->DataArray->size.load(std::memory_order_relaxed);
                    string kk = to_string(cur->Key)+"\t"+to_string(reinterpret_cast<Segment_pt*>(cur)->bound)+"\t"+to_string(ESIZE)+"\t"+
                    to_string(ESIZE*1.0/(reinterpret_cast<Segment_pt*>(cur)->DataArray->num_items))+"\n";
                    write_into_file("./segment.csv",kk.c_str());
                }
                cnt++;
                cur = cur->Next();
            }
            std::cout<<"SegmentNumber:"<<cnt<<std::endl;
        }

        //--------------------------skiplist member variable--------------------------//
        std::atomic<int> max_height_; 
        int MaxLevel;
		int gamma;
        int segment_max_size;
        double linearity;
        SNode* head_;
        SNode* tail_;
        skiplist();
        skiplist(int MaxLevel,int gamma):max_height_(1),MaxLevel(MaxLevel),gamma(gamma),segment_max_size(SEGMENT_MAX_SIZE),linearity(LINAERITY){
            vector<SNode*> left_index(MaxLevel+1,nullptr);//store the all level SNode of head_
            vector<SNode*> right_index(MaxLevel+1,nullptr);//store the all level SNode of tail_
            head_ = NewIndexArray(MaxLevel,0,0,left_index);
            tail_ = NewIndexArray(MaxLevel,UNINT_MAX,UNINT_MAX,right_index);

            int level_idx1 = RandLevel();
            vector<SNode*> middle_index(level_idx1+1,nullptr);//store the all level SNode of Idx1
            // SNode* Idx1 = NewIndexArray(level_idx1,1,UNINT_MAX,middle_index);
            NewIndexArray(level_idx1,1,UNINT_MAX,middle_index);

            //update the height of skiplist
            int max_height = 1;
            while (level_idx1 > max_height) {
                if (max_height_.compare_exchange_weak(max_height, level_idx1)) {
                // successfully updated it
                max_height = level_idx1;
                break;
                }
            }

            {//segment link
                reinterpret_cast<Segment_pt*>(left_index[0])->SetNext(middle_index[0]);
                reinterpret_cast<Segment_pt*>(middle_index[0])->SetNext(right_index[0]);
            }
            {//index link
                for(int i = 1;i<=level_idx1;i++){
                    reinterpret_cast<Index*>(left_index[i])->SetNext(middle_index[i]);
                    reinterpret_cast<Index*>(middle_index[i])->SetNext(right_index[i]);
                }
                for(int i = level_idx1+1;i<=MaxLevel;i++){
                    reinterpret_cast<Index*>(left_index[i])->SetNext(right_index[i]);
                }
            }
        }
};