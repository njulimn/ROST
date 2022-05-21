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
#define Gm 64
#define MAX_DEPTH 6
#define INSERT_ROUTE 0
#define SEGMENT_MAX_SIZE 2e5
#define WRITESEG 0
#define DELTA_INSERT 2e5
#define IsREAD_THREAD 1
#define StateType int

#define PLR_DATA_PREPROCESS 0 //control whether data of plr is preprocessed

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
    else{
        count = min(int(slot_*4),remain_size);
    }
    count = min(int(DELTA_INSERT),count);
    return count;
}

int ComputeSegmentMaxSize(int slot_,int esize,unsigned int key_space){
    int res = 0; 
    res = min((unsigned int)(2e5),key_space);
    return res;
}

template<class K, class V>
class skiplist {
    public:
        class Semaphore {
        public:
            explicit Semaphore(int count) : count_(count),work_(0) {
            }
            ~Semaphore(){}
            //reset count_ and notify blocked threads
            void SignalAll(const int seg_max_size,int esize,int slot_=DELTA_INSERT,int id=0){
                //seg_max_size:At present, the maximum elements that can be accommodated in the corresponding segment
                std::unique_lock<std::mutex> lock(mutex_);
                count_ = ComputeNextCount(seg_max_size,esize,slot_);
                cv_.notify_all();
            }
            //reset count_ as 1, nofity one blocked thread
            void SignalOne(int id=0){
                std::unique_lock<std::mutex> lock(mutex_);
                count_ = 1;
                cv_.notify_one();
            }

            int Wait(int read_thread = 0,int id=0) {
                std::unique_lock<std::mutex> lock(mutex_);
                --count_;
                while(count_ < 0){
                    cv_.wait(lock);
                    --count_;
                }

                if(read_thread){
                    ++count_;
                    if(count_ == 1) 
                        cv_.notify_one();
                }
                std::unique_lock<std::mutex> lock_(mutex_w);
                ++work_;
                return count_; 
            }

            //for work threads
            void AddWork(int id=0){
                std::unique_lock<std::mutex> lock(mutex_w);
                ++work_;
            }
            //for work threads V
            void SignalWork(int id=0){
                std::unique_lock<std::mutex> lock(mutex_w);
                --work_;
                if (work_ == 0){
                    cv_w.notify_one();
                }
            }
            //for work threads P
            void WaitWork(int id=0){
                std::unique_lock<std::mutex> lock(mutex_w);
                --work_;
                if(work_>0){
                    cv_w.wait(lock);
                }
                if(work_ == 0){
                    work_ = 1;
                }
            }

        private:
            std::mutex mutex_;//mutex for delta_insert
            std::condition_variable cv_;
            int count_;//Records the current number of elements that can be inserted

            std::mutex mutex_w;//mutex for semaphore 'work'
            std::condition_variable cv_w;
            int work_;//Records the current number of threads that are working
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

            bool isIndexNode(){
                return isIndex;
            }

            K Key;//the lower bound
            bool isIndex; //read only, true-> index node, false -> sgement node
            std::atomic<SNode*> next_; // next 指针用于构建链表，直接继承 
            SNode(K base,bool type):Key(base),isIndex(type){}  
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
                double v = this->slope * (static_cast<long double>(key)) + this->intercept;
                if(v > std::numeric_limits<int>::max() / 2) {
                    return this->num_items - 1;
                }
                if(v < 0 || static_cast<int>(v)<0){
                    return 0;
                }
                return std::min(this->num_items - 1, static_cast<int>(v));
            }

            //compute the range of subtree
            std::pair<unsigned int,unsigned int> computeRange(int pos){
                unsigned int start_ = (pos - this->intercept)/(this->slope * 1000.0)*1000.0;
                unsigned int stop_ = (pos + 1 - this->intercept)/(this->slope * 1000.0)*1000.0;
                if(pos == 0)
                    start_ = this->start;
                if(pos == this->num_items-1)
                    stop_ = this->stop;
                return {start_,stop_};
            }

            std::pair<bool,V> find_key_in_subtree(K key){
                subtree* n = this;
                int pos = n->predict(key);
                n->LockItem(pos);
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
                        int pos_next = next_subtree->predict(key);
                        next_subtree->LockItem(pos_next);
                        n->ReleaseItem(pos);
                        n = next_subtree;
                        pos = pos_next;
                    }
                }
                return {false,0};
            }

            void find_key_in_subtree(K key1,K key2,std::vector<std::pair<K,V>> &result){
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

            //get the lock of n->items[pos]
            inline void LockItem(int n){
                this->mlock->lock();
            }

            //release the lock of n->items[pos]
            inline void ReleaseItem(int n){
                this->mlock->unlock();
            }

            long long SpaceSize(){
                long long res = sizeof(subtree) + sizeof(mutex);
                const int bitmap_size_ = BITMAP_SIZE(num_items);
                res += (sizeof(bitmap_t) * bitmap_size_ * 2);
                res += (sizeof(Item) * num_items);
                for(int i = 0;i<num_items;i++){
                    if(BITMAP_GET(none_bitmap,i) == 0){
                        if(BITMAP_GET(child_bitmap,i) == 1){
                            res += ((items[i].comp.child)->SpaceSize());
                        }
                    }
                }
                return res;
            }

            int is_two = 0; // is special node for only two keys
            std::atomic<int> size;
            int num_items = 0; // size of items
            double slope = 0,intercept = 0;
            K start = UNINT_MAX,stop=UNINT_MAX;
            Item* items = nullptr;
            std::mutex *mlock = nullptr;
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
                       
            subtree* build_tree_none(K st,K ed){
                subtree* n = new_subtree(1);
                n->is_two = 0;
                n->size.store(0, std::memory_order_release);
                n->num_items = 1;
                n->slope = n->intercept = 0;
                n->start = st;
                n->stop = ed;
                n->items = new_items(1);
                memset(n->items,0,sizeof(Item));
                n->mlock = new mutex();
                n->none_bitmap = new_bitmap(1);
                n->none_bitmap[0] = 0;
                BITMAP_SET(n->none_bitmap, 0);
                n->child_bitmap = new_bitmap(1);
                n->child_bitmap[0] = 0;
                return n;
            }

            subtree* build_tree_two(K key1, int value1, K key2, int value2,K start=0,K stop=0,subtree* x = nullptr){
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
                    n->mlock = new mutex();
                    n->none_bitmap = new_bitmap(1);
                    n->child_bitmap = new_bitmap(1);
                    n->none_bitmap[0] = 0xff;
                    n->child_bitmap[0] = 0;
                }else{
                    n = tree_pool.top();
                    tree_pool.pop();
                }
                ReleasePool();

                if(!x){
                    n->start = start;
                    n->stop = stop;  
                }else{
                    n->start = x->start;
                    n->stop = x->stop;
                }
                
                const long double mid1_key = static_cast<long double>(key1);
                const long double mid2_key = static_cast<long double>(key2);

                const double mid1_target = n->num_items / 3;
                const double mid2_target = n->num_items * 2 / 3;

                n->slope = (mid2_target - mid1_target) / (mid2_key - mid1_key);
                RT_ASSERT(isinf(n->slope)==0);
                n->intercept = mid1_target - n->slope * mid1_key;
                RT_ASSERT(isfinite(n->slope));
                RT_ASSERT(isfinite(n->intercept));

                { // insert key1&value1
                    int pos = n->predict(key1);
                    RT_ASSERT(BITMAP_GET(n->none_bitmap, pos) == 1);
                    BITMAP_CLEAR(n->none_bitmap, pos);
                    n->items[pos].comp.data.key = key1;
                    n->items[pos].comp.data.value = value1;
                }
                { // insert key2&value2
                    int pos = n->predict(key2);
                    RT_ASSERT(BITMAP_GET(n->none_bitmap, pos) == 1);
                    BITMAP_CLEAR(n->none_bitmap, pos);
                    n->items[pos].comp.data.key = key2;
                    n->items[pos].comp.data.value = value2;
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

            inline void Lookup(K key1,K key2,std::vector<std::pair<K,V>> &result){
                this->DataArray->find_key_in_subtree(key1,key2,result);
            }
        
            void insert_subtree(K key,V value,int &path_size,long long &collision){
                subtree* n = DataArray;
                int pos = n->predict(key);
                StateType state_raw;
                n->mlock->lock();
                state_raw = BITMAP_GET(n->none_bitmap, pos) == 1?0:BITMAP_GET(n->child_bitmap, pos)+1;
                while(1){
                    path_size++;
                    if(state_raw == ItemState::Empty){
                        BITMAP_CLEAR(n->none_bitmap, pos);
                        n->items[pos].comp.data.key = key;
                        n->items[pos].comp.data.value = value;
                        n->size.fetch_add(1, std::memory_order_acquire);
                        n->mlock->unlock();
                        break;
                    }else if(state_raw == ItemState::Element){
                        //build model for DataArray
                        if(path_size == 1 &&  n->size.load(std::memory_order_acquire) == 1){
                            subtree* newtree =  build_tree_two(key, value, n->items[pos].comp.data.key, n->items[pos].comp.data.value,0,0,n);
                            memcpy(DataArray, newtree, sizeof(subtree));//shallow copy,so just free newtree
                            delete_subtree(newtree, 1);
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
                            path_size++;
                            n->mlock->unlock();
                            collision++;
                            break;
                        }
                    }else{//ItemState::Subtree
                        n->size.fetch_add(1, std::memory_order_acquire);
                        subtree* next_n = n->items[pos].comp.child;//n->items[pos].comp.child;
                        int next_pos = next_n->predict(key);
                        next_n->mlock->lock();
                        state_raw = BITMAP_GET(next_n->none_bitmap, next_pos) == 1?0:BITMAP_GET(next_n->child_bitmap, next_pos)+1;
                        n->mlock->unlock();
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
            }

            std::pair<subtree*,int> rebuild_tree(K *_keys,int* _values,int _size,K start,K stop,bool plr = false,double top_slope = 0,double top_intercept = 0){
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
                ret->start = start;
                ret->stop = stop;
                s.push({0, _size, 1, ret});
                int max_depth = 0;
                while(!s.empty()){
                    const int begin = s.front().begin;
                    const int end = s.front().end;
                    const int lvl = s.front().lvl;
                    max_depth = max(max_depth,lvl);
                    subtree* n = s.front().n;
                    s.pop();
                    RT_ASSERT(end - begin >= 2);
                    if(end-begin == 2){
                        //start and stop are set up in advance
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
                        n->size.store(size, std::memory_order_release);
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
                            RT_ASSERT(0 <= mid1_pos);
                            RT_ASSERT(mid1_pos < mid2_pos);
                            RT_ASSERT(mid2_pos < size - 1);
                            const long double mid1_key =
                                    (static_cast<long double>(keys[mid1_pos]) + static_cast<long double>(keys[mid1_pos + 1]) ) / 2;
                            const long double mid2_key =
                                    (static_cast<long double>(keys[mid2_pos]) + static_cast<long double>(keys[mid2_pos + 1]) ) / 2;

                            n->num_items = size * static_cast<int>(BUILD_GAP_CNT + 1);
                            
                            const double mid1_target = mid1_pos * static_cast<int>(BUILD_GAP_CNT + 1) + static_cast<int>(BUILD_GAP_CNT + 1) / 2;
                            const double mid2_target = mid2_pos * static_cast<int>(BUILD_GAP_CNT + 1) + static_cast<int>(BUILD_GAP_CNT + 1) / 2;
                            RT_ASSERT(abs(mid2_key - mid1_key)>1e-8);
                            n->slope = (mid2_target - mid1_target) / (mid2_key - mid1_key);
                            n->intercept = mid1_target - n->slope * mid1_key;
                            RT_ASSERT(isfinite(n->slope));
                            RT_ASSERT(isfinite(n->intercept));
                        }        
                        RT_ASSERT(n->slope >= 0);

                        n->items = new_items(n->num_items);
                        memset(n->items,0,n->num_items*sizeof(Item));
                        n->mlock = new mutex();
                        const int bitmap_size = BITMAP_SIZE(n->num_items);
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
                            RT_ASSERT(item_i < n->num_items);
                            if (next == offset + 1) {
                                RT_ASSERT(BITMAP_GET(n->none_bitmap,item_i) == 1);
                                RT_ASSERT(BITMAP_GET(n->child_bitmap,item_i) == 0);
                                BITMAP_CLEAR(n->none_bitmap, item_i);
                                n->items[item_i].comp.data.key = keys[offset];
                                n->items[item_i].comp.data.value = values[offset];
                            } else {
                                RT_ASSERT(BITMAP_GET(n->none_bitmap,item_i) == 1);
                                RT_ASSERT(BITMAP_GET(n->child_bitmap,item_i) == 0);
                                BITMAP_CLEAR(n->none_bitmap, item_i);
                                BITMAP_SET(n->child_bitmap, item_i);
                                n->items[item_i].comp.child = new_subtree(1);
                                std::pair<unsigned int,unsigned int> line = n->computeRange(item_i);
                                n->items[item_i].comp.child->start = min(line.first,keys[offset]);
                                n->items[item_i].comp.child->stop = line.second;
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
                RT_ASSERT(ret->start == start);
                RT_ASSERT(ret->stop == stop);
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
                    n->mlock = new mutex();
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
            
            //reset segment max size when build model or split
            void SetSegMax(int slot_,int esize,unsigned int key_space){
                SegmentMaxSize = ComputeSegmentMaxSize(slot_,esize,key_space);
            }

            long long SpaceSize(){
                long long space_sum = sizeof(Segment_pt);
                //DataArray
                space_sum += DataArray->SpaceSize();
                //tree_pool
                return space_sum;
            }

            K bound;//right guard,the upper bound
            int level;
            subtree *DataArray;
            Semaphore delta_inserts;//concurrency control
            int SegmentMaxSize;
            std::atomic<bool> lock_pool;
            std::stack<subtree*> tree_pool;
            
            Segment_pt(K base,K bnd,int lvl,int deta_insert_init = DELTA_INSERT):SNode(base,false),bound(bnd),level(lvl),
            delta_inserts(Semaphore(deta_insert_init)),lock_pool(false){}
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
            Index(K base):SNode(base,true){}
        };      

        //--------------------------new segment_pt--------------------------//

        SNode* NewSegment_pt(K base,K bound,int level,subtree* n){
            //data array = n,if no data,n = nullptr
            if(!n)
                return NewSegment_pt(base,bound,level);
            SNode *newseg = nullptr;
            int slots = (n->num_items);
            int ele_size = n->size.load(std::memory_order_relaxed);
            if(ele_size < 2){
                newseg = new Segment_pt(base,bound,level,1);//if the number of ele in subtree is smaller than 2,
                                                            //then set delta_insert = 2,else set delta_insert = DELTA_INSERT
            }
            else{
                int seg_max_size = ComputeSegmentMaxSize(slots,ele_size,bound - base);
                int count_set = ComputeNextCount(seg_max_size,ele_size,slots);
                newseg = new Segment_pt(base,bound,level,count_set);
            }
            reinterpret_cast<Segment_pt*>(newseg)->DataArray = n;
            reinterpret_cast<Segment_pt*>(newseg)->CreateMemPool();
            reinterpret_cast<Segment_pt*>(newseg)->SetNext(nullptr);
            reinterpret_cast<Segment_pt*>(newseg)->SetSegMax(slots,ele_size,bound - base);
            return newseg;
        }
        
        SNode* NewSegment_pt(K base,K bound,int level){
            SNode *newseg = new Segment_pt(base,bound,level,2);
            reinterpret_cast<Segment_pt*>(newseg)->CreateMemPool();
            //data array = build_tree_none
            reinterpret_cast<Segment_pt*>(newseg)->DataArray =  reinterpret_cast<Segment_pt*>(newseg)->build_tree_none(base,bound);
            reinterpret_cast<Segment_pt*>(newseg)->SetNext(nullptr);
            reinterpret_cast<Segment_pt*>(newseg)->SetSegMax(8,0,bound - base);
            return newseg;
        }
        

        //--------------------------new index--------------------------//
        
        //return the top index,SnodeArray[0] is segment layer,1-level are index layer
        SNode* NewIndexArray(int level,K base,K bound,vector<SNode*> &SnodeArray){
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
            for(int l = top_level;l>0;l--){
                //not sure key > curr's key
                RT_ASSERT(pred->isIndex);
                curr = reinterpret_cast<Index*>(pred)->Next();
                RT_ASSERT(curr!=nullptr);
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
                pred = reinterpret_cast<Index*>(pred)->Down();//pred will be updated as bottom segment node when l = 1
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
                    std::pair<int,V> res = locate->Lookup(key);
                    locate->delta_inserts.SignalWork();
                    return res;
                }
                locate->delta_inserts.SignalWork();
                locate = reinterpret_cast<Segment_pt*>(locate->Next());
            }
            return {-1,0};
        }
        
        //return the corresponding values whose keys satisfy [key1,key2)
        void Lookup(K key1,K key2,std::vector<std::pair<K,V>> &result){
            if(key2 < key1){
                return;
            }
            Segment_pt* locate = reinterpret_cast<Segment_pt*>(Scan(key1));
            while(true){
                locate->delta_inserts.Wait(IsREAD_THREAD);
                if(key1 < locate->bound){
                    locate->Lookup(key1,key2,result);
                    if(key2 <= locate->bound){
                        locate->delta_inserts.SignalWork();
                        return ;
                    } 
                }
                locate->delta_inserts.SignalWork();
                locate = reinterpret_cast<Segment_pt*>(locate->Next());
            }
            return;
        }

        //skiplist add a key 
        bool Add(K key,V value,int &path_depth,long long &scan_time,int &split_cnt,long long &collision){
            SNode* preds[MaxLevel+1];
            for(int i =0;i<=MaxLevel;i++){
                preds[i] = nullptr;
            }
            Segment_pt* locate = reinterpret_cast<Segment_pt*>(Scan(key,preds));
            while(true){    
                if(key < locate->bound){//skip when the bound of locate segment is smaller than key because of splitting
                    int state = locate->delta_inserts.Wait();
                    // locate->ExclusiveLock();
                    if(key < locate->bound){
                        return Add_key_in_Segment(key,value,preds,locate,state,path_depth,scan_time,split_cnt,collision);
                        //will finish operations on delta_inserts and work_threads
                    }
                    if(state == 0){//wake up one blocked thread to do split/rebuild 
                        locate->delta_inserts.SignalWork();
                        locate->delta_inserts.SignalOne();
                    }
                    else
                        locate->delta_inserts.SignalWork();//other case when state > 0
                                                            //can be optimistic that I inserted a repeated key, 
                                                            //and  don’t do other processes here.
                        
                }
                //skip 
                SNode* next = locate->Next();
                locate = reinterpret_cast<Segment_pt*>(next);
            }
            return false;	
        }

        //add a key-value into segment_pt
        bool Add_key_in_Segment(K key,V value,SNode** preds,Segment_pt* locate,int state,int &insert_depth,long long &scan_time,int &split_cnt,long long &collision){
            int path_size = 0;
            int flag = 0;//mark build model for tree_none
            if(state == 0){
                int ele_size = locate->DataArray->size.load(std::memory_order_acquire);
                if(ele_size <= 1){
                    locate->delta_inserts.WaitWork();//wait working threads before this,work_thread --count
                    flag = 1;
                }
            }
            locate->insert_subtree(key,value,path_size,collision);
            insert_depth = path_size;
            subtree *n = locate->DataArray;
            if(flag){
                //flag = 1,means state = 0 and the size of n <=2 after insert_subtree, DataArray was updated
                int new_size = n->size.load(std::memory_order_acquire);
                int new_slot_ = n->num_items;
                int key_space = n->stop - n->start;
                locate->SetSegMax(new_slot_,new_size, key_space);
                int new_segment_max_size = locate->SegmentMaxSize;
                locate->delta_inserts.SignalWork();
                locate->delta_inserts.SignalAll(new_segment_max_size,new_size,new_slot_);//wake up blocked threads
                return false;
            }
            const bool need_rebuild = ( state == 0);
            if(!need_rebuild){
                locate->delta_inserts.SignalWork();//working threads number --
                return false;
            }else{
                locate->delta_inserts.WaitWork();//wait working threads before this
                const int ESIZE = n->size.load(std::memory_order_acquire);
                int new_segment_max_size = locate->SegmentMaxSize;
                // int n_slots = n->num_items;
                const int REBUILD_SIZE = 64;
                const int SLOTS_ = 256;
                if(ESIZE < REBUILD_SIZE){
                    locate->delta_inserts.SignalWork();
                    locate->delta_inserts.SignalAll(new_segment_max_size,ESIZE,SLOTS_);//wake up blocked threads
                    return false;
                }
                const auto start_time = std::chrono::steady_clock::now();
                if(ESIZE >= new_segment_max_size ){//split
                    K* keys = new K[ESIZE];
                    V* values = new V[ESIZE];
                    unsigned int n_start = n->start,n_stop = n->stop;
                    locate->scan_and_destroy_subtree(n,keys,values,ESIZE);
                    locate->DataArray = nullptr;
                    SplitSegment(locate,preds,keys,values,ESIZE,n_start,n_stop);//SplitSegment already released the write lock of locate
                    int new_slot = locate->DataArray->num_items;
                    int new_size = locate->DataArray->size.load(std::memory_order_acquire);
                    new_segment_max_size = locate->SegmentMaxSize;
                    delete[] keys;
                    delete[] values;
                    split_cnt++;
                    const auto end_time = std::chrono::steady_clock::now();
                    const auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);
                    locate->delta_inserts.SignalWork();
                    locate->delta_inserts.SignalAll(new_segment_max_size,new_size,new_slot);//wake up blocked threads
                    scan_time+=duration.count();
                    return false;
                }else{
                    K* keys = new K[ESIZE];
                    V* values = new V[ESIZE];
                    unsigned int n_start = n->start,n_stop = n->stop;
                    locate->scan_and_destroy_subtree(n,keys,values,ESIZE);
                    locate->DataArray = nullptr;

                    std::pair<subtree*,int> res_ = locate->rebuild_tree(keys,values,ESIZE,n_start,n_stop);
                    locate->DataArray = res_.first;
                    int new_slot = locate->DataArray->num_items;
                    locate->SetSegMax(new_slot,ESIZE,n_stop - n_start);
                    delete[] keys;
                    delete[] values;
                    const auto end_time = std::chrono::steady_clock::now();
                    const auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);
                    // std::cout<<"rebuild:\t"<<ESIZE<<"\tkey_space:"<<n_start<<"---"<<n_stop<<std::endl;
                    new_segment_max_size = locate->SegmentMaxSize;

                    locate->delta_inserts.SignalWork();
                    locate->delta_inserts.SignalAll(new_segment_max_size,ESIZE,new_slot);//wake up blocked threads
                    scan_time+=duration.count();
                    return true;
                }
            }
            return false;
        }

        bool SplitSegment(Segment_pt *root,SNode** preds,K *_keys,int*  _values,int _size,K _start,K _stop){
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
                seg_res = plr->Process(static_cast<double>(_keys[i]) - static_cast<double>(_keys[st]), static_cast<double>(i-st));
                #else
                if(_stop - _start > 0x7fffffff)
                    seg_res = plr->Process(static_cast<double>(_keys[i]), static_cast<double>(i-st));
                else
                    seg_res = plr->Process(static_cast<double>(_keys[i]) - static_cast<double>(_keys[st]), static_cast<double>(i-st));
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
            unsigned int split_lower_bound = 0;//the right bound of the first segment node after split
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
                }

                //2.4 create segment ,index without data and link inner part
                SNode* next_seg = nullptr;//the bottom segment_pt
                vector<SNode*> SNodeArray(height_tmp+1,nullptr);
                if(i){
                    //2.4.1 create index array(stored in SNodeArray) with subtree
                    NewIndexArrayWithData(height_tmp,base,bound,new_subtree,SNodeArray);
                    next_seg = SNodeArray[0];
                    split_segments[i] = reinterpret_cast<Segment_pt*>(next_seg);//store segment in split_segments
                    reinterpret_cast<Segment_pt*>(next_seg)->delta_inserts.AddWork();//restrain split
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
                    root->SetSegMax(new_subtree->num_items,new_subtree->size.load(std::memory_order_acquire),bound - base);
                }
                  
            }

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
                    split_segments[i]->delta_inserts.SignalWork();
                }
                delete plr;
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

        long long SpaceSize(){
            long long space_sum  = 0;
            //skiplist class
            space_sum += sizeof(skiplist);
            
            SNode *x = head_;
            //get down to segment
            while(x->isIndexNode()){
                x = reinterpret_cast<Index*>(x)->Down();
            }
            //compute index and segment node size
            Segment_pt *n=nullptr;
            while(x->Key != UNINT_MAX){
                n = reinterpret_cast<Segment_pt*>(x);
                //index layer
                space_sum += (sizeof(Index) * n->level);
                //segment layer
                space_sum += (n->SpaceSize());
                x = x->Next();
            }
            //tail
            n = reinterpret_cast<Segment_pt*>(x);
            space_sum += (sizeof(Index) * n->level);
            space_sum += (n->SpaceSize());
            return space_sum;
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