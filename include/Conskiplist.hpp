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
#include <mutex>
#include <atomic>
#include <assert.h>
#include <random>

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
#define USEPLR 1
#define USEGREEDYPCCS 0
#define DEBUG 0
#define SHOWLINEAR 0
#define Gm 128
#define MAX_DEPTH 6
#define WRITESEG 0

// runtime assert
#define RT_ASSERT(expr) \
{ \
    if (!(expr)) { \
        fprintf(stderr, "RT_ASSERT Error at %s:%d, `%s`\n", __FILE__, __LINE__, #expr); \
        exit(0); \
    } \
}

using namespace std;

int collision = 0;
int rebuild_cnt = 0;
int split_cnt = 0;
int file_num = 0;
long long jmp_seg = 0;
long long jmp_subtree = 0;
long long scan_ = 0;
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

typedef struct node{
    unsigned key;
    int value;
}node;

template<class K, class V>
class skiplist {
    public:
        struct subtree;
        struct Item{
            union {
                node data;
                subtree* child;
            }comp;
        };
        struct subtree{
            int is_two = 0; // is special node for only two keys
            int build_size = 0; // tree size (include sub nodes) when node created
            int size = 0; // current tree size (include sub nodes)
            int fixed = 0; // fixed node will not trigger rebuild
            int num_inserts = 0, num_insert_to_data = 0;
            // int child_ptr = 0;
            int num_items = 0; // size of items
            double slope = 0,intercept = 0;
            unsigned int start = UNINT_MAX,stop=UNINT_MAX;
            Item* items = nullptr;
            bitmap_t* none_bitmap = nullptr; // 1 means None, 0 means Data or Child
            bitmap_t* child_bitmap = nullptr; // 1 means Child. will always be 0 when none_bitmap is 1
            
            inline int predict(unsigned int key){
                double v = this->slope * (static_cast<long double>(key)-static_cast<long double>(this->start)) + this->intercept;
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

        };
        struct Segment_pt{
            unsigned int anchor;//left guard
            unsigned int bound;//right guard
            subtree* DataArray;
            bool split;
            std::atomic<bool> lock;//false means segment is locked
            std::atomic<bool> spliting;//true means segment is locked
            std::atomic<Segment_pt*> next_;

            inline void InitLock(){
                lock.store(true, std::memory_order_release);
            }

            inline void InitSplit(){
                spliting.store(false, std::memory_order_release);
            }

            inline bool ReadLock(){
                return lock.load(std::memory_order_acquire);
            }

            inline bool ReadSplit(){
                return spliting.load(std::memory_order_acquire);
            }

            inline bool CASLock(bool free = true,bool block=false){
                return lock.compare_exchange_strong(free, block);
            }

            inline bool CASSplit(bool free = false,bool block=true){
                return spliting.compare_exchange_strong(free, block);
            }            

            inline void ReleaseLock(){
                lock.store(true, std::memory_order_release);
            }

            inline void ReleaseSplit(){
                spliting.store(false, std::memory_order_release);
            }

            Segment_pt* Next() {
                // Use an 'acquire load' so that we observe a fully initialized
                // version of the returned Node.
                return (next_.load(std::memory_order_acquire));
            }

            void SetNext(Segment_pt* x) {
                next_.store(x, std::memory_order_release);
            }

            bool CASNext(Segment_pt* expected, Segment_pt* x) {
                return next_.compare_exchange_strong(expected, x);
            }

            // No-barrier variants that can be safely used in a few locations.
            Segment_pt* NoBarrier_Next() {
                return next_.load(std::memory_order_relaxed);
            }
            
            void NoBarrier_SetNext(Segment_pt* x) {
                next_.store(x, std::memory_order_relaxed);
            }

        };

        struct Index{
            int level;
            unsigned key;
            std::atomic<Index*> *forward;
            Segment_pt *downward;
            //TODO:Index Append、Next、SetNext、newIndex
            inline void Append(Index *pred,Index *locate);
        };
        
        typedef struct {
            int begin;
            int end;
            int lvl; // top lvl = 1
            subtree* n;
        } Seg;
        
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

                // if (n->is_two) {
                // 	RT_ASSERT(n->build_size == 2);
                // 	RT_ASSERT(n->num_items == 32);
                // 	n->size = 2;
                // 	n->num_inserts = n->num_insert_to_data = 0;
                //     // n->num_items = 8;
                //     // n->none_bitmap[0] = 0xff;
                //     // n->child_bitmap[0] = 0;
                //     for(int _=0;_<4;_++){
                //         n->none_bitmap[_] = 0xff;
                //         n->child_bitmap[_] = 0;
                //     }
                // 	pending_two.push(n);
                // } else {
                    delete_items(n->items, n->num_items);
                    const int bitmap_size = BITMAP_SIZE(n->num_items);
                    delete_bitmap(n->none_bitmap, bitmap_size);
                    delete_bitmap(n->child_bitmap, bitmap_size);
                    delete_subtree(n, 1);
                // }
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
        bitmap_t* new_bitmap(int n)
        {
            bitmap_t* x = bitmap_allocator.allocate(n);
            // RT_ASSERT(x != NULL && x != (bitmap_t*)(-1));
            return x;
        }
        void delete_bitmap(bitmap_t* x, int n)
        {
            bitmap_allocator.deallocate(x, n);
        }
        
        inline int compare_(unsigned int k1,unsigned int  k2){
            if(k1<k2){
                return -1;
            }
            else if(k1 == k2){
                return 0;
            }else{
                return 1;
            }
        }

        void insert_subtree(Segment_pt* root,unsigned int key,int value,subtree* path[],int &path_size);
        
        std::pair<int,node*> find_subtree(subtree* n,unsigned int key){
            // subtree* n = this->DataArray;
            vector<subtree*> paths;
            vector<int> poss;
            while(1){
                paths.push_back(n);
                int pos = n->predict(key);
                poss.push_back(pos);
                if (BITMAP_GET(n->none_bitmap, pos) == 1) {
                    //checkTypeEqualNull
                    break;
                    // return {false,0};
                }
                else if(BITMAP_GET(n->child_bitmap, pos) == 0){
                    //checkTypeEqualData
                    return {n->items[pos].comp.data.key == key,&(n->items[pos].comp.data)};
                }
                else{
                    jmp_subtree++;
                    n = n->items[pos].comp.child;
                }
            }
            return {false,0};
        }

        std::pair<unsigned int,unsigned int> computeRange(int pos,subtree* x){
            unsigned int start = 0,stop = 0;
            start = (pos - x->intercept)/(x->slope * 1000.0)*1000.0;
            stop = (pos+1 - x->intercept)/(x->slope * 1000.0)*1000.0;
            return {start+x->start,stop+x->start};
        }

        inline int compute_gap_count(int size) {
            if (size >= 1000000) return 1;
            if (size >= 100000) return 2;
            return 5;
        }

        bool SplitSegment(Segment_pt *root,Index** preds,unsigned int *_keys,int*  _values,int _size,unsigned int _start,unsigned int _stop);

        subtree* build_tree_none();

        subtree* build_tree_two(unsigned int key1, int value1, unsigned int key2, int value2,unsigned int start=0,unsigned int stop=0,subtree* x = nullptr);

        subtree* rebuild_tree(unsigned int *_keys,int* _values,int _size,unsigned int start,unsigned int stop,bool plr = false,double top_slope = 0,double top_intercept = 0);

        std::pair<long double,int> scan_and_destroy_subtree(subtree* root,unsigned int*keys,int* values, bool destory = true);

        Segment_pt* NewSegment_pt(int level,unsigned int base,unsigned int bound,subtree* n){
            Segment_pt *newseg = new Segment_pt;
            newseg->anchor = base;
            newseg->bound = bound;
            newseg->DataArray = n;
            newseg->split = 0;
            newseg->InitLock();
            newseg->InitSplit();
            newseg->SetNext(nullptr);
            return newseg;
        }

        Segment_pt* NewSegment_pt_nodata(int level,unsigned int base,unsigned int bound){
            Segment_pt *newseg = new Segment_pt;
            newseg->anchor = base;
            newseg->bound = bound;
            newseg->DataArray = nullptr;
            newseg->split = 0;
            newseg->InitLock();
            newseg->InitSplit();
            newseg->SetNext(nullptr);
            return newseg;
        }
        
        Segment_pt* NewSegment_pt(int level,unsigned int base,unsigned int bound){
            Segment_pt *newseg = new Segment_pt;
            newseg->anchor = base;
            newseg->bound = bound;
            newseg->DataArray =  build_tree_none();
            newseg->split = 0;
            newseg->InitLock();
            newseg->InitSplit();
            newseg->SetNext(nullptr);
            return newseg;
        }
        
        /////////////////////////////////
        std::atomic<int> max_height_; 
        int MaxLevel;
		int gamma;
        int segment_max_size;
        double linearity;
        Segment_pt *head_;
        Segment_pt *tail_;
        skiplist();
        //TODO:skiplist 构造函数
        skiplist(int MaxLevel,int gamma):max_height_(1),MaxLevel(MaxLevel),gamma(gamma),segment_max_size(1e5),linearity(0.98){
            
        }

        inline int GetMaxHeight() const {
            return max_height_.load(std::memory_order_relaxed);
        }

        //for Lookup
        Segment_pt* Scan(K key);

        //for Add
        Segment_pt* Scan(K key, Index** preds);

        //skiplist lookup a key
        std::pair<int,int> Lookup(K key);

        //skiplist add a key 
        void Add(K key);

        int RandLevel(){
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
};
