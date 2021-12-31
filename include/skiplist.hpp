
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
#define LINAERITY 1.0//0.98
#define USEPLR 1
#define USEGREEDYPCCS 0
#define DEBUG 1
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
			// cerr<<"super class is called"<<endl;
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
		// cerr<<"in process:old state "<<this->state<<endl;
		// GreedyState newState = this->state;
		
		// cerr<<"new state(old)："<<newState;
		// cerr<<"GreedyState::Need2"<<GreedyState::Need2<<endl;
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
		// cerr<<"new state："<<newState;
		// this->state = newState;
		// cerr<<"plr state:"<<this->state<<endl;
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

class skiplist {
    public:
        std::atomic<int> max_height_; 
        int MaxLevel;
		int gamma;
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
            int level;
            unsigned int anchor;//left guard
            unsigned int bound;//right guard
            subtree* DataArray;
            // bool rebuild;
            // bool split;
            bool marked;
            std::atomic<bool> lock;//false means segment is locked
            std::atomic<bool> build;//true means 
            std::atomic<Segment_pt*> *next_;
            // std::atomic<Segment_pt*> next_[1];

            inline void InitLock(){
                lock.store(true, std::memory_order_release);
            }

            inline void InitBuildM(){
                build.store(false, std::memory_order_release);
            }

            inline bool ReadLock(){
                return lock.load(std::memory_order_acquire);
            }

            inline bool ReadBuildM(){
                return build.load(std::memory_order_acquire);
            }

            inline bool CASLock(bool free = true,bool block=false){
                return lock.compare_exchange_strong(free, block);
            }

            inline bool CASBuild(bool free = false,bool block=true){
                return build.compare_exchange_strong(free, block);
            }

            inline void RealseLock(){
                lock.store(true, std::memory_order_release);
            }

            inline void RealseBuild(){
                build.store(false, std::memory_order_release);
            }

            Segment_pt* Next(int n) {
                assert(n >= 0);
                // Use an 'acquire load' so that we observe a fully initialized
                // version of the returned Node.
                return (next_[n].load(std::memory_order_acquire));
            }

            void SetNext(int n, Segment_pt* x) {
                assert(n >= 0);
                next_[n].store(x, std::memory_order_release);
            }

            bool CASNext(int n, Segment_pt* expected, Segment_pt* x) {
                assert(n >= 0);
                return next_[n].compare_exchange_strong(expected, x);
            }

            // No-barrier variants that can be safely used in a few locations.
            Segment_pt* NoBarrier_Next(int n) {
                assert(n >= 0);
                return next_[n].load(std::memory_order_relaxed);
            }
            
            void NoBarrier_SetNext(int n, Segment_pt* x) {
                assert(n >= 0);
                next_[n].store(x, std::memory_order_relaxed);
            }

            long long Computespace(){
                long long space = 0;
                space+=sizeof(std::atomic<Segment_pt*>) * (level)+sizeof(Segment_pt);
                std::queue<subtree*> s;
                s.push(DataArray);
                while(!s.empty()){
                    subtree* n = s.front();
                    s.pop();
                    
                    space+=sizeof(subtree);
                    
                    if(n->size == 1){
                        space+=(sizeof(Item));
                        space+=(sizeof(bitmap_t) * 2);
                    }
                    else{
                        const int bitmap_size = BITMAP_SIZE(n->num_items);
                        //items
                        space+=(n->num_items * sizeof(Item));
                        //bitmap
                        space+=( bitmap_size * sizeof(bitmap_t) *2);
                        for(int i = 0;i<n->num_items;i++){
                            if(BITMAP_GET(n->child_bitmap,i)){
                                s.push(n->items[i].comp.child);
                            }
                        }
                    }
                }
                // space-=(DataArray->size*8);
                return space;
                
            }

        };

        struct State{
            unsigned int key;
            int value;
            Segment_pt** ppreds;
            Segment_pt** preds;
            Segment_pt** currs;
            Segment_pt** succs;

        };
        
        typedef struct {
            int begin;
            int end;
            int lvl; // top lvl = 1
            subtree* n;
        } Seg;
        
        typedef struct Splice{
            int height_ = 0;
            Segment_pt** prev_;
            Segment_pt** next_;
        }Splice;
        
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

        // void destory_pending()
        // {
        // 	while (!pending_two.empty()) {
        // 		subtree* n = pending_two.top(); pending_two.pop();
        //
        // 		delete_items(n->items, n->num_items);
        // 		const int bitmap_size = BITMAP_SIZE(n->num_items);
        // 		delete_bitmap(n->none_bitmap, bitmap_size);
        // 		delete_bitmap(n->child_bitmap, bitmap_size);
        // 		delete_subtree(n, 1);
        // 	}
        // }

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
        
        bool checkTree(subtree* n){
            int shouldbe = n->size;
            int real = 0;
            int prt = 0;
            for(int j = 0;j<n->num_items;j++){
                if(!BITMAP_GET(n->none_bitmap,j)){
                    if(BITMAP_GET(n->child_bitmap,j)){
                        if(checkTree(n->items[j].comp.child)){
                            real+=(n->items[j].comp.child->size);
                            prt++;
                        }
                        else
                            return false;
                    }
                    else{
                        real++;
                    }
                }
            }
            if(shouldbe!=real){
                cerr<<"rebuild "<<n<<"error"<<endl;
                return false;
            }
            return true;
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
            while(1){
                int pos = n->predict(key);
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

        bool SplitSegment(Segment_pt *root,State* proc_state,skiplist* list,unsigned int *_keys,int*  _values,int _size,long double pccs,int middle,unsigned int _start,unsigned int _stop);

        subtree* build_tree_none();

        subtree* build_tree_two(unsigned int key1, int value1, unsigned int key2, int value2,unsigned int start=0,unsigned int stop=0,subtree* x = nullptr);

        subtree* rebuild_tree(unsigned int *_keys,int* _values,int _size,unsigned int start,unsigned int stop,bool plr = false,double top_slope = 0,double top_intercept = 0);

        std::pair<long double,int> scan_and_destroy_subtree(subtree* root,unsigned int*keys,int* values, bool destory = true);

        // long long Computespace();

        // void show(){
        //     cerr<<this<<"\tstart:"<<this->DataArray->start<<"\tstop:"<<this->DataArray->stop<<endl;
        // }

        Segment_pt* NewSegment_pt(int level,unsigned int base,unsigned int bound,subtree* n){
            Segment_pt *newseg = new Segment_pt;
            newseg->level = level;
            newseg->anchor = base;
            newseg->bound = bound;
            newseg->DataArray = n;
            newseg->marked = 0;
            // newseg->rebuild = 0;
            // newseg->split = 0;
            newseg->InitLock();
            newseg->InitBuildM();
            newseg->next_ = new std::atomic<Segment_pt*>[level];
            return newseg;
        }

        Segment_pt* NewSegment_pt_nodata(int level,unsigned int base,unsigned int bound){
            Segment_pt *newseg = new Segment_pt;
            newseg->level = level;
            newseg->anchor = base;
            newseg->bound = bound;
            newseg->DataArray = nullptr;
            newseg->marked = 0;
            // newseg->rebuild = 0;
            // newseg->split = 0;
            newseg->InitLock();
            newseg->InitBuildM();
            newseg->next_ = new std::atomic<Segment_pt*>[level];
            return newseg;
        }
        
        Segment_pt* NewSegment_pt(int level,unsigned int base,unsigned int bound){
            Segment_pt *newseg = new Segment_pt;
            newseg->level = level;
            newseg->anchor = base;
            newseg->bound = bound;
            newseg->DataArray =  build_tree_none();
            newseg->marked = 0;
            // newseg->rebuild = 0;
            // newseg->split = 0;
            newseg->InitLock();
            newseg->InitBuildM();
            newseg->next_ = new std::atomic<Segment_pt*>[level];
            return newseg;
        }
        
        // void deleteSegment_pt(){
        //     destroy_subtree(DataArray);
        //     DataArray = NULL;
        //     // destory_pending();
        // }

        /////////////////////////////////
        
        // int segCnt;
        Segment_pt *head_;
        Segment_pt *tail_;
        skiplist();
        skiplist(int MaxLevel,int gamma):max_height_(1),MaxLevel(MaxLevel),gamma(gamma){
            head_ = NewSegment_pt(MaxLevel,0,0);
            head_->DataArray->start = 0;
            head_->DataArray->stop = 0;
            tail_ = NewSegment_pt(MaxLevel,UNINT_MAX,UNINT_MAX);
            tail_->DataArray->start = UNINT_MAX;
            tail_->DataArray->stop = UNINT_MAX;
            Segment_pt* newSeg = NewSegment_pt(this->RandLevel(),1,UNINT_MAX);//new Segment_pt(this->RandLevel());
            int max_height = 1;
            while (newSeg->level > max_height) {
                if (max_height_.compare_exchange_weak(max_height, newSeg->level)) {
                // successfully updated it
                max_height = newSeg->level;
                break;
                }
            }
            newSeg->DataArray->start = 1;
            newSeg->DataArray->stop = UNINT_MAX;
            for(int i = 0;i<newSeg->level;i++){
                head_->SetNext(i,newSeg);
                newSeg->SetNext(i,tail_);
                tail_->SetNext(i,nullptr);
            }  
            for(int i = newSeg->level;i<=MaxLevel;i++){
                head_->SetNext(i, tail_);
                tail_->SetNext(i,nullptr);
            }
            
        }

        inline bool KeyIsAfterNodeStart(unsigned int key, Segment_pt* n) const {
            // nullptr n is considered infinite
            assert(n != head_);
            return (n != nullptr) && (n->DataArray->start <= key);
        }

        inline bool KeyIsntAfterNode(unsigned int key, Segment_pt* n) const {
            // nullptr n is considered infinite
            assert(n != head_);
            return (n != nullptr) && (n->DataArray->stop > key);
        }

        inline bool KeyIsAfterNode(unsigned int key, Segment_pt* n) const {
            // nullptr n is considered infinite
            assert(n != head_);
            return (n != nullptr) && (n->DataArray->stop < key);
        }

        inline int GetMaxHeight() const {
            return max_height_.load(std::memory_order_relaxed);
        }

        Segment_pt* Scan(State* proc_state);

        std::pair<int,int> Seek(State* proc_state);

        std::pair<int,int> Lookup(State* proc_state);

		// std::pair<int,node*> Search(unsigned int key);

        int Add(State* proc_state,Segment_pt* locate);

		// bool Insert(unsigned int key,int value);

        bool Insert(State* proc_state);

		void insert_static(vector<Segment_pt*> &Update,Segment* seg,unsigned int st,unsigned int ed,node* input,int size,int level);

        bool LockAndValidateUntil(Segment_pt* from[],Segment_pt* to[],int level);

        bool LockAndValidateUntil(Segment_pt* from[],Segment_pt* to[],int from_level,int to_level);

        void UnlockUntil(Segment_pt* from[],int level);

        void UnlockUntil(Segment_pt* from[],int from_level,int to_level);

		void show();

        // void ShowNodeDis();
        // void Dump();

		// void ComputeSpace();
		        
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
        // int Delete(skiplist *list,int k);
        // int Insert(skiplist *list,int k,int v);

        void FindSpliceForLevel(unsigned int key,Segment_pt* before, Segment_pt* after,int level, Segment_pt** out_prev,Segment_pt** out_next,Segment_pt **locate) {
            while (true) {
                Segment_pt* next = before->Next(level);
                if (next != nullptr) {
                    PREFETCH(next->Next(level), 0, 1);
                }
                if (next != nullptr && level>0) {
                    PREFETCH(next->Next(level-1), 0, 1);
                }
                
                assert(before == head_ || next == nullptr ||KeyIsAfterNode(next->DataArray->stop, before));
                assert(before == head_ || KeyIsAfterNode(key, before));
                //key < next
                if (next == after || KeyIsntAfterNode(key, next)) {
                    // found it
                    *out_prev = before;
                    if(next != after && KeyIsAfterNodeStart(key,next)){
                        if(!(*locate))
                        (*locate) = next;
                        *out_next = next->Next(level);
                    }else{
                        *out_next = next;
                    }
                    return;
                }
                before = next;
            }
        }

        void RecomputeSpliceLevels(unsigned int key,Splice* splice,int recompute_level,Segment_pt **locate) {
            assert(recompute_level > 0);
            assert(recompute_level <= splice->height_);
            for (int i = recompute_level - 1; i >= 0; --i) {
                FindSpliceForLevel(key, splice->prev_[i + 1], splice->next_[i + 1], i,
                                &splice->prev_[i], &splice->next_[i],locate);
            }
        }
};


bool skiplist::SplitSegment(Segment_pt *root,State* proc_state,skiplist* list,unsigned int *_keys,int*  _values,int _size,long double pccs,int middle,unsigned int _start,unsigned int _stop){
    // cerr<<"split"<<endl;
#if USEPLR
    int st = 0,ed =0;//the rank of the first/last+1 ele of one segment
    GreedyPLR* plr = new GreedyPLR(Gm);
    vector<Segment*> seg;
    vector<int> segment_stIndex;
    vector<int> height_seq;
    int max_height_cur = list->GetMaxHeight();
    for (int i = 0; i < _size; i++) {
		Segment* seg_res = nullptr;
		seg_res = plr->Process(static_cast<double>(_keys[i])-static_cast<double>(_keys[st]), static_cast<double>(i-st));
		if(seg_res) {
            segment_stIndex.push_back(st);
            seg.push_back(seg_res);
            if(st != 0){
                int h = RandLevel();
                max_height_cur = max(max_height_cur,h);
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
            max_height_cur = max(max_height_cur,h);
        }
        
	}
    int seg_size = static_cast<int>(seg.size());
    cerr<<"split seg size:"<<seg_size<<endl;
    if(seg_size > 1){
        bool result = LockAndValidateUntil(proc_state->preds,proc_state->succs,root->level,max_height_cur);//root->level-1,
        if(!result){
            //memory release
            return false;
        }
    }
    //如果需要split，但是acquire失败，则进行rebuild
    // cerr<<"seg size:"<<seg.size()<<endl;
#if WRITESEG
    string kk = "\nseg size "+to_string(seg.size())+":\n";
    write_into_file("./split.txt",kk.c_str());
#endif
    vector<Segment_pt*> split_segments(seg_size);
    split_segments[0] = root;
    int max_seg_height = root->level;
    Segment_pt** preds = (Segment_pt**)(malloc(sizeof(Segment_pt*)*MaxLevel));
    memcpy(preds,proc_state->preds,sizeof(Segment_pt*)*MaxLevel);
    for(int i = 1;i<seg_size;i++){
       
        int height_tmp = height_seq[i-1];//RandLevel();
        unsigned int base = _keys[segment_stIndex[i]];
        unsigned int bound;
         if(i==1){
            root->bound = base;
        }
        if(i != seg_size){
            bound = _stop;
        }else{
            bound = _keys[segment_stIndex[i+1]];
        }
        Segment_pt* next_seg = NewSegment_pt_nodata(height_tmp,base,bound);
        split_segments[i] = next_seg;
        next_seg->CASBuild();
        // Acquire(next_seg);
        assert(next_seg->CASLock());
        max_seg_height = max(max_seg_height,height_tmp);
        for(int l = height_tmp-1;l>=0;l--){
            next_seg->NoBarrier_SetNext(l,preds[l]->NoBarrier_Next(l));
            preds[l]->NoBarrier_SetNext(l,next_seg);
            preds[l] = next_seg;
            /*
            if (proc_state->preds[l]->CASNext(l, proc_state->succs[l], next_seg)) {
                cerr<<" pre:"<<proc_state->preds[l]<<endl;
                proc_state->preds[l] = next_seg;
            }else{
                cerr<<"failed"<<endl;
                //delete next seg
                for(int j = i;j>=1;j--){
                    split_segments[j]->RealseLock();
                }
                return false;
            }
            */
        }

    }
    
    for(int i = 0;i<seg_size;i++){
        cerr<<split_segments[i]<<endl;
        for(int j = 0;j<split_segments[i]->level;j++){
            cerr<<"level "<<j<<"next:"<<split_segments[i]->NoBarrier_Next(j)<<endl;
        }
    }

    int max_height = list->max_height_.load(std::memory_order_relaxed);
    while (max_seg_height > max_height) {
        if (list->max_height_.compare_exchange_weak(max_height, max_seg_height)) {
            // successfully updated it
            max_height = max_seg_height;
            break;
        }
        // else retry, possibly exiting the loop because somebody else
        // increased it
    }

    for(int i=0;i<seg_size;i++){
        st = segment_stIndex[i];
        unsigned int stop = 0;// interval left value of one segment
        unsigned int start_local = 0;//interval right value of one segment
        if(st == 0){
            start_local = _start;//the parent's interval left
        }else{
            start_local = _keys[st];
        }
        if(i == seg_size-1){
            ed = _size;
            stop = _stop;//the parent's interval right
        }
        else{
            ed = segment_stIndex[i+1];
            stop = _keys[ed];//the interval right value is the interval left value of the next segment
        }
#if WRITESEG
        string kk = to_string(ed-st)+",";
        write_into_file("./split.txt",kk.c_str());
#endif
        RT_ASSERT(ed-st>0);
        subtree* new_subtree = nullptr;
        if(ed-st == 1){
            new_subtree = build_tree_none();
            new_subtree->start = start_local;
            new_subtree->stop = stop;
            BITMAP_CLEAR(new_subtree->none_bitmap,0);
            new_subtree->size++;
            new_subtree->num_inserts++;
            new_subtree->items[0].comp.data.key = _keys[st];
            new_subtree->items[0].comp.data.value = _values[st];
        }
        else{
            if(seg[i]->slope > 0)
                new_subtree = rebuild_tree(_keys+st,_values+st,ed-st,start_local,stop,true,seg[i]->slope,seg[i]->intercept);
            else
                new_subtree = rebuild_tree(_keys+st,_values+st,ed-st,start_local,stop);//,true,seg[i]->slope,seg[i]->intercept);
        }
        split_segments[i]->anchor = start_local;
        split_segments[i]->DataArray = new_subtree;
    }
    for(int i = 1;i<seg_size;i++){
        // Release(split_segments[i]);
        split_segments[i]->RealseBuild();
        split_segments[i]->RealseLock();
    }
    UnlockUntil(proc_state->preds,max_height_cur);
    cerr<<"split end"<<endl;
    return true;
#elif USEGREEDYPCCS
    long double sum_of_key = 0;//sigma(x-x0)
    long double sum_of_y = 0;//(static_cast<long double>( _root->size)-1)* _root->size/2;//simga(y)
    long double E_Y = 0;//( static_cast<long double>(_root->size)-1)/2.0;//E(Y)
    long double sum_of_keyMuly = 0;//sigma((x-x0)y)
    long double sum_of_key2 = 0;//sigma((x-x0)^2)
    long double sum_of_y2 = 0;//(static_cast<long double>(_root->size)-1) * _root->size * (2*_root->size - 1)/6;//sigma(y^2)
    long double E_key = 0;//sum_of_key/(static_cast<long double>(keylen));
    
    int st_idx = 0;
    int split_ = 0;
    long double  key_0 = _keys[st_idx];//x0
    unsigned int start_local = 0;
    unsigned int pre_stop= 0;
    sum_of_key = static_cast<long double>(_keys[1])-key_0;
    sum_of_keyMuly = static_cast<long double>(_keys[1])-key_0;
    sum_of_key2 = (static_cast<long double>(_keys[1])-key_0)*(static_cast<long double>(_keys[1])-key_0);
    // int split_cnt = 0;
    for(int i = 2; i < _size; i++){
        // if(i % 500 == 0 || i >= _size-2){
            long double ele_cnt = i-st_idx;
            sum_of_y = (ele_cnt* (ele_cnt-1))/2;//simga(y)
            E_Y = (ele_cnt-1)/2.0;
            sum_of_y2 = (ele_cnt-1) * ele_cnt * (2*ele_cnt - 1)/6;//sigma(y^2)
            RT_ASSERT(ele_cnt>0);
            E_key = sum_of_key/(ele_cnt);
            long double st_x = sum_of_key2 - sum_of_key*E_key;
            long double st_y = sum_of_y2 - sum_of_y*E_Y;
            RT_ASSERT(st_x>1e-6);
            RT_ASSERT(st_y>1e-6);
            long double PCCs = sum_of_keyMuly - sum_of_key*E_Y;
            PCCs/=sqrt(st_x);
            PCCs/=sqrt(st_y);

            if(abs(PCCs) < LINAERITY || i >= _size-2 ){//|| ele_cnt>1e5
                split_++;
                #if WRITESEG
                cerr<<"split pccs:"<<PCCs<<endl;
                #endif
                subtree* new_subtree = nullptr;
                if(st_idx == 0){
                    start_local = _start;
                }else{
                    start_local = _keys[st_idx];
                    RT_ASSERT(start_local == pre_stop);
                }
                if(i >= _size -2){
                    #if WRITESEG
                    string kk = to_string(_size-st_idx)+",";
                    write_into_file("./split.txt",kk.c_str());
                    #endif
                    new_subtree = rebuild_tree(_keys+st_idx,_values+st_idx,_size-st_idx,start_local,_stop);
                }else{
                    #if WRITESEG
                    string kk = to_string(i-st_idx)+",";
                    write_into_file("./split.txt",kk.c_str());
                    #endif
                    new_subtree = rebuild_tree(_keys+st_idx,_values+st_idx,i-st_idx,start_local,_keys[i]);
                    pre_stop = _keys[i];
                }
                if(st_idx == 0){
                    root->DataArray = new_subtree;
                    // for(int l = root->level;l>=1;l--){
                    //     preds[l] = root;
                    // }
                }
                else{
                    //Update 应保存 每层都在newSeg左边的segment
                    Segment_pt* newSeg = new Segment_pt(list->RandLevel(),new_subtree);
                    if(newSeg->level > list->max_height_){
                        // for(int l = newSeg->level;l>list->max_height_;l--){
                        //     Update[l] = list->head_;
                        // }
                        list->max_height_ = newSeg->level;
                    }
                    for(int l = newSeg->level;l>=1;l--){
                        preds[l]->forward[l] = newSeg;
                        newSeg->forward[l] = succs[l];
                        preds[l] = newSeg;
                    }
                    newSeg->fullyLinked = true;
                    list->segCnt++;
                }
                if(i >= _size -2)
                    break;
                st_idx = i;
                key_0 = _keys[st_idx];
                sum_of_key = static_cast<long double>(_keys[st_idx+1])-key_0;
                sum_of_keyMuly = (static_cast<long double>(_keys[st_idx+1])-key_0);
                sum_of_key2 = (static_cast<long double>(_keys[st_idx+1])-key_0)*(static_cast<long double>(_keys[st_idx+1])-key_0);
                i++;
                continue;
            }
        // }
        long double xi = static_cast<long double>(_keys[i])-key_0;
        sum_of_key+=(xi);
        sum_of_keyMuly+=(xi*static_cast<long double>(i-st_idx));
        sum_of_key2+=(xi*xi);
    }
    // if(st_idx != _size-1){
    //     // split_cnt++;
    //     split_++;
    //     subtree*new_subtree = rebuild_tree(_keys+st_idx,_values+st_idx,_size-st_idx,_keys[st_idx],_keys[_size-1]+1);
    //     if(st_idx == 0){
    //         root->DataArray = new_subtree;
    //         //for(int l = root->level;l>=1;l--){
    //         //    Update[l] = root;
    //         //}
    //     }
    //     else{
    //         //Update 应保存 每层都在newSeg左边的segment
    //         Segment_pt* newSeg = new Segment_pt(list->RandLevel(),new_subtree);
    //         if(newSeg->level > list->max_height_){
    //             for(int l = newSeg->level;l>list->max_height_;l--){
    //                 Update[l] = list->head_;
    //             }
    //             list->max_height_ = newSeg->level;
    //         }
    //         for(int l = newSeg->level;l>=1;l--){
    //             newSeg->forward[l] = Update[l]->forward[l];
    //             Update[l]->forward[l] = newSeg;
    //             Update[l] = newSeg;
    //         }
    //         list->segCnt++;
    //     }
    // }

    // cerr<<"split_:"<<split_<<endl;
#else
    int right_part = _size - middle;
    subtree* new_subtree1 = nullptr;
    subtree* new_subtree2 = nullptr;
    #if WRITESEG
    string kk = to_string(middle)+","+to_string(right_part);
    write_into_file("./split.txt",kk.c_str());
    #endif
   if(middle >= right_part){
        if(middle <= _size-2){
            new_subtree1 = rebuild_tree(_keys,_values,middle,_start,_keys[middle]);
            new_subtree2  = rebuild_tree(_keys+middle,_values+middle,right_part,_keys[middle],_stop);
        }
        else if(middle == _size -1){
            new_subtree1 = rebuild_tree(_keys,_values,middle,_start,_keys[middle]);
            new_subtree2  = build_tree_none();
            BITMAP_CLEAR(new_subtree2->none_bitmap,0);
            new_subtree2->size++;
            new_subtree2->num_inserts++;
            new_subtree2->items[0].comp.data.key =  _keys[middle];
            new_subtree2->items[0].comp.data.value = _values[middle];
            new_subtree2->start =  _keys[middle];
            new_subtree2->stop = _stop;
        }
        else{
            new_subtree1 = rebuild_tree(_keys,_values,middle,_start,_stop);
            root->DataArray = new_subtree1;
            // cerr<<"split end"<<endl;
            return ;
        }
#if DEBUG
        cerr<<"replace "<<root->DataArray<<"with "<<new_subtree1<<",insert a new top subtree "<<new_subtree2<<endl;
#endif        
   }
   else{
        if(middle >= 2){//right_part <= _size-2
            new_subtree1 = rebuild_tree(_keys,_values,middle,_start,_keys[middle]);
            new_subtree2  = rebuild_tree(_keys+middle,_values+middle,right_part,_keys[middle],_stop);
        }
        else if(middle == 1){
            new_subtree1 = build_tree_none();
            new_subtree1->start = _start;
            new_subtree1->stop = _keys[1];
            BITMAP_CLEAR(new_subtree1->none_bitmap,0);
            new_subtree1->size++;
            new_subtree1->num_inserts++;
            new_subtree1->items[0].comp.data.key = _keys[0];
            new_subtree1->items[0].comp.data.value = _values[0];

            new_subtree2  = rebuild_tree(_keys+1,_values+1,right_part,_keys[1],_stop);   
        }
        else{
            new_subtree1 = rebuild_tree(_keys,_values,right_part,_start,_stop);
            root->DataArray = new_subtree1;
#if DEBUG
            cerr<<"split end"<<endl;
#endif
            return ;
        }
#if DEBUG
        cerr<<"replace "<<root->DataArray<<"with "<<new_subtree1<<",insert a new top subtree "<<new_subtree2<<endl;root->DataArray = new_subtree1;
#endif
   }
    root->DataArray = new_subtree1;
    Segment_pt* newSeg = new Segment_pt(list->RandLevel(),new_subtree2);
    if(newSeg->level > list->max_height_){
        list->max_height_ = newSeg->level;
    }
    for(int l = newSeg->level;l>=1;l--){
        preds[l]->forward[l] = newSeg;
        newSeg->forward[l] = succs[l];
        preds[l] = newSeg;
    }
    newSeg->fullyLinked = true;
    list->segCnt++;
#endif  
    cerr<<"split end"<<endl;
    return true;
}

skiplist::subtree* skiplist::build_tree_none(){
    subtree* n = new_subtree(1);
    n->is_two = 0;
    n->build_size = 0;
    n->size = 0;
    n->fixed = 0;
    // n->child_ptr = 0;
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

skiplist::subtree* skiplist::build_tree_two(unsigned int key1, int value1, unsigned int key2, int value2,unsigned int start,unsigned int stop,subtree* x)
{
    if (key1 > key2) {
        std::swap(key1, key2);
        std::swap(value1, value2);
    }
    // RT_ASSERT(key1 < key2);
    // static_assert(BITMAP_WIDTH == 8);

    subtree* n = NULL;
        n = new_subtree(1);
        n->is_two = 1;
        n->build_size = 2;
        n->size = 2;
        n->fixed = 0;
        // n->child_ptr = 0;
        n->num_inserts = n->num_insert_to_data = 0;
        n->num_items = 32;
        n->items = new_items(n->num_items);
        memset(n->items,0,n->num_items);
        n->none_bitmap = new_bitmap(4);
        n->child_bitmap = new_bitmap(4);
        // n->none_bitmap[0] = 0xff;
        // n->child_bitmap[0] = 0;
        for(int i = 0;i<4;i++){
            n->none_bitmap[i] = 0xff;
            n->child_bitmap[i] = 0;
        }

    if(!x){
        n->start = start;
        n->stop = stop;  
    }else{
    //从size=1增加元素重新生成subtree,[start,stop)
        n->start = x->start;//max(x->start,2*key1-key2);
        n->stop = x->stop;//min(x->stop,2*key2-key1);
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
        int pos = n->predict(key1);
        // RT_ASSERT(BITMAP_GET(n->none_bitmap, pos) == 1);
        BITMAP_CLEAR(n->none_bitmap, pos);
        n->items[pos].comp.data.key = key1;
        n->items[pos].comp.data.value = value1;
    }
    { // insert key2&value2
        int pos = n->predict(key2);
        // RT_ASSERT(BITMAP_GET(n->none_bitmap, pos) == 1);
        BITMAP_CLEAR(n->none_bitmap, pos);
        n->items[pos].comp.data.key = key2;
        n->items[pos].comp.data.value = value2;
    }
    RT_ASSERT(n!=NULL);

    return n;
}

skiplist::subtree* skiplist::rebuild_tree(unsigned int *_keys,int*  _values,int _size,unsigned int start,unsigned int stop,bool plr,double top_slope,double top_intercept){
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
#if DEBUG
    cerr<<"rebuild tree new subtree=ret"<<ret<<endl;
#endif
    s.push({0, _size, 1, ret});
    while(!s.empty()){
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
                    (static_cast<long double>(keys[mid1_pos])-static_cast<long double>(n->start) + static_cast<long double>(keys[mid1_pos + 1])-static_cast<long double>(n->start) ) / 2;
            const long double mid2_key =
                    (static_cast<long double>(keys[mid2_pos]) -static_cast<long double>(n->start)+ static_cast<long double>(keys[mid2_pos + 1]) -static_cast<long double>(n->start)) / 2;

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
            
#if DEBUG           
            if(n->slope < 0){
                cerr<<"slope < 0"<<endl;
            }
#else
            RT_ASSERT(n->slope >= 0);
#endif
            //fix size?
            if (size > 1e5) {
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
            int item_i = n->predict(keys[0]);
            for (int offset = 0; offset < size; ) {
                //item_i 表示当前元素应插入的位置
                //next_i 表示下一个元素应插入的位置
                int next = offset + 1, next_i = -1;
                while (next < size) {     
                    next_i = n->predict( keys[next]);
                 
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
                    // n->child_ptr++;
                    n->items[item_i].comp.child = new_subtree(1);
                    
                    // n->items[item_i].comp.child->child_ptr = 0;
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
    RT_ASSERT(ret->start == start);
    RT_ASSERT(ret->stop == stop);
    return ret;

}

std::pair<long double,int> skiplist::scan_and_destroy_subtree(subtree* _root,unsigned int*keys,int* values, bool destory){
    typedef std::pair<int, subtree*> Seg; // <begin, subtree*>
    std::stack<Seg> s;
    s.push(Seg(0, _root));
    int keylen = _root->size;
    long double sum_of_key = 0;//sigma(x-x0)
    long double sum_of_y = (static_cast<long double>( _root->size)-1)* _root->size/2;//simga(y)
    long double E_Y = ( static_cast<long double>(_root->size)-1)/2.0;//E(Y)
    long double sum_of_keyMuly = 0;//sigma((x-x0)y)
    long double sum_of_key2 = 0;//sigma((x-x0)^2)
    long double sum_of_y2 = (static_cast<long double>(_root->size)-1) * _root->size * (2*_root->size - 1)/6;//sigma(y^2)
    long double  key_0 = -1;//x0
    int half_stone = _root->num_items/2;
    int mark = -1;
    while (!s.empty()) {
        int begin = s.top().first;
        subtree* n = s.top().second;
        const int SHOULD_END_POS = begin + n->size;
        s.pop();
        for (int i = 0; i < n->num_items; i ++) {
            if(mark == -1 && i>=half_stone){
                mark = begin;
            }
            if (BITMAP_GET(n->none_bitmap, i) == 0) {
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
            // if (n->is_two) {
            //     RT_ASSERT(n->build_size == 2);
            //     RT_ASSERT(n->num_items == 32);
            //     n->size = 2;
            //     n->num_inserts = n->num_insert_to_data = 0;
            //     // n->none_bitmap[0] = 0xff;
            //     // n->child_bitmap[0] = 0;
            //     for(int _=0;_<4;_++){
            //         n->none_bitmap[_] = 0xff;
            //         n->child_bitmap[_] = 0;
            //     }
            //     pending_two.push(n);
            // } 
            // else {
                delete_items(n->items, n->num_items);
                const int bitmap_size = BITMAP_SIZE(n->num_items);
                delete_bitmap(n->none_bitmap, bitmap_size);
                delete_bitmap(n->child_bitmap, bitmap_size);
                delete_subtree(n, 1);
            // }
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
        // cerr<<"PCCs"<<PCCs<<endl;
        // if(abs(PCCs) > LINAERITY){
        //     return {false,mark};
        // }else{
        //     // cerr<<"below threshold PCCs"<<PCCs<<endl;
        //     return {true,mark};
        // }
        return {PCCs,mark};
    }
    return {1,mark};
}

void skiplist::insert_subtree(Segment_pt* root,unsigned int key,int value,subtree* path[],int &path_size){
    int insert_to_data = 0;
    subtree* n = root->DataArray;
    while(1){
        path[path_size] = n;
        path_size++;
        // n->stop = max(n->stop,key);
        int pos = n->predict(key);
        if (BITMAP_GET(n->none_bitmap, pos) == 1) {
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
            if(path_size == 1 && n->size == 1){
                //从tree_one到tree_two，将tree_one作为参数传入后限制start、stop范围
                subtree* newtree =  build_tree_two(key, value, n->items[pos].comp.data.key, n->items[pos].comp.data.value,0,0,n);
                memcpy(root->DataArray, newtree, sizeof(subtree));
                path[0] = root->DataArray;
                delete_subtree(newtree, 1);
                break;
            }
            std::pair<unsigned int,unsigned int>line = computeRange(pos,n);
            BITMAP_SET(n->child_bitmap, pos);
            // n->child_ptr++;
            //产生非顶层subtree，此情况下star stop考虑其父tree
            if(n->is_two){
                n->is_two = 0;
            }
            n->size++;
            n->num_inserts++;
            n->items[pos].comp.child  = build_tree_two(key, value, n->items[pos].comp.data.key, n->items[pos].comp.data.value,min(line.first,min(key,n->items[pos].comp.data.key)),line.second);
            insert_to_data = 1;
            path[path_size] = n->items[pos].comp.child;
            path_size++;
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

}

/////////////////////////////////////////////

skiplist::Segment_pt* skiplist::Scan(skiplist::State* proc_state){
    Segment_pt* ppred = nullptr;
    Segment_pt* pred = head_;
    Segment_pt* curr = nullptr;
    Segment_pt* succ = nullptr;
    Segment_pt* locate = nullptr;
    int max_height = max_height_.load(memory_order_relaxed);
    for(int l = MaxLevel-1;l>=max_height;l--){
        proc_state->preds[l] = head_;
        proc_state->succs[l] = tail_;
    }
    for(int l = max_height-1;l>=0;l--){
        //not sure key > curr's anchor
        curr = pred->Next(l);
        succ = curr->Next(l);   
        while(succ && succ->anchor <= proc_state->key){
            ppred = pred;
            pred = curr;
            curr = succ;
            succ = succ->Next(l);
        }
        if(proc_state->key >= curr->bound){
            //key is bwtween curr and succ
            ppred = pred;
            pred = curr;//curr != tail_
        }else{
            //key is before or located in curr
            if(!locate && proc_state->key >= curr->anchor){
                //key is located in curr
                locate = curr;
                for(;l>=0;l--){
                    proc_state->preds[l] = curr;
                    proc_state->succs[l] = curr->Next(l);
                }
                return locate;
            }
            //key is rebore curr
            succ = curr;
        }
        if(pred->NoBarrier_Next(l) != succ){
            cerr<<"no";
        }
        proc_state->preds[l] = pred;
        proc_state->succs[l] = succ;
    }
    return locate;
}

std::pair<int,int> skiplist::Seek(State* proc_state){
    std::pair<int,node*> res =  find_subtree(proc_state->currs[0]->DataArray,proc_state->key);
    if(res.first){
        return {res.first,res.second->value};
    }else{
        return {0,0};
    }
}

bool skiplist::LockAndValidateUntil(Segment_pt* from[],Segment_pt* to[],int level){
    Segment_pt* pre_ = from[0];
    if(!pre_->CASLock()){
        return false;
    }
    for(int l = 0;l<level;l++){
        // Acquire(from[l]->lock);
        ///如果该segment被锁，当前process是应该挂起还是循环
        if(!from[l]){
            cerr<<"why empty"<<endl;
        }
        if((pre_ != from[l] || from[l] != head_) && !from[l]->CASLock())
            //避免重复锁
            return false;
        bool result = !(from[l]->marked) && from[l]->Next(l) == to[l];
        if(!result){
            //from[l]--from[0] get lock
            UnlockUntil(from,l+1);
            return false;
        }      
        pre_ = from[l]; 
    }
    return true;
}

bool skiplist::LockAndValidateUntil(Segment_pt* from[],Segment_pt* to[],int from_level,int to_level){
    if(from_level < to_level){
        Segment_pt* pre_ = from[from_level];
        if(!pre_->CASLock()){
            return false;
        }
        for(int l = from_level;l<to_level;l++){
            // Acquire(from[l]->lock);
            ///如果该segment被锁，当前process是应该挂起还是循环
            if(pre_ != from[l]  && !from[l]->CASLock())
                //避免重复锁
                return false;
            bool result = !(from[l]->marked) && from[l]->Next(l) == to[l];
            if(!result){
                //from[l]--from[from_level] get lock
                UnlockUntil(from,from_level,l+1);
                return false;
            }      
            pre_ = from[l]; 
        }
    }
    for(int l = 0;l<from_level;l++){
        bool result = !(from[l]->marked) && from[l]->Next(l) == to[l];
        if(!result)
            return false;
    }
    return true;
}

void skiplist::UnlockUntil(Segment_pt* from[],int level){
    //release[0,level)
    for(int l = level-1;l>=0;l--){
        // Release(from[l]->lock);
        from[l]->RealseLock();
    }
}

void skiplist::UnlockUntil(Segment_pt* from[],int from_level,int to_level){
    //release[from_level,to_level])
    for(int l = from_level;l<to_level;l++){
        from[l]->RealseLock();
    }
}

std::pair<int,int> skiplist::Lookup(State* proc_state){
    Segment_pt* locate = Scan(proc_state);
    if(locate == nullptr){
        return {0,0};
    }
    // if(proc_state->currs[0]->anchor > proc_state->key || proc_state->currs[0]->DataArray->stop >= proc_state->key ){
    //     return {0,0};
    // }
    std::pair<int,node*> res = find_subtree(locate->DataArray,proc_state->key);
    if(res.first){
        return {res.first,res.second->value};
    }else{
        return {0,0};
    }
}

// std::pair<int,node*> skiplist::Search(unsigned int key){
//     // Segment_pt* x = this->head_;
//     // // unsigned int pred;
//     // for(int i = this->max_height_;i>=1;i--){
//     //     while(key >= x->forward[i]->DataArray->stop){
//     //         x = x->forward[i];
//     //         jmp_seg++;
//     //     }
//     //     if(key >= x->forward[i]->DataArray->start){
//     //         x = x->forward[i];
//     //         std::pair<int,node*> res = find_subtree(x->DataArray,key);
//     //         //this->binarySearch(x,key);
//     //        return res;
//     //     }
//     //     jmp_seg++;
//     // }    
//     // return {0,0};
//     int max_height = this->GetMaxHeight();
//     Segment_pt* before = this->head_;
//     for(int level = max_height;level>=0;level--){
//         while (true) {
//             Segment_pt* next = before->Next(level);
//             if (next != nullptr) {
//                 PREFETCH(next->Next(level), 0, 1);
//             }
//             if (next != nullptr && level>0) {
//                 PREFETCH(next->Next(level-1), 0, 1);
//             }
//             // if(!(before == head_ || next == nullptr ||next->DataArray->start >= before->DataArray->stop)){
//             //     cerr<<"look"<<endl;
//             // }
//             assert(before == head_ || next == nullptr ||next->DataArray->start >= before->DataArray->stop);
//            
//             // if(!(before == head_ || key >= before->DataArray->stop)){
//             //     cerr<<"aa"<<endl;
//             // }
//             assert(before == head_ || key >= before->DataArray->stop);
//
//             //before-->next key may bwtween before and next ,may in next
//             if (next == nullptr || KeyIsntAfterNode(key, next)) {
//                 // found it
//                 if(next != nullptr && KeyIsAfterNodeStart(key,next)){
//                     // locate = next;
//                     return find_subtree(next->DataArray,key);
//                 }
//                 break;
//             }
//             before = next;
//         }
//     }
//     return {0,nullptr};
//
// }

// bool skiplist::Insert(unsigned int key,int value){
//     Segment_pt* prev[MaxLevel+1];
//     Segment_pt* next[MaxLevel+1];
//     Splice splice;
//     splice.prev_ = prev;
//     splice.next_ = next;
//     int maxlvl = GetMaxHeight();
//     splice.height_ = maxlvl;
//     for(int i=maxlvl;i<=MaxLevel;i++){
//         splice.prev_[i] = head_;
//         splice.next_[i] = nullptr;
//     }
//
//     Segment_pt* locate = nullptr;
//     while(!locate){
//         //可以优化
//         RecomputeSpliceLevels(key, &splice, maxlvl,&locate);
//     }
//     subtree* path[MAX_DEPTH];
//     int path_size = 0;
//     insert_subtree(locate,key,value,path,path_size);
//     //lock locate
//     for(int i = 0;i<path_size;i++){
//         subtree *n = path[i];
//         const bool need_rebuild =  path_size >= MAX_DEPTH ;
//         if(need_rebuild){
//             const int ESIZE = n->size;
//             unsigned int* keys = new unsigned int[ESIZE];
//             int* values = new int[ESIZE];
//             unsigned int n_start = n->start,n_stop = n->stop;
//             memset(keys,0,ESIZE);
//             memset(values,0,ESIZE);
//             std::pair<long double,int> res = scan_and_destroy_subtree(n,keys,values);
//             scan_+=ESIZE;
//             // cerr<<"lin:"<<res.first<<"\t"<<ESIZE<<endl;
//             if(!i && ((res.first < LINAERITY )|| ESIZE >1e6)){//|| ESIZE > 1e6
//                 split_cnt++;
// #if WRITESEG
//                 string kk = "\nheight:"+to_string(path_size)+"pccs:"+to_string(res.first)+"\tele number:"+to_string(ESIZE)+"\n";
//                 write_into_file("./split.txt",kk.c_str());
// #endif
//                 SplitSegment(locate,splice.prev_, splice.next_,this,keys,values,ESIZE,res.first,res.second,n_start,n_stop);
//             }
//             else{
//                 rebuild_cnt++;
//                 subtree* newsubtree = rebuild_tree(keys,values,ESIZE,n_start,n_stop);
//                 if(i){
//                     int pos = path[i-1]->predict(key);
//                     RT_ASSERT(BITMAP_GET(path[i-1]->none_bitmap, pos) == 0);
//                     RT_ASSERT(BITMAP_GET(path[i-1]->child_bitmap, pos) == 1);
//                     path[i-1]->items[pos].comp.child = newsubtree;
//                 }
//                 else{
//                     locate->DataArray = newsubtree;
//                 }
//             }    
//             delete[] keys;
//             delete[] values;            
//             return true;
//         }
//     }
//     //free locate
//     return true;
// }



bool skiplist::Insert(State* proc_state){
    while(true){
        Segment_pt* locate = Scan(proc_state);
        if(locate == head_ || locate == tail_ || !locate){
            for(int i = 0;i<10000;i++);
            continue;
        }
        if(locate->CASLock()){
            // Acquire(proc_state->currs[0].lock);
            Add(proc_state,locate);
            locate->RealseLock();
            // Release(proc_state->currs[0].lock);
            break;
        }else{
            for(int i = 0;i<10000;i++);
        }
        // int result = Add(proc_state);
        
        // if(result == FAILURE) continue;
        // else return result == ADDED?true :false;
    }
    return true;
}

int skiplist::Add(State* proc_state,Segment_pt* locate){
    subtree* path[MAX_DEPTH];
    int path_size = 0;
    insert_subtree(locate,proc_state->key,proc_state->value,
    path,path_size);
    int key = proc_state->key;
    for(int i = 0;i<path_size;i++){
        subtree *n = path[i];
        const bool need_rebuild =  path_size >= MAX_DEPTH ;
        if(need_rebuild){
            const int ESIZE = n->size;
            unsigned int* keys = new unsigned int[ESIZE];
            int* values = new int[ESIZE];
            unsigned int n_start = n->start,n_stop = n->stop;
            memset(keys,0,ESIZE);
            memset(values,0,ESIZE);
            if(!locate->CASBuild()){
                //to fix
                delete[] keys;
                delete[] values; 
                return true;
            }
            std::pair<long double,int> res = scan_and_destroy_subtree(n,keys,values);
            scan_+=ESIZE;
            cerr<<"lin:"<<res.first<<"\t"<<ESIZE<<endl;
            if(!i && ((res.first < LINAERITY )|| ESIZE >1e6)){//|| ESIZE > 1e6
                if(SplitSegment(locate,proc_state,this,keys,values,ESIZE,res.first,res.second,n_start,n_stop)){
                    locate->RealseBuild();
                    delete[] keys;
                    delete[] values;         
                    return true;
                }
                cerr<<"LockAndValidateUntil failed"<<endl;
            }
            cerr<<"rebuild"<<endl;
            subtree* newsubtree = rebuild_tree(keys,values,ESIZE,n_start,n_stop);
            if(i){
                int pos = path[i-1]->predict(key);
                RT_ASSERT(BITMAP_GET(path[i-1]->none_bitmap, pos) == 0);
                RT_ASSERT(BITMAP_GET(path[i-1]->child_bitmap, pos) == 1);
                path[i-1]->items[pos].comp.child = newsubtree;
            }
            else{
                locate->DataArray = newsubtree;
            }
            locate->RealseBuild();
            delete[] keys;
            delete[] values;            
            return true;
        }
    }
    return true;
}

// void skiplist::show(){
//     Segment_pt* x = head_;
//     int depth = 0;
//     unsigned int pre_stop = 0;
//     for(int i = 0;i<segCnt;i++){
//         x = x->forward[1];
//         // cerr<<"Segment_pt level:"<<x->level<<" start:"<<x->DataArray->start<<" stop:"<<x->DataArray->stop<<" slope:"<<x->DataArray->slope<<" intertectpt:"<<x->DataArray->intercept<<"\n";
//         // cerr<<"Segment_pt ele cnt"<<x->DataArray->size<<" start:"<<x->DataArray->start<<" stop:"<<x->DataArray->stop;
//         if(i && pre_stop != x->DataArray->start){
//             cerr<<"?"<<endl;
//         }
//         pre_stop = x->DataArray->stop;
//         typedef std::pair<int, subtree*> Seg; // <begin, subtree*>
//         std::stack<Seg> s;
//         s.push(Seg(1, x->DataArray));
//         int d = 0;
//         while (!s.empty()) {
//             Segment_pt::subtree* n = s.top().second;
//             int l = s.top().first;
//             d = max(d,l);
//             s.pop();
//             for (int j = 0; j < n->num_items; j ++) {
//                 if (BITMAP_GET(n->none_bitmap, j) == 0) {
//                     if (BITMAP_GET(n->child_bitmap, j) == 1) {
//                         s.push(Seg(l+1,n->items[j].comp.child));
//                     }
//                 }
//             }
//         }
//         // cerr<<"\tsegment max depth:"<<d<<endl;
//         depth+=d;
//     }
//     cerr<<"average depth:"<<depth/segCnt<<endl;
//
// }

// void skiplist::ComputeSpace(){
//     long long space = 0;
//     Segment_pt *x = this->head_;
//     space+=sizeof(this);
//     //head_
//     while (x && x->forward[1] != this->head_) {
//         space+=x->Computespace();
//         x = x->forward[1];
//         // cerr<<"segment predict error sum cost on per key:"<<x->DataArray->cost/x->DataArray->size<<endl;
//     }
//
//     cerr<<"space size:"<<space<<endl;
//     char space_size_data[] = "./space_size.csv";
//     string k = to_string(space)+"\n";
//     write_into_file(space_size_data,k.c_str());
// }






