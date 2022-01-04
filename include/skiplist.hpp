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

class skiplist {
    public:
        std::atomic<int> max_height_; 
        int MaxLevel;
		int gamma;
        int segment_max_size;
        double linearity;
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
            bool marked;
            bool split;
            std::atomic<bool> lock;//false means segment is locked
            // std::atomic<bool> build; 
            std::atomic<bool> spliting;//true means segment is locked
            std::atomic<Segment_pt*> *next_;
            // std::atomic<Segment_pt*> next_[1];

            inline void InitLock(){
                lock.store(true, std::memory_order_release);
            }

            // inline void InitBuildM(){
            //     build.store(false, std::memory_order_release);
            // }

            inline void InitSplit(){
                spliting.store(false, std::memory_order_release);
            }

            inline bool ReadLock(){
                return lock.load(std::memory_order_acquire);
            }
            //true means segment is spliting or rebuilding
            // inline bool ReadBuildM(){
            //     return build.load(std::memory_order_acquire);
            // }

            inline bool ReadSplit(){
                return spliting.load(std::memory_order_acquire);
            }

            inline bool CASLock(bool free = true,bool block=false){
                return lock.compare_exchange_strong(free, block);
            }

            // inline bool CASBuild(bool free = false,bool block=true){
            //     return build.compare_exchange_strong(free, block);
            // }

            inline bool CASSplit(bool free = false,bool block=true){
                return spliting.compare_exchange_strong(free, block);
            }            

            inline void ReleaseLock(){
                lock.store(true, std::memory_order_release);
            }

            // inline void ReleaseBuild(){
            //     build.store(false, std::memory_order_release);
            // }

            inline void ReleaseSplit(){
                spliting.store(false, std::memory_order_release);
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
            Segment_pt** preds;
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
            newseg->split = 0;
            // newseg->rebuild = 0;
            // newseg->split = 0;
            newseg->InitLock();
            // newseg->InitBuildM();
            newseg->InitSplit();
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
            newseg->split = 0;
            // newseg->rebuild = 0;
            // newseg->split = 0;
            newseg->InitLock();
            // newseg->InitBuildM();
            newseg->InitSplit();
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
            newseg->split = 0;
            // newseg->rebuild = 0;
            // newseg->split = 0;
            newseg->InitLock();
            // newseg->InitBuildM();
            newseg->InitSplit();
            newseg->next_ = new std::atomic<Segment_pt*>[level];
            return newseg;
        }
        
        // void deleteSegment_pt(){
        //     destroy_subtree(DataArray);
        //     DataArray = NULL;
        //     // destory_pending();
        // }

        /////////////////////////////////
        
        std::atomic<int> segCnt;
        Segment_pt *head_;
        Segment_pt *tail_;
        skiplist();
        skiplist(int MaxLevel,int gamma):max_height_(1),MaxLevel(MaxLevel),gamma(gamma),
        segment_max_size(1e5),linearity(0.98),segCnt(1){
            head_ = NewSegment_pt(MaxLevel,0,1);
            head_->DataArray->start = 0;
            head_->DataArray->stop = 1;
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

        inline int GetSegCnt() const{
            return segCnt.load(std::memory_order_relaxed);
        }

        Segment_pt* Scan(State* proc_state);

        // void FindStateForLevel(State* proc_state,int level);

        // std::pair<int,int> Seek(State* proc_state);

        std::pair<int,int> Lookup(State* proc_state);

		// std::pair<int,node*> Search(unsigned int key);

        int Add(State* proc_state,Segment_pt* locate);

		// bool Insert(unsigned int key,int value);

        bool Insert(State* proc_state);

		void insert_static(vector<Segment_pt*> &Update,Segment* seg,unsigned int st,unsigned int ed,node* input,int size,int level);

        bool LockAndValidateUntil(Segment_pt* from[],Segment_pt* to[],int level);

        bool LockAndValidateUntil(Segment_pt* from[],Segment_pt* to[],int from_level,int to_level);

        void UnlockUntil(Segment_pt* from[],int level);
        //release [from[from_level],from[to_level])
        void UnlockUntil(Segment_pt* from[],int from_level,int to_level);

        bool RecomputeAndLock(Segment_pt* from[],Segment_pt* to[],int level_from,int level_to,unsigned int key);

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

bool skiplist::RecomputeAndLock(Segment_pt* from[],Segment_pt* to[],int level_from,int level_to,unsigned int key){
    Segment_pt* pred = from[level_from];
    // Segment_pt* ppred = nullptr;
    Segment_pt* curr = nullptr;
    Segment_pt* succ = nullptr;
    for(int l = level_from;l>level_to;l--){
        curr = pred->Next(l);
        succ = pred->Next(l);
        while(key >= succ->anchor){
            // ppred = pred;
            pred = curr;
            curr = succ;
            succ = succ->Next(l);
        }
        //key < succ->anchor
        if(key >= curr->bound){
            //key is after curr, before succ
            if(!curr->CASLock()){
                //false;
                for(int j = l+1;j<=level_from;j++){
                    from[l]->ReleaseLock();
                }
                return false;
            }
            // ppred = pred;
            pred = curr;
        }else{
            succ = curr;
            if(!pred->CASLock()){
                //false;
                for(int j = l+1;j<=level_from;j++){
                    from[l]->ReleaseLock();
                }
                return false;
            }
        }
        if(succ == nullptr){
            succ = nullptr;
        }
        from[l] = pred;
        to[l] = succ;
    }
    return true;
}

/*
split root segment into several subsegments using greedy plr
before SplitSegment,both of the write lock of root and split lock of root are all obtained
At the end of SplitSegment function,is is needed to release the split lock firstly and release the write lock secondly.
*/
bool skiplist::SplitSegment(Segment_pt *root,State* proc_state,skiplist* list,unsigned int *_keys,int*  _values,int _size,long double pccs,int middle,unsigned int _start,unsigned int _stop){
    // cerr<<"split"<<endl;
    int st = 0,ed =0;//the rank of the first/last+1 ele of one segment
    GreedyPLR* plr = new GreedyPLR(Gm);
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
    int seg_size = static_cast<int>(seg.size());
    int to_level = root->level;
    //scan & lock
#if WRITESEG
    string kk = "\nseg size "+to_string(seg.size())+":\n";
    write_into_file("./split.txt",kk.c_str());
#endif
    vector<Segment_pt*> split_segments(seg_size);
    vector<Segment_pt*> split_inner_pred(MaxLevel,nullptr);
    vector<Segment_pt*> split_left(MaxLevel,nullptr);
    vector<Segment_pt*> split_right(MaxLevel,nullptr);
    int max_seg_height = to_level;
    Segment_pt** preds = (Segment_pt**)(malloc(sizeof(Segment_pt*)*MaxLevel));
    memcpy(preds,proc_state->preds,sizeof(Segment_pt*)*MaxLevel);
    root->split = 1;
    split_segments[0] = root;
    subtree *root_dataarray = nullptr;
    unsigned int split_lower_bound = 0;
    for(int i = 0;i<seg_size;i++){
        int height_tmp,ed,st=segment_stIndex[i];
        unsigned int base,bound;
        if(i == 0){
            height_tmp = max_seg_height;
            base = _start;
        }else{
            height_tmp = height_seq[i-1];
            base  = _keys[st];
        }
        if(i == seg_size-1){
            bound = _stop;
            ed = _size;
        }else{
            ed = segment_stIndex[i+1];
            bound = _keys[ed];
        }
        Segment_pt* next_seg = nullptr;
        if(i){
            next_seg = NewSegment_pt_nodata(height_tmp,base,bound);
            next_seg->split = 1;
            split_segments[i] = next_seg;
            // Acquire(next_seg);
            RT_ASSERT(next_seg->CASLock());
            RT_ASSERT(next_seg->CASSplit());
        }
        max_seg_height = max(max_seg_height,height_tmp);
        RT_ASSERT(ed-st>0);
        subtree* new_subtree = nullptr;
        if(ed-st == 1){
            new_subtree = build_tree_none();
            new_subtree->start = base;
            new_subtree->stop = bound;
            BITMAP_CLEAR(new_subtree->none_bitmap,0);
            new_subtree->size++;
            new_subtree->num_inserts++;
            new_subtree->items[0].comp.data.key = _keys[st];
            new_subtree->items[0].comp.data.value = _values[st];
        }
        else{
            if(seg[i]->slope > 0)
                new_subtree = rebuild_tree(_keys+st,_values+st,ed-st,base,bound,true,seg[i]->slope,seg[i]->intercept);
            else
                new_subtree = rebuild_tree(_keys+st,_values+st,ed-st,base,bound);//,true,seg[i]->slope,seg[i]->intercept);
        }
        if(i){
            //inner link
            split_segments[i]->DataArray = new_subtree;
            for(int l = height_tmp-1;l>=0;l--){
                next_seg->SetNext(l,nullptr);
                if(split_inner_pred[l]){
                    split_inner_pred[l]->SetNext(l,next_seg);
                    split_inner_pred[l] = next_seg;
                }else{
                    split_inner_pred[l] = next_seg;
                    split_left[l] = next_seg;
                }
                split_right[l] = next_seg;
            }
        }else{
            root_dataarray = new_subtree;
            split_lower_bound = bound;
        }
        
    }

    if(seg_size == 1){
        root->DataArray = root_dataarray;
        root->ReleaseSplit();
        root->ReleaseLock();
        return true;
    }
    else{
        #if DEBUG
        cerr<<"after link inner"<<endl;
        for(int i = 1;i<seg_size;i++){
            cerr<<"segment "<<split_segments[i]<<"["<<split_segments[i]->anchor
                <<","<<split_segments[i]->bound<<")"<<endl;
            for(int j = 0;j<split_segments[i]->level;j++){
                cerr<<j<<"->"<<split_segments[i]->Next(j)<<endl;
            }
            cerr<<endl;
        }
        #endif
        //link root
        if(split_right[0]){
            for(int i = 0;i<to_level;i++){
                if(split_right[i]){
                    split_right[i]->SetNext(i,root->Next(i));
                }
            }
            proc_state->succs[to_level-1] = root->Next(to_level-1);
            for(int i = 0;i<to_level;i++){
                if(split_left[i]){
                    root->SetNext(i,split_left[i]);
                }
            }
            root->bound = split_lower_bound;
            root->DataArray = root_dataarray;
            root->ReleaseSplit();
            root->ReleaseLock();
            for(int i = to_level;i<max_seg_height;i++){
                if(proc_state->succs[i]){
                    split_right[i]->SetNext(i,proc_state->succs[i]);
                }else{
                    split_right[i]->SetNext(i,proc_state->succs[i-1]);
                }
            }
            #if DEBUG
            cerr<<"after link root"<<endl;
            for(int i = 0;i<seg_size;i++){
                cerr<<"segment "<<split_segments[i]<<"["<<split_segments[i]->anchor
                    <<","<<split_segments[i]->bound<<")"<<endl;
                for(int j = 0;j<split_segments[i]->level;j++){
                    cerr<<j<<"->"<<split_segments[i]->Next(j)<<endl;
                }
                cerr<<endl;
            }
            #endif
        }

        for(int i = 1;i<seg_size;i++){
            split_segments[i]->ReleaseLock();
        }

        int max_height = list->max_height_.load(std::memory_order_relaxed);
        while (segment_max_height > max_height) {
            if (list->max_height_.compare_exchange_weak(max_height, max_seg_height)) {
                // successfully updated it
                max_height = segment_max_height;
                break;
            }
        }

        Segment_pt *pred = head_;
        Segment_pt *curr = nullptr;
        Segment_pt *succ = nullptr;
        for(int l = max_height-1;l>=to_level;l--){
            curr = pred->Next(l);
            RT_ASSERT(curr!=nullptr);
            succ = curr->Next(l);
            while(succ && succ->anchor <= split_lower_bound){
                pred = curr;
                curr = succ;
                succ = succ->Next(l);
            }
            //key < succ->anchor
            //not sure the location of key
            //there are 2 cases
            if(split_lower_bound >= curr->bound){
                //case2 : key is after curr ,but before succ
                pred = curr;//curr != tail_
            }else{
                //case1 : key is before
                #if DEBUG
                if(split_lower_bound >= curr->anchor){
                    cerr<<"split_lower_bound >= curr->anchor\n";
                }
                #endif
                RT_ASSERT(split_lower_bound < curr->anchor);
            }
            if(split_left[l]){
                while(!split_right[l]->CASLock()){
                    //wait 
                    //TODO:change the lock type
                    ;
                }
                while (! pred->CASLock())
                {
                    //wait
                    //TODO:change the lock type
                    ;
                }
                proc_state->preds[l] = pred;
                split_right[l]->SetNext(l,pred->Next(l));
                pred->SetNext(l,split_left[l]);
                // RT_ASSERT(split_right[l]->CASNext(l,nullptr,curr));
                // RT_ASSERT(pred->CASNext(l,curr,split_left[l]));
                pred->ReleaseLock();
                split_right[l]->ReleaseLock();
            }
        }
        #if DEBUG
        cerr<<endl;
        cerr<<"after link out,preds next"<<endl;
        for(int l = segment_max_height-1;l>=to_level;l--){
            cerr<<"level:"<<l<<" segment "<<proc_state->preds[l]<<"["<<proc_state->preds[l]->anchor
                <<","<<proc_state->preds[l]->bound<<")"<<"->"<<proc_state->preds[l]->Next(l)<<endl;
        }
        
        cerr<<endl;
        cerr<<"split next"<<endl;
        for(int i = 0;i<seg_size;i++){
            cerr<<"segment "<<split_segments[i]<<endl;
            for(int j = 0;j<split_segments[i]->level;j++){
                cerr<<j<<"->"<<split_segments[i]->Next(j)<<endl;
            }
            cerr<<endl;
        }
        
        int raw_segcnt = list->GetSegCnt();
        while(!(list->segCnt.compare_exchange_strong(raw_segcnt, raw_segcnt+seg_size-1))) {
            // successfully updated it
            raw_segcnt = list->GetSegCnt();
        }
        #endif
        for(int i = 1;i<seg_size;i++){
            split_segments[i]->ReleaseSplit();
        }
        #if DEBUG
        for(int i = 1;i<seg_size;i++){
            for(int j = 0;j<split_segments[i]->level;j++){
                if(split_segments[i]->Next(j) == nullptr){
                    cerr<<"no";
                }
            }
        }
        #endif
        return true;
    }
    // cerr<<"split end\n";
    // cerr<<seg_size<<endl;
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
    s.push({0, _size, 1, ret});
    while(!s.empty()){
        const int begin = s.front().begin;
        const int end = s.front().end;
        const int lvl = s.front().lvl;
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
            RT_ASSERT(n->slope >= 0);

            //TODO:fix size?
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
                    if (next_i == item_i) {
                        next ++;
                    } else {
                        break;
                    }
                }
                if (next == offset + 1) {
                    real_put_ele++;
                    RT_ASSERT(BITMAP_GET(n->none_bitmap,item_i) == 1);
                    RT_ASSERT(BITMAP_GET(n->child_bitmap,item_i) == 0);
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
                    s.push({begin + offset, begin + next, lvl + 1, n->items[item_i].comp.child});
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
                    begin += n->items[i].comp.child->size;
                }
            }
        }
        RT_ASSERT(SHOULD_END_POS == begin);
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
                for(int i = 0;i<path_size-1;i++){
                    path[i]->size--;
                    path[i]->num_inserts--;
                }
                // n->size++;
                // n->num_inserts++; 
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
            RT_ASSERT(n->items[pos].comp.child != nullptr);
            n = n->items[pos].comp.child;
        }
    }
    for(int i = 0;i<path_size;i++){
        path[i]->num_insert_to_data += insert_to_data;
    }

}

/////////////////////////////////////////////

// void skiplist::FindStateForLevel(skiplist::State* proc_state,int level){
//     Segment_pt* pred = proc_state->preds[level+1];
//     Segment_pt* curr = pred->Next(level);
//     Segment_pt* succ = nullptr;
//     Segment_pt* locate = nullptr;
//     succ = curr->Next(level);   
//     while(succ && succ->anchor <= proc_state->key){
//         pred = curr;
//         curr = succ;
//         succ = succ->Next(level);
//     }
//     if(proc_state->key >= curr->bound){
//         //key is bwtween curr and succ
//         pred = curr;//curr != tail_
//     }else{
//         //key is before or located in curr
//         if(!locate && proc_state->key >= curr->anchor){
//             //key is located in curr
//             locate = curr;
//             for(;level>=0;level--){
//                 proc_state->preds[level] = curr;
//                 proc_state->succs[level] = curr->Next(level);
//             }
//             return locate;
//         }
//         //key is rebore curr
//         succ = curr;
//     }
//     if(pred->NoBarrier_Next(level) != succ){
//         cerr<<"no";
//     }
//     proc_state->preds[level] = pred;
//     proc_state->succs[level] = succ;
//     return;
// }

skiplist::Segment_pt* skiplist::Scan(skiplist::State* proc_state){
    // Segment_pt* ppred = nullptr;
    Segment_pt* pred = head_;
    Segment_pt* curr = nullptr;
    Segment_pt* succ = nullptr;
    Segment_pt* locate = nullptr;
    int height_ = GetMaxHeight();
    for(int l = MaxLevel-1;l>=height_;l--){
        proc_state->succs[l] = tail_;
    }
    unsigned int key = proc_state->key;
    #if DEBUG
    vector<Segment_pt*> paths;
    vector<Segment_pt*> nexts;
    #endif
    for(int l = height_-1;l>=0;l--){
        //not sure key > curr's anchor
        #if DEBUG
        if(pred->Next(l) == nullptr){
            for(int j = pred->level-1;j>=0;j--){
                Segment_pt* mm = pred->Next(j);
                cerr<<"j:"<<mm;
                if(mm){
                    cerr<<mm->anchor<<"\t"<<mm->bound;
                }
                cerr<<endl;
            }
            curr = nullptr;
        }
        #endif
        curr = pred->Next(l);
        RT_ASSERT(curr!=nullptr);
        succ = curr->Next(l);   
        while(succ && succ->anchor <= key){
            // ppred = pred;
            pred = curr;
            curr = succ;
            succ = succ->Next(l);
        }
        //there are 3 case
        if(key >= curr->bound){
            //case 3:key is after curr but before succ
            // ppred = pred;
            pred = curr;//curr != tail_
            // curr = succ;
            proc_state->succs[l] = succ;
        }else{
            //key is before or located in curr
            if(!locate && key >= curr->anchor){
                //case 2:key is located in curr
                // locate = curr;
                return curr;
            }
            //case 1:key is after pred but before curr
            // succ = curr;
            proc_state->succs[l] = curr;
        }
        #if DEBUG
        paths.push_back(pred);
        nexts.push_back(curr);
        #endif
        
    }
    return locate;
}

bool skiplist::LockAndValidateUntil(Segment_pt* from[],Segment_pt* to[],int level){
    Segment_pt* pre_ = from[0];
    if(!pre_->CASLock()){
        return false;
    }
    for(int l = 0;l<level;l++){
        // Acquire(from[l]->lock);
        ///如果该segment被锁，当前process是应该挂起还是循环
        assert(from[l]!=nullptr);
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
            RT_ASSERT(from[l]!=nullptr);
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
        RT_ASSERT(from[l]!=nullptr);
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
        from[l]->ReleaseLock();
    }
}

void skiplist::UnlockUntil(Segment_pt* from[],int from_level,int to_level){
    //release[from_level,to_level])
    for(int l = from_level;l<to_level;l++){
        from[l]->ReleaseLock();
    }
}

std::pair<int,int> skiplist::Lookup(State* proc_state){   
    Segment_pt* locate = Scan(proc_state);
    if(locate == nullptr){
        string kk = to_string(proc_state->key)+"\tnullptr\t0"+"\n";
        write_into_file("./nofind.txt",kk.c_str());
        return {0,0};
    }else if(locate->ReadSplit()){
        string kk = to_string(proc_state->key)+"\tbuild lock\t0"+"\n";
        write_into_file("./nofind.txt",kk.c_str());
        return {0,1};//R-W conflicts
    }else if(locate->ReadLock() == false){
        string kk = to_string(proc_state->key)+"\t segment lock\t0"+"\n";
        write_into_file("./nofind.txt",kk.c_str());
        return {0,1};//R-W conflicts
    }
    // if(proc_state->currs[0]->anchor > proc_state->key || proc_state->currs[0]->DataArray->stop >= proc_state->key ){
    //     return {0,0};
    // }
    std::pair<int,node*> res = find_subtree(locate->DataArray,proc_state->key);
    if(res.first){
        // string kk = to_string(proc_state->key)+"\t"+to_string(locate->split)+"\t1\n";
        // write_into_file("./nofind.txt",kk.c_str());
        return {res.first,res.second->value};
    }else{
        string kk = to_string(proc_state->key)+"\t"+to_string(locate->split)+"\t0\n";
        write_into_file("./nofind.txt",kk.c_str());
        return {0,0};
    }
}

bool skiplist::Insert(State* proc_state){   
    while(true){
        Segment_pt* locate = Scan(proc_state);
        if(locate){
            if(locate->CASLock()){
                if(locate->anchor > proc_state->key || locate->bound <= proc_state->key){
                    locate->ReleaseLock();
                    continue;
                }
                Add(proc_state,locate);
                return true;
            }
        }
    }
    return true;
}

/*
before Add, the write lock of locate is obtained
At the end of Add, it is needed to release the write lock of locate
*/
int skiplist::Add(State* proc_state,Segment_pt* locate){
    subtree* path[MAX_DEPTH*2];
    int path_size = 0;
    insert_subtree(locate,proc_state->key,proc_state->value,path,path_size);
    // int key = proc_state->key;
    const bool need_rebuild = path_size >= MAX_DEPTH ;
    if(!need_rebuild){
        locate->ReleaseLock();
        return true;
    }
    else{
        if(! locate->CASSplit()){
            locate->ReleaseLock();
            return true;
        }
        subtree *n = locate->DataArray;
        const int ESIZE = n->size;
        unsigned int* keys = new unsigned int[ESIZE];
        int* values = new int[ESIZE];
        unsigned int n_start = n->start,n_stop = n->stop;
        memset(keys,0,ESIZE);
        memset(values,0,ESIZE);
        std::pair<long double,int> res = scan_and_destroy_subtree(n,keys,values);
        locate->DataArray = nullptr;
        scan_+=ESIZE;
        // cerr<<res.first<<endl;
        if((res.first < linearity )|| ESIZE >segment_max_size){
            if(SplitSegment(locate,proc_state,this,keys,values,ESIZE,res.first,res.second,n_start,n_stop)){
                //SplitSegment already released the write lock of locate
                delete[] keys;
                delete[] values;         
                return true;
            }
        }
        // cerr<<"rebuild"<<endl;
        subtree* newsubtree = rebuild_tree(keys,values,ESIZE,n_start,n_stop);
        locate->DataArray = newsubtree;
        locate->ReleaseSplit();
        locate->ReleaseLock();
        delete[] keys;
        delete[] values;         
        return true;
    }
    return true;
}