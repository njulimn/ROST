
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
#include <string>
#include <map>
#include <stack>
#include <cstring>
// #include "./fileOperation.hpp"

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


// runtime assert
#define RT_ASSERT(expr) \
{ \
    if (!(expr)) { \
        fprintf(stderr, "RT_ASSERT Error at %s:%d, `%s`\n", __FILE__, __LINE__, #expr); \
        exit(0); \
    } \
}

using namespace std;

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
			this->s0 = nullptr;
			this->s0 = k;
			this->s_last = nullptr;
			this->s_last = k;
			this->state = GreedyState::Need1;

			return current;
		}
		Point* s_u = k->upper_bound(this->gamma);
		Point* s_l = k->lower_bound(this->gamma);
		if (this->rho_upper->BeblowAssert(s_u)) {
			delete this->rho_upper;
			this->rho_upper = nullptr;
			this->rho_upper = new Line(this->sint, s_u);
		}

		if (this->rho_lower->AboveAssert(s_l)) {
			delete this->rho_lower;
			this->rho_lower = nullptr;
			this->rho_lower = new Line(this->sint, s_l);
		}
		return nullptr;

	}

	Segment* Process(double x, double y) {
		//delete this->s_last;
		Segment* res = nullptr;
		Point* newp = new Point(x, y);
		// cerr<<"in process:old state "<<this->state<<endl;
		GreedyState newState = this->state;
		
		// cerr<<"new state(old)："<<newState;
		// cerr<<"GreedyState::Need2"<<GreedyState::Need2<<endl;
		if (this->state == GreedyState::Need2) {
			this->s0 = nullptr;
			this->s_last = nullptr;
			newp->y = 0;
			this->s_last = newp;
			this->s0 = newp;
			newState = GreedyState::Need1;
			this->state = GreedyState::Need1;
		}
		else if (this->state == GreedyState::Need1) {
			//delete this->s1;
			this->s1 = nullptr;
			this->s_last = nullptr;
			this->s_last = newp;
			this->s1 = newp;
			this->setup(this->s0, this->s1);
			newState = GreedyState::Ready;
			this->state = GreedyState::Ready;
		}
		else if (this->state == GreedyState::Ready) {
			res = this->Process_pt(newp);
			if (res != nullptr) {
				newState = GreedyState::Need1;
				this->state = GreedyState::Need1;
			}
			else {
				newState = GreedyState::Ready;
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

typedef struct snode{
	unsigned int key;
	unsigned int level;
}snode;

class skiplist;

class Segment_pt{
    public:
		struct subtree;
		struct Item{
				union {
					node data;
					subtree* child;
				}comp;
		};
		struct subtree
		{
			int is_two = 0; // is special node for only two keys
			int build_size = 0; // tree size (include sub nodes) when node created
			int size = 0; // current tree size (include sub nodes)
			int fixed = 0; // fixed node will not trigger rebuild
			int num_inserts = 0, num_insert_to_data = 0;
			int child_ptr = 0;
			int num_items = 0; // size of items
			double slope = 0,intercept = 0;
			unsigned int start = UNINT_MAX,stop=UNINT_MAX;
			Item* items = nullptr;
			bitmap_t* none_bitmap = nullptr; // 1 means None, 0 means Data or Child
			bitmap_t* child_bitmap = nullptr; // 1 means Child. will always be 0 when none_bitmap is 1
			const double BUILD_LR_REMAIN = 0;
		};

		typedef struct {
            int begin;
            int end;
            int lvl; // top lvl = 1
            subtree* n;
		} Seg;

        vector<Segment_pt*> forward;//Segment_pt** forward;
        int level;
        subtree* DataArray;

        Segment_pt();
		Segment_pt(int level,subtree* n){
            this->level = level;
            vector<Segment_pt*> new_seg(level+1);//this->forward = (Segment_pt**)malloc(sizeof(Segment_pt*)*(level+1));
			this->forward = new_seg;
            for(int i=1;i<=level;i++){
                this->forward[i] = this;
            }
			std::vector<subtree*> trees;
            for (int _ = 0; _ < 1e3; _ ++) {
                subtree* tr = build_tree_two(0, 0, 1, 0,0,0,nullptr);
                trees.push_back(tr);
            }
            for (auto tr : trees) {
                destroy_subtree(tr);
            }

			DataArray = n;
        }
		Segment_pt(int level){
            this->level = level;
            vector<Segment_pt*> new_seg(level+1);//this->forward = (Segment_pt**)malloc(sizeof(Segment_pt*)*(level+1));
			this->forward = new_seg;
            for(int i=1;i<=level;i++){
                this->forward[i] = this;
            }
			if(level < 15){
				std::vector<subtree*> trees;
				for (int _ = 0; _ < 1e3; _ ++) {
					subtree* tr = build_tree_two(0, 0, 1, 0,0,0,nullptr);
					trees.push_back(tr);
				}
				for (auto tr : trees) {
					destroy_subtree(tr);
				}
			}
			DataArray = build_tree_none();
        }
		~Segment_pt(){
			destroy_subtree(DataArray);
			DataArray = NULL;
			destory_pending();
		}

		std::stack<subtree*> pending_two;

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

		void destroy_subtree(subtree* root);

		void destory_pending()
		{
			while (!pending_two.empty()) {
				subtree* n = pending_two.top(); pending_two.pop();

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

		void insert_subtree(unsigned int key,int value,Segment_pt** Update,skiplist* list);
		
		std::pair<int,Segment_pt::Item*> find_subtree(unsigned int key);

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

		inline int predict(subtree*root,unsigned int key){
			// cerr<<root->slope<<","<<root->intercept<<endl;
			double v = root->slope * (static_cast<long double>(key)-static_cast<long double>(root->start)) + root->intercept;
			// if (std::isnan(v)){
			// 	return root->num_items - 1;
			// }
			if(v > std::numeric_limits<int>::max() / 2) {
            	return root->num_items - 1;
			}
			if(v < 0 || static_cast<int>(v)<0){
				// cerr<<"prd:"<<v<<endl;
				return 0;
			}
			return std::min(root->num_items - 1, static_cast<int>(v));
		}

		bool SplitSegment(Segment_pt** Update,skiplist* list);

		/// build an empty tree
		subtree* build_tree_none();

		subtree* build_tree_two(unsigned int key1, int value1, unsigned int key2, int value2,unsigned int start=0,unsigned int stop=0,subtree* x = nullptr);

		subtree* rebuild_tree(unsigned int *_keys,int* _values,int _size,unsigned int start,unsigned int stop);

		void scan_and_destroy_subtree(subtree* root,unsigned int*keys,int* values, bool destory = true);

		void showArray(subtree* n){
			cerr<<"-";
			for (int i = 0; i < n->num_items; i ++) {
				if (BITMAP_GET(n->none_bitmap, i) == 0) {
					if (BITMAP_GET(n->child_bitmap, i) == 0) {
						cerr<<n->items[i].comp.data.key<<",";
					} else {
						showArray(n->items[i].comp.child);
					}
				}
			}
		}

		// void show(int w);

};

class skiplist {
    public:
        int skip_level;
        int segCnt;
        Segment_pt *header;
        int MaxLevel;
		int gamma;

        skiplist();
        skiplist(int MaxLevel,int gamma){
            Segment_pt *header = new Segment_pt(MaxLevel);//(snode *)malloc(sizeof(struct snode));
            this->header = header;
            this->skip_level = 1;
            this->segCnt = 0;
            this->MaxLevel = MaxLevel;
			this->gamma = gamma;
			// Segment_pt* newSeg = new Segment_pt(this->RandLevel(),0,0,0,0);
			// for(int i = newSeg->level;i>=1;i--){
			// 	this->header->forward[i] = newSeg;
			// 	newSeg->forward[i] = this->header;
			// }
        }

		void setup(vector<snode> input);

		node* binarySearch(Segment_pt* x,unsigned int key);
		std::pair<int,Segment_pt::Item*> Search(unsigned int key);

		void Insert(unsigned int key,int value);

		void insert_static(vector<Segment_pt*> &Update,Segment* seg,unsigned int st,unsigned int ed,node* input,int size,int level);

		void show();

        void ShowNodeDis();
        void Dump();

		void ComputeSpace();
		        
        int RandLevel(){
            int NewLevel = 1;
            int P = p*100;
            while(1){
                int t = rand() % 101;
                if(t<P)
                    NewLevel++;
                else
                    break;
            }
            return min(NewLevel,this->MaxLevel);
        }
        // int Delete(skiplist *list,int k);
        // int Insert(skiplist *list,int k,int v);
};






