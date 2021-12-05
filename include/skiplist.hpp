
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

int collsion = 0;
int rebuild_cnt = 0;
int split_cnt = 0;
int file_num = 0;
long long jmp_seg = 0;
long long jmp_subtree = 0;
long long scan_ = 0;

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
			// int child_ptr = 0;
			int num_items = 0; // size of items
			double slope = 0,intercept = 0;
			unsigned int start = UNINT_MAX,stop=UNINT_MAX;
			Item* items = nullptr;
			bitmap_t* none_bitmap = nullptr; // 1 means None, 0 means Data or Child
			bitmap_t* child_bitmap = nullptr; // 1 means Child. will always be 0 when none_bitmap is 1
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
		static std::stack<subtree*> pending_two;

        Segment_pt();
		Segment_pt(int level,subtree* n){
            this->level = level;
            vector<Segment_pt*> new_seg(level+1);//this->forward = (Segment_pt**)malloc(sizeof(Segment_pt*)*(level+1));
			this->forward = new_seg;
            for(int i=1;i<=level;i++){
                this->forward[i] = this;
            }
			std::vector<subtree*> trees;

			DataArray = n;
        }
		Segment_pt(int level,bool mempool = false){
            this->level = level;
            vector<Segment_pt*> new_seg(level+1);//this->forward = (Segment_pt**)malloc(sizeof(Segment_pt*)*(level+1));
			this->forward = new_seg;
            for(int i=1;i<=level;i++){
                this->forward[i] = this;
            }
			if(mempool){
				std::vector<subtree*> nodes;
				for (int _ = 0; _ < 1e5; _ ++) {
					subtree* node = build_tree_two(0,0,1,1);
					nodes.push_back(node);
				}
				for (auto node : nodes) {
					destroy_subtree(node);
				}
			}
			DataArray = build_tree_none();
        }
		~Segment_pt(){
			destroy_subtree(DataArray);
			DataArray = NULL;
			destory_pending();
		}

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

				if (n->is_two) {
					RT_ASSERT(n->build_size == 2);
					RT_ASSERT(n->num_items == 32);
					n->size = 2;
					n->num_inserts = n->num_insert_to_data = 0;
                    // n->num_items = 8;
                    // n->none_bitmap[0] = 0xff;
                    // n->child_bitmap[0] = 0;
                    for(int _=0;_<4;_++){
                        n->none_bitmap[_] = 0xff;
                        n->child_bitmap[_] = 0;
                    }
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
		
		std::pair<int,node*> find_subtree(unsigned int key);

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

		void SplitSegment(Segment_pt** Update,skiplist* list,unsigned int *_keys,int*  _values,int _size,long double pccs,int middle,unsigned int _start,unsigned int _stop);

		/// build an empty tree
		subtree* build_tree_none();

		subtree* build_tree_two(unsigned int key1, int value1, unsigned int key2, int value2,unsigned int start=0,unsigned int stop=0,subtree* x = nullptr);

		subtree* rebuild_tree(unsigned int *_keys,int* _values,int _size,unsigned int start,unsigned int stop,bool plr = false,double top_slope = 0,double top_intercept = 0);

		std::pair<long double,int> scan_and_destroy_subtree(subtree* root,unsigned int*keys,int* values, bool destory = true);

		long long Computespace();

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
            Segment_pt *header = new Segment_pt(MaxLevel,true);//(snode *)malloc(sizeof(struct snode));
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
		std::pair<int,node*> Search(unsigned int key);

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


std::stack<Segment_pt::subtree*> s;
std::stack<Segment_pt::subtree*> Segment_pt::pending_two = s;


void Segment_pt::SplitSegment(Segment_pt** Update,skiplist* list,unsigned int *_keys,int*  _values,int _size,long double pccs,int middle,unsigned int _start,unsigned int _stop){
    cerr<<"split"<<endl;
#if USEPLR
    int st = 0,ed =0;
    GreedyPLR* plr = new GreedyPLR(Gm);
    vector<Segment*> seg;
    vector<int> segment_stIndex;
    for (int i = 0; i < _size; i++) {
		Segment* seg_res = nullptr;
		seg_res = plr->Process(static_cast<double>(_keys[i])-static_cast<double>(_keys[st]), static_cast<double>(i-st));
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
    // cerr<<"seg size:"<<seg.size()<<endl;
#if WRITESEG
    string kk = "\nseg size "+to_string(seg.size())+":\n";
    write_into_file("./split.txt",kk.c_str());
#endif
    for(int i=0;i<static_cast<int>(seg.size());i++){
        st = segment_stIndex[i];
        unsigned int stop = 0;
        unsigned int start_local = 0;
        if(st == 0){
            start_local = _start;
        }else{
            start_local = _keys[st];
        }
        if(i == static_cast<int>(seg.size())-1){
            ed = _size;
            // stop = _keys[ed-1]+1;
            stop = _stop;
        }
        else{
            ed = segment_stIndex[i+1];
            stop = _keys[ed];
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
        if(i == 0){
            this->DataArray = new_subtree;
            for(int l = this->level;l>=1;l--){
                Update[l] = this;
            }
        }
        else{
            //Update 应保存 每层都在newSeg左边的segment
            Segment_pt* newSeg = new Segment_pt(list->RandLevel(),new_subtree);
            if(newSeg->level > list->skip_level){
                for(int l = newSeg->level;l>list->skip_level;l--){
                    Update[l] = list->header;
                }
                list->skip_level = newSeg->level;
            }
            for(int l = newSeg->level;l>=1;l--){
                newSeg->forward[l] = Update[l]->forward[l];
                Update[l]->forward[l] = newSeg;
                Update[l] = newSeg;
            }
            list->segCnt++;
        }
    }
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
                // cerr<<"split pccs:"<<PCCs<<endl;
                subtree* new_subtree = nullptr;
                if(st_idx == 0){
                    start_local = _start;
                }else{
                    start_local = _keys[st_idx];
                    RT_ASSERT(start_local == pre_stop);
                }
                if(i >= _size -2){
                    string kk = to_string(_size-st_idx)+",";
                    write_into_file("./split.txt",kk.c_str());
                    new_subtree = rebuild_tree(_keys+st_idx,_values+st_idx,_size-st_idx,start_local,_stop);
                }else{
                    string kk = to_string(i-st_idx)+",";
                    write_into_file("./split.txt",kk.c_str());
                    new_subtree = rebuild_tree(_keys+st_idx,_values+st_idx,i-st_idx,start_local,_keys[i]);
                    pre_stop = _keys[i];
                }
                if(st_idx == 0){
                    this->DataArray = new_subtree;
                    for(int l = this->level;l>=1;l--){
                        Update[l] = this;
                    }
                }
                else{
                    //Update 应保存 每层都在newSeg左边的segment
                    Segment_pt* newSeg = new Segment_pt(list->RandLevel(),new_subtree);
                    if(newSeg->level > list->skip_level){
                        for(int l = newSeg->level;l>list->skip_level;l--){
                            Update[l] = list->header;
                        }
                        list->skip_level = newSeg->level;
                    }
                    for(int l = newSeg->level;l>=1;l--){
                        newSeg->forward[l] = Update[l]->forward[l];
                        Update[l]->forward[l] = newSeg;
                        Update[l] = newSeg;
                    }
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
    //         this->DataArray = new_subtree;
    //         for(int l = this->level;l>=1;l--){
    //             Update[l] = this;
    //         }
    //     }
    //     else{
    //         //Update 应保存 每层都在newSeg左边的segment
    //         Segment_pt* newSeg = new Segment_pt(list->RandLevel(),new_subtree);
    //         if(newSeg->level > list->skip_level){
    //             for(int l = newSeg->level;l>list->skip_level;l--){
    //                 Update[l] = list->header;
    //             }
    //             list->skip_level = newSeg->level;
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
    string kk = to_string(middle)+","+to_string(right_part);
    write_into_file("./split.txt",kk.c_str());
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
            this->DataArray = new_subtree1;
            // cerr<<"split end"<<endl;
            return ;
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
            this->DataArray = new_subtree1;
#if DEBUG
            cerr<<"split end"<<endl;
#endif
            return ;
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

    list->segCnt++;
#endif  
    cerr<<"split end"<<endl;
}

Segment_pt::subtree* Segment_pt::build_tree_none(){
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
    } else {
        n = pending_two.top();
        // n->child_ptr = 0;
        // memset(n->child_bitmap,0,n->num_items);
        // memset(n->none_bitmap,0xff,n->num_items);
        pending_two.pop();
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

Segment_pt::subtree* Segment_pt::rebuild_tree(unsigned int *_keys,int*  _values,int _size,unsigned int start,unsigned int stop,bool plr,double top_slope,double top_intercept){
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

std::pair<long double,int> Segment_pt::scan_and_destroy_subtree(subtree* _root,unsigned int*keys,int* values, bool destory){
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
            if (n->is_two) {
                RT_ASSERT(n->build_size == 2);
                RT_ASSERT(n->num_items == 32);
                n->size = 2;
                n->num_inserts = n->num_insert_to_data = 0;
                // n->none_bitmap[0] = 0xff;
                // n->child_bitmap[0] = 0;
                for(int _=0;_<4;_++){
                    n->none_bitmap[_] = 0xff;
                    n->child_bitmap[_] = 0;
                }
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

void Segment_pt::insert_subtree(unsigned int key,int value,Segment_pt** Update,skiplist* list){
    // constexpr int MAX_DEPTH = 10;
    subtree* path[MAX_DEPTH];
    int path_size = 0;
    int insert_to_data = 0;
    subtree* n = this->DataArray;
    while(1){
        path[path_size++] = n;
        n->stop = max(n->stop,key);
        int pos = predict(n,key);
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
            // n->child_ptr++;
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
        // const int num_inserts = n->num_inserts;
        // const int num_insert_to_data = n->num_insert_to_data;
        //(n->size > 64 && (n->size*1.0/(n->stop - n->start) < 0.5))
        const bool need_rebuild =  path_size >= MAX_DEPTH ;//|| ( n->fixed == 0 && n->size >= n->build_size * 4 && n->size >= 64 && num_insert_to_data * 10 >= num_inserts);
        // int T;
        if(need_rebuild){
            const int ESIZE = n->size;
            unsigned int* keys = new unsigned int[ESIZE];
            int* values = new int[ESIZE];
            unsigned int n_start = n->start,n_stop = n->stop;
            memset(keys,0,ESIZE);
            memset(values,0,ESIZE);
            std::pair<long double,int> res = scan_and_destroy_subtree(n,keys,values);
            scan_+=ESIZE;
            // scan_and_destroy_subtree(n,keys,values);
            // cerr<<"lin:"<<res.first<<"\t"<<ESIZE<<endl;
            if(!i && ((res.first < LINAERITY && ESIZE>1e5 )|| ESIZE >1e6)){//|| ESIZE > 1e6
                split_cnt++;
#if WRITESEG
                string kk = "\nheight:"+to_string(path_size)+"pccs:"+to_string(res.first)+"\tele number:"+to_string(ESIZE)+"\n";
                write_into_file("./split.txt",kk.c_str());
#endif
                SplitSegment(Update,list,keys,values,ESIZE,res.first,res.second,n_start,n_stop);
            }
            else{
                rebuild_cnt++;
                subtree* newsubtree = rebuild_tree(keys,values,ESIZE,n_start,n_stop);
                if(i){
                    int pos = predict(path[i-1],key);
                    RT_ASSERT(BITMAP_GET(path[i-1]->none_bitmap, pos) == 0);
                    RT_ASSERT(BITMAP_GET(path[i-1]->child_bitmap, pos) == 1);
                    path[i-1]->items[pos].comp.child = newsubtree;
                }
                else{
                    this->DataArray = newsubtree;
                }
            }    
#if DEBUG
            cerr<<"(level:"<<i<<")replace it with"<<newsubtree<<" it's child ptr:"<<newsubtree->child_ptr<<endl;
#endif
            delete[] keys;
            delete[] values;            
            return;
        }
    }

}

std::pair<int,node*> Segment_pt::find_subtree(unsigned int key){
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
            return {n->items[pos].comp.data.key == key,&(n->items[pos].comp.data)};
        }
        else{
            jmp_subtree++;
            n = n->items[pos].comp.child;
        }
    }
    return {false,0};
}

long long Segment_pt::Computespace(){
    long long space = 0;
    space+=sizeof(Segment_pt);
    space+=(forward.capacity() * sizeof(Segment_pt*));
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

/////////////////////////////////////////////

std::pair<int,node*> skiplist::Search(unsigned int key){
    Segment_pt* x = this->header;
    // unsigned int pred;
    for(int i = this->skip_level;i>=1;i--){
        while(key >= x->forward[i]->DataArray->stop){
            x = x->forward[i];
            jmp_seg++;
        }
        if(key >= x->forward[i]->DataArray->start){
            x = x->forward[i];
            std::pair<int,node*> res = x->find_subtree(key);
            //this->binarySearch(x,key);
           return res;
        }
        jmp_seg++;
    }    
    return {0,0};

}

void skiplist::Insert(unsigned int key,int value){
    Segment_pt* x = this->header;
    Segment_pt* x_next = nullptr;
    Segment_pt** Update = nullptr;//store the predecessor in different level of locating segment
    Update = (Segment_pt**)malloc(sizeof(Segment_pt*)*(this->MaxLevel+1));
    // for(int i = 1;i<=this->MaxLevel;i++){
    //     Update[i] = this->header;
    // }
    Segment_pt* locate = nullptr;
    for(int i = this->skip_level;i>=1;i--){
        x_next = x->forward[i];
        while(key >= x_next->forward[i]->DataArray->start){
        //if segment's number of element is 1,x turn forward only when key >= stop
        //if segment's number of element > 1,x turn forward only when key > stop
            x = x_next;
            x_next = x_next->forward[i];
        }
        //key < x_next->for[i]->start
        if(key >= x_next->DataArray->stop){
            Update[i] = x_next;
            x = x_next;
            continue;
        }
        else
            Update[i] = x;
        //insure key < x_next stop 
        if(!locate && (key >= x_next->DataArray->start)){
            //(key > x_next->DataArray->start && x_next->DataArray->size == 1)||
            //if segment's number of element is 1,locate segment only when key > start
            //if segment's number of element >1,locate segment only when key>=start
            locate = x_next;
        }
    }

    if(locate){
        // std::pair<int,Segment_pt::Item*> find_ = locate->find_subtree(key);
        // if(find_.first){
        //     (find_.second)->comp.data.value = value;
        //     return;
        // }
        locate->insert_subtree(key,value,Update,this);
    }
    else{
        //create a new segment
        cerr<<"create "<<endl;
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
    int depth = 0;
    unsigned int pre_stop = 0;
    for(int i = 0;i<segCnt;i++){
        x = x->forward[1];
        // cerr<<"Segment_pt level:"<<x->level<<" start:"<<x->DataArray->start<<" stop:"<<x->DataArray->stop<<" slope:"<<x->DataArray->slope<<" intertectpt:"<<x->DataArray->intercept<<"\n";
        // cerr<<"Segment_pt ele cnt"<<x->DataArray->size<<" start:"<<x->DataArray->start<<" stop:"<<x->DataArray->stop;
        if(i && pre_stop != x->DataArray->start){
            cerr<<"?"<<endl;
        }
        pre_stop = x->DataArray->stop;
        typedef std::pair<int, Segment_pt::subtree*> Seg; // <begin, subtree*>
        std::stack<Seg> s;
        s.push(Seg(1, x->DataArray));
        int d = 0;
        while (!s.empty()) {
            Segment_pt::subtree* n = s.top().second;
            int l = s.top().first;
            d = max(d,l);
            s.pop();
            for (int j = 0; j < n->num_items; j ++) {
                if (BITMAP_GET(n->none_bitmap, j) == 0) {
                    if (BITMAP_GET(n->child_bitmap, j) == 1) {
                        s.push(Seg(l+1,n->items[j].comp.child));
                    }
                }
            }
        }
        // cerr<<"\tsegment max depth:"<<d<<endl;
        depth+=d;
    }
    cerr<<"average depth:"<<depth/segCnt<<endl;

}

void skiplist::ComputeSpace(){
    long long space = 0;
    Segment_pt *x = this->header;
    space+=sizeof(this);
    //header
    while (x && x->forward[1] != this->header) {
        space+=x->Computespace();
        x = x->forward[1];
        // cerr<<"segment predict error sum cost on per key:"<<x->DataArray->cost/x->DataArray->size<<endl;
    }

    cerr<<"space size:"<<space<<endl;
    char space_size_data[] = "./space_size.csv";
    string k = to_string(space)+"\n";
    write_into_file(space_size_data,k.c_str());
}






