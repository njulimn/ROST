
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
// #include "./fileOperation.hpp"

#define p 0.5
#define UNINT_MAX 0xffffffff
#define Epslion 1e-8

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
		if (!this->rho_lower->AboveAssert(k) && this->rho_upper->BeblowAssert(k)) {
			//重新开一个段
			Segment* current = this->CurrentSegment(k->x);
			delete this->s0;
			this->s0 = nullptr;
			this->s0 = k;
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
		this->s_last = nullptr;
		this->s_last = new Point(x, y);
		Segment* res = nullptr;
		// cerr<<"in process:old state "<<this->state<<endl;
		GreedyState newState = this->state;
		
		// cerr<<"new state(old)："<<newState;
		// cerr<<"GreedyState::Need2"<<GreedyState::Need2<<endl;
		if (this->state == GreedyState::Need2) {
			delete this->s0;
			this->s0 = nullptr;
			this->s0 = this->s_last;
			newState = GreedyState::Need1;
			this->state = GreedyState::Need1;
		}
		else if (this->state == GreedyState::Need1) {
			//delete this->s1;
			this->s1 = nullptr;
			this->s1 = this->s_last;
			this->setup(this->s0, this->s1);
			newState = GreedyState::Ready;
			this->state = GreedyState::Ready;
		}
		else if (this->state == GreedyState::Ready) {
			res = this->Process_pt(this->s_last);
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

class Segment_pt :public Segment{
    public:
        vector<Segment_pt*> forward;
        unsigned level;
        unsigned node_size;
        vector<node> nodes;//内部二分查找

        Segment_pt();
		Segment_pt(int level,unsigned int strt,unsigned int stp,double slp,double interc):Segment(strt,stp,slp,interc){
			// cerr<<"????????"<<endl;
            this->level = level;
            vector<Segment_pt*> new_seg(level+1);
            this->forward = new_seg;
            for(int i=1;i<=level;i++){
                this->forward[i] = this;
            }
			// cerr<<"@@@@@"<<endl;
        }

        Segment_pt(vector<node> input,int level,unsigned int strt,unsigned int stp,double slp,double interc):Segment(strt,stp,slp,interc)
		{
			// cerr<<"????????"<<endl;
            this->level = level;
            vector<Segment_pt*> new_seg(level+1);
            this->forward = new_seg;
            for(int i=1;i<=level;i++){
                this->forward[i] = this;
            }
			nodes.resize(input.size());
			nodes.assign(input.begin(),input.end());
			// cerr<<"new finish"<<endl;
        }

		void show(int w);

};

class skiplist {
    public:
        int level;
        int size;
        Segment_pt *header;
        int MaxLevel;
		int gamma;

        skiplist();
        skiplist(int MaxLevel,int gamma){
            int i;
            Segment_pt *header = new Segment_pt(MaxLevel,UNINT_MAX,UNINT_MAX,0,0);//(snode *)malloc(sizeof(struct snode));
            this->header = header;
            // for (i = 0; i <= MaxLevel; i++) {
            //     header->forward[i] = this->header;
            // }
            this->level = 1;
            this->size = 0;
            this->MaxLevel = MaxLevel;
			this->gamma = gamma;
        }

		void setup(vector<snode> input);

		node* binearySearch(Segment_pt* x,unsigned int key);
		node* Search(unsigned int key);

		void insert_static(vector<Segment_pt*> &Update,Segment* seg,unsigned int st,unsigned int ed,vector<node> input,int level);

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






