#pragma once
#include "../include/skiplist.hpp"
#include <fstream>
#include <string>
#include <algorithm>
#include <time.h>
// #include <unistd.h>

#define DEBUG 1
#define Gm 128

////////////////////////////////////

void verify_gamma(double gamma, vector<double> data, vector<Segment*>seg) {
	int j = 0;
	Segment* cur = seg[j];
	for (int i = 0; i < data.size(); i++)
	{
		double pred = 0;
		while (j < seg.size()) {
			if (data[i] >= cur->start) {
				if (data[i] < cur->stop) {
					break;
				}
				j++;
				cur = seg[j];
			}
			else {
				cerr << "data out of bound" << endl;
				break;
			}
		}
		if (j < seg.size()) {
			pred = cur->slope * data[i] + cur->intercept;
			if (fabs(pred - i) <= gamma) {
				cerr << "within error tolerence ";
			}
			else {
				cerr << "out of error tolerence ";
			}
			cerr << "data [" << i << "]" << data[i] << " pred:" << pred << endl;
		}
	}
}

void skiplist::setup(vector<snode> input){
    GreedyPLR* plr = new GreedyPLR(Gm);
    vector<Segment*> seg;
    for (int i = 0; i < input.size(); i++) {
		Segment* seg_res = nullptr;
		seg_res = plr->Process(input[i].key, i);
		if (seg_res) {
			seg.push_back(seg_res);
		}
	}
    Segment* seg_res = nullptr;
	seg_res = plr->finish();
	if (seg_res) {
		seg.push_back(seg_res);
	}
#if DEBUG
    cerr<<"seg size:"<<seg.size()<<endl;
	for (int i = 0; i < seg.size(); i++) {
		cerr << seg[i]->start << " " << seg[i]->stop << " " << seg[i]->slope << " " << seg[i]->intercept << endl;
	}
    cerr<<"insert node"<<endl;
#endif
    //insert segment into skiplist,考虑后面优化前面就整合
    vector<Segment_pt*> Update(this->MaxLevel+1);

#if DEBUG
    cerr<<"Updatevector"<<endl;
#endif

    for(int i = 0;i<=this->MaxLevel;i++){
        Update[i] = this->header;
    }
    unsigned int nodeLevel = 1;
    int st = 0,ed =0;
    for(int i=0,j=0;i<input.size()&&j<seg.size();i++){

        // cerr<<"seg number"<<j<<"\t"<<"input[i].key>=seg[j]->start:"<<(input[i].key>=seg[j]->start)<<"\t"<<"input[i].key < seg[j]->stop"<<(input[i].key < seg[j]->stop)<<endl;

        if(i == input.size()-1)
        {
            nodeLevel = max(input[i].level,nodeLevel);
            ed = i;
        }
        if(j<seg.size() && input[i].key>=seg[j]->start && input[i].key < seg[j]->stop && i<input.size()-1){
            nodeLevel = max(input[i].level,nodeLevel);
            ed = i;
        }
        else{
            vector<node> inputPart(ed-st+1);
            for(int l=0;l<ed-st+1;l++){
                inputPart[l].key = input[l+st].key;
                inputPart[l].value = inputPart[l].key;
            }
            // cerr<<"create new segment...";
            Segment_pt* newSeg = new Segment_pt(inputPart,nodeLevel,st,ed,seg[j]->slope,seg[j]->intercept);
            // cerr<<"new succeed..";
            nodeLevel = input[i].level;
            //insert
            for(int l = 1;l<=nodeLevel;l++){
                newSeg->forward[l] = Update[l]->forward[l];
                Update[l]->forward[l] = newSeg;
                Update[l] = newSeg;
            }
            j++;
            st = ed = i;
            // cerr<<" insert segment success"<<endl;
        }
        // cerr<<"value i"<<i<<" nodeLevel:"<<nodeLevel<<endl;
    }

}

node* skiplist::Search(unsigned int key){
    Segment_pt* x = this->header;
    bool locate = false;
    for(int i = this->level;i>=1;i--){
        while (1)
        {
            unsigned int pred = x->forward[i]->slope*key+x->forward[i]->intercept;
            if(pred < x->forward[i]->start){
                return nullptr;
            }
            else if(pred > x->forward[i]->stop){
                x = x->forward[i];
            }
            else{
                locate = true;
                break;
            }
        }
        if(locate) break;        
    }
    if(locate){
        x = x->forward[1];
        //bineary search
        int l=0,r=x->nodes.size()-1,mid;
        while(l<=r){
            mid = (l+r)>>1;
            if(x->nodes[mid].key == key){
                return &(x->nodes[mid]);
            }
            else if(x->nodes[mid].key > key){
                r = mid-1;
            }
            else{
                l = mid+1;
            }
        }
    }
    
    return nullptr;

}

void skiplist::ComputeSpace(){
    unsigned int space = 0;
    Segment_pt *x = this->header;
    while (x && x->forward[1] != this->header) {
        // cerr<<"node:"<<sizeof(*(x->forward[1]))<<endl;
        space+=(sizeof(*(x->forward[1])));
        space+=((x->forward[1]->forward).capacity() * sizeof(Segment_pt*));
        space+=((x->forward[1]->nodes).capacity() * sizeof(node));
        x = x->forward[1];
    }
    cerr<<"space size:"<<space<<endl;
}

void skiplist::ShowNodeDis(){
    Segment_pt *x = this->header;
    vector<int> Count(this->level+1,0);
    while (x && x->forward[1] != this->header) {
        x->forward[1]->show();
        x = x->forward[1];
    }
}
