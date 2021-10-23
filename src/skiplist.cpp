#pragma once
#include "../include/skiplist.hpp"
#include <fstream>
#include <string>
#include <algorithm>
#include <time.h>
// #include <unistd.h>

#define DEBUG 1
#define Gm 128


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
/////////////////////////////////////////////



////////////////////////////////////

void verify_gamma(double gamma, vector<unsigned int> data, vector<Segment*>seg) {
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

void skiplist::insert_static(vector<Segment_pt*> &Update,Segment* seg,unsigned int st,unsigned int ed,vector<node> input,int level){
    // cerr<<"create new segment...";
    Segment_pt* newSeg = new Segment_pt(input,level,st,ed,seg->slope,seg->intercept);
    // cerr<<"new succeed..";
    // nodeLevel = input[i].level;
    //insert
    if(level > this->level){
        for(int l = level;l>this->level;l--){
            Update[l] = this->header;
        }
        this->level = level;
    }
    for(int l = 1;l<=level;l++){
        newSeg->forward[l] = Update[l]->forward[l];
        Update[l]->forward[l] = newSeg;
        Update[l] = newSeg;
    }
}

void skiplist::setup(vector<snode> input){
    //init
    vector<Segment_pt*> Update(this->MaxLevel+1);
    for(int i = 0;i<=this->MaxLevel;i++){
        Update[i] = this->header;
    }
    int st = 0,ed =0;

    GreedyPLR* plr = new GreedyPLR(Gm);
    vector<Segment*> seg;
    for (int i = 0; i < input.size(); i++) {
		Segment* seg_res = nullptr;
		seg_res = plr->Process(input[i].key, i);
		if(seg_res) {
            seg.push_back(seg_res);
            vector<node> inputPart(ed-st+1);
            for(int l=0;l<ed-st+1;l++){
                inputPart[l].key = input[l+st].key;
                inputPart[l].value = inputPart[l].key;
            }
            this->insert_static(Update,seg_res,input[st].key,input[ed].key,inputPart,input[st].level);
            st = ed = i;
		}
        else{
            ed = i;
        }
	}
    Segment* seg_res = nullptr;
	seg_res = plr->finish();
	if (seg_res) {
        seg.push_back(seg_res);
		vector<node> inputPart(ed-st+1);
        for(int l=0;l<ed-st+1;l++){
            inputPart[l].key = input[l+st].key;
            inputPart[l].value = inputPart[l].key;
        }
        this->insert_static(Update,seg_res,input[st].key,UNINT_MAX,inputPart,input[st].level);
	}
#if DEBUG
    cerr<<"seg size:"<<seg.size()<<endl;
	for (int i = 0; i < seg.size(); i++) {
		cerr << seg[i]->start << " " << seg[i]->stop << " " << seg[i]->slope << " " << seg[i]->intercept << endl;
	}
    cerr<<"insert node"<<endl;
#endif
    // unsigned int nodeLevel = 1;

}

node* skiplist::binearySearch(Segment_pt* x,unsigned int key){
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
    return nullptr;
}

node* skiplist::Search(unsigned int key){
    Segment_pt* x = this->header;
    unsigned int pred;
    bool locate = false;
    for(int i = this->level;i>=1;i--){
        while (1)
        {
            pred = x->forward[i]->slope*key+x->forward[i]->intercept;
            if(pred < x->forward[i]->start){
                return nullptr;
            }
            else if(pred > x->forward[i]->stop){
                x = x->forward[i];
            }
            else{
                locate = true;
                x = x->forward[i];
                break;
            }
        }
        if(locate) break;        
    }
    if(locate){
        node* res = nullptr;
        res = this->binearySearch(x,key);
        if(res){
            return res;
        }
    }
    // else{
    //     //  cerr<<"locate key"<<key<<" fasle"<<endl;
    // }
    
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
        x->forward[1]->show(1);
        x = x->forward[1];
    }
}

void Segment_pt::show(int w){
    char filePath[] = "./log/exp_log.txt";
    string s1 ="Segment_pt level:"+to_string(this->level)+" start:"+to_string(this->start)+" stop:"+to_string(this->stop)+" slope:"+to_string(this->slope)+" intertectpt:"+to_string(this->intercept)+"\n";
    string s2 = "node count:"+to_string(nodes.size())+" node range:"+to_string(nodes[0].key)+" "+to_string(nodes[nodes.size()-1].key)+"\n";
    if(w){
        write_into_file(filePath,s1.c_str());
        write_into_file(filePath,s2.c_str());
    }
    else{
        cerr<<s1<<endl;
        cerr<<s2<<endl;
    }
    

}
