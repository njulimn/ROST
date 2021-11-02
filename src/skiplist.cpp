#pragma once
#include "../include/skiplist.hpp"
#include <fstream>
#include <string>
#include <algorithm>
#include <time.h>
#include <time.h>
// #include <unistd.h>

#define DEBUG 0 
// #define Gm 128

// map<Segment_pt*,node*> Seg2Data;

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

void skiplist::insert_static(vector<Segment_pt*> &Update,Segment* seg,unsigned int st,unsigned int ed,node* input,int size,int level){
    Segment_pt* newSeg = new Segment_pt(size,input,level,st,ed,seg->slope,seg->intercept);
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
    // Seg2Data[newSeg] = input;
    this->size++;
}

void skiplist::setup(vector<snode> input){
    //init
    vector<Segment_pt*> Update(this->MaxLevel+1);
    for(int i = 0;i<=this->MaxLevel;i++){
        Update[i] = this->header;
    }
    int st = 0,ed =0;

    srand((int)time(0));

    GreedyPLR* plr = new GreedyPLR(this->gamma);
    vector<Segment*> seg;
    vector<int> segment_stIndex;
    for (int i = 0; i < input.size(); i++) {
		Segment* seg_res = nullptr;
		seg_res = plr->Process(input[i].key, i-st);
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
		// vector<node> inputPart(ed-st+1);
        // for(int l=0;l<ed-st+1;l++){
        //     inputPart[l].key = input[l+st].key;
        //     inputPart[l].value = inputPart[l].key;
        // }
        // this->insert_static(Update,seg_res,input[st].key,UNINT_MAX,inputPart,input[st].level);
	}
    cerr<<"seg size:"<<seg.size()<<endl;
    this->MaxLevel = (log(seg.size())/log(2));

    //segment_data = (int**)malloc(sizeof(int*)*seg.size());//申请seg.size个空间存储每段数据的数组首地址

    for(int i=0;i<seg.size();i++){
        st = segment_stIndex[i];
        if(i == seg.size()-1)
            ed = input.size()-1;
        else
            ed = segment_stIndex[i+1]-1;
        node* data_segement = nullptr;
        data_segement = (node*)malloc(sizeof(node)*(ed-st+1));
        // vector<node> inputPart(ed-st+1);
        for(int l=0;l<ed-st+1;l++){
            data_segement[l].key = input[l+st].key;
            data_segement[l].value = data_segement[l].key;
        }
        int level = this->RandLevel();
        this->insert_static(Update,seg[i],input[st].key,input[ed].key,data_segement,ed-st+1,level);
    }

#if DEBUG
    
	for (int i = 0; i < seg.size(); i++) {
		cerr << seg[i]->start << " " << seg[i]->stop << " " << seg[i]->slope << " " << seg[i]->intercept << endl;
	}
    cerr<<"insert node"<<endl;
#endif
    // unsigned int nodeLevel = 1;

}

node* skiplist::binarySearch(Segment_pt* x,unsigned int key){
    int pred = x->slope*key+x->intercept;
    int l=max((pred-this->gamma),0);
    // node* data = Seg2Data[x];
    int r=min(x->node_size-1,pred+this->gamma),mid;
    while(l<=r){
        mid = l+(r-l)/2;
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
    // unsigned int pred;
    for(int i = this->level;i>=1;i--){
        while(key > x->forward[i]->stop){
            x = x->forward[i];
        }
        if(key >= x->forward[i]->start){
            x = x->forward[i];
            node* res = nullptr;
            res = this->binarySearch(x,key);
            if(res){
                return res;
            }
        }
    }    
    return nullptr;

}

void skiplist::ComputeSpace(){
    unsigned int space = 0;
    Segment_pt *x = this->header;
    space+=sizeof(this);
    while (x && x->forward[1] != this->header) {
        // cerr<<"node:"<<sizeof(*(x->forward[1]))<<endl;
        space+=(sizeof(*(x->forward[1])));
        space+=((this->level+1) * sizeof(Segment_pt*));
        space+=((x->forward[1]->node_size) * sizeof(node));
        x = x->forward[1];
    }
    // space+=( (sizeof(Segment_pt*)+sizeof(node*) )*Seg2Data.size() );

    cerr<<"space size:"<<space<<endl;
}

void skiplist::ShowNodeDis(){
    Segment_pt *x = this->header;
    char filePath[] = "./log/exp_log_dis.txt";
    char filePath2[] = "./log/exp_log_plr.csv";
    // char filePath2[] = "./log/exp_log_dis2.txt";
    vector<int> Count(this->level+1,0);
    map<int,int> nodemap;
    vector<vector<int>> jump(this->level+1,vector<int>(1));
    int i = 1;
    while (x && x->forward[1] != this->header) {
        x->forward[1]->show(1);
        if(!Count[x->forward[1]->level]){
            jump[x->forward[1]->level][0] = i;
        }
        else{
            jump[x->forward[1]->level].push_back(i);
        }
        Count[x->forward[1]->level]++;
        if(!nodemap.count(x->forward[1]->node_size)){
            nodemap[x->forward[1]->node_size] = 1;
        }
        else{
            nodemap[x->forward[1]->node_size]++;
        }
        i++;
        x = x->forward[1];
    }
    
    string s1 = "NodeDistribution---------------------\nthe number of segment:"+to_string (this->size)+"\n";
    write_into_file(filePath,s1.c_str());

    map<int,int>::iterator it;
    for(it = nodemap.begin();it!=nodemap.end();it++){
        string sit = to_string((*it).first)+","+to_string((*it).second)+"\n";
        write_into_file(filePath2,sit.c_str());
    }

    for(int i=1;i<=this->level;i++){
        string s2 = "level "+to_string (i)+",segment count "+to_string(Count[i])+"\n";
        write_into_file(filePath,s2.c_str());
        for(int j = 0;j<jump[i].size();j++){
            string s3 = to_string(jump[i][j])+"\t";
            write_into_file(filePath,s3.c_str());
        }
        write_into_file(filePath,"\n");
    }
}

void Segment_pt::show(int w){
    char filePath[] = "./log/exp_log.txt";
    string s1 ="Segment_pt level:"+to_string(this->level)+" start:"+to_string(this->start)+" stop:"+to_string(this->stop)+" slope:"+to_string(this->slope)+" intertectpt:"+to_string(this->intercept)+"\n";
    string s2 = "node count:"+to_string(this->node_size)+" node range:"+to_string(this->nodes[0].key)+" "+to_string(this->nodes[this->node_size-1].key)+"\n";
    if(w==1){
        write_into_file(filePath,s1.c_str());
        write_into_file(filePath,s2.c_str());
    }
    else if(w==0){
        cerr<<s1;
        cerr<<s2<<endl;
    }
    else if(w == 2){
        string s3 ="Segment_pt level:"+to_string(this->level)+" start:"+to_string(this->start)+" stop:"+to_string(this->stop)+" node count:"+to_string(this->node_size)+"\n";
        cerr<<s3;
    }
    

}
