#pragma once
#include "../include/read.hpp"
#include <fstream>
#include <string>
#include <algorithm>
#include <time.h>
// #include <unistd.h>

#define DEBUG 1

#define Gm 5

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
    // cerr<<"seg size:"<<seg.size()<<endl;
	for (int i = 0; i < seg.size(); i++) {
		cerr << seg[i]->start << " " << seg[i]->stop << " " << seg[i]->slope << " " << seg[i]->intercept << endl;
	}
    cerr<<"insert node"<<endl;
#endif
    //insert segment into skiplist,考虑后面优化前面就整合
    vector<Segment_pt*> Update(this->MaxLevel+1);
    
    cerr<<"Updatevector"<<endl;

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
            Segment_pt* newSeg = new Segment_pt(inputPart,nodeLevel,seg[j]->start,seg[j]->stop,seg[j]->slope,seg[j]->intercept);
            // cerr<<"new succeed..";
            nodeLevel = input[i].level;
            //insert
            for(int l = 1;l<=nodeLevel;l++){
                newSeg->forward[l] = Update[l]->forward[l];
                Update[l]->forward[l] = newSeg;
            }
            for(int l = 1;l<=nodeLevel;l++){
                Update[l] = newSeg;
            }
            j++;
            st = ed = i;
            // cerr<<" insert segment success"<<endl;
        }
        // cerr<<"value i"<<i<<" nodeLevel:"<<nodeLevel<<endl;
    }

}

void skiplist::ShowNodeDis(){
    Segment_pt *x = this->header;
    vector<int> Count(this->level+1,0);
    while (x && x->forward[1] != this->header) {
        x->forward[1]->show();
        x = x->forward[1];
    }
}

int main()
{
    // node x;
    // x.key = 1;
    // x.value = 1;
    // vector<node> m;
    // m.push_back(x);
    // Segment_pt* ma = new Segment_pt(m,12,1,5,1,0);

    srand((int)time(0));
    // vector<unsigned int> data_x;
    vector<snode> data;
    unsigned int MaxLevel = 0;
    MaxLevel = readFromCSV(data);
    // cerr<<"Max Level:"<<MaxLevel<<endl;
    skiplist* list = new skiplist(MaxLevel);
    list->setup(data);
    list->ShowNodeDis();
    
}