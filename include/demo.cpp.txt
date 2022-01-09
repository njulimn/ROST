#include <cassert>
#include <iostream>
#include <atomic>
#include <cassert>

using namespace std;


struct Node{ //指针统一用node类型，
    int64_t Key; // next 指针用于构建链表，直接继承
    bool isIndex; //只读变量，指示是否是index node, true-> index node, false sgement node
    std::atomic<Node*> next_; // next 指针用于构建链表，直接继承     
    void AppendNext(Node* next) //链表添加后置元素，无论index sgemeng统一用append next
    {
        //TODO: add cas append
        next_.store(next);
    }
};

struct Index: Node //继承自node
{
    Index(){isIndex = true;};  // 设置指示变量
    std::atomic<Node*> downward_; // 指向底层segment 
    void AppendBelow(Node* below) // 这里node 有可能是index，也有可能是segment，统一这样处理了 
    {
        downward_.store(below);
    }
};

struct Segment_pt: Node // 继承自node
{
    Segment_pt(){isIndex = false;};  // 设置指示变量
};

int main()
{
    Node* head = new Index(); // 新建head节点    
    Node* node1 = new Index(); 
    
    //set head->node1
    head->AppendNext(node1);

    //set head->node1->node2
    Node* node2 = new Index(); 
    node1->AppendNext(node2);    

    //set head->node1->node2
    //            |
    //           node3
    Node* node3 = new Index(); 
    if(node1->isIndex == true) // 为防止错误，先检查node1类型
    {
        reinterpret_cast<Index*>(node1)->AppendBelow(node3);
    }

    //set head->node1->node2
    //            |
    //           node3
    //            |
    //           seg1
    Node* seg_1 =  new Segment_pt();
    if(node3->isIndex == true) // 为防止错误，先检查node1类型
    {
       reinterpret_cast<Index*>(node3)->AppendBelow(seg_1);
    }
    return 0;
}