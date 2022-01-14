#include <iostream>
#include <cstdio>
#include <stdlib.h>
#include <cmath>
#include <vector>
#include <algorithm>
#include <ctime>
#include <queue>
#include <string>
#include <parlay/sequence.h>
#include <parlay/parallel.h>
#include <parlay/primitives.h>
#include <parlay/utilities.h>
#include "get_time.h"


#define ARRAYSIZE 1000000000
#define MINFREQ 1
#define MAXFREQ 100000
//#define TEST

using namespace std;
using namespace parlay;

using data_type = long;


data_type total = 0;

timer t_minheap;
timer t_findindex;
timer t_merge;
timer t_add;

struct NODE{
  int index;
  data_type value;

};

void printArray(sequence<data_type> &arr, int size, int start=0){
	for(int i=start;i<start+size;++i){
		cout<<arr[i]<<endl;
	}
	cout<<endl;
}

void initRandomArray(sequence<data_type> &arr, int size, data_type min, data_type max){
	//srand(seed);
	parallel_for(0, size, [&](int i){
		arr[i]=min+hash32(i)%(max-min);
	});
	sort_inplace(arr.cut(0,size));
}

data_type findMinHeap(sequence<data_type> &NodeArray, int &progLeaf, int &progInter, int &interSize){
  t_minheap.start();
	int midSize=0;
	data_type fa[4];
	if(progLeaf<ARRAYSIZE){
		fa[midSize++]=NodeArray[progLeaf];
	}
	if(progLeaf<ARRAYSIZE-1){
		fa[midSize++]=NodeArray[progLeaf+1];
	}
	if(progInter<interSize){
		fa[midSize++]=NodeArray[ARRAYSIZE+progInter];
	}
	if(progInter<interSize-1){
		fa[midSize++]=NodeArray[ARRAYSIZE+progInter+1];
	}
	sort(fa,fa+midSize);
	data_type minFreq = fa[0]+fa[1];	
  t_minheap.stop();
	return minFreq;
}

int buildHeapLayer(data_type Freq, sequence<data_type> &NodeArray, sequence<int> &leftSeq, sequence<int> &rightSeq, int &progLeaf, int &progInter, int &interSize){
	t_findindex.start();
	int leafIndex = lower_bound(NodeArray.begin()+progLeaf,NodeArray.begin()+ARRAYSIZE,Freq+1)-NodeArray.begin();
	int interIndex = lower_bound(NodeArray.begin()+ARRAYSIZE+progInter,NodeArray.begin()+ARRAYSIZE+interSize,Freq+1)-NodeArray.begin();
  t_findindex.stop();

  t_merge.start();
  sequence<NODE> leafSeq(max(1,leafIndex-progLeaf));
	sequence<NODE> interSeq(max(1,interIndex-ARRAYSIZE-progInter));	
	parallel_for(0, leafIndex-progLeaf, [&](int i){
		leafSeq[i].index = i+progLeaf;
    leafSeq[i].value = NodeArray[i+progLeaf];
	});
	parallel_for(0, interIndex-progInter-ARRAYSIZE, [&](int i){
		interSeq[i].index = i+ARRAYSIZE+progInter;
    interSeq[i].value = NodeArray[ARRAYSIZE+progInter+i];
	});
	auto mid = merge(leafSeq.cut(0,leafIndex-progLeaf),interSeq.cut(0,interIndex-ARRAYSIZE-progInter),[&](NODE a, NODE b){return a.value<b.value;});
	progLeaf=leafIndex;
  progInter=interIndex-ARRAYSIZE;
	
	if(mid.size()%2==1){
		if(mid[mid.size()-1].index<ARRAYSIZE)progLeaf--;
		else progInter--;
	}
  t_merge.stop();
  #ifdef TEST
  cout<<"mid"<<endl;
	printArray(mid,mid.size());
  cout<<"***"<<endl;
  #endif
  t_add.start();
  int startPoint = interSize+ARRAYSIZE;
	int addSize = mid.size()/2;
	
	parallel_for(0, addSize, [&](int i){
		leftSeq[interSize+i]=mid[2*i].index;
		rightSeq[interSize+i]=mid[2*i+1].index;
		NodeArray[startPoint+i]=mid[2*i].value+mid[2*i+1].value;
	});
  interSize+=addSize;
	t_add.stop();
	return addSize;
}

bool accuracyTest(sequence<data_type> &NodeArray){
	data_type ans = 0;
  int progInter=0, progLeaf=0;
  int interSize=0;
  vector<data_type> SeqArray(ARRAYSIZE-1);
  while(progInter!=ARRAYSIZE-2){
    int flag=0;
    data_type k;
    data_type mn = numeric_limits<data_type>::max();
    if(progLeaf<ARRAYSIZE-1){
        k = NodeArray[progLeaf]+NodeArray[progLeaf+1];
        if(mn>k){
          mn=k;
          flag=1;
        }
    }
    if(progLeaf<ARRAYSIZE&&progInter<interSize){
      k = NodeArray[progLeaf]+SeqArray[progInter];
      if(mn>k){
        mn=k;
        flag=2;
      }
    }
    if(progInter<interSize-1){
      k = SeqArray[progInter]+SeqArray[progInter+1];
      if(mn>k){
        mn=k;
        flag=3;
      }
    }
    SeqArray[interSize++]=mn;
    ans+=mn;
    if(flag==1)progLeaf+=2;
    if(flag==2){
      progLeaf++;progInter++;
    }
    if(flag==3)progInter+=2;
  }
	
	if(ans==total)return true;
	return false;
}

int main(){
  int ROUND = 5;
  t_minheap.reset();
  t_findindex.reset();
  t_merge.reset();
  t_add.reset();
  timer t;
  t.reset();
  timer t_seq;
  t_seq.reset();

  for(int i=ROUND;i>0;--i){
    sequence<data_type> NodeArray(2*ARRAYSIZE-1);
    sequence<int> leftSeq(ARRAYSIZE);
    sequence<int> rightSeq(ARRAYSIZE);
    int progLeaf = 0, progInter = 0;
    int interSize = 0;
    int layer = 0;
    printf("init\n");
    initRandomArray(NodeArray,ARRAYSIZE,MINFREQ,MAXFREQ);
    //printArray(NodeArray,ARRAYSIZE);
    printf("start\n");
    //timer t_timer;
    t.start();
    while(progInter!=ARRAYSIZE-2){
      data_type minFreq = findMinHeap(NodeArray, progLeaf, progInter, interSize);  
      buildHeapLayer(minFreq, NodeArray, leftSeq, rightSeq, progLeaf, progInter, interSize);
      layer++;
      //t_timer.next("Round "+to_string(layer)+"--AddSize "+to_string(addSize));
      #ifdef TEST
      printArray(NodeArray,interSize,ARRAYSIZE);
      printf("array: ");
      for(int i = 0; i < 2*ARRAYSIZE; i++) {
        printf("%d ", NodeArray[i]);
      }
      puts("");
      #endif
    }
    t.stop();
    #ifdef TEST
    printArray(leftSeq,ARRAYSIZE-1);
    printArray(rightSeq,ARRAYSIZE-1);
    #endif
    total = reduce(NodeArray.cut(ARRAYSIZE,2*ARRAYSIZE-1));
    t_seq.start();
    bool result=accuracyTest(NodeArray);
    t_seq.stop();
    
    if(!result){
      cout<<"WRONG ANSWER"<<endl;
      break;
    }
  }

  cout<<"total:"<<total/ROUND<<endl;
  cout<<"Find min heap time: "<<t_minheap.get_total()/ROUND<<endl;
  cout<<"Find index time: "<<t_findindex.get_total()/ROUND<<endl;
  cout<<"Merge Sequence time: "<<t_merge.get_total()/ROUND<<endl;
  cout<<"Add up time: "<<t_add.get_total()/ROUND<<endl;
  cout<<"parallel time: "<<t.get_total()/ROUND<<endl;

  cout<<"sequential time: "<<t_seq.get_total()/ROUND<<endl;

	return 0;
}
