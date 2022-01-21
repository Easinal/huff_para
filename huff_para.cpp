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

using data_type = long;

constexpr size_t ARRAYSIZE = 1e9+10;
constexpr data_type MINFREQ = 1;
constexpr data_type MAXFREQ = 100000;
constexpr size_t GRANULARITY = 1024;
constexpr int b_size = 10000;
 
//#define TEST

using namespace std;
using namespace parlay;

data_type total = 0;

timer t_minheap;
timer t_findindex;
timer t_merge;
timer t_bisearch;

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

void initRandomArray(sequence<data_type> &arr, int size, data_type min, data_type max, int seed){
	//srand(seed);
	parallel_for(0, size, [&](int i){
		arr[i]=min+hash32(i+seed*ARRAYSIZE)%(max-min);
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

void paraMerge(sequence<data_type> &arr, int left1, int right1, int left2, int right2, int left3){
  int size1 = right1-left1;
  int size2 = right2-left2;
  int b_num = (size1+size2+b_size-1)/b_size;
  int *up = new int[b_num+1];
  int *down = new int[b_num+1];
  up[0]=left1,down[0]=left2;
  up[b_num]=right1,down[b_num]=right2;
  t_bisearch.start();
  parallel_for(1,b_num,[&](int k){
    int size = k*b_size;
    int left = left1, right = min(right1,left1+size);

    int mid=left, i, j;
    while(left<=right){

      i = (left+right)/2;
      data_type v = ((i==right1)?numeric_limits<data_type>::max():arr[i]); 
      j = lower_bound(arr.begin()+left2,arr.begin()+right2,v)-arr.begin();
      if((i+j-left1-left2)>=size){
        mid = i;
        right = i-1;
      }else{
        left = i+1;
      }
    }
    up[k]=mid,down[k]=left1+left2+size-mid;
  });
  t_bisearch.stop();
  t_merge.start();
  parallel_for(0,b_num,[&](int i){
      int progress = i*b_size/2+left3;
      int ending = min((i+1)*b_size/2+left3,left3+(size1+size2)/2);
      int prog1 = up[i];
      int prog2 = down[i];
      int end1 = up[i+1];
      int end2 = down[i+1];
      while(progress!=ending){
      int flag=0;
      data_type k;
      data_type mn = numeric_limits<data_type>::max();
      if(prog1<end1-1){
          k = arr[prog1]+arr[prog1+1];
          if(mn>k){
            mn=k;
            flag=1;
          }
      }
      if(prog1<end1&&prog2<end2){
        k = arr[prog1]+arr[prog2];
        if(mn>k){
          mn=k;
          flag=2;
        }
      }
      if(prog2<end2-1){
        k = arr[prog2]+arr[prog2+1];
        if(mn>k){
          mn=k;
          flag=3;
        }
      }
      arr[progress++]=mn;
      if(flag==1)prog1+=2;
      if(flag==2){
        prog1++;prog2++;
      }
      if(flag==3)prog2+=2;
    }
  });
  t_merge.stop();
  //cout<<"MERGE TIME: "<<t_merge.get_total()<<endl;
  delete[] up;
  delete[] down;
  return;
}

int buildHeapLayer(data_type Freq, sequence<data_type> &NodeArray, sequence<int> &leftSeq, sequence<int> &rightSeq, int &progLeaf, int &progInter, int &interSize, bool flag){
  timer t_detail;
	t_findindex.start();
	int leafIndex = lower_bound(NodeArray.begin()+progLeaf,NodeArray.begin()+ARRAYSIZE,Freq+1)-NodeArray.begin();
	int interIndex = lower_bound(NodeArray.begin()+ARRAYSIZE+progInter,NodeArray.begin()+ARRAYSIZE+interSize,Freq+1)-NodeArray.begin();
  t_findindex.stop();
  
  if(flag)t_detail.next("Find Index:");

  if((leafIndex-progLeaf+interIndex-ARRAYSIZE-progInter)%2==1){
    if(interIndex==progInter+ARRAYSIZE||NodeArray[leafIndex-1]>NodeArray[interIndex-1]){
      leafIndex--;
    }else{
      interIndex--;  
    }
  }
  t_findindex.stop();
  paraMerge(NodeArray,progLeaf,leafIndex,progInter+ARRAYSIZE,interIndex,interSize+ARRAYSIZE);
  int addSize = (leafIndex+interIndex-ARRAYSIZE-progInter-progLeaf)/2;	 
  progLeaf=leafIndex;
  progInter=interIndex-ARRAYSIZE;
  interSize+=addSize;
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
	 
  for(int i=0;i<ARRAYSIZE;++i){
    if(NodeArray[i+ARRAYSIZE]!=SeqArray[i]){
      for(int j=i-5;j<i+5;++j){
        cout<<"j= "<<j<<" Actual Answer: "<<NodeArray[j+ARRAYSIZE]<<" True Answer: "<<SeqArray[j]<<endl;
      }
      cout<<"i= "<<i<<" Actual Answer: "<<NodeArray[i+ARRAYSIZE]<<" True Answer: "<<SeqArray[i]<<endl;
      return false;
    }
  }

  return false;
}

int runParallel(sequence<data_type> &NodeArray, sequence<int> &leftSeq, sequence<int> &rightSeq, timer &t, bool printDetail = false){
  int progLeaf = 0, progInter = 0;
  int interSize = 0;
  int layer = 0;
  timer t_timer;
  t.start();
  while(progInter!=ARRAYSIZE-2){
    bool flag = false;
    data_type minFreq = findMinHeap(NodeArray, progLeaf, progInter, interSize);  
    int addSize = buildHeapLayer(minFreq, NodeArray, leftSeq, rightSeq, progLeaf, progInter, interSize, flag);
    layer++;
    if(printDetail){
      t_timer.next("Round "+to_string(layer)+"--AddSize "+to_string(addSize));
    }
  }
  t.stop();
  return 0;
}

int main(){
  int ROUND = 5;
  sequence<data_type> NodeArray(2*ARRAYSIZE-1);
  sequence<int> leftSeq(ARRAYSIZE);
  sequence<int> rightSeq(ARRAYSIZE);
  printf("init\n");
  initRandomArray(NodeArray,ARRAYSIZE,MINFREQ,MAXFREQ,0);
  //printArray(NodeArray,ARRAYSIZE);
  printf("start\n");
  timer t;


  runParallel(NodeArray, leftSeq, rightSeq,t,true);
  t_minheap.reset();
  t_findindex.reset();
  t_merge.reset();
  t_bisearch.reset();
  t.reset();
  timer t_seq;
  t_seq.reset();
  
  for(int i=0;i<ROUND;++i){
    runParallel(NodeArray,leftSeq, rightSeq,t); 
  }

  total = reduce(NodeArray.cut(ARRAYSIZE,2*ARRAYSIZE-1));
  t_seq.start();
  bool result=accuracyTest(NodeArray);
  t_seq.stop();
  if(!result){
    cout<<"WRONG ANSWER"<<endl;
    return 0;
  }
 
  cout<<"total:"<<total/ROUND<<endl;
  cout<<"Find min heap time: "<<t_minheap.get_total()/ROUND<<endl;
  cout<<"Find index time: "<<t_findindex.get_total()/ROUND<<endl;
  cout<<"Bisearch time: "<<t_bisearch.get_total()/ROUND<<endl;
  cout<<"Merge Sequence time: "<<t_merge.get_total()/ROUND<<endl;
  cout<<"parallel time: "<<t.get_total()/ROUND<<endl;

  cout<<"sequential time: "<<t_seq.get_total()<<endl;

	return 0;
}
