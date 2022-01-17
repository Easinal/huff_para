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


#define ARRAYSIZE 100000000
#define MINFREQ 1
#define MAXFREQ 100
#define GRANULARITY 1024
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

void initRandomArray(sequence<data_type> &arr, int size, data_type min, data_type max, int seed = 0){
	//srand(seed);
	parallel_for(0, size, [&](int i){
		arr[i]=hash32(i+seed*size);
    //arr[i]=MINFREQ+rand()%(MAXFREQ-MINFREQ);
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

int buildHeapLayer(data_type Freq, sequence<data_type> &NodeArray, sequence<int> &leftSeq, sequence<int> &rightSeq, int &progLeaf, int &progInter, int &interSize, bool flag){
  timer t_detail;
	t_findindex.start();
	int leafIndex = lower_bound(NodeArray.begin()+progLeaf,NodeArray.begin()+ARRAYSIZE,Freq+1)-NodeArray.begin();
	int interIndex = lower_bound(NodeArray.begin()+ARRAYSIZE+progInter,NodeArray.begin()+ARRAYSIZE+interSize,Freq+1)-NodeArray.begin();
  t_findindex.stop();
  
  if(flag)t_detail.next("Find Index:");

  t_merge.start();
  auto mid = merge(NodeArray.cut(progLeaf,leafIndex),NodeArray.cut(ARRAYSIZE+progInter,interIndex),[&](int a, int     b){return a<b;});
  if((leafIndex-progLeaf+interIndex-ARRAYSIZE-progInter)%2==0){
    progLeaf=leafIndex;
    progInter=interIndex-ARRAYSIZE;
  }
  else if(NodeArray[leafIndex-1]>NodeArray[interIndex-1]){
    progLeaf=leafIndex-1;
    progInter=interIndex-ARRAYSIZE;
  }else{
    progLeaf=leafIndex;
    progInter=interIndex-1-ARRAYSIZE;
  }
  
  t_merge.stop();
   
  if(flag)t_detail.next("Merge:");

  t_add.start();
  int startPoint = interSize+ARRAYSIZE;
	int addSize = mid.size()/2;
	
	parallel_for(0, addSize, [&](int i){
		NodeArray[startPoint+i]=mid[i<<1]+mid[i<<1|1];
	}, GRANULARITY);
  interSize+=addSize;
	t_add.stop();
	 
  if(flag)t_detail.next("Add:");

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

void paraMerge(sequence<data_type> &arr1, sequence<data_type> &arr2, sequence<data_type> &arr3, int left1, int right1, int left2, int right2, int left3, int dep=0){
  int addSize;
  if((right1+right2-left1-left2)%2==1){
  //if((left3==ARRAYSIZE-1)){
  // cout<<"ARR1: left1: "<<left1<<" right1: "<<right1<<endl;
  // printArray(arr1,right1-left1,left1);
  // cout<<"ARR2: left2: "<<left2<<" right2: "<<right2<<endl;
  // printArray(arr2,right2-left2,left2);
  // cout<<"ARR3: left3: "<<left3<<endl;
  // printArray(arr3,left3);
  }
  if((left1==right1-1)&&(left2==right2-1)){
    arr3[left3]=arr1[left1]+arr2[left2];
    return;
  }
  if(left1==right1){
    addSize=(right2-left2)/2;
    parallel_for(0,addSize,[&](int i){
      arr3[left3+i]=arr2[left2+2*i]+arr2[left2+2*i+1];
    });
    return;
  }
  if(left2==right2){
    addSize=(right1-left1)/2;
    parallel_for(0,addSize,[&](int i){
      arr3[left3+i]=arr1[left1+2*i]+arr1[left1+2*i+1];
    });
    return;
  }
  int mid = (left1+right1-1)/2;
  int index = lower_bound(arr2.begin()+left2,arr2.begin()+right2,arr1[mid])-arr2.begin();
  if(mid==left1 &&index==left2){
    if(mid+1==right1||arr1[mid+1]>arr2[index]){
      index++;
    }else{
      mid++;
    }
  }
  if(mid==left1 &&index==right2){
    int gap = right1-left1;
    if(gap%2==0){
      mid--;
    }else{
      index--;
    }
  }
  //cout<<"pivot: "<<mid<<"--"<<arr1[mid]<<endl;
  if((index+mid-left1-left2)%2==0){
    paraMerge(arr1,arr2,arr3,left1,mid,left2,index,left3,dep+1);
    paraMerge(arr1,arr2,arr3,mid,right1,index,right2,left3+(index+mid-left1-left2)/2,dep+1);
    return;
  }
  else{
    paraMerge(arr1,arr2,arr3,left1,mid+1,left2,index,left3,dep+1);
    paraMerge(arr1,arr2,arr3,mid+1,right1,index,right2,left3+(index+mid-left1-left2+1)/2,dep+1);
  }
  /*
  if((index+mid-left1-left2)%2==0){
    auto A = [&]() {
      paraMerge(arr1,arr2,arr3,left1,mid,left2,index,left3);
    };
    auto B = [&]() {
      paraMerge(arr1,arr2,arr3,mid,right1,index,right2,left3+(index+mid-left1-left2)/2);
    };
    par_do(A, B);
  }
  else{
    auto A = [&]() {
      paraMerge(arr1,arr2,arr3,left1,mid-1,left2,index,left3);
    };
    auto B = [&]() {
      paraMerge(arr1,arr2,arr3,mid-1,right1,index,right2,left3+(index+mid-left1-left2-1)/2);
    };
    par_do(A, B);
  }
  */
  return;
}

int main(){
  int ROUND = 1;
  timer t_timer;
  t_timer.reset();
  sequence<data_type> NodeArray(2*ARRAYSIZE);
  sequence<data_type> leftSeq(ARRAYSIZE);
  sequence<data_type> rightSeq(ARRAYSIZE);
  printf("init\n");
  initRandomArray(leftSeq,ARRAYSIZE,  MINFREQ,MAXFREQ,1);
  initRandomArray(rightSeq,ARRAYSIZE,  MINFREQ,MAXFREQ,2);
  //printArray(leftSeq,ARRAYSIZE);
  //printArray(rightSeq,ARRAYSIZE);
  printf("start\n");
  
  auto mid = merge(leftSeq,rightSeq,[&](data_type a, data_type  b){return a<b;});
  int startPoint = ARRAYSIZE;
  parallel_for(0, ARRAYSIZE, [&](int i){
    NodeArray[startPoint+i]=mid[i<<1]+mid[i<<1|1];
  }, GRANULARITY);
  // printArray(mid,2*ARRAYSIZE);

  for(int i=0;i<ROUND;++i){
    t_timer.reset();
    t_timer.start();
    paraMerge(leftSeq,rightSeq,NodeArray,0,ARRAYSIZE,0,ARRAYSIZE,0);
    t_timer.stop();
    //printArray(NodeArray,2*ARRAYSIZE);
   }
  for(int i=0;i<ARRAYSIZE;++i){
    if(NodeArray[i]!=NodeArray[ARRAYSIZE+i]){
      cout<<"FALSE"<<endl;
      return 0;
    }
  }
  cout<<"TRUE"<<endl;
 return 0;
}
