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

void paraMerge(sequence<data_type> &arr, int left1, int right1, int left2, int right2, int left3){
  int addSize;
  if((left1==right1-1)&&(left2==right2-1)){
    arr[left3]=arr[left1]+arr[left2];
    return;
  }
  if(left1==right1){
    addSize=(right2-left2)/2;
    parallel_for(0,addSize,[&](int i){
      arr[left3+i]=arr[left2+2*i]+arr[left2+2*i+1];
    });
    return;
  }
  if(left2==right2){
    addSize=(right1-left1)/2;
    parallel_for(0,addSize,[&](int i){
      arr[left3+i]=arr[left1+2*i]+arr[left1+2*i+1];
    });
    return;
  }
  int mid = (left1+right1-1)/2;
  int index = lower_bound(arr.begin()+left2,arr.begin()+right2,arr[mid])-arr.begin();
  if(mid==left1 &&index==left2){
    if(mid+1==right1||arr[mid+1]>arr[index]){
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
  if((index+mid-left1-left2)%2==0){
    auto A = [&]() {
      paraMerge(arr,left1,mid,left2,index,left3);
    };
    auto B = [&]() {
      paraMerge(arr,mid,right1,index,right2,left3+(index+mid-left1-left2)/2);
    };
    par_do(A, B);
  }
  else{
    auto A = [&]() {
      paraMerge(arr,left1,mid+1,left2,index,left3);
    };
    auto B = [&]() {
      paraMerge(arr,mid+1,right1,index,right2,left3+(index+mid-left1-left2+1)/2);
    };
    par_do(A, B);
  }
  return;
}

int buildHeapLayer(data_type Freq, sequence<data_type> &NodeArray, sequence<int> &leftSeq, sequence<int> &rightSeq, int &progLeaf, int &progInter, int &interSize, bool flag){
  timer t_detail;
	t_findindex.start();
	int leafIndex = lower_bound(NodeArray.begin()+progLeaf,NodeArray.begin()+ARRAYSIZE,Freq+1)-NodeArray.begin();
	int interIndex = lower_bound(NodeArray.begin()+ARRAYSIZE+progInter,NodeArray.begin()+ARRAYSIZE+interSize,Freq+1)-NodeArray.begin();
  t_findindex.stop();
  
  if(flag)t_detail.next("Find Index:");

  t_merge.start();
  if((leafIndex-progLeaf+interIndex-ARRAYSIZE-progInter)%2==1){
    if(NodeArray[leafIndex-1]>NodeArray[interIndex-1]){
      leafIndex--;
    }else{
      interIndex--;  
    }
  }
  
  paraMerge(NodeArray,progLeaf,leafIndex,progInter+ARRAYSIZE,interIndex,interSize+ARRAYSIZE);
  int addSize = (leafIndex+interIndex-ARRAYSIZE-progInter-progLeaf)/2;	 
  progLeaf=leafIndex;
  progInter=interIndex-ARRAYSIZE;
  
  t_merge.stop();
  return addSize;

  if(flag)t_detail.next("Merge:");

  t_add.start();
  auto mid = merge(NodeArray.cut(progLeaf,leafIndex),NodeArray.cut(ARRAYSIZE+progInter,interIndex),[&](int a, int     b){return a<b;});
  int startPoint = interSize+ARRAYSIZE;
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

  for(int i=0;i<ROUND;++i){
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
    timer t_timer;
    t.start();
    while(progInter!=ARRAYSIZE-2){
      bool flag = false;
      //if(layer==15 || layer==16 || layer==17)flag=true;
      data_type minFreq = findMinHeap(NodeArray, progLeaf, progInter, interSize);  
      int addSize = buildHeapLayer(minFreq, NodeArray, leftSeq, rightSeq, progLeaf, progInter, interSize, flag);
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
