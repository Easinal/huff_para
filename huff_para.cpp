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

using data_type = unsigned long long;

constexpr size_t ARRAYSIZE = 1e9+10;
constexpr data_type MINFREQ = 1;
constexpr data_type MAXFREQ = 100000;
constexpr size_t GRANULARITY = 1024;
constexpr size_t b_size = 10000;
 
//#define TEST

using namespace std;
using namespace parlay;

data_type total = 0;

timer t_minheap;
timer t_findindex;
timer t_merge;
timer t_bisearch;

struct NODE{
  size_t index;
  data_type value;

};

void printArray(sequence<data_type> &arr, size_t size, size_t start=0){
	for(size_t i=start;i<start+size;++i){
		cout<<arr[i]<<endl;
	}
	cout<<endl;
}

void exponential_distribution(sequence<data_type> &arr, size_t &size, double lamda) {
  //constexpr int MAX = 10000;
  //sequence<double> tmp(MAX+1);
  //sequence<long long> nums(MAX+1);
  
  parallel_for(0,size,[&](size_t i) {
    arr[i] = 10 * exp(-i*1.0/size);
  });
  /*
  double sum = 0;
  for(int i = 0; i < MAX; i++) {
    tmp[i] = lamda*exp(-lamda*i);
    sum+=tmp[i];
  }
  parallel_for(0,MAX,[&](size_t i){
    nums[i] = round(tmp[i]/sum*size);
  });
  nums[MAX] = 0;
  size = scan_inplace(nums);
  arr = sequence<data_type>(2*size);
  parallel_for(0,MAX,[&](size_t i) {
    parallel_for(nums[i],nums[i+1],[&](size_t j) {
      arr[j] = i+1;
    });
  });
  */
  sort_inplace(arr.cut(0,size));
}

void initRandomArray(sequence<data_type> &arr, size_t &size, data_type min, data_type max, size_t seed){
	//srand(seed);
  arr.resize(2*size);
	if(max==min){
    parallel_for(0,size,[&](size_t i){
      arr[i]=min;
    });
    return;
  }
  parallel_for(0, size, [&](size_t i){
		arr[i]=min+hash64(i+seed*ARRAYSIZE)%(max-min);
	});
	sort_inplace(arr.cut(0,size));
}
void zipfian_distribution(sequence<data_type> &arr, size_t &size){
  size_t MAX = 10000;
  sequence<double> tmp(MAX+1);
  sequence<long long> nums(MAX+1);
  double sum = 0;
  for(size_t i = 0; i < MAX; i++) {
    tmp[i] = 1.0/(i+1);
    sum+=tmp[i];
  }
  parallel_for(0,MAX,[&](size_t i){
    nums[i] = round(tmp[i]/sum*size);
  });
  nums[MAX] = 0;
  size = scan_inplace(nums);
  arr = sequence<data_type>(2*size);

  parallel_for(0,MAX,[&](size_t i) {
    parallel_for(nums[i],nums[i+1],[&](size_t j) {
      arr[j] = i+1;
    });
  });
 sort_inplace(arr.cut(0,size));
}

data_type findMinHeap(sequence<data_type> &NodeArray, size_t size, size_t &progLeaf, size_t &progInter, size_t &interSize){
  t_minheap.start();
	size_t midSize=0;
	data_type fa[4];
	if(progLeaf<size){
		fa[midSize++]=NodeArray[progLeaf];
	}
	if(progLeaf+1<size){
		fa[midSize++]=NodeArray[progLeaf+1];
	}
	if(progInter<interSize){
		fa[midSize++]=NodeArray[size+progInter];
	}
	if(progInter+1<interSize){
		fa[midSize++]=NodeArray[size+progInter+1];
	}
	sort(fa,fa+midSize);
  assert(fa[0] < ULLONG_MAX-fa[1]);
	data_type minFreq = fa[0]+fa[1];	
  t_minheap.stop();
	return minFreq;
}

void paraMerge(sequence<data_type> &arr, size_t left1, size_t right1, size_t left2, size_t right2, size_t left3){
  size_t size1 = right1-left1;
  size_t size2 = right2-left2;
  size_t b_num = (size1+size2+b_size-1)/b_size;
  size_t *up = new size_t[b_num+1];
  size_t *down = new size_t[b_num+1];
  up[0]=left1,down[0]=left2;
  up[b_num]=right1,down[b_num]=right2;
  t_bisearch.start();
  parallel_for(1,b_num,[&](size_t k){
    size_t size = k*b_size;
    size_t left = left1, right = min(right1,left1+size);

    size_t mid=left, i, j;
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
  parallel_for(0,b_num,[&](size_t i){
      size_t progress = i*b_size/2+left3;
      size_t ending = min((i+1)*b_size/2+left3,left3+(size1+size2)/2);
      size_t prog1 = up[i];
      size_t prog2 = down[i];
      size_t end1 = up[i+1];
      size_t end2 = down[i+1];
      while(progress!=ending){
      size_t flag=0;
      data_type k;
      data_type mn = numeric_limits<data_type>::max();
      if(prog1+1<end1){
          assert(arr[prog1]<ULLONG_MAX-arr[prog1+1]);
          k = arr[prog1]+arr[prog1+1];
          if(mn>k){
            mn=k;
            flag=1;
          }
      }
      if(prog1<end1&&prog2<end2){
        assert(arr[prog1]<ULLONG_MAX-arr[prog2]);
        k = arr[prog1]+arr[prog2];
        if(mn>k){
          mn=k;
          flag=2;
        }
      }
      if(prog2+1<end2){
        assert(arr[prog2]<ULLONG_MAX-arr[prog2+1]);
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

size_t buildHeapLayer(data_type Freq, sequence<data_type> &NodeArray, sequence<size_t> &leftSeq, sequence<size_t> &rightSeq, size_t size, size_t &progLeaf, size_t &progInter, size_t &interSize, bool flag){
  timer t_detail;
	t_findindex.start();
	size_t leafIndex = lower_bound(NodeArray.begin()+progLeaf,NodeArray.begin()+size,Freq+1)-NodeArray.begin();
	size_t interIndex = lower_bound(NodeArray.begin()+size+progInter,NodeArray.begin()+size+interSize,Freq+1)-NodeArray.begin();
  t_findindex.stop();
  
  if(flag)t_detail.next("Find Index:");

  if((leafIndex-progLeaf+interIndex-size-progInter)%2==1){
    if(interIndex==progInter+size||NodeArray[leafIndex-1]>NodeArray[interIndex-1]){
      leafIndex--;
    }else{
      interIndex--;  
    }
  }
  t_findindex.stop();
  paraMerge(NodeArray,progLeaf,leafIndex,progInter+size,interIndex,interSize+size);
  size_t addSize = (leafIndex+interIndex-size-progInter-progLeaf)/2;	 
  progLeaf=leafIndex;
  progInter=interIndex-size;
  interSize+=addSize;
  return addSize;
}

bool accuracyTest(sequence<data_type> &NodeArray, size_t size, timer & t_seq){
	data_type ans = 0;
  size_t progInter=0, progLeaf=0;
  size_t interSize=0;
  vector<data_type> SeqArray(size-1);
  t_seq.start();
  while(progInter!=size-2){
    int flag=0;
    data_type k;
    data_type mn = numeric_limits<data_type>::max();
    if(progLeaf+1<size){
        assert(NodeArray[progLeaf]<ULLONG_MAX-NodeArray[progLeaf+1]);
        k = NodeArray[progLeaf]+NodeArray[progLeaf+1];
        if(mn>k){
          mn=k;
          flag=1;
        }
    }
    if(progLeaf<size&&progInter<interSize){
        assert(NodeArray[progLeaf]<ULLONG_MAX-SeqArray[progInter]);
        k = NodeArray[progLeaf]+SeqArray[progInter];
        if(mn>k){
          mn=k;
          flag=2;
        }
    }
    if(progInter+1<interSize){
        assert(SeqArray[progInter]<ULLONG_MAX-SeqArray[progInter+1]);
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
  t_seq.stop();
 	if(ans==total)return true;
	 
  for(size_t i=0;i<size;++i){
    if(NodeArray[i+size]!=SeqArray[i]){
      for(size_t j=i-5;j<i+5;++j){
        cout<<"j= "<<j<<" Actual Answer: "<<NodeArray[j+size]<<" True Answer: "<<SeqArray[j]<<endl;
      }
      cout<<"i= "<<i<<" Actual Answer: "<<NodeArray[i+size]<<" True Answer: "<<SeqArray[i]<<endl;
      return false;
    }
  }

  return false;
}

size_t runParallel(sequence<data_type> &NodeArray, sequence<size_t> &leftSeq, sequence<size_t> &rightSeq, size_t size, timer &t, bool printDetail = false){
  size_t progLeaf = 0, progInter = 0;
  size_t interSize = 0;
  size_t layer = 0;
  timer t_timer;
  t.start();
  while(progInter!=size-2){
    bool flag = false;
    data_type minFreq = findMinHeap(NodeArray, size, progLeaf, progInter, interSize);  
    size_t addSize = buildHeapLayer(minFreq, NodeArray, leftSeq, rightSeq, size, progLeaf, progInter, interSize, flag);
    layer++;
    if(printDetail){
      t_timer.next("Round "+to_string(layer)+"--AddSize "+to_string(addSize));
    }
  }
  t.stop();
  return layer;
}

int run_final(int ROUND, int type, size_t max_freq, size_t size, double lambda=0.001){
  sequence<data_type> NodeArray;
  sequence<size_t> leftSeq;
  sequence<size_t> rightSeq;
  //printf("init\n");
  if(type==0){
    cout<<endl<<"Uniform: size = "<<size<<" max freq = "<<max_freq<<endl;
    initRandomArray(NodeArray,size,MINFREQ,max_freq,0);
  }
  if(type==1){
    exponential_distribution(NodeArray, size, lambda);
    cout<<endl<<"Exponential: size = "<<size<<" lambda = "<<lambda<<endl;
  }
  if(type==2){
     zipfian_distribution(NodeArray, size);
     cout<<endl<<"Zipfian:size = "<<size<<endl;
  }
  leftSeq.resize(size);
  rightSeq.resize(size);
  //printArray(NodeArray,ARRAYSIZE);
  //printf("start\n");
  timer t;


  int k = runParallel(NodeArray, leftSeq, rightSeq, size, t, false);
  t_minheap.reset();
  t_findindex.reset();
  t_merge.reset();
  t_bisearch.reset();
  t.reset();
  timer t_seq;
  t_seq.reset();
  
  cout<<"ROUND: "<<k<<endl;
  for(int i=0;i<ROUND;++i){
    runParallel(NodeArray,leftSeq, rightSeq, size, t);
  }
  cout<<"root: "<<NodeArray[2*size-2]<<endl;

  total = reduce(NodeArray.cut(size,2*size-1));

  bool result=accuracyTest(NodeArray,size,t_seq);

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

int main(){
  int ROUND = 5;
  
  vector<size_t> uniform_size{100000,300000,1000000,3000000,10000000,30000000,100000000,300000000,1000000000,3000000000};
  vector<size_t> freq_maximum{1,10,100,1000,10000,100000,(size_t)1e6,(size_t)1e7,(size_t)1e8,(size_t)1e9,(size_t)1e10};
  vector<double> lambda{0.001,0.005,0.01,0.05,0.1,0.5,1};
  //run_final(1,0,freq_maximum.back(),uniform_size.back());
  /*
  for(size_t i=0;i<uniform_size.size();++i){
    for(size_t k=0;k<freq_maximum.size();++k){
      run_final(ROUND, 0, freq_maximum[k],uniform_size[i]);
    }
  }
  */
  for(size_t i=0;i<lambda.size();++i){
    for(size_t k=0;k<uniform_size.size();++k){
      run_final(ROUND, 1,0,uniform_size[k],lambda[i]);
    }
  }
  /*
  for(size_t k=0;k<uniform_size.size();++k){
    run_final(ROUND, 2,0,uniform_size[k]);
  }
  */
  //run_final(ROUND, 1, (int)1e8, uniform_size.back()); 
  return 0;
}
