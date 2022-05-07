#include <cstdio>
#include <cstdlib>

__global__ void count(int *key, int *bucket){
  int i = threadIdx.x;
  atomicAdd(&bucket[key[i]], 1);
}

__global__ void bucketsort(int *key, int *bucket){
  int i = threadIdx.x;
  int j = 0;

  for (int k=0; k<i; k++) {
    j += bucket[k];
  }

  for (; bucket[i]>0; bucket[i]--) {
      key[j++] = i;
    }
}

int main() {
  const int n = 50;
  int range = 5;
  int *key, *bucket;
  cudaMallocManaged(&key, n*sizeof(int));
  cudaMallocManaged(&bucket, range*sizeof(int));
  for (int i=0; i<n; i++) {
    key[i] = rand() % range;
    printf("%d ",key[i]);
  }
  printf("\n");

  for (int i=0; i<range; i++) {
    bucket[i] = 0;
  }

  count<<<1,n>>>(key, bucket);
  bucketsort<<<1,range>>>(key, bucket);
  cudaDeviceSynchronize();
  
  for (int i=0; i<n; i++) {
    printf("%d ",key[i]);
  }
  printf("\n");
  cudaFree(key);
  cudaFree(bucket);
}
