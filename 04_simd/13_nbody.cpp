#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <immintrin.h>

int main() {
  const int N = 8;
  float x[N], y[N], m[N], fx[N], fy[N];
  for(int i=0; i<N; i++) {
    m[i] = drand48();
    fx[i] = fy[i] = 0;
  }
  for(int i=0; i<N; i++) {
  __m256 xvec = _mm256_load_ps(x);
  __m256 yvec = _mm256_load_ps(y);
  __m256 mvec = _mm256_load_ps(m);
  __m256 fxvec = _mm256_load_ps(fx); 
  __m256 fyvec = _mm256_load_ps(fy); 
  //broadcast
  __m256 xivec = _mm256_broadcast_ps(x[i]);
  __m256 yivec = _mm256_broadcast_ps(y[i]);
  //computing
  __m256 rxvec = _mm256_sub_ps(xivec,xvec);
  __m256 ryvec = _mm256_sub_ps(yivec,yvec);
  //rx*rx+ry*ry
  __m256 sumvec = _mm256_add_ps(_mm256_mul_ps(rxvec,rxvec),_mm256_mul_ps(ryvec,ryvec));
  //1/r
  __m256 rinvvec = _mm256_rsqrt_ps(sumvec);
  
  __m256 rinv3vec = _mm256_mul_ps(rinvvec,_mm256_mul_ps(rinvvec,rinvvec));
  __m256 fxvec = _mm256_sub_ps(fxvec,_mm256_mul_ps(rxvec,_mm256_mul(mvec,rinv3vec)));
  _mm256_store_ps(fx, fxvec);
/*
    for(int j=0; j<N; j++) {
      if(i != j) {
        float rx = x[i] - x[j];
        float ry = y[i] - y[j];
        float r = std::sqrt(rx * rx + ry * ry);
        fx[i] -= rx * m[j] / (r * r * r);
        fy[i] -= ry * m[j] / (r * r * r);
      }
    }
*/
    printf("%d %g %g\n",i,fx[i],fy[i]);
  }
}

