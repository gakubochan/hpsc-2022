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
  //array to vector
  __m256 xvec = _mm256_load_ps(x);
  __m256 yvec = _mm256_load_ps(y);
  __m256 mvec = _mm256_load_ps(m);
  __m256 zerovec = _mm256_setzero_ps();
  
  for(int i=0; i<N; i++) {
 　　　　　　__m256 xjvec = _mm256_set1_ps(x[i]);
    __m256 yjvec = _mm256_set1_ps(y[i]);
    
    //rx = x[i]-x[j]
    __m256 rxvec = _mm256_sub_ps(xvec, xjvec);
　　　　　　　　//ry = y[i]-y[j]
    __m256 ryvec = _mm256_sub_ps(yvec, yjvec);
    
    //r^2=rx^2+ry^2
    __m256 r2vec = zerovec;
    r2vec = _mm256_fmadd_ps(rxvec, rxvec, r2vec);
    r2vec = _mm256_fmadd_ps(ryvec, ryvec, r2vec);
　　　　　　　　
    //1/rを求める
    __m256 rinvvec = _mm256_rsqrt_ps(r2vec);
    
    //r2vecの0より大きい要素を持つインデックスには1,それ以外は0
    __m256 mask = _mm256_cmp_ps(r2vec, zerovec, _CMP_NEQ_OQ);
    //i=jに対応するインデックスには０が入っている
    rinvvec = _mm256_blendv_ps(zerovec, rinvvec, mask);
　　　　　　　　//1/r^3を求める
    __m256 rinv2vec = _mm256_mul_ps(rinvvec, rinvvec);
    __m256 rinv3vec = _mm256_mul_ps(rinv2vec, rinv2vec);
    //rx*m[j]/(r^3)
    __m256 fxjvec = _mm256_mul_ps(rxvec, _mm256_mul_ps(mvec, rinv3vec));
    //ry*m[j]/(r^3)
    __m256 fyjvec = _mm256_mul_ps(ryvec, _mm256_mul_ps(mvec, rinv3vec));
    
    //reduction 
    __m256 fxivec = _mm256_permute2f128_ps(fxjvec, fxjvec, 1);
    fxivec = _mm256_add_ps(fxivec, fxjvec);
    fxivec = _mm256_hadd_ps(fxivec, fxivec);
    fxivec = _mm256_hadd_ps(fxivec, fxivec);

    __m256 fyivec = _mm256_permute2f128_ps(fyjvec, fyjvec, 1);
    fyivec = _mm256_add_ps(fyivec, fyjvec);
    fyivec = _mm256_hadd_ps(fyivec, fyivec);
    fyivec = _mm256_hadd_ps(fyivec, fyivec);
    
    //vector to array
    _mm256_store_ps(fx, fxivec);
    _mm256_store_ps(fy, fyivec);

    printf("%d %g %g\n",i,fx[i],fy[i]);
  }
}

