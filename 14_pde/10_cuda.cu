#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <vector>
#include <cmath>
#include <cstdio>
#include <cstdlib>

#define BLOCK_X 16
#define BLOCK_Y 16

//using namespace std;

__device__ float pow2(float x)
{
	return x*x;
}


__global__ void cavity(float *u,float *v,float *b,float *p,float *un,float *vn,float *pn,int ny,int nx,int nit,int dx,int dy,float dt,int rho,float nu)
{
	int i = blockDim.x*blockIdx.x+threadIdx.x+1;
	int j = blockDim.y*blockIdx.y+threadIdx.y+1;

	if((j < ny-1) && (i < nx-1)){	

 		b[i*ny+j] = rho*(1 / dt *\
				((u[(i+1)*ny+j] - u[(i-1)*ny+j])/(2*dx)+(v[(i+1)*ny+j] - v[(i-1)*ny+j])/(2*dy))-\
				pow2((u[(i+1)*ny+j] - u[(i-1)*ny+j])/(2*dx)) - 2 * ((u[i*ny+j+1] - u[i*ny+j-1])/(2*dy) *\
				(v[(i+1)*ny+j] - v[(i-1)*ny+j])/(2*dx)) - ((v[i*ny+j+1] - v[i*ny+j-1]) / pow2(2*dy)));
		
	 
		for(int it = 0; it < nit; it++){
			//pn=p;
			pn[i*ny+j] = p[i*ny+j];
			/*
			for(int j=0;j<ny;j++){
				for(int i=0;i<nx;i++){
					pn[i*ny+j] = p[i*ny+j];
				}
			}		
			*/


			p[i*ny+j] = (pow2(dy) * (pn[(i+1)*ny+j] + pn[(i-1)*ny+j]) +\
					pow2(dx) * (pn[i*ny+j+1] + pn[i*ny+j-1]) -\
					b[i*ny+j] * pow2(dx) * pow2(dy))\
					/ (2 * (pow2(dx) + pow2(dy)));			
			//p[:, -1] = p[:, -2] & p[:, 0] = p[:, 1]
			p[(ny-1)*ny+j] = p[(ny-2)*ny+j];
			p[j] = p[ny+j];
			//p[0, :] = p[1, :] & p[-1, :] = 0
			p[ny*i] = p[ny*i+1];
			p[ny*i+nx-1] = 0;	
		}

		//un = u;
		un[i*ny+j] = u[i*ny+j];
		/*
		for(int j=0;j<ny;j++){
			for(int i=0;i<nx;i++){
				un[i*ny+j] = u[i*ny+j];
			}
		}
		*/

		//vn = v;
		vn[i*ny+j] = v[i*ny+j];
		/*
		for(int j=0;j<ny;j++){
			for(int i=0;i<nx;i++){
				vn[i*ny+j] = v[i*ny+j];
			}
		}
		*/


		u[i*ny+j] = un[i*ny+j] - un[i*ny+j] * dt / dx * (un[i*ny+j] - un[(i-1)*ny+j])\
				- un[i*ny+j] * dt / dy * (un[i*ny+j] - un[i*ny+j-1])\
				- dt / (2 * rho * dx) * (p[(i+1)*ny+j] - p[(i-1)*ny+j])\
				+ nu * dt / pow2(dx) * (un[(i+1)*ny+j] - 2 * un[i*ny+j] + un[(i-1)*ny+j])\
				+ nu * dt / pow2(dy) * (un[i*ny+j+1] - 2 * un[i*ny+j] + un[i*ny+j-1]);
		v[i*ny+j] = vn[i*ny+j] - vn[i*ny+j] * dt / dx * (vn[i*ny+j] - vn[(i-1)*ny+j])\
                	        - vn[i*ny+j] * dt / dy * (vn[i*ny+j] - vn[i*ny+j-1])\
                        	- dt / (2 * rho * dx) * (p[i*ny+j+1] - p[i*ny+j-1])\
                        	+ nu * dt / pow2(dx) * (vn[(i+1)*ny+j] - 2 * vn[i*ny+j] + vn[(i-1)*ny+j])\
                        	+ nu * dt / pow2(dy) * (vn[i*ny+j+1] - 2 * vn[i*ny+j] + vn[i*ny+j-1]);
		//u[0, :]  = 0 & u[-1, :] = 1 & v[0, :]  = 0 & v[-1, :] = 0
		u[i*ny] = 0;
		u[i*ny + nx-1] = 1;
		v[i*ny] = 0;
		v[i*ny + nx-1] = 0;

		//u[:, 0]  = 0 & u[:, -1] = 0 & v[:, 0]  = 0 & v[:, -1] = 0
		u[j] = 0;
        	u[i*(ny-1) + j] = 0;
        	v[j] = 0;
        	v[i*(ny-1) + j] = 0;

 		return;
	}	
}


long time_diff_us(struct timeval st, struct timeval et)
{
    return (et.tv_sec-st.tv_sec)*1000000+(et.tv_usec-st.tv_usec);
}

int main(){
	int nx = 41;
	int ny = 41;
	int nt = 500;
	int nit = 50;
	int dx = 2 / (nx-1);
	int dy = 2 / (ny-1);
	float dt = .01;
	int rho = 1;
	float nu = .02;
	size_t size = ny*nx*sizeof(float);
	float *u,*v,*p,*b,*un,*vn,*pn;
	
	struct timeval st;
	struct timeval et;
	long us;
	
	float *du,*dv,*dp,*db,*dun,*dvn,*dpn;
	
	//allocate region
	u = (float *)malloc(size);
        v = (float *)malloc(size);
        b = (float *)malloc(size);
	p = (float *)malloc(size);
	un = (float *)malloc(size);
	vn = (float *)malloc(size);
	pn = (float *)malloc(size);
	
	gettimeofday(&st, NULL); /* get start time */

	cudaMalloc((void **)&du,size);
	cudaMalloc((void **)&dv,size);
	cudaMalloc((void **)&db,size);
	cudaMalloc((void **)&dp,size);
	cudaMalloc((void **)&dun,size);
	cudaMalloc((void **)&dvn,size);
	cudaMalloc((void **)&dpn,size);
	
	cudaMemcpy(du, u, size,cudaMemcpyHostToDevice);
        cudaMemcpy(dv, v, size,cudaMemcpyHostToDevice);
        cudaMemcpy(db, b, size,cudaMemcpyHostToDevice);
	cudaMemcpy(dp, p, size,cudaMemcpyHostToDevice);
        cudaMemcpy(dun, un, size,cudaMemcpyHostToDevice);
        cudaMemcpy(dvn, vn, size,cudaMemcpyHostToDevice);
	cudaMemcpy(dpn, pn, size,cudaMemcpyHostToDevice);

	dim3 grid(((nx-2)+BLOCK_X-1)/BLOCK_X,((ny-2)+BLOCK_Y-1)/BLOCK_Y, 1);
	dim3 block(BLOCK_X, BLOCK_Y, 1);

	for(int n = 0; n < nt; n++){
		cavity<<<grid,block>>>(du,dv,db,dp,dun,dvn,dpn,ny,nx,nit,dx,dy,dt,rho,nu);	
		cudaDeviceSynchronize();
	}

	cudaMemcpy(u, du, size,cudaMemcpyDeviceToHost);
	cudaMemcpy(v, dv, size,cudaMemcpyDeviceToHost);
	cudaMemcpy(p, dp, size,cudaMemcpyDeviceToHost);
	
	cudaFree(du);
	cudaFree(dv);
	cudaFree(db);
	cudaFree(dp);
	cudaFree(dun);
	cudaFree(dvn);
	cudaFree(dpn);

	gettimeofday(&et, NULL); /* get end time */

        us = time_diff_us(st, et);
        printf("Matmul took %ld us\n",us);
}
