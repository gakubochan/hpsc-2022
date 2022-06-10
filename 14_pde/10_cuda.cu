#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <vector>
#include <cmath>
#include <cstdio>
#include <cstdlib>

//#define BLOCK_X 16
//#define BLOCK_Y 16

using namespace std;

__global__ void cavity(float *u,float *v,float *b,float *p,float *un,float *vn,float *pn,int ny,int nx,int nit,int dx,int dy,float dt,int rho,float nu)
{
	//int nn = ny*nx;
	int i = blockDim.x*blockIdx.x+threadIdx.x+1;
	int j = blockDim.y*blockIdx.y+threadIdx.y+1;

	/*
	float u_j_i+1_ = u[(i+1)*ny+j];
	float u_j_i-1_ = u[(i-1)*ny+j];
	float v_j_i+1_ = v[(i+1)*ny+j];
        float v_j_i-1_ = v[(i-1)*ny+j];
	float u_j+1_i_ = u[i*ny+j+1];
        float u_j-1_i_ = u[i*ny+j-1];
        float v_j+1_i_ = v[i*ny+j+1];
        float v_j-1_i_ = v[i*ny+j-1];
	*/

 	b[i*ny+j] = rho*(1 / dt *\
				((u[(i+1)*ny+j]/*u_j_i+1_*/ - u[(i-1)*ny+j]/*u_j_i-1_*/)/(2*dx)+(v[(i+1)*ny+j]/*v_j_i+1_*/ - v[(i-1)*ny+j]/*v_j_i-1_*/)/(2*dy))-\
				std::pow((u[(i+1)*ny+j]/*u_j_i+1_*/ - u[(i-1)*ny+j]/*u_j_i-1_*/)/(2*dx),2.0) - 2 * ((u[i*ny+j+1]/*u_j+1_i_*/ - u[i*ny+j-1]/*u_j-1_i_*/)/(2*dy) *\
				(v[(i+1)*ny+j]/*v_j_i+1_*/ - v[(i-1)*ny+j]/*v_j_i-1_*/)/(2*dx)) - ((v[i*ny+j+1]/*v_j+1_i_*/ - v[i*ny+j-1]/*v_j-1_i_*/) / std:pow(2*dy,2.0)));
		
	 
	for(int it = 0; it < nit; it++){
		//pn=p;
		for(int j=0;j<ny;j++){
			for(int i=0;i<nx;i++){
				pn[i*ny+j] = p[i*ny+j];
			}
		}
		/*
		float pn_j_i+1_ = pn[(i+1)*ny+j];
        	float pn_j_i-1_ = pn[(i-1)*ny+j];
		float pn_j+1_i_ = pn[i*ny+j+1];
        	float pn_j-1_i_ = pn[i*ny+j-1];
		*/
		p[i*ny+j] = (std::pow(dy,2.0) * (pn[(i+1)*ny+j]/*pn_j_i+1_*/ + pn[(i-1)*ny+j]/*pn_j_i-1_*/) +\
                                        std::pow(dx,2.0) * (pn[i*ny+j+1]/*pn_j+1_i_*/ + pn[i*ny+j-1]/*pn_j-1_i_*/) -\
                                        b[i*ny+j] * std::pow(dx,2.0) * std::pow(dy,2.0))\
                                        / (2 * (std::pow(dx,2.0) + std::pow(dy,2.0)));			
		//p[:, -1] = p[:, -2] & p[:, 0] = p[:, 1]
		p[(ny-1)*ny+j] = p[(ny-2)*ny+j];
		p[j] = p[ny+j];
		//p[0, :] = p[1, :] & p[-1, :] = 0
		p[ny*i] = p[ny*i+1];
		p[ny*i+nx-1] = 0;	
	}

	//un = u;
	for(int j=0;j<ny;j++){
		for(int i=0;i<nx;i++){
			un[i*ny+j] = u[i*ny+j];
		}
	}
	//vn = v;
	for(int j=0;j<ny;j++){
		for(int i=0;i<nx;i++){
			vn[i*ny+j] = v[i*ny+j];
		}
	}
	/*
	float un_j_i_ = un[i*ny+j];
	float un_j_i+1_ = un[(i+1)*ny+j];
	float un_j_i-1_ = un[(i-1)*ny+j];
	float un_j+1_i_ = un[i*ny+j+1];
	float un_j-1_i_ = un[i*ny+j-1];

	float vn_j_i_ = vn[i*ny+j];
	float vn_j_i+1_ = vn[(i+1)*ny+j];
        float vn_j_i-1_ = vn[(i-1)*ny+j];
        float vn_j+1_i_ = vn[i*ny+j+1];
        float vn_j-1_i_ = vn[i*ny+j-1];

	float p_j_i+1_ = p[(i+1)*ny+j];
        float p_j_i-1_ = p[(i-1)*ny+j];
        float p_j+1_i_ = p[i*ny+j+1];
        float p_j-1_i_ = p[i*ny+j-1];
	*/	

	u[i*ny+j] = un[i*ny+j]/*un_j_i_*/ - un[i*ny+j]/*un_j_i_*/ * dt / dx * (un[i*ny+j]/*un_j_i_*/ - un[(i-1)*ny+j]/*un_j_i-1_*/)\
			- un[i*ny+j]/*un_j_i_*/ * dt / dy * (un[i*ny+j]/*un_j_i_*/ - un[i*ny+j-1]/*un_j-1_i_*/)\
			- dt / (2 * rho * dx) * (p[(i+1)*ny+j]/*p_j_i+1_*/ - p[(i-1)*ny+j]/*p_j_i-1_*/)\
			+ nu * dt / std::pow(dx,2.0) * (un[(i+1)*ny+j]/*un_j_i+1_*/ - 2 * un[i*ny+j]/*un_j_i_*/ + un[(i-1)*ny+j]/*un_j_i-1_*/)\
			+ nu * dt / std::pow(dy,2.0) * (un[i*ny+j+1]/*un_j+1_i_*/ - 2 * un[i*ny+j]/*un_j_i_*/ + un[i*ny+j-1]/*un_j-1_i_*/);
	v[i*ny+j] = vn[i*ny+j]/*vn_j_i_*/ - vn[i*ny+j]/*vn_j_i_*/ * dt / dx * (vn[i*ny+j]/*vn_j_i_*/ - vn[(i-1)*ny+j]/*vn_j_i-1_*/)\
                        - vn[i*ny+j]/*vn_j_i_*/ * dt / dy * (vn[i*ny+j]/*vn_j_i_*/ - vn[i*ny+j-1]/*vn_j-1_i_*/)\
                        - dt / (2 * rho * dx) * (p[i*ny+j+1]/*p_j+1_i_*/ - p[i*ny+j-1]/*p_j-1_i_*/)\
                        + nu * dt / std::pow(dx,2.0) * (vn[(i+1)*ny+j]/*vn_j_i+1_*/ - 2 * vn[i*ny+j]/*vn_j_i_*/ + vn[(i-1)*ny+j]/*vn_j_i-1_*/)\
                        + nu * dt / std::pow(dy,2.0) * (vn[i*ny+j+1]/*vn_j+1_i_*/ - 2 * vn[i*ny+j]/*vn_j_i_*/ + vn[i*ny+j-1]/*vn_j-1_i_*/);
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
	float *u,*v,*p,*b,*un,*vn,*pn;
	
	struct timeval st;
	struct timeval et;
	long us;
	
	float *du,*dv,*dp,*db,*dun,*dvn,*dpn;
	
	//allocate region
	/*
	u = (float *)malloc(2*sizeof(float));
	u[0] = (float *)malloc(ny*nx*sizeof(float));
	u[1] = (float *)malloc(ny*nx*sizeof(float));
	v = (float *)malloc(2*sizeof(float));
	v[0] = (float *)malloc(ny*nx*sizeof(float));
	v[1] = (float *)malloc(ny*nx*sizeof(float));
	b = (float *)malloc(2*sizeof(float));
	b[0] = (float *)malloc(ny*nx*sizeof(float));
	b[1] = (float *)malloc(ny*nx*sizeof(float));
	p = (float *)malloc(2*sizeof(float));
	p[0] = (float *)malloc(ny*nx*sizeof(float));
	p[1] = (float *)malloc(ny*nx*sizeof(float));
	*/
	u = (float *)malloc(ny*nx*sizeof(float));
        v = (float *)malloc(ny*nx*sizeof(float));
        b = (float *)malloc(ny*nx*sizeof(float));
	p = (float *)malloc(ny*nx*sizeof(float));
	un = (float *)malloc(ny*nx*sizeof(float));
	vn = (float *)malloc(ny*nx*sizeof(float));
	pn = (float *)malloc(ny*nx*sizeof(float));
	
	gettimeofday(&st, NULL); /* get start time */

	cudaMalloc((void **)&du,2*ny*nx*sizeof(float));
	cudaMalloc((void **)&dv,2*ny*nx*sizeof(float));
	cudaMalloc((void **)&db,2*ny*nx*sizeof(float));
	cudaMalloc((void **)&dp,2*ny*nx*sizeof(float));
	cudaMalloc((void **)&dun,ny*nx*sizeof(float));
	cudaMalloc((void **)&dvn,ny*nx*sizeof(float));
	cudaMalloc((void **)&dpn,ny*nx*sizeof(float));
	
	cudaMemcpy(du, u, 2*ny*nx*sizeof(float),cudaMemcpyHostToDevice);
        cudaMemcpy(dv, v, 2*ny*nx*sizeof(float),cudaMemcpyHostToDevice);
        cudaMemcpy(db, b, 2*ny*nx*sizeof(float),cudaMemcpyHostToDevice);
	cudaMemcpy(dp, p, 2*ny*nx*sizeof(float),cudaMemcpyHostToDevice);
        cudaMemcpy(dun, un, 2*ny*nx*sizeof(float),cudaMemcpyHostToDevice);
        cudaMemcpy(dvn, vn, 2*ny*nx*sizeof(float),cudaMemcpyHostToDevice);
	cudaMemcpy(dpn, pn, 2*ny*nx*sizeof(float),cudaMemcpyHostToDevice);

	dim3 grid(((nx-2)+BLOCK_X-1)/BLOCK_X,((ny-2)+BLOCK_Y-1)/BLOCK_Y, 1);
	dim3 block(BLOCK_X, BLOCK_Y, 1);

	for(int n = 0; n < nt; n++){
		cavity<<grid,block/*ny*nx/16,16*/>>(du,dv,db,dp,dun,dvn,dpn,ny,nx,nit,dx,dy,dt,rho,nu);	
		cudaDeviceSynchronize();
		//current = 1 - current;
		//next = 1 - next;
	}
	cudaMemcpy(u, du, 2*ny*nx*sizeof(float),cudaMemcpyDeviceToHost);
	cudaMemcpy(v, dv, 2*ny*nx*sizeof(float),cudaMemcpyDeviceToHost);
	cudaMemcpy(p, dp, 2*ny*nx*sizeof(float),cudaMemcpyDeviceToHost);
	
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
