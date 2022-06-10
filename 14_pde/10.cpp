#define _GLIBCXX_DEBUG

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <vector>
#include <cmath>


using namespace std;

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
	
	struct timeval st;
	struct timeval et;
	long us;

	float u[ny][nx];
	float v[ny][nx];
	float p[ny][nx];
	float b[ny][nx];
        float un[ny][nx];
        float vn[ny][nx];
        float pn[ny][nx];
	
	gettimeofday(&st, NULL); /* get start time */

	int n;

	for(n = 0; n < nt; n++){
		for(int j =1; j < ny-1; j++){
			for(int i =1; i < nx-1; i++){
				b[j][i] = rho * (1 / dt *\
					((u[j][i+1] - u[j][i-1]) / (2 * dx) + (v[j+1][i] - v[j-1][i]) / (2 * dy)) -\
		     			std::pow((u[j][i+1] - u[j][i-1]) / (2 * dx), 2.0) - 2 * ((u[j+1][i] - u[j-1][i]) / (2 * dy) *\
		     			(v[j][i+1] - v[j][i-1]) / (2 * dx)) - ((v[j+1][i] - v[j-1][i]) / std::pow(2 * dy, 2.0)));
			}
		}
		
		for(int it = 0; it < nit; it++){
			//pn = p;
			for(int j =0; j < ny; j++){
                        	for(int i =0; i < nx; i++){
					pn[j][i] = p[j][i];
				}
			}
			for(int j =1; j < ny-1; j++){
	            		for(int i =1; i < nx-1; i++){
					p[j][i] = (std::pow(dy,2.0) * (pn[j][i+1] + pn[j][i-1]) +\
                           		std::pow(dx,2.0) * (pn[j+1][i] + pn[j-1][i]) -\
                           		b[j][i] * std::pow(dx,2.0) * std::pow(dy,2.0))\
                          		/ (2 * (std::pow(dx,2.0) + std::pow(dy,2.0)));						
				}
			}
			//p[:, -1] = p[:, -2] & p[:, 0] = p[:, 1]
			for(int j =0; j < ny; j++){
				p[j][ny-1] = p[j][ny-2];
				p[j][0] = p[j][1];
			}
			//p[0, :] = p[1, :] & p[-1, :] = 0
			for(int i =0; i < nx; i++){
				p[0][i] = p[1][i];
				p[nx-1][i] = 0;
			}   
		}
		//un = u;
		for(int j =0; j < ny; j++){
			for(int i =0; i < nx; i++){
				un[j][i] = u[j][i];
			}
		}
		//vn = v;
		for(int j =0; j < ny; j++){
			for(int i =0; i < nx; i++){
				vn[j][i] = v[j][i];
			}
		}
		for(int j =1; j < ny-1; j++){
            		for(int i =1; i < nx-1; i++){ 
				u[j][i] = un[j][i] - un[j][i] * dt / dx * (un[j][i] - un[j][i - 1])\
                               - un[j][i] * dt / dy * (un[j][i] - un[j - 1][i])\
                               - dt / (2 * rho * dx) * (p[j][i+1] - p[j][i-1])\
                               + nu * dt / std::pow(dx,2.0) * (un[j][i+1] - 2 * un[j][i] + un[j][i-1])\
                               + nu * dt / std::pow(dy,2.0) * (un[j+1][i] - 2 * un[j][i] + un[j-1][i]);
            			v[j][i] = vn[j][i] - vn[j][i] * dt / dx * (vn[j][i] - vn[j][i - 1])\
                               - vn[j][i] * dt / dy * (vn[j][i] - vn[j - 1][i])\
                               - dt / (2 * rho * dx) * (p[j+1][i] - p[j-1][i])\
                               + nu * dt / std::pow(dx,2.0) * (vn[j][i+1] - 2 * vn[j][i] + vn[j][i-1])\
                               + nu * dt / std::pow(dy,2.0) * (vn[j+1][i] - 2 * vn[j][i] + vn[j-1][i]);				
			}
		}
		//u[0, :]  = 0 & u[-1, :] = 1 & v[0, :]  = 0 & v[-1, :] = 0
		for(int i =0; i < nx; i++){
            		u[0][i] = 0;
			u[nx-1][i] = 1;
            		v[0][i] = 0;
			v[nx-1][i] = 0;
        	}
		//u[:, 0]  = 0 & u[:, -1] = 0 & v[:, 0]  = 0 & v[:, -1] = 0
		for(int j =0; j < ny; j++){
			u[j][0] = 0;
			u[j][ny-1] = 0;
			v[j][0] = 0;
			v[j][ny-1] = 0;
		}
	}
	gettimeofday(&et, NULL); /* get end time */
	us = time_diff_us(st, et);
        printf("Matmul took %ld us\n",us);
}
