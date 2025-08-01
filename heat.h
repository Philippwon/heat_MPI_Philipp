/*
 * heat.h
 *
 * Global definitions for the iterative solver
 */

#ifndef JACOBI_H_INCLUDED
#define JACOBI_H_INCLUDED

#include <stdio.h>

// configuration

typedef struct
{
    float posx;
    float posy;
    float range;
    float temp;
}
heatsrc_t;

typedef struct
{
    unsigned maxiter;       // maximum number of iterations
    unsigned act_res;
    unsigned max_res;       // spatial resolution
    unsigned initial_res;
    unsigned res_step_size;
    int algorithm;          // 0=>Jacobi, 1=>Gauss

    unsigned visres;        // visualization resolution
  
    double *u, *uhelp;
    double *uvis;

    unsigned   numsrcs;     // number of heat sources
    heatsrc_t *heatsrcs;
}
algoparam_t;


// function declarations

// misc.c
int initialize( algoparam_t *param );
// int initialize_MPI( algoparam_t *param, int width, int height, double x_stepsize, double y_stepsize, int x_coord, int y_coord, int count_node_width, int count_node_height, int node_pixel_width, int node_pixel_height);
int initialize_MPI( algoparam_t *param, int width, int height, double x_stepsize, double y_stepsize, int x_coord, int y_coord, int count_node_width, int count_node_height, int node_pixel_width, int node_pixel_height, int visres_x, int visres_y);
// void write_image_MPI( FILE * f, double *u, unsigned sizex, unsigned sizey, int size_x_global, int x_pos_start, int y_pos_start);
void write_image_MPI(FILE * f, double *u, unsigned sizex, unsigned sizey, int size_x_global, int size_y_global, int x_pos_start, int y_pos_start);
int finalize( algoparam_t *param );
void write_image( FILE * f, double *u,
		  unsigned sizex, unsigned sizey );
int coarsen(double *uold, unsigned oldx, unsigned oldy ,
	    double *unew, unsigned newx, unsigned newy );

// Gauss-Seidel: relax_gauss.c
double residual_gauss( double *u,
		       unsigned sizex, unsigned sizey );
void relax_gauss( double *u, 
		  unsigned sizex, unsigned sizey  );

// Jacobi: relax_jacobi.c
double residual_jacobi( double *u, double *utmp,
			unsigned sizex, unsigned sizey );
void relax_jacobi( double *u, double *utmp,
		   unsigned sizex, unsigned sizey ); 


#endif // JACOBI_H_INCLUDED
