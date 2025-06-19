/*
 * misc.c
 *
 * Helper functions for
 * - initialization
 * - finalization,
 * - writing out a picture
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include <stdlib.h>
#include <string.h>

#include "heat.h"

/*
 * Initialize the iterative solver
 * - allocate memory for matrices
 * - set boundary conditions according to configuration
 */
int initialize( algoparam_t *param )
{
    int i, j;
    double dist;

    // total number of points (including border)
    const int np = param->act_res + 2;
  
    //
    // allocate memory
    //
    (param->u)     = (double*) aligned_alloc(64,sizeof(double) * np*np);//(double*)calloc( sizeof(double),np*np );
    (param->uhelp) = (double*) aligned_alloc(64,sizeof(double) * np*np);
    (param->uvis)  = (double*) aligned_alloc(64,sizeof(double) * (param->visres+2)*(param->visres+2));

	memset(param->u,0,sizeof(double) * np*np);
	memset(param->uhelp,0,sizeof(double) * np*np);
	memset(param->uvis,0,sizeof(double) * (param->visres+2)*(param->visres+2));

    if( !(param->u) || !(param->uhelp) || !(param->uvis) )
    {
	fprintf(stderr, "Error: Cannot allocate memory\n");
	return 0;
    }

    for( i=0; i<param->numsrcs; i++ )
    {
	/* top row */
	for( j=0; j<np; j++ )
	{
	    dist = sqrt( pow((double)j/(double)(np-1) - 
			     param->heatsrcs[i].posx, 2)+
			 pow(param->heatsrcs[i].posy, 2));
	  
	    if( dist <= param->heatsrcs[i].range )
	    {
		(param->u)[j] +=
		    (param->heatsrcs[i].range-dist) /
		    param->heatsrcs[i].range *
		    param->heatsrcs[i].temp;
	    }
	}
      
	/* bottom row */
	for( j=0; j<np; j++ )
	{
	    dist = sqrt( pow((double)j/(double)(np-1) - 
			     param->heatsrcs[i].posx, 2)+
			 pow(1-param->heatsrcs[i].posy, 2));
	  
	    if( dist <= param->heatsrcs[i].range )
	    {
		(param->u)[(np-1)*np+j]+=
		    (param->heatsrcs[i].range-dist) / 
		    param->heatsrcs[i].range * 
		    param->heatsrcs[i].temp;
	    }
	}
      
	/* leftmost column */
	for( j=1; j<np-1; j++ )
	{
	    dist = sqrt( pow(param->heatsrcs[i].posx, 2)+
			 pow((double)j/(double)(np-1) - 
			     param->heatsrcs[i].posy, 2)); 
	  
	    if( dist <= param->heatsrcs[i].range )
	    {
		(param->u)[ j*np ]+=
		    (param->heatsrcs[i].range-dist) / 
		    param->heatsrcs[i].range *
		    param->heatsrcs[i].temp;
	    }
	}
      
	/* rightmost column */
	for( j=1; j<np-1; j++ )
	{
	    dist = sqrt( pow(1-param->heatsrcs[i].posx, 2)+
			 pow((double)j/(double)(np-1) - 
			     param->heatsrcs[i].posy, 2)); 
	  
	    if( dist <= param->heatsrcs[i].range )
	    {
		(param->u)[ j*np+(np-1) ]+=
		    (param->heatsrcs[i].range-dist) /
		    param->heatsrcs[i].range *
		    param->heatsrcs[i].temp;
	    }
	}
    }

    return 1;
}




// normal width and height refer to the number of pixels the current process is responsible for
// node height and node width refer to the number of nodes used to distribute the work e.g 3x5 nodes refers to 3 lines of 5 nodes
// x and y_coord refer to the position of the current node inside this grid starting from (0,0) in the top left corner up to (h-1,w-1) in bottom right corner
// node_pixel_width/width refer to the width/height of all but the last node in a line/column. it is used to calculate the exakt ccordinates e.g: (0.784,0.342) of the upper left corner of a nodes grid
int initialize_MPI( algoparam_t *param, int width, int height, double x_stepsize, double y_stepsize, int x_coord, int y_coord, int count_node_width, int count_node_height, int node_pixel_width, int node_pixel_height)
{
    int i, j;
    double dist;

	int width_in = width+2;
	int height_in = height+2;

  

    (param->u)     = (double*) aligned_alloc(64,sizeof(double) * width_in*height_in);//(double*)calloc( sizeof(double),np*np );
    (param->uhelp) = (double*) aligned_alloc(64,sizeof(double) * width_in*height_in);
    (param->uvis)  = (double*) aligned_alloc(64,sizeof(double) * (param->visres+2)*(param->visres+2));

	memset(param->u,0,sizeof(double) * width_in*height_in);
	memset(param->uhelp,0,sizeof(double) * width_in*height_in);
	memset(param->uvis,0,sizeof(double) * (param->visres+2)*(param->visres+2));

    if( !(param->u) || !(param->uhelp) || !(param->uvis) )
	{
		fprintf(stderr, "Error: Cannot allocate memory\n");
		return 0;
    }



	double x_start = x_coord * node_pixel_width * x_stepsize;
	double y_start = y_coord * node_pixel_height * y_stepsize;

	// check if this nodes top edge is part of the top edge of the entire screen
	if(y_coord == 0){
		for( i=0; i<param->numsrcs; i++ )
		{
			for( j=0; j<width_in; j++)
			{
				dist = sqrt( pow(((double)j*x_stepsize + x_start) - param->heatsrcs[i].posx, 2) + pow(param->heatsrcs[i].posy, 2));
				if( dist <= param->heatsrcs[i].range )
				{
					(param->u)[j] += (param->heatsrcs[i].range-dist) / param->heatsrcs[i].range * param->heatsrcs[i].temp;
				}
			}

		}

	}

	// check if this nodes bottom edge is part of the bottom edge of the entire screen
	if(y_coord == count_node_height - 1){
		for( i=0; i<param->numsrcs; i++ )
		{
			for( j=0; j<width_in; j++)
			{
				dist = sqrt( pow(((double)j * x_stepsize + x_start) - param->heatsrcs[i].posx, 2) + pow(1 - param->heatsrcs[i].posy, 2));
				if( dist <= param->heatsrcs[i].range )
				{
					(param->u)[width_in * (height_in-1) + j] += (param->heatsrcs[i].range-dist) / param->heatsrcs[i].range * param->heatsrcs[i].temp;
				}
			}
		}
	}


	// check if this nodes left edge is part of the left edge of the entire screen
	if(x_coord == 0){
		for( i=0; i<param->numsrcs; i++ )
		{
			for( j=0; j<height_in; j++)
			{
				dist = sqrt( pow(param->heatsrcs[i].posx, 2) + pow(((double)j * y_stepsize + y_start) - param->heatsrcs[i].posy, 2));
				if( dist <= param->heatsrcs[i].range )
				{
					(param->u)[j * width_in ] += (param->heatsrcs[i].range-dist) / param->heatsrcs[i].range * param->heatsrcs[i].temp;
				}
			}
		}
	}

	// check if this nodes right edge is part of the right edge of the entire screen
	if(x_coord == count_node_width - 1){
		for( i=0; i<param->numsrcs; i++ )
		{
			for( j=0; j<height_in; j++)
			{
				dist = sqrt( pow(1 - param->heatsrcs[i].posx, 2) + pow(((double)j * y_stepsize + y_start) - param->heatsrcs[i].posy, 2));
				if( dist <= param->heatsrcs[i].range )
				{
					(param->u)[j * width_in + width_in-1] += (param->heatsrcs[i].range-dist) / param->heatsrcs[i].range * param->heatsrcs[i].temp;
				}
			}
		
		
		}

	}

	char txt_files[100];

	snprintf(txt_files, sizeof(txt_files), "data/output%d.txt", (count_node_width*y_coord+x_coord));

	FILE *pfile = fopen(txt_files, "w");

	if(pfile == NULL){
		fprintf(stderr, "couldnt open write File");
	}

	for( int j = 0; j < height_in; j++){
		for( int i = 0; i < width_in; i++)
			{
				fprintf(pfile, "%f ; ", (param->u)[j * (width_in) + i]);
			}
			fprintf(pfile, "\n");
	}
	fprintf(pfile, "\n\n\n");

	fclose(pfile);
    return 1;
}







































/*
 * free used memory
 */
int finalize( algoparam_t *param )
{
    if( param->u ) {
	free(param->u);
	param->u = 0;
    }

    if( param->uhelp ) {
	free(param->uhelp);
	param->uhelp = 0;
    }

    if( param->uvis ) {
	free(param->uvis);
	param->uvis = 0;
    }

    return 1;
}


/*
 * write the given temperature u matrix to rgb values
 * and write the resulting image to file f
 */
void write_image( FILE * f, double *u,
		  unsigned sizex, unsigned sizey ) 
{
    // RGB table
    unsigned char r[1024], g[1024], b[1024];
    int i, j, k;
  
    double min, max;

    j=1023;

    // prepare RGB table
    for( i=0; i<256; i++ )
    {
	r[j]=255; g[j]=i; b[j]=0;
	j--;
    }
    for( i=0; i<256; i++ )
    {
	r[j]=255-i; g[j]=255; b[j]=0;
	j--;
    }
    for( i=0; i<256; i++ )
    {
	r[j]=0; g[j]=255; b[j]=i;
	j--;
    }
    for( i=0; i<256; i++ )
    {
	r[j]=0; g[j]=255-i; b[j]=255;
	j--;
    }


    min=DBL_MAX;
    max=-DBL_MAX;

    // find minimum and maximum 
    for( i=0; i<sizey; i++ )
    {
	for( j=0; j<sizex; j++ )
	{
	    if( u[i*sizex+j]>max )
		max=u[i*sizex+j];
	    if( u[i*sizex+j]<min )
		min=u[i*sizex+j];
	}
    }
  

    fprintf(f, "P3\n");
    fprintf(f, "%u %u\n", sizex, sizey);
    fprintf(f, "%u\n", 255);

    for( i=0; i<sizey; i++ )
    {
	for( j=0; j<sizex; j++ )
	{
	    k=(int)(1024.0*(u[i*sizex+j]-min)/(max-min));
		if (k==1024) k=1023;

	    fprintf(f, "%d %d %d  ", r[k], g[k], b[k]);
	}
	fprintf(f, "\n");
    }
}


int coarsen( double *uold, unsigned oldx, unsigned oldy ,
	     double *unew, unsigned newx, unsigned newy )
{
    int i, j;

    int stepx;
    int stepy;
    int stopx = newx;
    int stopy = newy;

    if (oldx>newx)
		stepx=oldx/newx;
    else {
		stepx=1;
		stopx=oldx;
    }

    if (oldy>newy)
		stepy=oldy/newy;
    else {
		stepy=1;
		stopy=oldy;
    }

    // NOTE: this only takes the top-left corner,
    // and doesnt' do any real coarsening 
    for( i=0; i<stopy-1; i++ )
    {
	for( j=0; j<stopx-1; j++ )
        {
	    unew[i*newx+j]=uold[i*oldx*stepy+j*stepx];
        }
    }

  return 1;
}


