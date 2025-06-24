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
#include <mpi.h>

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
int initialize_MPI( algoparam_t *param, int width, int height, double x_stepsize, double y_stepsize, int x_coord, int y_coord, int count_node_width,
	 int count_node_height, int node_pixel_width, int node_pixel_height, int visres_x, int visres_y)
{
    int i, j;
    double dist;

	int width_in = width+2;
	int height_in = height+2;



    (param->u)     = (double*) aligned_alloc(64,sizeof(double) * width_in*height_in);//(double*)calloc( sizeof(double),np*np );
    (param->uhelp) = (double*) aligned_alloc(64,sizeof(double) * width_in*height_in);
    // (param->uvis)  = (double*) aligned_alloc(64,sizeof(double) * (param->visres+2)*(param->visres+2));
	(param->uvis)  = (double*) aligned_alloc(64,sizeof(double) * (visres_x)*(visres_y));

	memset(param->u,0,sizeof(double) * width_in*height_in);
	memset(param->uhelp,0,sizeof(double) * width_in*height_in);

	// memset(param->uvis,0,sizeof(double) * (param->visres+2)*(param->visres+2));
	memset(param->uvis,0,sizeof(double) * (visres_x)*(visres_y));

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








/*
 * write the given temperature u matrix to rgb values
 * and write the resulting image to file f
 */
void write_image_MPI_old( FILE * f, double *u, unsigned sizex, unsigned sizey, int size_x_global, int x_pos_start, int y_pos_start) 
{
    // RGB table
    unsigned char r[1024], g[1024], b[1024];
    int i, j, k;
  
    double min, max;
	double global_min, global_max;

    j=1023;

    // prepare RGB table
	{
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

	MPI_Allreduce(&max,&global_max,1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
	MPI_Allreduce(&min,&global_min,1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);

	MPI_Status status;
	MPI_File fhw;
	int ierr;
	ierr = MPI_File_open(MPI_COMM_WORLD,"datafile",MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &fhw);
	MPI_File_set_view(fhw, 0, MPI_UNSIGNED_CHAR, MPI_UNSIGNED_CHAR, "native", MPI_INFO_NULL);

	int offset;

	fprintf(stderr, " wir kommen hierhin und min ist: %f\n",global_min);

    // fprintf(f, "P3\n");
    // fprintf(f, "%u %u\n", sizex, sizey);
    // fprintf(f, "%u\n", 255);

	// if (ierr != MPI_SUCCESS) {
	// 	char err_string[MPI_MAX_ERROR_STRING];
	// 	int err_length;
	// 	MPI_Error_string(ierr, err_string, &err_length);
	// 	fprintf(stderr, "[%d] Fehler bei MPI_File_open: %s\n", size_x_global * y_pos_start +  x_pos_start , err_string);
	// 	MPI_Abort(MPI_COMM_WORLD, ierr);  // optional: beende das Programm sauber
	// }

	// if (rank == 0) {
    //     char header[64];
    //     int header_len = snprintf(header, sizeof(header), "P3\n%d %d\n255\n", size_x_global, size_y_global);
    //     MPI_File_write_at(fh, 0, header, header_len, MPI_CHAR, &status);
    // }

    for( i=0; i<sizey; i++ )
    {
		for( j=0; j<sizex; j++ )
		{	
			offset = y_pos_start*size_x_global+ x_pos_start*3 + i*size_x_global + j*3;

			offset = (y_pos_start + i) * size_x_global * 3 + (x_pos_start + j) * 3;


			k=(int)(1024.0*(u[i*sizex+j]-global_min)/(global_max-global_min));
			if (k==1024) k=1023;
			// fprintf(f, "%d %d %d  ", r[k], g[k], b[k]);


			// if(offset < 1000)
			// 	fprintf(stderr, " wir kommen hierhin und printen \"%c\" der offset ist %d \n",dome1, offset);

			// MPI_File_write_at(fhw, offset, &r[k], 1, MPI_UNSIGNED_CHAR, &status);
			// MPI_File_write_at(fhw, offset+1, &g[k], 1, MPI_UNSIGNED_CHAR, &status);
			// MPI_File_write_at(fhw, offset+2, &b[k], 1, MPI_UNSIGNED_CHAR, &status);

			char pixelbuf[16];
            int pixel_len = snprintf(pixelbuf, sizeof(pixelbuf), "%3d %3d %3d ", r[k], g[k], b[k]);

            MPI_File_write_at(fhw, offset, pixelbuf, pixel_len, MPI_CHAR, &status);

		}
		// fprintf(f, "\n");
    }
	fprintf(stderr, "[%d] finished writing\n",  x_pos_start );

	MPI_File_close(&fhw);

}









void write_image_MPI(FILE * f, double *u, unsigned sizex, unsigned sizey, int size_x_global, int size_y_global, int x_pos_start, int y_pos_start)
{
    unsigned char r[1024], g[1024], b[1024];
    int i, j, k;
    double min = DBL_MAX, max = -DBL_MAX;
    double global_min, global_max;

    // Create RGB gradient
    int jx = 1023;
    for (i = 0; i < 256; i++) { r[jx]=255; g[jx]=i; b[jx]=0; jx--; }
    for (i = 0; i < 256; i++) { r[jx]=255-i; g[jx]=255; b[jx]=0; jx--; }
    for (i = 0; i < 256; i++) { r[jx]=0; g[jx]=255; b[jx]=i; jx--; }
    for (i = 0; i < 256; i++) { r[jx]=0; g[jx]=255-i; b[jx]=255; jx--; }

    // Find local min/max
    for (i = 0; i < sizey; i++) {
        for (j = 0; j < sizex; j++) {
            double val = u[i*sizex + j];
            if (val > max) max = val;
            if (val < min) min = val;
        }
    }

    // Global min/max
    MPI_Allreduce(&min, &global_min, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(&max, &global_max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

    MPI_File fh;
    MPI_Status status;
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    MPI_File_open(MPI_COMM_WORLD, "data.ppm", MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);

    // Header schreiben (nur Rank 0)
    if (rank == 0) {
        char header[64];
        int header_len = snprintf(header, sizeof(header), "P3\n%d %d\n255\n", size_x_global, size_y_global);
        MPI_File_write_at(fh, 0, header, header_len, MPI_CHAR, &status);
    }

    // Synchronisieren, damit Header fertig ist
    MPI_Barrier(MPI_COMM_WORLD);

    // Offset berechnen: wie viele Zeichen hat der Header?
    MPI_Offset header_offset = 0;
    if (rank == 0) {
        header_offset = snprintf(NULL, 0, "P3\n%d %d\n255\n", size_x_global, size_y_global);
    }
    MPI_Bcast(&header_offset, 1, MPI_OFFSET, 0, MPI_COMM_WORLD);

    // Jeder Prozess schreibt seine Zeilen
    for (i = 0; i < sizey; i++) {
        int global_y = y_pos_start + i;
        MPI_Offset line_offset = header_offset;

        // Jeder vorherige Pixel erzeugt ca. "xxx yyy zzz " → max 12 Byte
        line_offset += (global_y * size_x_global) * 12;


        // Berechne Position innerhalb dieser Zeile
        line_offset += x_pos_start * 12;

        // Zeile schreiben
        for (j = 0; j < sizex; j++) {
            int global_x = x_pos_start + j;
            int index = i * sizex + j;
            double val = u[index];

            k = (int)(1024.0 * (val - global_min) / (global_max - global_min));
            if (k > 1023) k = 1023;

            char pixelbuf[16];
            int pixel_len = snprintf(pixelbuf, sizeof(pixelbuf), "%3d %3d %3d ", r[k], g[k], b[k]);

            MPI_File_write_at(fh, line_offset, pixelbuf, pixel_len, MPI_CHAR, &status);
            line_offset += pixel_len;
        }

        // Nur Rank 0 jedes Prozesses schreibt am Ende der Zeile ein \n
        if (x_pos_start + sizex == size_x_global) {
            MPI_File_write_at(fh, line_offset, "\n", 1, MPI_CHAR, &status);
        }
    }

    MPI_File_close(&fh);
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





// oldx is the number of pixels for the entire width of the old screen 
// newx is the number of pixels for the entire width of the resolution screen (including 2 for padding)
// int coarsen_MPI( double *uold, unsigned oldx, unsigned oldy, double *unew, unsigned newx, unsigned newy, int x_coord, int y_coord, int count_proc_row, int count_proc_col )
// {

// 	double x_stepsize = (oldx-1) / (newx-1);
// 	double y_stepsize = (oldy-1) / (newy-1);

	


// 	int upper_left_corner_x = x_coord > 0 ? ((int)((oldx-2)/count_proc_row)) * x_coord + 1 : 0;
// 	int upper_left_corner_y = y_coord > 0 ? ((int)((oldy-2)/count_proc_col)) * y_coord + 1 : 0;




// 	double x_pos = round(upper_left_corner_x/x_stepsize + 0.5) * x_stepsize;
// 	double y_pos = round(upper_left_corner_y/y_stepsize + 0.5) * y_stepsize;



// 	while
	



	




// }


