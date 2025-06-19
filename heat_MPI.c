/*
 * heat.h
 *
 * Iterative solver for heat distribution
 */

#include "heat.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>

#include "input.h"
#include "timing.h"

// void usage(char *s) {
// 	fprintf(stderr, "Usage: %s <input file> [result file]\n\n", s);
// }


void usage(char *s) {
	fprintf(stderr, "Usage: %s <input file> [result file] <width> <height>\n\n", s);
}



int main(int argc, char *argv[]) {

	int size, rank;
	MPI_Init (&argc, &argv);
	MPI_Comm_size (MPI_COMM_WORLD, &size);
	MPI_Comm_rank (MPI_COMM_WORLD, &rank);
	MPI_Comm row_comm, col_comm;

	unsigned iter;
	FILE *infile, *resfile;
	char *resfilename;

	// algorithmic parameters
	algoparam_t param;
	int np,i;

	double runtime, flop;
	double residual;
	double time[1000];
	double floprate[1000];
	int resolution[1000];
	int experiment=0;








	// check arguments
	if (argc < 5) {
		usage(argv[0]);
		return 1;
	}
	// check input file
	if (!(infile = fopen(argv[1], "r"))) {
		fprintf(stderr, "\nError: Cannot open \"%s\" for reading.\n\n", argv[1]);
		usage(argv[0]);
		return 1;
	}
	// check result file
	resfilename = (argc >= 3) ? argv[2] : "heat.ppm";
	if (!(resfile = fopen(resfilename, "w"))) {
		fprintf(stderr, "\nError: Cannot open \"%s\" for writing.\n\n", resfilename);
		usage(argv[0]);
		return 1;
	}
	// check input
	if (!read_input(infile, &param)) {
		fprintf(stderr, "\nError: Error parsing input file.\n\n");
		usage(argv[0]);
		return 1;
	}
	fprintf(stderr, "\n################################################################################\n");
	print_params(&param);
	fprintf(stderr, "\nwidth ist \"%s\" height ist \"%s\".\n\n", argv[3], argv[4]);










	// set the visualization resolution
	param.visres = 1024;

	param.u = 0;
	param.uhelp = 0;
	param.uvis = 0;

	param.act_res = param.initial_res;


	int width = atoi(argv[3]);
	int height = atoi(argv[4]);

	fprintf(stderr, "\nwidth ist \"%d\" height ist \"%d\" mein rank ist: \"%d\" daher ist meine posx: \"%d\" .\n\n", width, height, rank, rank % width);


	int pos_x = rank % width;
	int pos_y = rank / width;

	fprintf(stderr, "\npos_x ist \"%d\" pos_y ist \"%d\".\n", pos_x, pos_y);


	int height_pixels = pos_y == height - 1 ? param.act_res - pos_y * (param.act_res / height) : (param.act_res / height);
	int width_pixels = pos_x == width - 1 ? param.act_res - pos_x * (param.act_res / width) : (param.act_res / width);

	fprintf(stderr, "process: %d  	my coordinates (x,y) are: (%d,%d)\nquadrant width = %d and height = %d\n\n",rank, pos_x, pos_y, width_pixels, height_pixels);


	double x_stepsize = 1.0 / param.act_res;
	double y_stepsize = 1.0 / param.act_res;

	int normal_height_pixels = param.act_res / height;
	int normal_width_pixels = param.act_res / width;


	// loop over different resolutions
	while (1) {

		// free allocated memory of previous experiment
		// if (param.u != 0)
		// 	finalize(&param);

		// if (!initialize(&param)) {
		// 	fprintf(stderr, "Error in Jacobi initialization.\n\n");

		// 	usage(argv[0]);
		// }
		// fprintf(stderr, "Resolution: %5u\r", param.act_res);


		// node height and node width refer to the number of nodes used to distribute the work e.g 3x5 nodes refers to 3 lines of 5 nodes
		// x and y_coord refer to the position of the current node inside this grid starting from (0,0) in the top left corner up to (h-1,w-1) in bottom right corner
		// node_pixel_width/width refer to the width/height of all but the last node in a line/column. it is used to calculate the exakt ccordinates e.g: (0.784,0.342) of the upper left corner of a nodes grid
		// int initialize_MPI( algoparam_t *param, int width, int height, double x_stepsize, double y_stepsize, int x_coord, int y_coord, int count_node_width, int count_node_height, int node_pixel_width, int node_pixel_height)

		if (param.u != 0)
			finalize(&param);

		if (!initialize_MPI(&param, width_pixels, height_pixels, x_stepsize, y_stepsize, pos_x, pos_y, width, height, normal_width_pixels, normal_height_pixels)) {
			fprintf(stderr, "Error in Jacobi initialization for process: %d.\n\n", rank);
			usage(argv[0]);
		}

		fprintf(stderr, "process: %d height, width : (%d,%d)   coordinates (y,x) (%d,%d) \n", rank, height_pixels, width_pixels, pos_y, pos_x);





		int pos_x = rank % width;
		int pos_y = rank / width;

		//comm, color (has to be the same), key (determines the rank), newcomm
		MPI_Comm_split(MPI_COMM_WORLD, pos_y, pos_x, &row_comm);
		MPI_Comm_split(MPI_COMM_WORLD, pos_x, pos_y, &col_comm);

		fprintf(stderr, "HUUHUHUHUH");

		double * left_halo_zone = (double*) aligned_alloc(64,sizeof(double) * height_pixels);
		double * right_halo_zone = (double*) aligned_alloc(64,sizeof(double) * height_pixels);

		double * send_left_halo = (double*) aligned_alloc(64,sizeof(double) * height_pixels);
		double * send_right_halo = (double*) aligned_alloc(64,sizeof(double) * height_pixels);

		double * upper_halo_zone = (double*) aligned_alloc(64,sizeof(double) * width_pixels);
		double * lower_halo_zone = (double*) aligned_alloc(64,sizeof(double) * width_pixels);




		
		// full size (param.act_res are only the inner points)
		np = param.act_res + 2;

		// starting time
		runtime = wtime();
		residual = 999999999;

		iter = 0;

		if(param.algorithm == 0){
			memcpy(param.uhelp, param.u, sizeof(double) * (width_pixels + 2) * (height_pixels + 2));// memcpy(param.uhelp, param.u, sizeof(double) * np * np);
		}

		fprintf(stderr, "HUUHUHUHUH 2");

		while (1) {
			fprintf(stderr, "Process: %d   iteration:%d\n",rank, iter);

			switch (param.algorithm) {

			case 0: // JACOBI
				residual = iter % 2 == 0 ? residual_jacobi(param.u, param.uhelp, width_pixels + 2, height_pixels + 2) : residual_jacobi(param.uhelp, param.u, width_pixels + 2, height_pixels + 2); // residual = iter % 2 == 0 ? residual_jacobi(param.u, param.uhelp, np, np) : residual_jacobi(param.uhelp, param.u, np, np);
				break;

			case 1: // GAUSS
				residual = residual_gauss(param.u, width_pixels + 2, height_pixels + 2); // residual = residual_gauss(param.u, np, np);
				break;
			}

			iter++;

			// solution good enough ?
			// if (residual < 0.000005)
			// 	break;

			// max. iteration reached ? (no limit with maxiter=0)
			if (param.maxiter > 0 && iter >= param.maxiter)
				break;



			// HERE THE COMMUNICATION PHASE BETWEEN THE DIFFERENT PORCESSES NEEDS TO HAPPEN
			// WE HAVE TO ALSO LOOK AT THE ITERATION SINCE OUR CORRECT VALUES MIGHT BE IN DIFFERENT SPOTS
			// WE SHOULD ALSO LOOK OUT FOR EDGE CASES LIKE ONLY 1 ROW OR 1 COLUMN
			double * current_data = iter % 2 == 1 ? param.uhelp : param.u;

			int upper_partner = pos_y > 0 ? pos_y - 1 : -1;
			int lower_partner = pos_y < height - 1 ? pos_y + 1 : -1;

			int left_partner = pos_x > 0 ? pos_x - 1 : -1;
			int right_partner = pos_x < width - 1 ? pos_x + 1 : -1;

			fprintf(stderr, "before communication    Process: %d   iteration:%d\n",rank, iter);

			// VERTICAL EXCHANGE (UP,DOWN)
			if(pos_x % 2 == 0){ // 0,2,4,6
				if(upper_partner != -1) // above
					MPI_Sendrecv(current_data + (width_pixels+2) + 1, width_pixels, MPI_DOUBLE, upper_partner, 0, 
								 upper_halo_zone, width_pixels, MPI_DOUBLE, upper_partner, 0, col_comm, MPI_STATUS_IGNORE);
				if(lower_partner != -1)
					MPI_Sendrecv(current_data + (width_pixels+2) * (height_pixels) + 1, width_pixels, MPI_DOUBLE, lower_partner, 0, 
								 lower_halo_zone, width_pixels, MPI_DOUBLE, lower_partner, 0, col_comm, MPI_STATUS_IGNORE);			 
			}
			else { // 1,3,5,7
				if(lower_partner != -1)// below
					MPI_Sendrecv(current_data + (width_pixels+2) * (height_pixels) + 1, width_pixels, MPI_DOUBLE, lower_partner, 0, 
								 lower_halo_zone, width_pixels, MPI_DOUBLE, lower_partner, 0, col_comm, MPI_STATUS_IGNORE);
				if(upper_partner != -1) // unnecessary because 1,3,5 and so on will always have an upper partner but is more understandable
					MPI_Sendrecv(current_data + (width_pixels+2) + 1, width_pixels, MPI_DOUBLE, upper_partner, 0, 
								 upper_halo_zone, width_pixels, MPI_DOUBLE, upper_partner, 0, col_comm, MPI_STATUS_IGNORE);
			}


			fprintf(stderr, "after vertical communication    Process: %d   iteration:%d\n",rank, iter);
			for(int i = 1; i < height_pixels+1; i++){
				*(send_left_halo+i) = *(current_data + (width_pixels+2)*i +1);
				*(send_right_halo+i) = *(current_data + (width_pixels+2)*i +(width_pixels));
			}


			// HORIZONTAL EXCHANGE (LEFT,RIGHT)
			if(pos_y % 2 == 0){ // 0,2,4,6
				if(left_partner != -1)
					MPI_Sendrecv(send_left_halo, height_pixels, MPI_DOUBLE, left_partner, 0,
								 left_halo_zone, height_pixels, MPI_DOUBLE, left_partner, 0, row_comm, MPI_STATUS_IGNORE);
				if(right_partner != -1)
					MPI_Sendrecv(send_right_halo, height_pixels, MPI_DOUBLE, right_partner, 0,
								 right_halo_zone, height_pixels, MPI_DOUBLE, right_partner, 0, row_comm, MPI_STATUS_IGNORE);
			}
			else {
				if(right_partner != -1)
					MPI_Sendrecv(send_right_halo, height_pixels, MPI_DOUBLE, right_partner, 0,
								 right_halo_zone, height_pixels, MPI_DOUBLE, right_partner, 0, row_comm, MPI_STATUS_IGNORE);
				if(left_partner != -1)
					MPI_Sendrecv(send_left_halo, height_pixels, MPI_DOUBLE, left_partner, 0,
								 left_halo_zone, height_pixels, MPI_DOUBLE, left_partner, 0, row_comm, MPI_STATUS_IGNORE);
			}



			if(pos_y > 0){ // there exists an upper partner
				memcpy(current_data + 1 , upper_halo_zone, sizeof(double) * width_pixels);
			}
			if(pos_y < height -1){ // there exists an lower partner
				memcpy(current_data + (height_pixels+1) * (width_pixels+2) + 1, lower_halo_zone,  sizeof(double) * width_pixels);
			}

			if(pos_x > 0){ // there exists a left partner
				for(int i = 1; i < height_pixels+1;i++){
					*(current_data + (width_pixels+2)*i) = *(left_halo_zone+i);
				}
			}

			if(pos_x < width -1){ // there exists a right partner
				for(int i = 1; i < height_pixels+1;i++){
					*(current_data + (width_pixels+2)*i + (width_pixels)) = *(right_halo_zone+i);
				}
			}


			fprintf(stderr, "end of iteration    Process: %d   iteration:%d\n",rank, iter);
			MPI_Barrier(MPI_COMM_WORLD);
		}



		if(param.algorithm == 0 && iter % 2 == 1 ){
			// memcpy(param.u, param.uhelp, sizeof(double) * np * np);
			memcpy(param.u, param.uhelp, sizeof(double) * (width_pixels + 2) * (height_pixels + 2));
		}





		char txt_files[100];
		snprintf(txt_files, sizeof(txt_files), "data/final_output%d.txt", (width*pos_y+pos_x));
		FILE *pfile = fopen(txt_files, "w");
		if(pfile == NULL){
			fprintf(stderr, "couldnt open write File");
		}

		for( int j = 0; j < height_pixels+2; j++){
			for( int i = 0; i < width_pixels+2; i++)
				{
					fprintf(pfile, "%f ; ", (param.u)[j * (width_pixels+2) + i]);
				}
				fprintf(pfile, "\n");
		}
		fprintf(pfile, "\n\n\n");

		fclose(pfile);












		// Flop count after <i> iterations
		flop = iter * 7.0 * param.act_res * param.act_res;
		// stopping time
		runtime = wtime() - runtime;

		fprintf(stderr, "Resolution: %5u, ", param.act_res);
		fprintf(stderr, "Time: %04.3f ", runtime);
		fprintf(stderr, "(%3.3f GFlop => %6.2f MFlop/s, ", flop / 1000000000.0, flop / runtime / 1000000);
		fprintf(stderr, "residual %f, %d iterations)\n", residual, iter);

		// for plot...
		time[experiment]=runtime;
		floprate[experiment]=flop / runtime / 1000000;
		resolution[experiment]=param.act_res;
		experiment++;

		if (param.act_res + param.res_step_size > param.max_res)
			break;
		param.act_res += param.res_step_size;
		
	}








	for (i=0;i<experiment; i++){
		printf("%5d; %5.3f; %5.3f\n", resolution[i], time[i], floprate[i]);

	}

	// coarsen(param.u, np, np, param.uvis, param.visres + 2, param.visres + 2);

	write_image(resfile, param.uvis, param.visres + 2, param.visres + 2);

	finalize(&param);
	MPI_Finalize ();

	return 0;
}
