/*
 * relax_jacobi.c
 *
 * Jacobi Relaxation
 *
 */

#include "heat.h"

/*
 * Residual (length of error vector)
 * between current solution and next after a Jacobi step
 */
double residual_jacobi(double *u, double *utmp, unsigned sizex, unsigned sizey) {
	unsigned i, j;
	double diff, sum = 0.0;

	const unsigned xBlock = 512;
	const unsigned yBlock = 512;

	for (i = 1; i < sizey - 1; i+=yBlock) {
		for (j = 1; j < sizex - 1; j+=xBlock) {
			const unsigned xRemBlock = xBlock <= ((sizex - 1) - j) ? xBlock : ((sizex - 1) - j);
			const unsigned yRemBlock = yBlock <= ((sizey - 1) - i) ? yBlock : ((sizey - 1) - i);
			for (unsigned ii = i; ii < i + yRemBlock; ii++) {
				for (unsigned jj = j; jj < j + xRemBlock; jj++) {
					utmp[ii * sizex + jj] = 0.25 * (u[ii * sizex + (jj - 1)] +  // left
												u[ii * sizex + (jj + 1)] +  // right
												u[(ii - 1) * sizex + jj] +  // top
												u[(ii + 1) * sizex + jj]); // bottom

					diff = utmp[ii * sizex + jj] - u[ii * sizex + jj];
					sum += diff * diff;
				}
			}
		}
	}

	return sum;
}

/*
 * One Jacobi iteration step
 */
void relax_jacobi(double *u, double *utmp, unsigned sizex, unsigned sizey) {
	int i, j;

	for (i = 1; i < sizey - 1; i++) {
		for (j = 1; j < sizex - 1; j++) {
			utmp[i * sizex + j] = 0.25 * (u[i * sizex + (j - 1)] +  // left
						u[i * sizex + (j + 1)] +  // right
						u[(i - 1) * sizex + j] +  // top
						u[(i + 1) * sizex + j]); // bottom
		}
	}
}
