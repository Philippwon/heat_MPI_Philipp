/*
 * relax_gauss.c
 *
 * Gauss-Seidel Relaxation
 *
 */

#include "heat.h"
#include "math.h"

/*
 * Residual (length of error vector)
 * between current solution and next after a Gauss-Seidel step
 *
 * Temporary array utmp needed to not change current solution
 *
 * Flop count in inner body is 7
 */

double residual_gauss(double *u, unsigned sizex, unsigned sizey) {
	unsigned i, j;
	double ucurrent, diff, sum = 0.0;

	const unsigned xBlock = 4;
	const unsigned yBlock = 4;

	for (i = 1; i < sizey - 1; i+=yBlock) {
		for (j = 1; j < sizex - 1; j+=xBlock) {
			const unsigned xRemBlock = xBlock <= ((sizex - 1) - j) ? xBlock : ((sizex - 1) - j);
			const unsigned yRemBlock = yBlock <= ((sizey - 1) - i) ? yBlock : ((sizey - 1) - i);
			for (unsigned ii = i; ii < i + yRemBlock; ii++) {
				for (unsigned jj = j; jj < j + xRemBlock; jj++) {
					ucurrent = u[ii * sizex + jj];
					u[ii * sizex + jj] = 0.25 * (u[ii * sizex + (jj - 1)] +  // new left
											u[ii * sizex + (jj + 1)] +  // right
											u[(ii - 1) * sizex + jj] +  // new top
											u[(ii + 1) * sizex + jj]); // bottom

					diff = u[ii * sizex + jj] - ucurrent;
					sum += diff * diff;
				}
			}
		}
	}

	return sum;
}

/*
 * One Gauss-Seidel iteration step
 *
 * Flop count in inner body is 4
 */
void relax_gauss(double *u, unsigned sizex, unsigned sizey) {
	unsigned i, j;

	for (i = 1; i < sizey - 1; i++) {
		for (j = 1; j < sizex - 1; j++) {
			u[i * sizex + j] = 0.25 * (u[i * sizex + (j - 1)] + u[i * sizex + (j + 1)] + u[(i - 1) * sizex + j] + u[(i + 1) * sizex + j]);
		}
	}
}
