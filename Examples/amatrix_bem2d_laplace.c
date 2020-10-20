#include <stdio.h>
#include <unistd.h>

#include "basic.h"
#include "laplacebem2d.h"


/****************************************************
 * This examples sets up single layer potential operator(SLP),
 * double layer potential operator(DLP) as well as the mass matrix M
 * as dense matrices in order to solve the interior Dirichlet
 * problem for the Laplace equation.
 ****************************************************/

int
main(int argc, char **argv)
{
  pstopwatch sw;
  pcurve2d gr;
  pbem2d bem_slp;
  basisfunctionbem2d basis;
  pamatrix  V;
  uint q_reg, edges;
  real t, size, norm;

  /* Init the H2Lib, should be called before any other function. */
  init_h2lib(&argc, &argv);

  /****************************************************
   * Set up basic parameters
   ****************************************************/

  /* Number of quadrature points */
  q_reg = 4;
  /* Basis functions that should be used. */
  basis = BASIS_CONSTANT_BEM2D;

  /* absolute norm of the residuum for CG-method */
  /* eps_solve = 1.0e-4; */

  /* maximum number of CG-steps that should be performed. */
  /* maxiter = 250; */

  /* Stopwatch for measuring the time. */
  sw = new_stopwatch();

  /****************************************************
   * Create geometry
   ****************************************************/

  /* Create mesh of 2D curve. */
  edges = argc > 1 ? atoi(argv[1]) : 64;
  gr = new_circle_curve2d(edges, 1.0);

  printf("Created geometry with %d vertices and %d edges\n",
	 gr->vertices, gr->edges);

  /* print_curve2d(gr); */

  /****************************************************
   * Set up bem objects
   ****************************************************/

  /* Create a new BEM-object, that can compute entries of SLP operator. */
  bem_slp = new_slp_laplace_bem2d(gr, q_reg, basis);

  /****************************************************
   * Assemble Dense matrix SLP
   ****************************************************/

  /* printf("Assemble dense matrix V:\n"); */

  /* Create amatrix structure. */
  V = new_amatrix(gr->vertices, gr->vertices);

  start_stopwatch(sw);

  /* Assemble entries of V. */
  bem_slp->nearfield(NULL, NULL, bem_slp, false, V);

  t = stop_stopwatch(sw);
  /* Get the total memory footprint for V. */
  size = getsize_amatrix(V) / 1024.0 / 1024.0;

  printf("  %.2f s\n", t);
  printf("  %.3f MB\n", size);
  printf("rows = %d, cols = %d\n", V->rows, V->cols);

  char outFile[100];
  sprintf(outFile, "h2lib_bem2d_circle_laplace_%d.csv", gr->vertices);
  freopen(outFile, "w", stdout);
  print_amatrix(V);
  fclose(stdout);
  printf("Matrix successfully written into %s\n", outFile);

  char geomFile[100];
  sprintf(geomFile, "h2lib_bem2d_circle_laplace_%d.geom", bem_slp->gr->vertices);
  freopen(geomFile, "w", stdout);
  print_curve2d(bem_slp->gr);
  fclose(stdout);
  printf("Geometry information successfully written into %s\n", geomFile);

  /****************************************************
   * cleanup
   ****************************************************/

  del_amatrix(V);
  del_bem2d(bem_slp);
  del_curve2d(gr);
  del_stopwatch(sw);

  /* Uninit the H2Lib. */
  uninit_h2lib();

  return 0;
}
