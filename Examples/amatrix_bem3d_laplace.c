#include <stdio.h>

#include "basic.h"
#include "krylovsolvers.h"
#include "laplacebem3d.h"

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
  pmacrosurface3d mg;
  psurface3d gr;
  pbem3d    bem_slp, bem_dlp;
  uint      q_reg, q_sing;
  basisfunctionbem3d basis;
  pamatrix  V, KM;
  pavector  gd, b, x;
  real      eps_solve;
  uint      maxiter;
  real      t, size, norm;

  /* Init the H2Lib, should be called before any other function. */
  init_h2lib(&argc, &argv);

  /****************************************************
   * Set up basic parameters
   ****************************************************/

  /* Number of quadrature points for regular integrals. */
  q_reg = 4;
  /* Number of quadrature points for singular integrals. */
  q_sing = q_reg + 2;
  /* Basis functions that should be used. */
  basis = BASIS_CONSTANT_BEM3D;

  /* absolute norm of the residuum for CG-method */
  eps_solve = 1.0e-4;

  /* maximum number of CG-steps that should be performed. */
  maxiter = 250;

  /* Stopwatch for measuring the time. */
  sw = new_stopwatch();

  /****************************************************
   * Create geometry
   ****************************************************/

  /* Create abstract geometry of a sphere. */
  mg = new_sphere_macrosurface3d();
  /* Mesh the abstract geometry with 32 levels of refinement. */
  gr = build_from_macrosurface3d_surface3d(mg, argc > 1 ? atoi(argv[1]) : 4);
  printf("Created geometry with %d vertices, %d edges and %d triangles\n",
	 gr->vertices, gr->edges, gr->triangles);

  /****************************************************
   * Set up bem objects
   ****************************************************/

  /* Create a new BEM-object, that can compute entries of SLP operator. */
  bem_slp = new_slp_laplace_bem3d(gr, q_reg, q_sing, basis, basis);
  /* Create a new BEM-object, that can compute entries of DLP operator
   * and 0.5*I. */
  bem_dlp = new_dlp_laplace_bem3d(gr, q_reg, q_sing, basis, basis, 0.5);

  /****************************************************
   * Assemble Dense matrix SLP
   ****************************************************/

  printf("Assemble dense matrix V:\n");

  /* Create amatrix structure. */
  V = new_amatrix(gr->triangles, gr->triangles);

  start_stopwatch(sw);
  /* Assemble entries of V. */
  assemble_bem3d_amatrix(bem_slp, V);
  t = stop_stopwatch(sw);
  /* Get the total memory footprint for V. */
  /* size = getsize_amatrix(V) / 1024.0 / 1024.0; */

  printf("  %.2f s\n", t);
  printf("rows = %d, cols = %d\n", V->rows, V->cols);

  char outFile[100];
  sprintf(outFile, "h2lib_bem3d_laplace_%d.csv", gr->triangles);
  freopen(outFile, "w", stdout);
  print_amatrix(V);
  fclose(stdout);

  /****************************************************
   * cleanup
   ****************************************************/

  del_amatrix(V);
  del_bem3d(bem_slp);
  del_bem3d(bem_dlp);
  del_macrosurface3d(mg);
  del_surface3d(gr);
  del_stopwatch(sw);

  /* Uninit the H2Lib. */
  uninit_h2lib();

  return 0;
}
