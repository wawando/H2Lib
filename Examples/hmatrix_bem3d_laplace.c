#include <stdio.h>

#include "basic.h"
#include "krylovsolvers.h"
#include "laplacebem3d.h"

/****************************************************
 * This examples sets up single layer potential operator(SLP),
 * double layer potential operator(DLP) as well as the mass matrix M
 * as H-matrices in order to solve the interior Dirichlet
 * problem for the Laplace equation.
 ****************************************************/

int
main(int argc, char **argv)
{
  pstopwatch sw;
  pmacrosurface3d mg;
  psurface3d gr;
  pbem3d    bem_slp;
  uint      q_reg, q_sing;
  basisfunctionbem3d basis;
  pcluster  root;
  uint      clf;
  pblock    broot;
  real      eta;
  phmatrix  V;
  uint      m;
  real      eps_solve;
  uint      maxiter;
  real      t, size, norm;

  /* Init the H2Lib, should be called before any other function. */
  init_h2lib(&argc, &argv);

  /****************************************************
   * Set up basic parameters
   ****************************************************/

  /* Number of quadrature points for regular integrals. */
  q_reg = 2;
  /* Number of quadrature points for singular integrals. */
  q_sing = q_reg + 2;
  /* Basis functions that should be used. */
  basis = BASIS_CONSTANT_BEM3D;

  /* Number of interpolation points */
  m = 4;

  /* Minimal leaf size for cluster tree construction. */
  clf = 2 * m * m * m;
  /* Parameter 'eta' within the admissibilty condition. */
  eta = 1.4;

  /* absolute norm of the residuum for CG-method */
  eps_solve = 1.0e-10;

  /* maximum number of CG-steps that should be performed. */
  maxiter = 500;

  /* Stopwatch for measuring the time. */
  sw = new_stopwatch();

  /****************************************************
   * Create geometry
   ****************************************************/

  /* Create abstract geometry of a sphere. */
  mg = new_sphere_macrosurface3d();
  /* Mesh the abstract geometry with variable levels of refinement. */
  gr = build_from_macrosurface3d_surface3d(mg, argc > 1 ? atoi(argv[1]) : 4);
  printf("Created geometry with %d vertices, %d edges and %d triangles\n",
	 gr->vertices, gr->edges, gr->triangles);

  /****************************************************
   * Set up basis data structures for H-matrix approximations
   ****************************************************/

  /* Create a new BEM-object, that can compute entries of SLP operator. */
  bem_slp = new_slp_laplace_bem3d(gr, q_reg, q_sing, basis, basis);

  /* Create cluster tree. */
  root = build_bem3d_cluster(bem_slp, clf, basis);
  
  /* Create block tree. */
  broot = build_nonstrict_block(root, root, &eta, admissible_2_cluster);

  /* Set up interpolation approximation scheme for H-matrix V. */
  setup_hmatrix_aprx_inter_row_bem3d(bem_slp, root, root, broot, m);

  /* /\**************************************************** */
  /*  * Assemble H-matrix SLP */
  /*  ****************************************************\/ */

  printf("Assemble H-matrix V:\n");

  /* /\* Create H-matrix structure from block tree. *\/ */
  V = build_from_block_hmatrix(broot, m * m * m);

  start_stopwatch(sw);
  /* /\* Assemble near- and farfield entries of V. *\/ */
  assemble_bem3d_hmatrix(bem_slp, broot, V);
  t = stop_stopwatch(sw);
  /* /\* Get the total memory footprint for V. *\/ */
  size = getsize_hmatrix(V) / 1024.0 / 1024.0;

  printf("  %.2f s\n", t);
  printf("  %.3f MB\n", size);

  printf("%d x %d blocks\n", V->rsons, V->csons);

  /****************************************************
   * cleanup
   ****************************************************/

  del_hmatrix(V);
  del_block(broot);
  /* Permutation array for Dofs was automatically created by
   * 'build_bem3d_cluster', has to be free before the cluster tree. */
  freemem(root->idx);
  del_cluster(root);
  del_bem3d(bem_slp);
  del_macrosurface3d(mg);
  del_surface3d(gr);
  del_stopwatch(sw);

  /* Uninit the H2Lib. */
  uninit_h2lib();

  return 0;
}
