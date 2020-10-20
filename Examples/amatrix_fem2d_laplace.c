/*----------------------------------------------------------------------------*/
/* Working with linear p1 elements and hierarchical matrices for laplace                                            */
/*----------------------------------------------------------------------------*/
#include <stdio.h>

#include "parameters.h"		/* Read parameters interactively */

#include "tri2d.h"		/* 2-dimensional mesh */
#include "tri2dp1.h"		/* discretisation with linear p1 elements based on tri2d.h */

#include "tet3d.h"		/* 3-dimensional mesh */
#include "tet3dp1.h"		/* discretisation with linear p1 elements based on tet3d.h */

#include "ddcluster.h"		/* Domain decomposition clustering */
#include "hmatrix.h"		/* Hierarchical matrices */
#include "harith.h"		/* Arithmetic for hierarchical matrices */

#include "truncation.h"		/* Auxiliary function for truncation */
#include "hcoarsen.h"		/* Coarsening of hierarchical matrices */
#include "matrixnorms.h"	/* Computing difference norms for two matrices */

real
norm2lu_sparsematrix(pchmatrix LU, pcsparsematrix sp)
{
  avector   tmp1, tmp2;
  uint      rows = LU->rc->size;
  uint      cols = LU->cc->size;
  pavector  x, y;
  real      norm;
  uint      j;

  assert(sp->rows == rows);
  assert(sp->cols == cols);

  x = init_avector(&tmp1, cols);
  y = init_avector(&tmp2, rows);

  random_avector(x);
  norm = norm2_avector(x);
  scale_avector(1.0 / norm, x);

  for (j = 0; j < NORM_STEPS; j++) {
    // printf("norm = %g \n", sqrt( norm));
    clear_avector(y);
    mvm_sparsematrix_avector(1.0, false, sp, x, y);
    triangularsolve_hmatrix_avector(true, true, false, LU, y);
    triangularsolve_hmatrix_avector(false, false, false, LU, y);
    add_avector(-1.0, y, x);
    copy_avector(x, y);
    triangularsolve_hmatrix_avector(false, false, true, LU, y);
    triangularsolve_hmatrix_avector(true, true, true, LU, y);
    mvm_sparsematrix_avector(-1.0, true, sp, y, x);
    norm = norm2_avector(x);
    scale_avector(1.0 / norm, x);
  }

  uninit_avector(y);
  uninit_avector(x);

  return REAL_SQRT(norm);
}

real
norm2chol_sparsematrix(pchmatrix LU, pcsparsematrix sp)
{
  avector   tmp1, tmp2;
  uint      rows = LU->rc->size;
  uint      cols = LU->cc->size;
  pavector  x, y;
  real      norm;
  uint      j;

  assert(sp->rows == rows);
  assert(sp->cols == cols);

  x = init_avector(&tmp1, cols);
  y = init_avector(&tmp2, rows);

  random_avector(x);
  norm = norm2_avector(x);
  scale_avector(1.0 / norm, x);

  for (j = 0; j < NORM_STEPS; j++) {
    // printf("norm = %g \n", sqrt( norm));
    clear_avector(y);
    mvm_sparsematrix_avector(1.0, false, sp, x, y);
    triangularsolve_hmatrix_avector(true, false, false, LU, y);
    triangularsolve_hmatrix_avector(true, false, true, LU, y);
    add_avector(-1.0, y, x);
    copy_avector(x, y);
    triangularsolve_hmatrix_avector(true, false, false, LU, y);
    triangularsolve_hmatrix_avector(true, false, true, LU, y);
    mvm_sparsematrix_avector(-1.0, true, sp, y, x);
    norm = norm2_avector(x);
    scale_avector(1.0 / norm, x);
  }

  uninit_avector(y);
  uninit_avector(x);

  return REAL_SQRT(norm);
}

void
print_dof_tri2d(pctri2d t2) {
  const real(*x)[2] = (const real(*)[2]) t2->x;
  uint ndof = 0;
  uint i;
  for(i = 0; i < t2->vertices; i++) {
    if(t2->xb[i] == 0) ndof++;
  }
  printf("%u 2\n", ndof);
  for(i = 0; i < t2->vertices; i++) {
    if(t2->xb[i] == 0) {
      printf("%.10e %.10e\n", x[i][0], x[i][1]);
    }
  }
}


int
main(int argc, char **argv)
{

  uint      L;			/* Number of grid refinements */
  uint      i, j;		/* Auxiliary variable for loops */
  ptri2d   *gr_2d;		/* 2d mesh hierarchy */
  ptri2dp1  p1_2d;		/* Linear p1 basis functions in 2d */
  psparsematrix sp;		/* Sparsematrix object */
  pamatrix V;

  /* First initialise the library */
  init_h2lib(&argc, &argv);

  /* L = askforint("Refinement?\n", "h2lib_L", 4); */
  L = argc > 1 ? atoi(argv[1]) : 4;

  /* Build geometry and discretisation for laplace equation */
  printf("========================================\n"
         "  Create and fill fem2d sparsematrix\n");
  /* Mesh hierarchy */
  gr_2d = (ptri2d *) allocmem((size_t) sizeof(ptri2d) * (L + 1));
  /* gr_2d[0] = new_unitsquare_tri2d();	/\* Set domain *\/ */
  gr_2d[0] = new_unitcircle_tri2d();
  //gr_2d[0] = new_lshape_tri2d();
  for (i = 0; i < L; i++) {	/* Mesh refinements */
    gr_2d[i + 1] = refine_tri2d(gr_2d[i], NULL);
  }
  check_tri2d(gr_2d[L]);	/* Check mesh for inconsistencies */
  printf("Created geometry with %u vertices, %u edges, %u triangles\n", gr_2d[L]->vertices, gr_2d[L]->edges, gr_2d[L]->triangles);

  p1_2d = new_tri2dp1(gr_2d[L]);	/* Build discretisation */
  sp = build_tri2dp1_sparsematrix(p1_2d);	/* Build corresponding sparsematrix */
  assemble_tri2dp1_laplace_sparsematrix(p1_2d, sp, 0);	/* Fill the sparsematrix */
  /* Convert sparsematrix to amatrix for printing*/
  V = new_amatrix(sp->rows, sp->cols);
  init_zero_amatrix(V, sp->rows, sp->cols);
  add_sparsematrix_amatrix(1.0, false, sp, V);

  printf("rows = %d, cols = %d\n", V->rows, V->cols);

  char outFile[100];
  sprintf(outFile, "h2lib_fem2d_circle_laplace_%d.csv", V->rows);
  freopen(outFile, "w", stdout);
  print_amatrix(V);
  fclose(stdout);

  char geomFile[100];
  sprintf(geomFile, "h2lib_fem2d_circle_laplace_%d.geom", V->rows);
  freopen(geomFile, "w", stdout);
  print_dof_tri2d(gr_2d[L]);
  fclose(stdout);

  del_tri2dp1(p1_2d);
  for (i = 0; i <= L; i++) {
    j = L - i;
    del_tri2d(gr_2d[j]);
  }
  freemem(gr_2d);

  /* Cleaning up */
  del_sparsematrix(sp);
  del_amatrix(V);

  uninit_h2lib();
  return 0;
}
