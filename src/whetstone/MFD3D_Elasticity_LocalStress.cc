/*
  WhetStone, Version 2.2
  Release name: naka-to.

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  The mimetic finite difference method for elasticity: local stress approach
*/

#include <cmath>
#include <iterator>
#include <vector>

#include "Mesh.hh"
#include "Point.hh"
#include "errors.hh"

#include "DenseMatrix.hh"
#include "MFD3D_Elasticity.hh"
#include "Tensor.hh"

namespace Amanzi {
namespace WhetStone {

/* ******************************************************************
* Stiffness matrix: a wrapper for other low-level routines
****************************************************************** */
int MFD3D_Elasticity::StiffnessMatrix_LocalStress(
   int v, const std::vector<Tensor>& T, DenseMatrix& A)
{
  Entity_ID_List cells;
  mesh_->node_get_cells(v, AmanziMesh::Parallel_type::ALL, &cells);
  int ncells = cells.size();
  int nrows = 2 * ncells;

  A.Reshape(nrows, nrows);
  A.PutScalar(0.0);

  for (int n = 0; n < nrows; n++) {
    A(n, n) = T[0](0, 0);
  }
  return WHETSTONE_ELEMENTAL_MATRIX_OK;
}


/* ******************************************************************
* Mass matrix in the space of stresses.
****************************************************************** */
int MFD3D_Elasticity::MassMatrix_StressStress_(
   int v, const std::vector<Tensor>& T, DenseMatrix& M)
{
  Entity_ID_List cells, faces, corner_faces;

  mesh_->node_get_faces(v, AmanziMesh::Parallel_type::ALL, &faces);
  mesh_->node_get_cells(v, AmanziMesh::Parallel_type::ALL, &cells);
  int ncells = cells.size();

  M.Reshape(d_ * ncells, d_ * ncells);
  M.PutScalar(0.0);

  DenseMatrix Mcorner(d_ * d_, d_ * d_);

  for (int n = 0; n < ncells; ++n) {
    int c = cells[n];
    mesh_->node_get_cell_faces(v, c, Parallel_type::ALL, &corner_faces);
    AMANZI_ASSERT(corner_faces.size() == d_);

    // generate corner matrix

    // assemble corner matrix
    for (int i = 0; i < d_; i++) {
      int k = std::distance(faces.begin(), std::find(faces.begin(), faces.end(), corner_faces[i]));

      for (int j = i; j < d_; ++j) {
        int l = std::distance(faces.begin(), std::find(faces.begin(), faces.end(), corner_faces[j]));
        M(k, l) += Mcorner(i, j);
        M(l, k) = M(k, l);
      }
    }
  }
}

}  // namespace WhetStone
}  // namespace Amanzi

