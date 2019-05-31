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
  DenseMatrix M, D, S;
  LocalStressMatrices_(v, T, M, D, S);

  DenseMatrix DT(D.NumCols(), D.NumRows()), ST(S.NumCols(), S.NumRows());
  DT.Transpose(D);
  ST.Transpose(S);
  M.Inverse();

  auto Q11 = DT * M * D;
  auto Q12 = DT * M * S;
  auto Q22 = ST * M * S;
  DenseMatrix Q21(Q12.NumCols(), Q12.NumRows());
  Q21.Transpose(Q12);
  Q22.Inverse();

  A = Q11 - Q12 * Q22 * Q21;
  
  return WHETSTONE_ELEMENTAL_MATRIX_OK;
}


/* ******************************************************************
* Mass matrix in the space of stresses.
****************************************************************** */
void MFD3D_Elasticity::LocalStressMatrices_(
   int v, const std::vector<Tensor>& T,
   DenseMatrix& M, DenseMatrix& D, DenseMatrix& S)
{
  Entity_ID_List cells, faces, cfaces, vcfaces;
  std::vector<int> cdirs;

  mesh_->node_get_faces(v, AmanziMesh::Parallel_type::ALL, &faces);
  mesh_->node_get_cells(v, AmanziMesh::Parallel_type::ALL, &cells);
  int nfaces = faces.size();
  int ncells = cells.size();

  int nk = d_ * (d_ - 1) / 2;
  M.Reshape(d_ * nfaces, d_ * nfaces);
  D.Reshape(d_ * nfaces, d_ * ncells);
  S.Reshape(d_ * nfaces, nk);

  M.PutScalar(0.0);
  D.PutScalar(0.0);
  S.PutScalar(0.0);

  DenseMatrix Mcorner, N(d_, d_), R(d_, d_);
  DenseVector Dcorner(d_);

  // skew-symmetric matrix
  int m(0);
  std::vector<DenseMatrix> Skew(nk);
  for (int i = 0; i < d_; ++i) {
    for (int j = i + 1; j < d_; ++j) {
      Skew[m].Reshape(d_, d_);
      Skew[m].PutScalar(0.0);

      Skew[m](i, j) = 1.0;
      Skew[m](j, i) =-1.0;
      m++;
    }
  }

  for (int n = 0; n < ncells; ++n) {
    int c = cells[n];
    const AmanziGeometry::Point& xc = mesh_->cell_centroid(c);

    mesh_->cell_get_faces_and_dirs(c, &cfaces, &cdirs);
    mesh_->node_get_cell_faces(v, c, Parallel_type::ALL, &vcfaces);
    int nvc = vcfaces.size();
    AMANZI_ASSERT(nvc == d_);

    // generate corner matrix for one component
    for (int i = 0; i < nvc; i++) {
      int f = vcfaces[i];
      int k = std::distance(cfaces.begin(), std::find(cfaces.begin(), cfaces.end(), f));

      double area = mesh_->face_area(f);
      const AmanziGeometry::Point& xf = mesh_->face_centroid(f);
      AmanziGeometry::Point conormal = T[n] * mesh_->face_normal(f);

      double facet_area = area / 2;  // FIXME
      for (int j = 0; j < d_; ++j) {
        N(i, j) = conormal[j] / area;
        R(i, j) = cdirs[k] * (xf[j] - xc[j]) * facet_area;
      }

      Dcorner(i) = cdirs[k] * facet_area;
    }

    N.Inverse();
    Mcorner = R * N;

    // assemble mass matrices
    for (int i = 0; i < nvc; i++) {
      int k = std::distance(faces.begin(), std::find(faces.begin(), faces.end(), vcfaces[i]));

      for (int j = 0; j < nvc; ++j) {
        int l = std::distance(faces.begin(), std::find(faces.begin(), faces.end(), vcfaces[j]));

        for (int s = 0; s < d_; ++s) { 
          M(d_ * k + s, d_ * l + s) += Mcorner(i, j);
        }
      }
    }

    // assemble divergence matrices
    for (int i = 0; i < nvc; i++) {
      int k = std::distance(faces.begin(), std::find(faces.begin(), faces.end(), vcfaces[i]));
      for (int s = 0; s < d_; ++s) { 
        D(d_ * k + s, d_ * n + s) = Dcorner(i);
      }
    }

    // assemble rotation matrices
    for (int m = 0; m < nk; ++m) {
      auto tmp = R * Skew[m];

      for (int i = 0; i < nvc; i++) {
        int k = std::distance(faces.begin(), std::find(faces.begin(), faces.end(), vcfaces[i]));
        for (int s = 0; s < d_; ++s) { 
          S(d_ * k + s, m) += tmp(i, s);
        }
      }
    }
  }
}

}  // namespace WhetStone
}  // namespace Amanzi

