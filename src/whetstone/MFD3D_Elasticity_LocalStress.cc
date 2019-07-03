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
   int v, const std::vector<Tensor>& T, DenseMatrix& A, DenseMatrix& B)
{
  DenseMatrix M, D, S;
  LocalStressMatrices_(v, T, M, D, S);

  DenseMatrix DT(D.NumCols(), D.NumRows()), ST(S.NumCols(), S.NumRows());
  DT.Transpose(D);
  ST.Transpose(S);
  M.InverseMoorePenrose();

  auto Q11 = DT * M * D;
  auto Q12 = DT * M * S;
  auto Q22 = ST * M * S;
  auto Q21 = ST * M * D;
  Q22.InverseMoorePenrose();

  A = Q11 - Q12 * Q22 * Q21;
  B = (DT - Q12 * Q22 * ST) * M;
  
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
  AMANZI_ASSERT(ncells > 0);

  int nd = d_ * (d_ + 1) / 2;  // symmetric tensors
  int nk = d_ * (d_ - 1) / 2;  // skew-symmetric tensors
  M.Reshape(d_ * nfaces, d_ * nfaces);
  D.Reshape(d_ * nfaces, d_ * ncells);
  S.Reshape(d_ * nfaces, nk);

  M.PutScalar(0.0);
  D.PutScalar(0.0);
  S.PutScalar(0.0);

  // basis of stresses (symmetric stresses + non-symmetric)
  std::vector<Tensor> vE;

  for (int i = 0; i < d_; ++i) {
    Tensor E(d_, 2);
    E(i, i) = 1.0;
    vE.push_back(E);
  }

  for (int i = 0; i < d_; ++i) {
    for (int j = i + 1; j < d_; ++j) {
      Tensor E(d_, 2);
      E(i, j) = E(j, i) = 1.0;
      vE.push_back(E);
    }
  }

  for (int i = 0; i < d_; ++i) {
    for (int j = i + 1; j < d_; ++j) {
      Tensor E(d_, 2);
      E(i, j) = 1.0;
      E(j, i) =-1.0;
      vE.push_back(E);
    }
  }

  // node
  AmanziGeometry::Point xv(d_);
  mesh_->node_get_coordinates(v, &xv);

  // main loop over cells around the given node
  for (int n = 0; n < ncells; ++n) {
    int c = cells[n];
    const AmanziGeometry::Point& xc = mesh_->cell_centroid(c);

    mesh_->cell_get_faces_and_dirs(c, &cfaces, &cdirs);
    mesh_->node_get_cell_faces(v, c, Parallel_type::ALL, &vcfaces);
    int nvc = vcfaces.size();
    AMANZI_ASSERT(nvc == d_);

    // -- generate auxiliary corner matrix for one component
    int mx = nvc * d_, nx = d_ * d_;
    DenseMatrix Mcorner(mx, nx), Ncorner(mx, nx), Rcorner(mx, nx);
    DenseVector Dcorner(nvc);

    for (int i = 0; i < nvc; i++) {
      int f = vcfaces[i];
      int m = std::distance(cfaces.begin(), std::find(cfaces.begin(), cfaces.end(), f));

      double area = mesh_->face_area(f);
      const AmanziGeometry::Point& xf = mesh_->face_centroid(f);
      const AmanziGeometry::Point& normal = mesh_->face_normal(f);

      double facet_area = area / 2;  // FIXME
      Dcorner(i) = cdirs[m] * facet_area;

      for (int j = 0; j < d_ * d_; ++j) {
        auto conormal = (T[n] * vE[j]) * (normal / area);
        // auto dx = vE[j] * ((xf + xv) / 2 - xc);
        auto dx = vE[j] * (xf - xc);

        for (int k = 0; k < d_; ++k) {
          Ncorner(nvc * k + i, j) = conormal[k];
          Rcorner(nvc * k + i, j) = cdirs[m] * dx[k] * facet_area;
        }
      }
    }

    /*
    auto R = Rcorner.SubMatrix(0, mx, 0, nd);
    auto N = Ncorner.SubMatrix(0, mx, 0, nd);

    DenseMatrix NT(R.NumCols(), N.NumRows());
    NT.Transpose(N);
    
    auto NN = NT * N;
    NN.Inverse();
    Mcorner = R * NN * NT;
    StabilityScalar_(N, Mcorner);
    */

    Ncorner.InverseMoorePenrose();
    Mcorner = Rcorner * Ncorner;

    // assemble mass matrices
    for (int i = 0; i < nvc; i++) {
      int k = std::distance(faces.begin(), std::find(faces.begin(), faces.end(), vcfaces[i]));

      for (int j = 0; j < nvc; ++j) {
        int l = std::distance(faces.begin(), std::find(faces.begin(), faces.end(), vcfaces[j]));

        for (int s = 0; s < d_; ++s) { 
          for (int r = 0; r < d_; ++r) { 
            M(d_ * k + s, d_ * l + r) += Mcorner(nvc * s + i, nvc * r + j);
          }
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
      for (int i = 0; i < nvc; i++) {
        int k = std::distance(faces.begin(), std::find(faces.begin(), faces.end(), vcfaces[i]));
        for (int s = 0; s < d_; ++s) { 
          S(d_ * k + s, m) += Rcorner(nvc * s + i, nd + m);
        }
      }
    }
  }
}

}  // namespace WhetStone
}  // namespace Amanzi

