/*
  WhetStone, Version 2.2
  Release name: naka-to.

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  Mimetic schemes for generalized polyhedra, i.e. 3D only.
*/

#include "Mesh.hh"

#include "MFD3D_GeneralizedElectromagnetics.hh"
#include "SurfaceCoordinateSystem.hh"
#include "WhetStoneDefs.hh"

namespace Amanzi {
namespace WhetStone {

/* ******************************************************************
* Consistency condition for inner product on a generized polyhedron.
****************************************************************** */
int MFD3D_GeneralizedElectromagnetics::L2consistency(
    int c, const Tensor& T, DenseMatrix& N, DenseMatrix& Mc, bool symmetry)
{
  Entity_ID_List fnodes, edges, fedges, faces;
  std::vector<int> edirs, fdirs, map;

  mesh_->cell_get_faces_and_dirs(c, &faces, &fdirs);
  int nfaces = faces.size();

  mesh_->cell_get_edges(c, &edges);
  int nedges = edges.size();

  N.Reshape(nedges, d_);
  Mc.Reshape(nedges, nedges);

  AmanziGeometry::Point xf(d_), v1(d_), v2(d_), v3(d_), tau(d_);
  AmanziGeometry::Point vv[3];

  // To calculate matrix R, we re-use matrix N
  N.PutScalar(0.0);

  const AmanziGeometry::Point& xc = mesh_->cell_centroid(c);
  double volume = mesh_->cell_volume(c);

  for (int n = 0; n < nfaces; ++n) {
    int f = faces[n];

    // calculate linear map E_int = A * E_ext, where A = inv(B1 + a ee^T) * B2
    // constant a is selected to get exact lifting for constant functions
    // -- extract geometry assuming this order: v0, e0, v1, e1, ..., vN.
    mesh_->face_get_nodes(f, &fnodes);
    mesh_->face_get_edges_and_dirs(f, &fedges, &edirs, true);
    int nfedges = fedges.size();

    DenseVector dint, dext;
    std::vector<double> tri_area;
    std::vector<AmanziGeometry::Point> tri_normal, tri_center;
    std::vector<std::vector<AmanziGeometry::Point> > rt_tri, nd_tri;

    CurvedFaceGeometry_(fnodes, fdirs[n], xf, dint, dext, tri_normal, tri_center, tri_area);
    CurvedFaceBasisFunctions_(fnodes, edirs, xf, dint, dext, tri_normal, tri_area, rt_tri, nd_tri);
    // DenseMatrix A = CurvedFaceLifting_(xf, fnodes, edirs, tri_area, dint, dext);
    DenseMatrix A = CurvedFaceLiftingConstant_(xf, fnodes, edirs, dint, dext);
    
    // integrate over triangles
    AmanziGeometry::Point p0(d_), p1(d_), a(d_);
    std::vector<AmanziGeometry::Point> xe(3);
    DenseVector rext(nfedges), rint(nfedges), rtmp(nfedges);
    
    mesh_->face_to_cell_edge_map(f, c, &map);

    // -- exterior unit normals
    for (int k = 0; k < d_; ++k) {
      AmanziGeometry::Point u0(d_);
      u0[k] = 1.0;

      rint.PutScalar(0.0);

      for (int i = 0; i < nfedges; ++i) {
        int j = (i + 1) % nfedges;
        mesh_->node_get_coordinates(fnodes[i], &p0);
        mesh_->node_get_coordinates(fnodes[j], &p1);

        xe[0] = (p0 + p1) / 2;
        xe[1] = (p1 + xf) / 2;
        xe[2] = (p0 + xf) / 2;

        Tensor R90 = RotationMatrix90(tri_normal[i], false);

        double tmp0(0.0), tmp1(0.0), tmp2(0.0);
        for (int l = 0; l < 3; ++l) {
          a = R90 * (tri_normal[i] ^ (u0 ^ (xe[l] - xc)));
          tmp0 += a * rt_tri[i][l];  
          tmp1 += a * rt_tri[i][l + 3];
          tmp2 += a * rt_tri[i][l + 6];  
        }

        rext(i)  = tmp0 * tri_area[i] / 6;  
        rint(j) += tmp1 * tri_area[i] / 6;  
        rint(i) += tmp2 * tri_area[i] / 6;  
      }

      A.Multiply(rint, rtmp, true);

      // sum up into the matrix R 
      for (int i = 0; i < nfedges; ++i) {
        N(map[i], k) += (rext(i) + rtmp(i)) * fdirs[n];
      }
    }
  }

  // calculate Mc = R (R^T N)^{-1} R^T 
  Tensor Tinv(T);
  Tinv.Inverse();

  for (int i = 0; i < nedges; i++) {
    for (int k = 0; k < d_; ++k) v1[k] = N(i, k);
    v2 = Tinv * v1;

    for (int j = i; j < nedges; j++) {
      for (int k = 0; k < d_; ++k) v3[k] = N(j, k);
      Mc(i, j) = (v2 * v3) / volume;
    }
  }

  // Rows of matrix N are simply tangents. Since N goes to the
  // Gramm-Schmidt orthogonalizetion procedure, we can skip scaling
  // tensorial factor T.
  for (int i = 0; i < nedges; i++) {
    int e = edges[i];
    const AmanziGeometry::Point& tau = mesh_->edge_vector(e);
    double len = mesh_->edge_length(e);
    for (int k = 0; k < d_; ++k) N(i, k) = tau[k] / len;
  }

  return WHETSTONE_ELEMENTAL_MATRIX_OK;
}


/* ******************************************************************
* Consistency condition for a stiffness matrix.
****************************************************************** */
int MFD3D_GeneralizedElectromagnetics::H1consistency(
    int c, const Tensor& T, DenseMatrix& N, DenseMatrix& Ac)
{
  Entity_ID_List fnodes, edges, fedges, faces;
  std::vector<int> fdirs, edirs, map;

  mesh_->cell_get_faces_and_dirs(c, &faces, &fdirs);
  int nfaces = faces.size();

  mesh_->cell_get_edges(c, &edges);
  int nedges = edges.size();

  int nd = 6;  // order_ * (order_ + 2) * (order_ + 3) / 2;
  N.Reshape(nedges, nd);
  Ac.Reshape(nedges, nedges);

  // To calculate matrix R, we re-use matrix N
  AmanziGeometry::Point p0(d_), p1(d_), xf(d_), v1(d_), v2(d_), v3(d_);

  const AmanziGeometry::Point& xc = mesh_->cell_centroid(c);
  double volume = mesh_->cell_volume(c);

  N.PutScalar(0.0);

  for (int n = 0; n < nfaces; ++n) {
    int f = faces[n];

    mesh_->face_get_nodes(f, &fnodes);
    mesh_->face_get_edges_and_dirs(f, &fedges, &edirs);
    int nfedges = fedges.size();

    mesh_->face_to_cell_edge_map(f, c, &map);

    DenseVector dint, dext;
    std::vector<double> tri_area;
    std::vector<AmanziGeometry::Point> tri_normal, tri_center;
    std::vector<std::vector<AmanziGeometry::Point> > rt_tri, nd_tri;

    CurvedFaceGeometry_(fnodes, fdirs[n], xf, dint, dext, tri_normal, tri_center, tri_area);
    CurvedFaceBasisFunctions_(fnodes, edirs, xf, dint, dext, tri_normal, tri_area, rt_tri, nd_tri);
    DenseMatrix A = CurvedFaceLiftingLinear_(xf, fnodes, edirs, dint, dext, tri_normal, tri_center, tri_area);

    // -- exterior
    for (int i = 0; i < nfedges; ++i) {
      int e = fedges[i];
      const AmanziGeometry::Point& xe = mesh_->edge_centroid(e);
      v1 = (xe - tri_center[i]) * (2 * dext(i) * edirs[i] * fdirs[n]);
      for (int k = 0; k < d_; ++k) N(map[i], k) += v1[k];
    }

    // -- interior
    DenseMatrix Rint(nfedges, d_);
    Rint.PutScalar(0.0);

    for (int i = 0; i < nfedges; ++i) {
      int j = (i + 1) % nfedges;
      mesh_->node_get_coordinates(fnodes[i], &p0);
      mesh_->node_get_coordinates(fnodes[j], &p1);

      v1 = ((xf + p0) / 2 - tri_center[i]) * (2 * dint(i) * fdirs[n]);
      v2 = ((xf + p1) / 2 - tri_center[i]) * (2 * dint(j) * fdirs[n]);

      for (int k = 0; k < d_; ++k) {
        Rint(i, k) -= v1[k];
        Rint(j, k) += v2[k];
      }
    }

    A.Transpose();
    auto Rtmp = A * Rint;

    for (int i = 0; i < nfedges; ++i) {
      for (int k = 0; k < d_; ++k) N(map[i], k) += Rtmp(i, k);
    }
  }
  
  // calculate Ac = R (R^T N)^{+} R^T
  for (int i = 0; i < nedges; i++) {
    for (int k = 0; k < d_; ++k) v1[k] = N(i, k);
    v2 = T * v1;

    for (int j = i; j < nedges; j++) {
      for (int k = 0; k < d_; ++k) v3[k] = N(j, k);
      Ac(i, j) = (v2 * v3) / (4 * volume);
    }
  }

  // Matrix N(:, 1:3) are simply tangents
  for (int i = 0; i < nedges; i++) {
    int e = edges[i];
    const AmanziGeometry::Point& xe = mesh_->edge_centroid(e);
    AmanziGeometry::Point tau = mesh_->edge_vector(e);
    double len = mesh_->edge_length(e);

    tau /= len;
    v1 = xe - xc;
    v3 = tau ^ v1;

    for (int k = 0; k < d_; ++k) N(i, k) = tau[k];
    for (int k = 0; k < d_; ++k) N(i, d_ + k) = v3[k];
  }

  return WHETSTONE_ELEMENTAL_MATRIX_OK;
}


/* ******************************************************************
* Mass matrix for genelized polyhedron
****************************************************************** */
int MFD3D_GeneralizedElectromagnetics::MassMatrix(int c, const Tensor& T, DenseMatrix& M)
{
  DenseMatrix N;

  int ok = L2consistency(c, T, N, M, true);
  if (ok) return ok;

  StabilityScalar_(N, M);
  return WHETSTONE_ELEMENTAL_MATRIX_OK;
}


/* ******************************************************************
* Stiffness matrix is calculated by a hybridization algorithm.
****************************************************************** */
int MFD3D_GeneralizedElectromagnetics::StiffnessMatrix(
    int c, const Tensor& T, DenseMatrix& A)
{
  DenseMatrix N;

  int ok = H1consistency(c, T, N, A);
  if (ok) return ok;

  WhetStone::AddGradient(mesh_, c, N);

  StabilityScalar_(N, A);
  return WHETSTONE_ELEMENTAL_MATRIX_OK;
}


/* ******************************************************************
* Geometry of curved face is based on the piecewise-linear facet model.
****************************************************************** */
void MFD3D_GeneralizedElectromagnetics::CurvedFaceGeometry_(
    const Entity_ID_List& nodes, int fdir,
    AmanziGeometry::Point& xf, DenseVector& dint, DenseVector& dext, 
    std::vector<AmanziGeometry::Point>& tri_normal,
    std::vector<AmanziGeometry::Point>& tri_center,
    std::vector<double>& tri_area)
{
  tri_normal.clear();
  tri_center.clear();
  tri_area.clear();

  int nnodes = nodes.size();
  dint.Reshape(nnodes);
  dext.Reshape(nnodes);

  AmanziGeometry::Point p0(d_), p1(d_), v1(d_), v2(d_), normal(d_);

  // calculate geometric center that defines face geometry by mesh assumption
  xf.set(0.0);
  for (int i = 0; i < nnodes; ++i) {
    mesh_->node_get_coordinates(nodes[i], &p0);
    xf += p0;
  }
  xf /= nnodes;

  for (int i = 0; i < nnodes; ++i) {
    int j = (i + 1) % nnodes; 
    mesh_->node_get_coordinates(nodes[i], &p0);
    v1 = xf - p0;
    dint(i) = norm(v1);

    mesh_->node_get_coordinates(nodes[j], &p1);
    v2 = p1 - p0;
    dext(i) = norm(v2);

    normal = (v2 ^ v1) / 2;
    double area = norm(normal);
    normal /= fdir * area;

    tri_normal.push_back(normal);
    tri_center.push_back((p0 + p1 + xf) / 3); 
    tri_area.push_back(area);
  }
}


/* ******************************************************************
* Oriented Raviart Thomas basis functions at quadrature points.
* The basis functions are ordered counter-clockwise.
****************************************************************** */
void MFD3D_GeneralizedElectromagnetics::CurvedFaceBasisFunctions_(
    const Entity_ID_List& nodes, const std::vector<int>& dirs,
    const AmanziGeometry::Point& xf, const DenseVector& dint, const DenseVector& dext, 
    const std::vector<AmanziGeometry::Point>& tri_normal,
    const std::vector<double>& tri_area,
    std::vector<std::vector<AmanziGeometry::Point> >& rt_tri,
    std::vector<std::vector<AmanziGeometry::Point> >& nd_tri)
{
  int nedges = dext.NumRows();
  rt_tri.resize(nedges);
  nd_tri.resize(nedges);

  double area, scale;
  AmanziGeometry::Point p0(d_), p1(d_), xe0(d_), xe1(d_), xe2(d_);

  for (int i = 0; i < nedges; ++i) {
    int j = (i + 1) % nedges;
    mesh_->node_get_coordinates(nodes[i], &p0);
    mesh_->node_get_coordinates(nodes[j], &p1);

    xe0 = (p0 + p1) / 2;
    xe1 = (p1 + xf) / 2;
    xe2 = (p0 + xf) / 2;

    // order: { phi0(x_k) }, { phi1(x_k) }, { phi2(x_k) }
    rt_tri[i].resize(9);

    area = 2 * tri_area[i];
    scale = dirs[i] * dext(i) / area;
    rt_tri[i][0] = (xe0 - xf) * scale;
    rt_tri[i][1] = (xe1 - xf) * scale;
    rt_tri[i][2] = (xe2 - xf) * scale;

    scale = dint(j) / area;
    rt_tri[i][3] = (xe0 - p0) * scale;
    rt_tri[i][4] = (xe1 - p0) * scale;
    rt_tri[i][5] = (xe2 - p0) * scale;

    scale = -dint(i) / area;
    rt_tri[i][6] = (xe0 - p1) * scale;
    rt_tri[i][7] = (xe1 - p1) * scale;
    rt_tri[i][8] = (xe2 - p1) * scale;

    // order: { phi0(x_k) }, { phi1(x_k) }, { phi2(x_k) }
    nd_tri[i].resize(6);

    Tensor R90 = RotationMatrix90(tri_normal[i], true);
    area = 2 * tri_area[i];
    scale = dirs[i] * dext(i) / area;
    nd_tri[i][0] = (R90 * (p0 - xf)) * scale;
    nd_tri[i][1] = (R90 * (p1 - xf)) * scale;

    scale = dint(j) / area;
    nd_tri[i][2] = (R90 * (p1 - p0)) * scale;
    nd_tri[i][3] = (R90 * (xf - p0)) * scale;

    scale = -dint(i) / area;
    nd_tri[i][4] = (R90 * (p0 - p1)) * scale;
    nd_tri[i][5] = (R90 * (xf - p1)) * scale;
  }
}


/* ******************************************************************
* Lifting of boundary data that preserves constant functions.
****************************************************************** */
DenseMatrix MFD3D_GeneralizedElectromagnetics::CurvedFaceLifting_(
    const AmanziGeometry::Point& xf, const Entity_ID_List& nodes,
    const std::vector<int>& dirs, const std::vector<double>& tri_area,
    const DenseVector& dint, const DenseVector& dext)
{
  int nedges = dext.NumRows();

  // form constraints equation
  DenseMatrix B1(nedges, nedges), B2(nedges, nedges);
  B1.PutScalar(0.0);

  double area(0.0);
  for (int i = 0; i < nedges; ++i)
    area += tri_area[i];
    
  for (int i = 0; i < nedges; ++i) {
    int j = (i + 1) % nedges; 
    B1(i, i) =-dint(i);
    B1(i, j) = dint(j);

    double weight = tri_area[i] / area;
    for (int k = 0; k < nedges; ++k) {
      B2(i, k) = dext(k) * dirs[k] * weight;
    }
    B2(i, i) -= dext(i) * dirs[i];
  }

  // form blocks of quadratic cost functional
  AmanziGeometry::Point p0(d_), p1(d_);
  DenseMatrix L1(nedges, d_), L2(nedges, d_);
  DenseMatrix L1T(d_, nedges), L2T(d_, nedges);

  for (int i = 0; i < nedges; ++i) {
    int j = (i + 1) % nedges;
    mesh_->node_get_coordinates(nodes[i], &p0);
    mesh_->node_get_coordinates(nodes[j], &p1);

    for (int k = 0; k < d_; ++k) {
      L1(i, k) = (xf[k] - p0[k]) / dint(i);
      L2(i, k) = (p1[k] - p0[k]) * dirs[i] / dext(i); 
    }
  }

  L1T.Transpose(L1);
  L2T.Transpose(L2);

  auto S = L1T * L1 + L2T * L2;
  S.InverseMoorePenrose();

  auto M1 = L1 * S * L1T;
  auto M2 = L1 * S * L2T;

  M1 *= -1.0;
  M2 *= -1.0;
  for (int i = 0; i < nedges; ++i) M1(i, i) += 1.0;

  // form the lifting matrix
  DenseMatrix A(nedges, nedges), B1T(B1), RT(nedges, nedges);

  B1T.Transpose();
  M1.Inverse();

  auto Q = B1 * M1 * B1T;
  Q.InverseMoorePenrose();

  auto R = B1 * M1;
  RT.Transpose(R);

  A = RT * Q * (R * M2 + B2) - M1 * M2;
  return A;
}


/* ******************************************************************
* Lifting of boundary data that preserves constant functions.
****************************************************************** */
DenseMatrix MFD3D_GeneralizedElectromagnetics::CurvedFaceLiftingConstant_(
    const AmanziGeometry::Point& xf, const Entity_ID_List& nodes,
    const std::vector<int>& dirs,
    const DenseVector& dint, const DenseVector& dext)
{
  int nedges = dext.NumRows();

  // form blocks of quadratic cost functional
  AmanziGeometry::Point p0(d_), p1(d_);
  DenseMatrix L1(nedges, d_), L2(nedges, d_);
  DenseMatrix L1T(d_, nedges), L2T(d_, nedges);

  for (int i = 0; i < nedges; ++i) {
    int j = (i + 1) % nedges;
    mesh_->node_get_coordinates(nodes[i], &p0);
    mesh_->node_get_coordinates(nodes[j], &p1);

    for (int k = 0; k < d_; ++k) {
      L1(i, k) = (xf[k] - p0[k]) / dint(i);
      L2(i, k) = (p1[k] - p0[k]) * dirs[i] / dext(i); 
    }
  }

  L1T.Transpose(L1);
  L2T.Transpose(L2);

  auto S = L1T * L1 + L2T * L2;
  S.InverseMoorePenrose();

  auto M1 = L1 * S * L1T;
  auto M2 = L1 * S * L2T;

  M1 *= -1.0;
  for (int i = 0; i < nedges; ++i) M1(i, i) += 1.0;

  M1.Inverse();
  auto A = M1 * M2;
  return A;
}


/* ******************************************************************
* Lifting of boundary data that preserves Linear functions.
****************************************************************** */
DenseMatrix MFD3D_GeneralizedElectromagnetics::CurvedFaceLiftingLinear_(
    const AmanziGeometry::Point& xf, const Entity_ID_List& nodes,
    const std::vector<int>& dirs,
    const DenseVector& dint, const DenseVector& dext,
    const std::vector<AmanziGeometry::Point>& tri_normal,
    const std::vector<AmanziGeometry::Point>& tri_center,
    const std::vector<double>& tri_area)
{
  int nedges = dext.NumRows();

  // form blocks of quadratic cost functional
  AmanziGeometry::Point p0(d_), p1(d_), v0(d_), v1(d_), v2(d_), v3(d_), v4(d_), v5(d_), v6(d_);
  DenseMatrix L1(nedges, 2 * d_), L2(nedges + 3, 2 * d_);
  DenseMatrix L1T(2 * d_, nedges), L2T(2 * d_, nedges + 3);

  auto normal = tri_normal[0];
  SurfaceCoordinateSystem coordsys(xf, normal);
  const auto& tau = *coordsys.tau();

  L2.PutScalar(0.0);

  double area(0.0);
  for (int i = 0; i < nedges; ++i) area += tri_area[i];

  for (int i = 0; i < nedges; ++i) {
    int j = (i + 1) % nedges;
    mesh_->node_get_coordinates(nodes[i], &p0);
    mesh_->node_get_coordinates(nodes[j], &p1);

    v0 = xf - p0;
    v1 = p1 - p0;

    v2 = ((xf + p0) ^ v0) / 2;
    v3 = ((p1 + p0) ^ v1) / 2;

    for (int k = 0; k < d_; ++k) {
      L1(i, k) = v0[k] / dint(i);
      L2(i, k) = v1[k] * dirs[i] / dext(i); 

      L1(i, d_ + k) = v2[k] / dint(i);
      L2(i, d_ + k) = v3[k] * dirs[i] / dext(i); 

      v6 = AmanziGeometry::Point(1.0, 0.0, 0.0);
      v4 = (tri_center[i] - normal) ^ v6;
      v5 = (tri_center[i] - normal) ^ (v6 ^ tri_center[i]);

      L2(nedges, k) += v4[k] * tri_area[i] / area;
      L2(nedges, d_ + k) += v5[k] * tri_area[i] / area;

      v6 = AmanziGeometry::Point(0.0, 1.0, 0.0);
      v4 = (tri_center[i] - normal) ^ v6;
      v5 = (tri_center[i] - normal) ^ (v6 ^ tri_center[i]);

      L2(nedges + 1, k) += v4[k] * tri_area[i] / area;
      L2(nedges + 1, d_ + k) += v5[k] * tri_area[i] / area;

      v6 = AmanziGeometry::Point(0.0, 0.0, 1.0);
      v4 = (tri_center[i] - normal) ^ v6;
      v5 = (tri_center[i] - normal) ^ (v6 ^ tri_center[i]);

      L2(nedges + 2, k) += v4[k] * tri_area[i] / area;
      L2(nedges + 2, d_ + k) += v5[k] * tri_area[i] / area;
    }
  }

  L1T.Transpose(L1);
  L2T.Transpose(L2);

  auto S = L1T * L1 + L2T * L2;
  S.InverseMoorePenrose();

  auto M1 = L1 * S * L1T;
  auto M2 = L1 * S * L2T;

  M1 *= -1.0;
  for (int i = 0; i < nedges + 3; ++i) M1(i, i) += 1.0;

  M1.InverseMoorePenrose();
  auto A = M1 * M2;
  return A;
}

}  // namespace WhetStone
}  // namespace Amanzi
