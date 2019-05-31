/*
  WhetStone, Version 2.2
  Release name: naka-to.

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  Mimetic schemes for generalized polyhedra.
*/

#ifndef AMANZI_MFD3D_GENERALIZED_ELECTROMAGNETICS_HH_
#define AMANZI_MFD3D_GENERALIZED_ELECTROMAGNETICS_HH_

#include "Teuchos_RCP.hpp"

#include "Mesh.hh"

#include "BilinearFormFactory.hh"
#include "DenseMatrix.hh"
#include "MFD3D.hh"
#include "Tensor.hh"

namespace Amanzi {
namespace WhetStone {

class MFD3D_GeneralizedElectromagnetics : public MFD3D { 
 public:
  MFD3D_GeneralizedElectromagnetics(const Teuchos::ParameterList& plist,
                                     const Teuchos::RCP<const AmanziMesh::Mesh>& mesh)
    : MFD3D(mesh) {};
  ~MFD3D_GeneralizedElectromagnetics() {};

  // required member functions
  // -- schema for this element
  virtual std::vector<SchemaItem> schema() const override {
    return std::vector<SchemaItem>(1, std::make_tuple(AmanziMesh::EDGE, DOF_Type::SCALAR, d_));
  }

  // -- mass matrices
  virtual int L2consistency(int c, const Tensor& T, DenseMatrix& N, DenseMatrix& Mc, bool symmetry) override;
  virtual int MassMatrix(int c, const Tensor& T, DenseMatrix& M) override; 

  // -- stiffness matrices
  virtual int H1consistency(int c, const Tensor& T, DenseMatrix& N, DenseMatrix& Ac) override;
  virtual int StiffnessMatrix(int c, const Tensor& K, DenseMatrix& A) override;

 private:
  void CurvedFaceGeometry_(
      const Entity_ID_List& nodes, int fdir,
      AmanziGeometry::Point& xf, DenseVector& dint, DenseVector& dext, 
      std::vector<AmanziGeometry::Point>& tri_normal,
      std::vector<AmanziGeometry::Point>& tri_center,
      std::vector<double>& tri_area);

  void CurvedFaceBasisFunctions_(
      const Entity_ID_List& nodes, const std::vector<int>& dirs,
      const AmanziGeometry::Point& xf, const DenseVector& dint, const DenseVector& dext, 
      const std::vector<AmanziGeometry::Point>& tri_normal,
      const std::vector<double>& tri_area,
      std::vector<std::vector<AmanziGeometry::Point> >& rt_tri,
      std::vector<std::vector<AmanziGeometry::Point> >& nd_tri);

  DenseMatrix CurvedFaceLifting_(
      const AmanziGeometry::Point& xf, const Entity_ID_List& nodes,
      const std::vector<int>& dirs, const std::vector<double>& tri_area,
      const DenseVector& dint, const DenseVector& dext);

  DenseMatrix CurvedFaceLiftingConstant_(
      const AmanziGeometry::Point& xf, const Entity_ID_List& nodes, const std::vector<int>& dirs,
      const DenseVector& dint, const DenseVector& dext);

  DenseMatrix CurvedFaceLiftingLinear_(
      const AmanziGeometry::Point& xf, const Entity_ID_List& nodes,
      const std::vector<int>& dirs,
      const DenseVector& dint, const DenseVector& dext,
      const std::vector<AmanziGeometry::Point>& tri_normal,
      const std::vector<AmanziGeometry::Point>& tri_center,
      const std::vector<double>& tri_area);

 private:
  static RegisteredFactory<MFD3D_GeneralizedElectromagnetics> factory_;
};

}  // namespace WhetStone
}  // namespace Amanzi

#endif

