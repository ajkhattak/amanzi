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
  virtual int H1consistency(int c, const Tensor& K, DenseMatrix& N, DenseMatrix& Ac) override {
    Errors::Message msg("H1 consistency is not implemented for generalized electromagnetics scheme.");
    Exceptions::amanzi_throw(msg);
    return 0;
  }
  virtual int StiffnessMatrix(int c, const Tensor& K, DenseMatrix& A) override;

 private:
  void CurvedFaceGeometry_(
      const Entity_ID_List& nodes, AmanziGeometry::Point& xf,
      DenseVector& dint, DenseVector& dext,
      std::vector<AmanziGeometry::Point>& tri_normal,
      std::vector<AmanziGeometry::Point>& tri_center,
      std::vector<double>& tri_area);

  void CurvedFaceRaviartThomas_(
      const AmanziGeometry::Point& xf, const Entity_ID_List& nodes,
      const std::vector<int>& dirs, const std::vector<double>& tri_area,
      const DenseVector& dint, const DenseVector& dext,
      std::vector<std::vector<AmanziGeometry::Point> >& rt_tri);

  DenseMatrix CurvedFaceLifting_(
      const AmanziGeometry::Point& xf, const Entity_ID_List& nodes,
      const std::vector<int>& dirs, const std::vector<double>& tri_area,
      const DenseVector& dint, const DenseVector& dext);

 private:
  static RegisteredFactory<MFD3D_GeneralizedElectromagnetics> factory_;
};

}  // namespace WhetStone
}  // namespace Amanzi

#endif

