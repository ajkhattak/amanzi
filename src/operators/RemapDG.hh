/*
  Operators

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  The helper advection-based base class for various remap methods. It
  provides support of time integration and calculation of various static
  and dynamic geometric quantities. The actual time-step loop could be
  implemented differently by an application.

  The integration is performed on the pseudo-time interval from 0 to 1. 
  The remap velocity u is constant, but since the integration is performed
  in the reference coordinate system associated with mesh0, the transformed
  velocity v is the time-dependent quantity. The call it the co-velocity,
  v = C^t u where C is the matrix of co-factors for the Jacobian matrix J.
  Recall that C = det(J) J^{-T}. 

  Input parameter list describes operators, limiters, and mesh maps,
  see native spec for more detail. 
*/

#ifndef AMANZI_OPERATOR_REMAP_DG_HH_
#define AMANZI_OPERATOR_REMAP_DG_HH_

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <vector>

// TPLs
#include "Epetra_MultiVector.h"
#include "Teuchos_RCP.hpp"

// Amanzi
#include "Explicit_TI_RK.hh"
#include "MeshMapsFactory.hh"
#include "VectorPolynomial.hh"
#include "WhetStoneDefs.hh"

// Amanzi::Operators
#include "LimiterCell.hh"
#include "OperatorDefs.hh"
#include "PDE_Abstract.hh"
#include "PDE_AdvectionRiemann.hh"
#include "PDE_Reaction.hh"

namespace Amanzi {

class RemapDG : public Explicit_TI::fnBase<CompositeVector> {
 public:
  RemapDG(const Teuchos::RCP<const AmanziMesh::Mesh> mesh0,
          const Teuchos::RCP<AmanziMesh::Mesh> mesh1,
          Teuchos::ParameterList& plist);
  ~RemapDG() {};

  // main members required by the time integration class
  // u is the conservative quantity
  virtual void FunctionalTimeDerivative(double t, const CompositeVector& u, CompositeVector& f) override;

  // initialization routines
  void InitializeOperators();
  void InitializeFaceVelocity();
  void InitializeJacobianMatrix();
  void InitializeConsistentJacobianDeterminant(const Amanzi::WhetStone::DG_Modal& dg); 

  // dynamic geometric quantities
  virtual void DynamicJacobianMatrix(
      int c, double t, const WhetStone::MatrixPolynomial& J, WhetStone::MatrixPolynomial& Jt);
  virtual void DynamicFaceVelocity(double t);
  virtual void DynamicCellVelocity(double t, bool consistent_det);

  // limiters
  void ApplyLimiter(const Amanzi::WhetStone::DG_Modal& dg, CompositeVector& u); 

 protected:
  Teuchos::RCP<const AmanziMesh::Mesh> mesh0_;
  Teuchos::RCP<AmanziMesh::Mesh> mesh1_;
  int ncells_owned_, ncells_wghost_, nfaces_owned_, nfaces_wghost_;
  int dim_;

  const Teuchos::ParameterList plist_;
  std::shared_ptr<WhetStone::MeshMaps> maps_;
  bool consistent_jac_;

  // operators
  int order_;
  Teuchos::RCP<Operators::PDE_Abstract> op_adv_;
  Teuchos::RCP<Operators::PDE_AdvectionRiemann> op_flux_;
  Teuchos::RCP<Operators::PDE_Reaction>op_reac_;
  
  int bc_type_;
  std::string smoothness_;

  // geometric data
  std::vector<WhetStone::VectorPolynomial> uc_;
  std::vector<WhetStone::MatrixPolynomial> J_;
  Teuchos::RCP<std::vector<WhetStone::VectorPolynomial> > jac_;

  Teuchos::RCP<std::vector<WhetStone::VectorPolynomial> > velc_;
  Teuchos::RCP<std::vector<WhetStone::Polynomial> > velf_;
  std::vector<WhetStone::VectorPolynomial> velf_vec_;

  // new determinant
  std::vector<WhetStone::Polynomial> jac0_, jac1_;

  // phisycal field to be remap
  Teuchos::RCP<CompositeVector> field_;

  // statistics
  double sharp_;
};

} // namespace Amanzi

#endif
