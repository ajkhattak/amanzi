/*
  WhetStone, Version 2.2
  Release name: naka-to.

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  Discontinuous Galerkin modal method. Efficient implementation
  requires to cache various data for all mesh cells.
*/

#ifndef AMANZI_WHETSTONE_DG_MODAL_HH_
#define AMANZI_WHETSTONE_DG_MODAL_HH_

#include <memory>

#include "Teuchos_RCP.hpp"

// Amanzi
#include "Mesh.hh"
#include "Point.hh"

// WhetStone
#include "Basis.hh"
#include "BilinearForm.hh"
#include "BilinearFormFactory.hh"
#include "DenseMatrix.hh"
#include "DenseVector.hh"
#include "NumericalIntegration.hh"
#include "Polynomial.hh"
#include "PolynomialOnMesh.hh"
#include "Tensor.hh"
#include "VectorObjects.hh"
#include "WhetStoneDefs.hh"
#include "WhetStoneFunction.hh"

namespace Amanzi {
namespace WhetStone {

class Polynomial;

class DG_Modal : public BilinearForm {
 public:
  DG_Modal(const Teuchos::ParameterList& plist,
           const Teuchos::RCP<const AmanziMesh::Mesh>& mesh);
  ~DG_Modal() {};

  // basic member functions
  // -- schema
  virtual std::vector<SchemaItem> schema() const override;

  // -- mass matrices
  virtual int MassMatrix(int c, const Tensor& K, DenseMatrix& M) override;
  virtual int MassMatrix(int c, const VectorPolynomial& K, DenseMatrix& M) override {
    int ok;
    if (K.size() == 1) {
      ok = MassMatrixPoly_(c, K[0], M);
    } else {
      ok = MassMatrixPiecewisePoly_(c, K, M);
    }
    return ok;
  }
  int MassMatrix(int c, const Tensor& K, PolynomialOnMesh& integrals, DenseMatrix& M);

  // -- inverse mass matrices
  virtual int MassMatrixInverse(int c, const Tensor& K, DenseMatrix& W) override {
    int ok = MassMatrix(c, K, W);
    W.Inverse();
    return ok;
  }

  virtual int MassMatrixInverse(int c, const VectorPolynomial& K, DenseMatrix& W) override {
    int ok = MassMatrix(c, K, W);
    W.Inverse();
    return ok;
  }

  // -- stiffness matrices. General coefficient requires to specify the quadrature order
  virtual int StiffnessMatrix(int c, const Tensor& K, DenseMatrix& A) override;
  virtual int StiffnessMatrix(int c, const WhetStoneFunction& K, DenseMatrix& A, int order) override;

  // -- advection matrices
  virtual int AdvectionMatrix(int c, const VectorPolynomial& uc,
                              DenseMatrix& A, bool grad_on_test) override {
    int ok;
    if (uc.size() == d_) {
      ok = AdvectionMatrixPoly_(c, uc, A, grad_on_test);
    } else {
      ok = AdvectionMatrixPiecewisePoly_(c, uc, A, grad_on_test);
    }
   return ok;
  }

  // -- flux matrices
  //    returns point flux value (u.n) in the last parameter
  int FluxMatrix(int f, const Polynomial& uf, DenseMatrix& A, bool upwind, bool jump_on_test, double* flux);
  int FluxMatrixRusanov(int f, const VectorPolynomial& uc1, const VectorPolynomial& uc2,
                        const Polynomial& uf, DenseMatrix& A);
  int FluxMatrixGaussPoints(int f, const Polynomial& uf, DenseMatrix& A, bool upwind, bool jump_on_test);

  // -- interface matrices: jumps and penalty
  int FaceMatrixJump(int f, const Tensor& K1, const Tensor& K2, DenseMatrix& A);
  int FaceMatrixJump(int f, const WhetStoneFunction& K1, const WhetStoneFunction& K2, DenseMatrix& A, int order);
  int FaceMatrixPenalty(int f, double Kf, DenseMatrix& A);

  // interfaces that are not used
  virtual int L2consistency(int c, const Tensor& T, DenseMatrix& N, DenseMatrix& Mc, bool symmetry) override { return 0; }
  virtual int L2consistencyInverse(int c, const Tensor& T, DenseMatrix& R, DenseMatrix& Wc, bool symmetry) override { return 0; }
  virtual int H1consistency(int c, const Tensor& T, DenseMatrix& N, DenseMatrix& Mc) override { return 0; }

  // miscalleneous
  // -- order of polynomials in each cell
  void set_order(int order) { order_ = order; }
  int order() const { return order_; }

  // -- access
  const Basis<AmanziMesh::Mesh>& cell_basis(int c) const { return *basis_[c]; }
  Polynomial& monomial_integrals(int c) { return monomial_integrals_[c]; }

 private:
  int MassMatrixPoly_(int c, const Polynomial& K, DenseMatrix& M);
  int MassMatrixPiecewisePoly_(int c, const VectorPolynomial& K, DenseMatrix& M);

  int AdvectionMatrixPoly_(int c, const VectorPolynomial& uc, DenseMatrix& A, bool grad_on_test);
  int AdvectionMatrixPiecewisePoly_(int c, const VectorPolynomial& uc, DenseMatrix& A, bool grad_on_test);

 private:
  NumericalIntegration<AmanziMesh::Mesh> numi_;

  std::vector<Polynomial> monomial_integrals_;  // integrals of non-normalized monomials
  std::vector<std::shared_ptr<Basis<AmanziMesh::Mesh> > > basis_;

  static RegisteredFactory<DG_Modal> factory_;
};

}  // namespace WhetStone
}  // namespace Amanzi

#endif

