/*
  Operators

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  Base class for remap methods.
*/

#ifndef AMANZI_OPERATOR_REMAP_DG_TESTS_HH_
#define AMANZI_OPERATOR_REMAP_DG_TESTS_HH_

// Amanzi::Operators
#include "OperatorDefs.hh"
#include "RemapDG.hh"

#include "MeshDeformation.hh"

namespace Amanzi {

template<class AnalyticDG>
class RemapDG_Tests : public Operators::RemapDG {
 public:
  RemapDG_Tests(const Teuchos::RCP<const AmanziMesh::Mesh> mesh0,
                const Teuchos::RCP<AmanziMesh::Mesh> mesh1,
                Teuchos::ParameterList& plist) 
    : RemapDG(mesh0, mesh1, plist),
      tprint_(0.0),
      l2norm_(-1.0),
      dt_output_(0.1) {};
  ~RemapDG_Tests() {};

  // CFL condition
  double StabilityCondition();

  // output
  void CollectStatistics(double t, const CompositeVector& u);
  virtual double global_time(double t) { return t; }
  void set_dt_output(double dt) { dt_output_ = dt; }

 protected:
  std::vector<WhetStone::Polynomial> det0_, det1_;

  // statistics
  double tprint_, dt_output_, l2norm_;
};


/* *****************************************************************
* Initialization of the consistent jacobian determinant
***************************************************************** */
template<class AnalyticDG>
double RemapDG_Tests<AnalyticDG>::StabilityCondition()
{
  double dt(1e+99), alpha(0.2), tmp;

  for (int f = 0; f < nfaces_wghost_; ++f) {
    double area = mesh0_->face_area(f);
    const AmanziGeometry::Point& xf = mesh0_->face_centroid(f);
    velf_vec_[f].Value(xf).Norm2(&tmp);
    dt = std::min(dt, area / tmp);
  }

  return dt * alpha / (2 * order_ + 1);
}


/* *****************************************************************
* Print statistics using conservative field u
***************************************************************** */
template<class AnalyticDG>
void RemapDG_Tests<AnalyticDG>::CollectStatistics(double t, const CompositeVector& u)
{
  double tglob = global_time(t);
  if (tglob >= tprint_) {
    op_reac_->UpdateMatrices(t);
    auto& matrices = op_reac_->local_op()->matrices;
    for (int n = 0; n < matrices.size(); ++n) matrices[n].Inverse();

    auto& rhs = *op_reac_->global_operator()->rhs();
    op_reac_->global_operator()->Apply(u, rhs);
    rhs.Dot(u, &l2norm_);

    Epetra_MultiVector& xc = *rhs.ViewComponent("cell");
    int nk = xc.NumVectors();
    double xmax[nk], xmin[nk], lmax(-1.0), lmin(-1.0), lavg(-1.0);
    xc.MaxValue(xmax);
    xc.MinValue(xmin);

    if (limiter() != Teuchos::null) {
      const auto& lim = *limiter()->limiter();
      lim.MaxValue(&lmax);
      lim.MinValue(&lmin);
      lim.MeanValue(&lavg);
    }

    if (mesh0_->get_comm()->MyPID() == 0) {
      printf("t=%8.5f  L2=%9.5g  nfnc=%5d  sharp=%5.1f%%  limiter: %6.3f %6.3f %6.3f  umax/umin: %9.5g %9.5g\n",
             tglob, l2norm_, nfun_, sharp_, lmax, lmin, lavg, xmax[0], xmin[0]);
    }

    tprint_ += dt_output_;
    sharp_ = 0.0;
  } 
}

} // namespace Amanzi

#endif
