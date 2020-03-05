/*
  This is the energy component of the Amanzi code. 

  Copyright 2010-2012 held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include <cmath>
#include <iostream>
#include <string>
#include <vector>

// TPLs
#include "AmanziComm.hh"
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "UnitTest++.h"

// Amanzi
#include "GMVMesh.hh"
#include "LeastSquare.hh"
#include "Mesh.hh"
#include "MeshFactory.hh"
#include "secondary_variable_field_evaluator.hh"
#include "State.hh"

// Energy
#include "Analytic01.hh"
#include "EnergyOnePhase_PK.hh"

using namespace Amanzi;
using namespace Amanzi::AmanziMesh;
using namespace Amanzi::AmanziGeometry;
using namespace Amanzi::Energy;

namespace Amanzi {

class TestEnthalpyEvaluator : public SecondaryVariableFieldEvaluator {
 public:
  explicit TestEnthalpyEvaluator(Teuchos::ParameterList& plist) :
      SecondaryVariableFieldEvaluator(plist) {
    my_key_ = "enthalpy";
    temperature_key_ = "temperature";
    dependencies_.insert(temperature_key_);
  };
  TestEnthalpyEvaluator(const TestEnthalpyEvaluator& other) :
     SecondaryVariableFieldEvaluator(other) {};

  virtual Teuchos::RCP<FieldEvaluator> Clone() const {
    return Teuchos::rcp(new TestEnthalpyEvaluator(*this));
  }

  // Required methods from SecondaryVariableFieldEvaluator
  virtual void EvaluateField_(
          const Teuchos::Ptr<State>& S,
          const Teuchos::Ptr<CompositeVector>& result) {
    const Epetra_MultiVector& temp_c = *S->GetFieldData("temperature")->ViewComponent("cell");
    Epetra_MultiVector& result_c = *result->ViewComponent("cell");

    int ncomp = result->size("cell", false);
    for (int i = 0; i != ncomp; ++i) {
      result_c[0][i] = std::pow(temp_c[0][i], 3.0);
    }
  }

  virtual void EvaluateFieldPartialDerivative_(
          const Teuchos::Ptr<State>& S,
          Key wrt_key, const Teuchos::Ptr<CompositeVector>& result) {
    const Epetra_MultiVector& temp_c = *S->GetFieldData("temperature")->ViewComponent("cell");
    Epetra_MultiVector& result_c = *result->ViewComponent("cell");

    int ncomp = result->size("cell", false);
    for (int i = 0; i != ncomp; ++i) {
      result_c[0][i] = 3.0 * std::pow(temp_c[0][i], 2.0);
    }
  }

 protected:
  Key temperature_key_;
};

}  // namespace Amanzi


TEST(ENERGY_CONVERGENCE) {
  auto comm = Amanzi::getDefaultComm();
  int MyPID = comm->MyPID();
  if (MyPID == 0) std::cout <<"Convergence analysis on three random meshes" << std::endl;

  std::string xmlFileName = "test/energy_convergence.xml";
  Teuchos::RCP<Teuchos::ParameterList> plist = Teuchos::getParametersFromXmlFile(xmlFileName);

  // convergence estimate: Use n=3 for the full test.
  int nmeshes = plist->get<int>("number of meshes", 1);
  std::vector<double> h, error;

  int nx(20);
  double dt(0.02);
  for (int n = 0; n < nmeshes; n++, nx *= 2) {
    Teuchos::ParameterList region_list = plist->get<Teuchos::ParameterList>("regions");
    Teuchos::RCP<Amanzi::AmanziGeometry::GeometricModel> gm =
        Teuchos::rcp(new Amanzi::AmanziGeometry::GeometricModel(2, region_list, *comm));
    
    Preference pref;
    pref.clear();
    pref.push_back(Framework::MSTK);

    MeshFactory meshfactory(comm,gm);
    meshfactory.set_preference(pref);
    Teuchos::RCP<const Mesh> mesh;
    if (n == 0) {
      mesh = meshfactory.create(1.0, 0.0, 2.0, 1.0, 20, 10);
      // mesh = meshfactory.create("test/random_mesh1.exo");
    } else if (n == 1) {
      mesh = meshfactory.create(1.0, 0.0, 2.0, 1.0, 40, 10);
      // mesh = meshfactory.create("test/random_mesh2.exo");
    } else if (n == 2) {
      mesh = meshfactory.create(1.0, 0.0, 2.0, 1.0, 80, 10);
      // mesh = meshfactory.create("test/random_mesh3.exo");
    }

    // create a simple state and populate it
    Teuchos::ParameterList state_list = plist->get<Teuchos::ParameterList>("state");
    Teuchos::RCP<State> S = Teuchos::rcp(new State(state_list));
    S->RegisterDomainMesh(Teuchos::rcp_const_cast<Mesh>(mesh));

    Teuchos::ParameterList pk_tree = plist->sublist("PKs").sublist("energy");
    Teuchos::RCP<TreeVector> soln = Teuchos::rcp(new TreeVector());
    Teuchos::RCP<EnergyOnePhase_PK> EPK = Teuchos::rcp(new EnergyOnePhase_PK(pk_tree, plist, S, soln));

    // overwrite enthalpy with a different model
    Teuchos::ParameterList ev_list;
    S->RequireField("enthalpy")->SetMesh(mesh)->SetGhosted()->AddComponent("cell", AmanziMesh::CELL, 1);
    auto enthalpy = Teuchos::rcp(new TestEnthalpyEvaluator(ev_list));
    S->SetFieldEvaluator("enthalpy", enthalpy);

    EPK->Setup(S.ptr());
    S->Setup();
    S->InitializeFields();
    S->InitializeEvaluators();

    EPK->Initialize(S.ptr());
    S->CheckAllFieldsInitialized();

    // constant time stepping 
    int itrs(0);
    double t(0.0), t1(0.5), dt_next;
    while (t < t1) {
      // swap conserved quntity (no backup, we check dt_next instead)
      const CompositeVector& e = *S->GetFieldData("energy");
      CompositeVector& e_prev = *S->GetFieldData("prev_energy", "thermal");
      e_prev = e;

      if (itrs == 0) {
        Teuchos::RCP<TreeVector> udot = Teuchos::rcp(new TreeVector(*soln));
        udot->PutScalar(0.0);
        EPK->bdf1_dae()->SetInitialState(t, soln, udot);
        EPK->UpdatePreconditioner(t, soln, dt);
      }

      EPK->bdf1_dae()->TimeStep(dt, dt_next, soln);
      CHECK(dt_next >= dt);
      EPK->bdf1_dae()->CommitSolution(dt, soln);
      EPK->temperature_eval()->SetFieldAsChanged(S.ptr());

      t += dt;
      itrs++;
    }

    EPK->CommitStep(0.0, 1.0, S);

    // calculate errors
    Teuchos::RCP<const CompositeVector> temp = S->GetFieldData("temperature");
    Analytic01 ana(temp, mesh);

    double l2_norm, l2_err, inf_err;  // error checks
    ana.ComputeCellError(*temp->ViewComponent("cell"), t1, l2_norm, l2_err, inf_err);

    h.push_back(1.0 / nx);
    error.push_back(l2_err);

    printf("mesh=%d bdf1_steps=%d  L2_temp_err=%7.3e L2_temp=%7.3e\n", n, itrs, l2_err, l2_norm);
    CHECK(l2_err < 8e-1);

    // save solution
    GMV::open_data_file(*mesh, (std::string)"energy.gmv");
    GMV::start_data();
    GMV::write_cell_data(*temp->ViewComponent("cell"), 0, "temperature");
    GMV::close_data_file();
  }

  // check convergence rate
  double l2_rate = Amanzi::Utils::bestLSfit(h, error);
  printf("convergence rate: %10.2f\n", l2_rate);
  CHECK(l2_rate > 0.84);
}


TEST(ENERGY_PRECONDITIONER) {
  Comm_ptr_type comm = Amanzi::getDefaultComm();
  int MyPID = comm->MyPID();
  if (MyPID == 0) std::cout << "\nImpact of advection term in the preconditioner" << std::endl;

  std::string xmlFileName = "test/energy_convergence.xml";
  Teuchos::RCP<Teuchos::ParameterList> plist = Teuchos::getParametersFromXmlFile(xmlFileName);

  // preconditioner with (loop=0) and without (loop=1) enthalpy term.
  int num_itrs[2];
  for (int loop = 0; loop < 2; loop++) {
    Teuchos::ParameterList region_list = plist->get<Teuchos::ParameterList>("regions");
    auto gm = Teuchos::rcp(new Amanzi::AmanziGeometry::GeometricModel(2, region_list, *comm));
    
    Preference pref;
    pref.clear();
    pref.push_back(Framework::MSTK);

    MeshFactory meshfactory(comm,gm);
    meshfactory.set_preference(pref);
    Teuchos::RCP<const Mesh> mesh;
    mesh = meshfactory.create(1.0, 0.0, 2.0, 1.0, 30, 30);
    // mesh = meshfactory.create("test/random_mesh1.exo");

    // create a simple state and populate it
    Teuchos::ParameterList state_list = plist->get<Teuchos::ParameterList>("state");
    Teuchos::RCP<State> S = Teuchos::rcp(new State(state_list));
    S->RegisterDomainMesh(Teuchos::rcp_const_cast<Mesh>(mesh));

    Teuchos::ParameterList pk_tree = plist->sublist("PKs").sublist("energy");
    Teuchos::RCP<TreeVector> soln = Teuchos::rcp(new TreeVector());
    Teuchos::RCP<EnergyOnePhase_PK> EPK = Teuchos::rcp(new EnergyOnePhase_PK(pk_tree, plist, S, soln));

    // overwrite enthalpy with a different model
    Teuchos::ParameterList ev_list;
    S->RequireField("enthalpy")->SetMesh(mesh)->SetGhosted()->AddComponent("cell", AmanziMesh::CELL, 1);
    auto enthalpy = Teuchos::rcp(new TestEnthalpyEvaluator(ev_list));
    S->SetFieldEvaluator("enthalpy", enthalpy);

    EPK->Setup(S.ptr());
    S->Setup();
    S->InitializeFields();
    S->InitializeEvaluators();

    EPK->Initialize(S.ptr());
    S->CheckAllFieldsInitialized();

    // constant time stepping 
    int itrs(0);
    double t(0.0), t1(0.5), dt(0.02), dt_next;
    while (t < t1) {
      // swap conserved quntity (no backup, we check dt_next instead)
      const CompositeVector& e = *S->GetFieldData("energy");
      CompositeVector& e_prev = *S->GetFieldData("prev_energy", "thermal");
      e_prev = e;

      if (itrs == 0) {
        Teuchos::RCP<TreeVector> udot = Teuchos::rcp(new TreeVector(*soln));
        udot->PutScalar(0.0);
        EPK->bdf1_dae()->SetInitialState(t, soln, udot);
        EPK->UpdatePreconditioner(t, soln, dt);
      }

      EPK->bdf1_dae()->TimeStep(dt, dt_next, soln);
      CHECK(dt_next >= dt);
      EPK->bdf1_dae()->CommitSolution(dt, soln);
      EPK->temperature_eval()->SetFieldAsChanged(S.ptr());

      t += dt;
      itrs++;
    }

    EPK->CommitStep(0.0, 1.0, S);
    num_itrs[loop] = EPK->bdf1_dae()->number_nonlinear_steps();
    printf("number of nonlinear steps: %d\n", num_itrs[loop]);
    plist->sublist("PKs").sublist("energy").sublist("operators")
        .set<bool>("include enthalpy in preconditioner", false);
  }
  CHECK(num_itrs[1] > num_itrs[0]);
}

