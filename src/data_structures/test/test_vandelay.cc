/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
   Amanzi

   License: see $AMANZI_DIR/COPYRIGHT
   Author: Ethan Coon

   Unit tests for the composite vector.
   ------------------------------------------------------------------------- */

#include <vector>

#include "UnitTest++.h"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_RCP.hpp"
#include "Epetra_Vector.h"

// #include "Mesh_MSTK.hh"
#include "MeshFactory.hh"
#include "CompositeVector.hh"

using namespace Amanzi;
using namespace Amanzi::AmanziMesh;

struct test_cv_vandelay {
  Comm_ptr_type comm;
  Teuchos::RCP<Mesh> mesh;
  Teuchos::RCP<CompositeVectorSpace> x_space;
  Teuchos::RCP<CompositeVector> x;

  test_cv_vandelay() {
    comm = getDefaultComm();
    Preference pref;
    pref.clear();
    pref.push_back(Framework::MSTK);

    MeshFactory meshfactory(comm);
    meshfactory.set_preference(pref);

    mesh = meshfactory.create(0.0, 0.0, 0.0, 4.0, 4.0, 4.0, 2, 2, 2);
    // mesh = Teuchos::rcp(new Mesh_MSTK(0.0, 0.0, 0.0, 4.0, 4.0, 4.0, 2, 2, 2, comm, NULL));

    std::vector<Entity_kind> locations(2);
    locations[0] = CELL;
    locations[1] = FACE;

    std::vector<Entity_kind> locations_v(2);
    locations_v[0] = CELL;
    locations_v[1] = BOUNDARY_FACE;

    std::vector<std::string> names(2);
    names[0] = "cell";
    names[1] = "face";

    std::vector<std::string> names_v(2);
    names_v[0] = "cell";
    names_v[1] = "boundary_face";

    std::vector<int> num_dofs(2);
    num_dofs[0] = 2;
    num_dofs[1] = 1;

    x_space = Teuchos::rcp(new CompositeVectorSpace());
    x_space->SetMesh(mesh)->SetGhosted()
        ->SetComponents(names, locations, num_dofs);
    x = Teuchos::rcp(new CompositeVector(*x_space));
  }
  ~test_cv_vandelay() {  }
};


SUITE(VANDELAY_COMPOSITE_VECTOR) {
  // test the vector's putscalar
  TEST_FIXTURE(test_cv_vandelay, CVVandelay) {
    std::cout << "X has " << x->NumComponents() << " components" << std::endl;
    x->PutScalar(2.0);

    {
      auto v_c = x->ViewComponent<AmanziDefaultHost>("cell", false);
      CHECK_CLOSE(2., v_c(0,0), 0.00001);
      CHECK_CLOSE(2., v_c(0,1), 0.00001);

      auto v_f = x->ViewComponent<AmanziDefaultHost>("face", 0, false);
      CHECK_CLOSE(2.0, v_f(0), 0.00001);
    }

    {
      auto v_bf = x->ViewComponent<AmanziDefaultHost>("boundary_face", 0, false);
      CHECK_CLOSE(2.0, v_bf(0), 0.00001);
    }
  }
}

