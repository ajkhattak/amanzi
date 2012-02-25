#include <cstdlib>
#include <cmath>
#include <iostream>
#include <vector>

#include "UnitTest++.h"

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"

#include "Mesh_simple.hh"
#include "MeshAudit.hh"
#include "Point.hh"

#include "State.hpp"
#include "Transport_PK.hpp"


double f_step(const Amanzi::AmanziGeometry::Point& x, double t) { 
  if (x[0] <= 1 + t) return 1;
  return 0;
}

double f_smooth(const Amanzi::AmanziGeometry::Point& x, double t) { 
  return 0.5 - atan(50*(x[0]-5-t)) / M_PI;
}

double f_cubic(const Amanzi::AmanziGeometry::Point& x, double t) {
  if( x[0] < 1 + t ) return 1;
  if( x[0] > 3 + t ) return 0;
  double z = (x[0]-1-t) / 2;
  return 2*z*z*z - 3*z*z + 1;
}


/* **************************************************************** */
TEST(CONVERGENCE_ANALYSIS_DONOR) {
  using namespace Teuchos;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziTransport;
  using namespace Amanzi::AmanziGeometry;

  std::cout << "TEST: convergence analysis of the donor scheme" << endl;
  Epetra_SerialComm  *comm = new Epetra_SerialComm();

  // read parameter list
  ParameterList parameter_list;
  string xmlFileName = "test/transport_convergence.xml";
  updateParametersFromXmlFile(xmlFileName, &parameter_list);

  for (int nx=20; nx<321; nx*=2 ) {
    // create an MSTK mesh framework 
    ParameterList region_list = parameter_list.get<Teuchos::ParameterList>("Regions");
    GeometricModelPtr gm = new GeometricModel(3, region_list, (Epetra_MpiComm *)comm);
    RCP<Mesh> mesh = rcp(new Mesh_simple(0.0,0.0,0.0, 5.0,1.0,1.0, nx, 2, 2, (const Epetra_MpiComm *)comm, gm));

    // create a transport state with one component 
    int num_components = 1;
    State mpc_state(num_components, mesh);
    RCP<Transport_State> TS = rcp(new Transport_State(mpc_state));
 
    Point u(1.0, 0.0, 0.0);
    TS->analytic_darcy_flux(u);
    TS->analytic_total_component_concentration(f_cubic);
    TS->analytic_porosity(1.0);
    TS->analytic_water_saturation(1.0);
    TS->analytic_water_density(1.0);

    ParameterList transport_list =  parameter_list.get<Teuchos::ParameterList>("Transport");
    Transport_PK TPK(transport_list, TS);
    TPK.InitPK();
    TPK.set_standalone_mode(true);
    TPK.spatial_disc_order = TPK.temporal_disc_order = 1;
    if (nx == 20) TPK.print_statistics();
    TPK.verbosity_level = 0;
 
    // advance the state
    int i, k, iter = 0;
    double T = 0.0, T1 = 1.0;

    RCP<Transport_State> TS_next = TPK.get_transport_state_next();
    RCP<Epetra_MultiVector> tcc = TS->get_total_component_concentration();
    RCP<Epetra_MultiVector> tcc_next = TS_next->get_total_component_concentration();

    while (T < T1) {
      double dT = std::min(TPK.calculate_transport_dT(), T1 - T);
      TPK.advance(dT);
      T += dT;

      *tcc = *tcc_next;
      iter++;
    }

    // calculate L1 and L2 errors
    double L1, L2;
    TS->error_total_component_concentration(f_cubic, T, &L1, &L2);
    printf("nx=%3d  L1 error=%7.5f  L2 error=%7.5f  dT=%7.4f\n", nx, L1, L2, T1 / iter);

    delete gm;
  }

  delete comm;
}


/* **************************************************************** */
TEST(CONVERGENCE_ANALYSIS_2ND) {
  using namespace Teuchos;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziTransport;
  using namespace Amanzi::AmanziGeometry;

  cout << "Test: Convergence analysis, 2nd order scheme" << endl;
  Epetra_SerialComm  *comm = new Epetra_SerialComm();

  /* read parameter list */
  ParameterList parameter_list;
  string xmlFileName = "test/transport_convergence.xml";
  updateParametersFromXmlFile(xmlFileName, &parameter_list);
  
  /* create an MSTK mesh framework */
  ParameterList region_list = parameter_list.get<Teuchos::ParameterList>("Regions");
  GeometricModelPtr gm = new GeometricModel(3, region_list, (Epetra_MpiComm *)comm);
 
  for (int nx=10; nx<81; nx*=2 ) {
    RCP<Mesh> mesh = rcp(new Mesh_simple(0.0,0.0,0.0, 5.0,1.0,1.0, nx, 2, 1, (const Epetra_MpiComm *)comm, gm)); 

    // create a transport states with one component
    int num_components = 1;
    State mpc_state(num_components, mesh);
    RCP<Transport_State> TS = rcp(new Transport_State(mpc_state));

    Point u(1.0, 0.0, 0.0);
    TS->analytic_darcy_flux(u);
    TS->analytic_total_component_concentration(f_cubic);
    TS->analytic_porosity(1.0);
    TS->analytic_water_saturation(1.0);
    TS->analytic_water_density(1.0);

    ParameterList transport_list =  parameter_list.get<Teuchos::ParameterList>("Transport");
    Transport_PK TPK(transport_list, TS);
    TPK.InitPK();
    if (nx == 10) TPK.print_statistics();
    TPK.verbosity_level = 0;
    TPK.spatial_disc_order = TPK.temporal_disc_order = 2;

    // advance the state
    int i, k, iter = 0;
    double T = 0.0, T1 = 2.0;

    RCP<Transport_State> TS_next = TPK.get_transport_state_next();
    RCP<Epetra_MultiVector> tcc = TS->get_total_component_concentration();
    RCP<Epetra_MultiVector> tcc_next = TS_next->get_total_component_concentration();

    double dT, dT0;
    if (nx==10) dT0 = TPK.calculate_transport_dT();
    else dT0 /= 2;

    while (T < T1) {
      dT = std::min(TPK.calculate_transport_dT(), T1 - T);
      dT = std::min(dT, dT0);

      TPK.advance(dT);
      T += dT;
      if (TPK.internal_tests) {
        TPK.check_tracer_bounds(*tcc_next, 0, 0.0, 1.0, 1e-12);
      }

      *tcc = *tcc_next;
      iter++;
    }
    //for (int k=0; k<nx; k++) cout << (*tcc_next)[0][k] << endl;

    double L1, L2;  // L1 and L2 errors
    TS->error_total_component_concentration(f_cubic, T, &L1, &L2);
    printf("nx=%3d  L1 error=%10.8f  L2 error=%10.8f  dT=%7.4f\n", nx, L1, L2, T1 / iter);
  }

  delete gm;
  delete comm;
}



