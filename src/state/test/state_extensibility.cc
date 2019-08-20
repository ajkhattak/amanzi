/*
 State

 Tests for state as a container of data

 NOTE: this test passes if it compiles!
*/

// TPLs
#include "UnitTest++.h"

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_ParameterXMLFileReader.hpp"
#include "Teuchos_RCP.hpp"

#include "AmanziComm.hh"
#include "MeshFactory.hh"

#include "MeshFactory.hh"
#include "State.hh"
#include "errors.hh"

struct MyPoint {
  double a;
  double b;
};

using MyPointList = std::vector<MyPoint>;

bool inline UserInitialize(Teuchos::ParameterList &plist, MyPointList &t,
                           const Amanzi::Key &fieldname,
                           const std::vector<std::string> &subfieldnames) {
  std::cout << "found it!" << std::endl;
  return true;
}

void UserWriteVis(const Amanzi::Visualization &vis,
                  const Amanzi::Key &fieldname,
                  const std::vector<std::string> &subfieldnames,
                  const MyPointList &vec) {}

void UserWriteCheckpoint(const Amanzi::Checkpoint &chkp,
                         const Amanzi::Key &fieldname, const MyPointList &vec) {
}
void UserReadCheckpoint(const Amanzi::Checkpoint &chkp,
                        const Amanzi::Key &fieldname, MyPointList &vec) {}

TEST(STATE_EXTENSIBILITY_CREATION) {
  using namespace Amanzi;

  auto comm = Amanzi::getDefaultComm();
  Teuchos::ParameterList region_list;
  Teuchos::RCP<Amanzi::AmanziGeometry::GeometricModel> gm = Teuchos::rcp(
      new Amanzi::AmanziGeometry::GeometricModel(3, region_list, *comm));

  Amanzi::AmanziMesh::MeshFactory meshfactory(comm, gm);
  auto m = meshfactory.create(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 8, 1, 1);

  std::string xmlFileName = "test/state_extensibility.xml";
  Teuchos::ParameterXMLFileReader xmlreader(xmlFileName);
  auto plist = Teuchos::parameterList(xmlreader.getParameters());

  State s(Teuchos::sublist(plist, "state"));
  s.RegisterDomainMesh(m);
  s.Require<MyPointList>("my_points", "", "my_points");
  s.GetRecordW("my_points", "my_points").set_io_vis();
  s.Setup();
  s.InitializeFields();

  Visualization vis(plist->sublist("visualization"));
  vis.set_mesh(m);
  vis.CreateFiles();
  WriteVis(vis, s);

  Checkpoint chkp(plist->sublist("checkpoint"), comm);
  WriteCheckpoint(chkp, comm, s, 0.0);
}
