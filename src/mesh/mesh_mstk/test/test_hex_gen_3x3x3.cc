#include <UnitTest++.h>

#include <iostream>

#include "../Mesh_MSTK.hh"


#include "Epetra_Map.h"
#include "Epetra_MpiComm.h"

#include "MeshAudit.hh"

// Test generation of hex mesh in serial

TEST(MSTK_HEX_GEN_3x3x3)
{

  int i, j, k, err, nc, nf, nv;
  Amanzi::AmanziMesh::Set_ID faces[6], nodes[8];
  int facedirs[6];

  int NV = 64;
  int NF = 108;
  int NC = 27;

  Teuchos::RCP<Epetra_MpiComm> comm(new Epetra_MpiComm(MPI_COMM_WORLD));

  // Generate a mesh consisting of 3x3x3 elements 

  Teuchos::RCP<Amanzi::AmanziMesh::Mesh> mesh(new Amanzi::AmanziMesh::Mesh_MSTK(0.0,0.0,0.0,1.0,1.0,1.0,3,3,3,comm.get()));

  nv = mesh->num_entities(Amanzi::AmanziMesh::NODE,Amanzi::AmanziMesh::OWNED);
  CHECK_EQUAL(NV,nv);
  
  nf = mesh->num_entities(Amanzi::AmanziMesh::FACE,Amanzi::AmanziMesh::OWNED);
  CHECK_EQUAL(NF,nf);
  
  nc = mesh->num_entities(Amanzi::AmanziMesh::CELL,Amanzi::AmanziMesh::OWNED);
  CHECK_EQUAL(NC,nc);


  Amanzi::AmanziMesh::Entity_ID_List  c2f(6);
  std::vector<int> c2fdirs(6);
  Epetra_Map cell_map(mesh->cell_map(false));
  Epetra_Map face_map(mesh->face_map(false));
  for (int c=cell_map.MinLID(); c<=cell_map.MaxLID(); c++)
    {
      CHECK_EQUAL(cell_map.GID(c),mesh->GID(c,Amanzi::AmanziMesh::CELL));
      mesh->cell_get_faces_and_dirs(c, &c2f, &c2fdirs, true);
      for (int j=0; j<6; j++)
	{
	  int f = face_map.LID(c2f[j]);
	  CHECK( f == c2f[j] );
	}

    }
  
  std::ofstream outfile("test/mstk_hex_gen_3x3x3.out");
  Amanzi::MeshAudit auditor(mesh,outfile);
  auditor.Verify();


}

