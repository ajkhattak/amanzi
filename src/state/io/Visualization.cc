/*
  State

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Markus Berndt
           Ethan Coon (ecoon@lanl.gov)

  Visualization of data.
*/

// TPLs
#include "Teuchos_VerboseObjectParameterListHelpers.hpp"

#include "dbc.hh"
#include "errors.hh"
#include "exceptions.hh"

#include "OutputXDMF.hh"
#if ENABLE_Silo
#include "OutputSilo.hh"
#endif

#include "Visualization.hh"

namespace Amanzi {

// -----------------------------------------------------------------------------
// Constructor
// -----------------------------------------------------------------------------
Visualization::Visualization(Teuchos::ParameterList &plist) : IOEvent(plist) {
  // ReadParameters_();

  // // set the line prefix for output
  // this->setLinePrefix("Amanzi::Visualization  ");
  // // make sure that the line prefix is printed
  // this->getOStream()->setShowLinePrefix(true);

  // // Read the sublist for verbosity settings.
  // Teuchos::readVerboseObjectSublist(&plist_, this);
}

// -----------------------------------------------------------------------------
// Constructor for a disabled Vis.
// -----------------------------------------------------------------------------
Visualization::Visualization() : IOEvent() {}

// -----------------------------------------------------------------------------
// Set up control from parameter list.
// -----------------------------------------------------------------------------
void Visualization::ReadParameters_() {
  // Teuchos::Array<std::string> no_regions(0);
  // Teuchos::ParameterList &tmp = plist_.sublist("write regions");

  // regions_.clear();
  // for (Teuchos::ParameterList::ConstIterator it = tmp.begin(); it != tmp.end();
  //      ++it) {
  //   regions_[it->first] =
  //       tmp.get<Teuchos::Array<std::string>>(it->first, no_regions);
  // }
  // write_partition_ = plist_.get<bool>("write partitions", false);

  // dynamic_mesh_ = plist_.get<bool>("dynamic mesh", false);
}

// -----------------------------------------------------------------------------
// Write a field with region information
// -----------------------------------------------------------------------------
void Visualization::WriteRegions() {
  // if (regions_.size() > 0) {
  //   for (std::map<std::string, Teuchos::Array<std::string>>::const_iterator it =
  //            regions_.begin();
  //        it != regions_.end(); ++it) {
  //     // first make an Epetra_Vector to hold the region information
  //     Epetra_Vector reg(mesh_->cell_map(false), true);

  //     // loop over the regions and initialize the reg array
  //     double reg_index = 1.0;
  //     for (Teuchos::Array<std::string>::const_iterator
  //              reg_it = (it->second).begin();
  //          reg_it != (it->second).end(); ++reg_it, reg_index += 1.0) {
  //       // only do something if the user provided a valid region name
  //       // for a region that consists of cells
  //       if (mesh_->valid_set_name(*reg_it, AmanziMesh::CELL)) {
  //         AmanziMesh::Entity_ID_List ids;
  //         mesh_->get_set_entities(*reg_it, AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED,
  //                                 &ids);

  //         for (AmanziMesh::Entity_ID_List::const_iterator rit = ids.begin();
  //              rit != ids.end(); ++rit) {
  //           reg[*rit] = reg_index;
  //         }
  //       }
  //     }

  //     Write(it->first, reg);
  //   }
  // }
}

// -----------------------------------------------------------------------------
// Write a field with region information
// -----------------------------------------------------------------------------
void Visualization::WritePartition() {
  // if (write_partition_) {
  //   // first make an Epetra_Vector to hold the partitions information
  //   Epetra_Vector reg(mesh_->cell_map(false), false);
  //   // loop over the regions and initialize the reg array
  //   double part_index = static_cast<double>(mesh_->get_comm()->MyPID());
  //   reg.putScalar(part_index);

  //   Write("partition", reg);
  // }
}

// -----------------------------------------------------------------------------
// Writing to files
// -----------------------------------------------------------------------------
void Visualization::CreateFiles() {
//   AMANZI_ASSERT(mesh_ != Teuchos::null);

//   std::string file_format = plist_.get<std::string>("file format", "XDMF");

//   if (file_format == "XDMF" || file_format == "xdmf") {
//     //visualization_output_ = Teuchos::rcp(new OutputXDMF(plist_, mesh_, true, dynamic_mesh_));
// #if ENABLE_Silo
//   } else if (file_format == "Silo" || file_format == "SILO" ||
//              file_format == "silo") {
//     //visualization_output_ = Teuchos::rcp(new OutputSilo(plist_, mesh_, true, dynamic_mesh_));
// #endif
//   } else {
//     Errors::Message msg("Visualization: Unknown file format: \"" + file_format +
//                         "\"");
//     throw(msg);
//   }
}

void Visualization::CreateTimestep(const double &time, const int &cycle) {
  //visualization_output_->InitializeCycle(time, cycle);
}

void Visualization::FinalizeTimestep() const {
  //visualization_output_->FinalizeCycle();
}

} // namespace Amanzi
