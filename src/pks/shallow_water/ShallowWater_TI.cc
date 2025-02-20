/*
 Shallow water PK
 
 Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL.
 Amanzi is released under the three-clause BSD License.
 The terms of use and "as is" disclaimer for this license are
 provided in the top-level COPYRIGHT file.
 
 Author: Svetlana Tokareva (tokareva@lanl.gov)
 */

#include "ShallowWater_PK.hh"

namespace Amanzi {
namespace ShallowWater {
    
void ShallowWater_PK::FunctionalTimeDerivative(double t, const Epetra_Vector& component,
                                               Epetra_Vector& f_component) {};

} // namespace ShallowWater
} // namespace Amanzi
