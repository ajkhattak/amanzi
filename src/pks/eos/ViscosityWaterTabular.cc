/*
  EOS
   
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  EOS for liquid water for temperature between 0.5 and 800 C and 
  pressure between 634 Pa and 110 MPa
*/

#include "ViscosityWaterTabular.hh"

namespace Amanzi {
namespace AmanziEOS {

ViscosityWaterTabular::ViscosityWaterTabular(Teuchos::ParameterList& eos_plist)
  : ViscosityBase(eos_plist)
{
  table_ = Teuchos::rcp(new LookupTable(eos_plist_));
}

}  // namespace AmanziEOS
}  // namespace Amanzi
