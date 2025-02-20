/*
  EOS
   
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon (ecoon@lanl.gov)

  EOS for a combination of air and vapor pressure. Mass density is not
  available, not because it can't be calculated, but because it depends 
  upon omega. It's not really needed, and if it were, would not fit the 
*/

#include "EOSDensityFactory.hh"
#include "EOS_DensityVaporInGas.hh"

namespace Amanzi {
namespace AmanziEOS {

/* *******************************************************************
* Constructor takes a parameter list with EOS parameters.
******************************************************************* */
EOS_DensityVaporInGas::EOS_DensityVaporInGas(Teuchos::ParameterList& eos_plist)
  : EOS_Density(eos_plist) {
  InitializeFromPlist_();
}


void EOS_DensityVaporInGas::InitializeFromPlist_()
{
  Teuchos::ParameterList gas_plist = eos_plist_.sublist("gas EOS parameters");
  EOSDensityFactory eos_factory;
  gas_eos_ = eos_factory.CreateEOS(gas_plist);
}

}  // namespace AmanziEOS
}  // namespace Amanzi
