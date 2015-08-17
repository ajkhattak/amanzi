/*
  This is the input component of the Amanzi code. 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include <algorithm>
#include <sstream>
#include <string>

//TPLs
#include <boost/algorithm/string.hpp>
#include <xercesc/dom/DOM.hpp>

// Amanzi's
#include "errors.hh"
#include "exceptions.hh"
#include "dbc.hh"

#include "InputConverterU.hh"
#include "InputConverterU_Defs.hh"

namespace Amanzi {
namespace AmanziInput {

XERCES_CPP_NAMESPACE_USE

/* ******************************************************************
* Create operators sublist
****************************************************************** */
Teuchos::ParameterList InputConverterU::TranslateTimeIntegrator_(
    const std::string& err_options, const std::string& nonlinear_solver,
    bool modify_correction, const std::string& unstr_controls)
{
  Teuchos::ParameterList out_list;

  Teuchos::OSTab tab = vo_->getOSTab();
  if (vo_->getVerbLevel() >= Teuchos::VERB_HIGH)
    *vo_->os() << "Translating time integrator" << std::endl;

  MemoryManager mm;
  DOMNodeList* node_list;
  DOMNode* node;
  DOMElement* element;

  // error control options
  std::vector<std::string> tmp = CharToStrings_(err_options.c_str());
  out_list.set<Teuchos::Array<std::string> >("error control options", tmp);

  // linear solver
  out_list.set<std::string>("linear solver", TI_SOLVER);
  out_list.set<std::string>("preconditioner", TI_PRECONDITIONER);
  out_list.set<std::string>("preconditioner enhancement", "none");

  // pressure-lambda constraints
  Teuchos::ParameterList& plamb = out_list.sublist("pressure-lambda constraints");
  plamb.set<std::string>("method", "projection");
  plamb.set<bool>("inflow krel correction", true);
  plamb.set<std::string>("linear solver", TI_PLAMBDA_SOLVER);

  // time stepping method
  out_list.set<std::string>("time integration method", "BDF1");
  Teuchos::ParameterList& bdf1 = out_list.sublist("BDF1");

  // use standard timestep controller type
  bdf1.set<std::string>("timestep controller type", TI_TIMESTEP_CONTROLLER);
  Teuchos::ParameterList& controller = bdf1.sublist("timestep controller standard parameters");
  controller.set<int>("max iterations", TI_MAX_ITERATIONS)
      .set<int>("min iterations", TI_MIN_ITERATIONS)
      .set<double>("time step increase factor", TI_TS_INCREASE_FACTOR)
      .set<double>("time step reduction factor", TI_TS_REDUCTION_FACTOR)
      .set<double>("max time step", MAXIMUM_TIMESTEP)
      .set<double>("min time step", MINIMUM_TIMESTEP);

  // nonlinear solver
  Teuchos::ParameterList* solver;

  if (nonlinear_solver == std::string("newton") ||
      nonlinear_solver == std::string("newton-picard")) {
    bdf1.set<std::string>("solver type", "Newton");
    Teuchos::ParameterList& test = bdf1.sublist("Newton parameters");
    solver = &test;
    solver->set<double>("nonlinear tolerance", NONLINEAR_TOLERANCE);
    solver->set<double>("diverged tolerance", NKA_DIVERG_TOL);
    solver->set<double>("max du growth factor", INC_DIVERG_FACTOR);
    solver->set<int>("max divergent iterations", MAX_DIVERG_ITERATIONS);
    solver->set<int>("limit iterations", NKA_LIMIT_ITERATIONS);
    solver->set<bool>("modify correction", true);
  }
  else if (nonlinear_solver == "nka") {
    bdf1.set<std::string>("solver type", "nka");
    Teuchos::ParameterList& test = bdf1.sublist("nka parameters");
    solver = &test;
    solver->set<double>("nonlinear tolerance", NONLINEAR_TOLERANCE);
    solver->set<double>("diverged tolerance", NKA_DIVERG_TOL);
    solver->set<double>("max du growth factor", INC_DIVERG_FACTOR);
    solver->set<int>("max divergent iterations", MAX_DIVERG_ITERATIONS);
    solver->set<int>("max nka vectors", NKA_NUM_VECTORS);
    solver->set<int>("limit iterations", NKA_LIMIT_ITERATIONS);
    solver->set<bool>("modify correction", modify_correction);
  }
  else if (nonlinear_solver == std::string("jfnk")) {
    bdf1.set<std::string>("solver type", "JFNK");
    Teuchos::ParameterList& test = bdf1.sublist("JFNK parameters");
    solver = &test;
    solver->set<double>("typical solution value", 1.0);

    Teuchos::ParameterList& tmp = solver->sublist("nonlinear solver");     
    tmp.set<std::string>("solver type", "Newton");
    Teuchos::ParameterList& newton = tmp.sublist("Newton parameters");     
    newton.set<double>("diverged tolerance", NKA_DIVERG_TOL);
    newton.set<double>("max du growth factor", INC_DIVERG_FACTOR);
    newton.set<int>("max divergent iterations", MAX_DIVERG_ITERATIONS);
    newton.set<int>("max nka vectors", NKA_NUM_VECTORS);
    newton.set<int>("limit iterations", NKA_LIMIT_ITERATIONS);

    Teuchos::ParameterList& jfmat = solver->sublist("JF matrix parameters");
    jfmat.set<double>("finite difference epsilon", 1.0e-8);
    jfmat.set<std::string>("method for epsilon", "Knoll-Keyes");

    Teuchos::ParameterList& linop = solver->sublist("linear operator");
    linop.set<std::string>("iterative method", "gmres");
    Teuchos::ParameterList& gmres = linop.sublist("gmres parameters");
    gmres.set<double>("error tolerance", 1e-7);
    gmres.set<int>("maximum number of iterations", 100);
    std::vector<std::string> criteria;
    criteria.push_back("relative rhs");
    criteria.push_back("relative residual");
    gmres.set<Teuchos::Array<std::string> >("convergence criteria", criteria);
  } 
  else {
    Errors::Message msg;
    msg << "In the definition of \"unstr_nonlinear_solver\" you must specify either"
        << " 'nka', 'newton', 'jfnk', or 'newton-picard'.\n";
    Exceptions::amanzi_throw(msg);
  }

  // remaining BDF1 parameters
  bdf1.set<int>("max preconditioner lag iterations", TI_MAX_PC_LAG);
  bdf1.set<bool>("extrapolate initial guess", true);
  bdf1.set<double>("restart tolerance relaxation factor", TI_TOL_RELAX_FACTOR);
  bdf1.set<double>("restart tolerance relaxation factor damping", TI_TOL_RELAX_FACTOR_DAMPING);

  bool flag;
  node = GetUniqueElementByTagsString_(unstr_controls + ", max_iterations", flag);
  if (flag) controller.set<int>("max iterations",
      strtol(mm.transcode(node->getTextContent()), NULL, 10));

  node = GetUniqueElementByTagsString_(unstr_controls + ", min_iterations", flag);
  if (flag) controller.set<int>("min iterations",
      strtol(mm.transcode(node->getTextContent()), NULL, 10));

  node = GetUniqueElementByTagsString_(unstr_controls + ", limit_iterations", flag);
  if (flag) solver->set<int>("limit iterations",
      strtol(mm.transcode(node->getTextContent()), NULL, 10));

  node = GetUniqueElementByTagsString_(unstr_controls + ", nonlinear_tolerance", flag); 
  if (flag) solver->set<double>("nonlinear tolerance",
      strtod(mm.transcode(node->getTextContent()), NULL));

  node = GetUniqueElementByTagsString_(unstr_controls + ", time_step_reduction_factor", flag); 
  if (flag) controller.set<double>("time step reduction factor",
      strtod(mm.transcode(node->getTextContent()), NULL));

  node = GetUniqueElementByTagsString_(unstr_controls + ", time_step_increase_factor", flag); 
  if (flag) controller.set<double>("time step increase factor",
      strtod(mm.transcode(node->getTextContent()), NULL));

  node = GetUniqueElementByTagsString_(unstr_controls + ", max_preconditioner_lag_iterations", flag); 
  if (flag) bdf1.set<int>("max preconditioner lag iterations",
      strtol(mm.transcode(node->getTextContent()), NULL, 10));

  node = GetUniqueElementByTagsString_(unstr_controls + ", max_divergent_iterations", flag); 
  if (flag) solver->set<int>("max divergent iterations",
      strtol(mm.transcode(node->getTextContent()), NULL, 10));

  node = GetUniqueElementByTagsString_(unstr_controls + ", nonlinear_iteration_damping_factor", flag); 
  if (flag) bdf1.set<double>("nonlinear iteration damping factor",
      strtod(mm.transcode(node->getTextContent()), NULL));

  node = GetUniqueElementByTagsString_(
      unstr_controls + ", nonlinear_iteration_initial_guess_extrapolation_order", flag); 
  if (flag) bdf1.set<int>("nonlinear iteration initial guess extrapolation order",
      strtol(mm.transcode(node->getTextContent()), NULL, 10));

  node = GetUniqueElementByTagsString_(unstr_controls + ", restart_tolerance_relaxation_factor", flag); 
  if (flag) bdf1.set<double>("restart tolerance relaxation factor",
      strtod(mm.transcode(node->getTextContent()), NULL));

  node = GetUniqueElementByTagsString_(
      unstr_controls + ", restart_tolerance_relaxation_factor_damping", flag); 
  if (flag) bdf1.set<double>("restart tolerance relaxation factor damping",
      strtod(mm.transcode(node->getTextContent()), NULL));

  node = GetUniqueElementByTagsString_(
      unstr_controls + ", nonlinear_iteration_divergence_factor", flag); 
  if (flag) solver->set<double>("max du growth factor",
       strtod(mm.transcode(node->getTextContent()), NULL));

  node = GetUniqueElementByTagsString_(unstr_controls + ", error_control_options", flag); 
  if (flag) out_list.set<Teuchos::Array<std::string> >("error control options",
      CharToStrings_(mm.transcode(node->getTextContent())));

  node = GetUniqueElementByTagsString_(unstr_controls + ", preconditioner", flag); 
  if (flag) {
    std::string text = mm.transcode(node->getTextContent());
    if (text == "hypre_amg") text = "Hypre AMG";
    if (text == "trilinos_ml") text = "Trilinos ML";
    if (text == "block_ilu") text = "Block ILU";
    out_list.set<std::string>("preconditioner", text);
  }

  // special cases
  if (flow_single_phase_) {
    node = GetUniqueElementByTagsString_(unstr_controls + ", time_step_increase_factor", flag); 
    if (flag) controller.set<double>("time step increase factor",
        strtol(mm.transcode(node->getTextContent()), NULL, 10));
  }

  // initialization
  node = GetUniqueElementByTagsString_(unstr_controls + ", unstr_initialization", flag); 
  if (flag) {
    Teuchos::ParameterList& init = out_list.sublist("initialization");
    init = TranslateInitialization_(unstr_controls);
  }

  // overwrite parameters for special solvers
  if (nonlinear_solver == "newton" || 
      nonlinear_solver == "newton-picard") {
    bdf1.set<int>("max preconditioner lag iterations", 0);
    bdf1.set<bool>("extrapolate initial guess", false);	    
    out_list.set<std::string>("linear solver", "GMRES for Newton");
    out_list.set<std::string>("preconditioner enhancement", "GMRES for Newton");
  }

  bdf1.sublist("VerboseObject") = verb_list_.sublist("VerboseObject");
  return out_list;
}


/* ******************************************************************
* Translate initializa sublist for time integrator
****************************************************************** */
Teuchos::ParameterList InputConverterU::TranslateInitialization_(
    const std::string& unstr_controls)
{
  Teuchos::ParameterList out_list;

  MemoryManager mm;
  DOMNode* node;

  // set defaults
  out_list.set<std::string>("method", "saturated solver");
  out_list.set<std::string>("linear solver", TI_SOLVER);

  // overwite defaults using numerical controls
  bool flag;
  std::string method;
  std::string controls(unstr_controls + ", unstr_initialization");

  node = GetUniqueElementByTagsString_(controls + ", method", flag); 
  if (flag) {
    method = GetTextContentS_(node, "picard, darcy_solver");
    if (method == "darcy_solver") method = "saturated solver";
    out_list.set<std::string>("method", method);
  }

  node = GetUniqueElementByTagsString_(controls + ", clipping_saturation", flag); 
  if (flag) out_list.set<double>("clipping saturation value",
      strtod(mm.transcode(node->getTextContent()), NULL));

  node = GetUniqueElementByTagsString_(controls + ", clipping_pressure", flag); 
  if (flag) out_list.set<double>("clipping pressure value",
      strtod(mm.transcode(node->getTextContent()), NULL));

  node = GetUniqueElementByTagsString_(controls + ", linear_solver", flag); 
  if (flag) {
    std::string text = mm.transcode(node->getTextContent());
    if (text == "aztecoo") text = "AztecOO";
    out_list.set<std::string>("linear solver", text);
  }

  if (method == "picard") {
    Teuchos::ParameterList& pic_list = out_list.sublist("picard parameters");
    pic_list.set<std::string>("linear solver", PICARD_SOLVER);
    pic_list.set<double>("convergence tolerance", PICARD_TOLERANCE);
    pic_list.set<int>("maximum number of iterations", PICARD_MAX_ITERATIONS);

    node = GetUniqueElementByTagsString_(controls + ", convergence_tolerance", flag); 
    if (flag) pic_list.set<double>("convergence tolerance", 
        strtod(mm.transcode(node->getTextContent()), NULL));

    node = GetUniqueElementByTagsString_(controls + ", max_iterations", flag); 
    pic_list.set<int>("maximum number of iterations", 
        strtol(mm.transcode(node->getTextContent()), NULL, 10));
  }

  return out_list;
}


/* ******************************************************************
* Create operators sublist
****************************************************************** */
Teuchos::ParameterList InputConverterU::TranslateDiffusionOperator_(
    const std::string& disc_method, const std::string& pc_method,
    const std::string& nonlinear_solver, const std::string& extensions)
{
  Teuchos::ParameterList out_list;
  Teuchos::ParameterList tmp_list;

  std::string tmp = boost::replace_all_copy(disc_method, "-", ": ");
  replace(tmp.begin(), tmp.end(), '_', ' ');
  if (tmp == "mfd: two point flux approximation") tmp = "mfd: two-point flux approximation";

  tmp_list.set<std::string>("discretization primary", boost::to_lower_copy(tmp));
  tmp_list.set<std::string>("discretization secondary", "mfd: optimized for sparsity");

  if (disc_method != "fv: default") {
    Teuchos::Array<std::string> stensil(2);
    stensil[0] = "face";
    stensil[1] = "cell";
    tmp_list.set<Teuchos::Array<std::string> >("schema", stensil);

    if (pc_method != "linearized_operator") stensil.remove(1);
    tmp_list.set<Teuchos::Array<std::string> >("preconditioner schema", stensil);
    tmp_list.set<bool>("gravity", true);
  } else {
    Teuchos::Array<std::string> stensil(1);
    stensil[0] = "cell";
    tmp_list.set<Teuchos::Array<std::string> >("schema", stensil);

    tmp_list.set<Teuchos::Array<std::string> >("preconditioner schema", stensil);
    tmp_list.set<bool>("gravity", true);
  }

  // create two operators for matrix and preconditioner
  out_list.sublist("diffusion operator").sublist("matrix") = tmp_list;
  out_list.sublist("diffusion operator").sublist("preconditioner") = tmp_list;

  // extensions
  if (extensions == "vapor matrix") {
    Teuchos::ParameterList& vapor = out_list.sublist("diffusion operator").sublist("vapor matrix");
    vapor = tmp_list;
    vapor.set<std::string>("nonlinear coefficient", "standard: cell");
    vapor.set<bool>("exclude primary terms", false);
    vapor.set<bool>("scaled constraint equation", false);
    vapor.set<bool>("gravity", "false");
    vapor.set<std::string>("newton correction", "none");
  }

  // fixing miscalleneous scenarious
  if (pc_method == "linearized_operator") {
    out_list.sublist("diffusion operator").sublist("preconditioner")
        .set<std::string>("newton correction", "approximate jacobian");
  }

  if (nonlinear_solver == "Newton") {
    Teuchos::ParameterList& pc_list = 
        out_list.sublist("diffusion operator").sublist("preconditioner");
    pc_list.set<std::string>("newton correction", "true jacobian");

    Teuchos::ParameterList& slist = pc_list.sublist("linear operator");
    slist.set<std::string>("iterative method", "gmres");
    Teuchos::ParameterList& gmres_list = slist.sublist("gmres parameters");
    gmres_list.set<double>("error tolerance", NONLINEAR_TOLERANCE * 1e-2);
    gmres_list.set<int>("maximum number of iterations", 50);

    std::vector<std::string> criteria;
    criteria.push_back("relative rhs");
    criteria.push_back("relative residual");
    gmres_list.set<Teuchos::Array<std::string> >("convergence criteria", criteria);

    gmres_list.sublist("VerboseObject").set<std::string>("Verbosity Level", "low");
  }

  return out_list;
}

}  // namespace AmanziInput
}  // namespace Amanzi

