/*
  This is the input component of the Amanzi code. 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Erin Barker (original version)
           Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include <string>

// TPLs
#include <boost/algorithm/string.hpp>

#include "Teuchos_ParameterList.hpp"

#include "InputConverterU.hh"
#include "XMLParameterListWriter.hh"

namespace Amanzi {
namespace AmanziInput {

XERCES_CPP_NAMESPACE_USE

/* ******************************************************************
* Empty
****************************************************************** */
Teuchos::ParameterList InputConverterU::Translate(int rank, int num_proc)
{
  rank_ = rank;
  num_proc_ = num_proc;
  Teuchos::ParameterList out_list;
  
  // grab verbosity early
  verb_list_ = TranslateVerbosity_();
  Teuchos::ParameterList tmp_list(verb_list_);
  vo_ = new VerboseObject("InputConverter", tmp_list);
  Teuchos::OSTab tab = vo_->getOSTab();

  // parsing of miscalleneous lists
  ParseSolutes_();
  ParseConstants_();

  out_list.set<bool>("Native Unstructured Input", "true");

  out_list.sublist("Miscalleneous") = TranslateMisc_();  
  out_list.sublist("Mesh") = TranslateMesh_();
  out_list.sublist("Domain").set<int>("Spatial Dimension", dim_);
  out_list.sublist("Regions") = TranslateRegions_();

  const Teuchos::ParameterList& tmp = TranslateOutput_();
  for (Teuchos::ParameterList::ConstIterator it = tmp.begin(); it != tmp.end(); ++it)
    out_list.sublist(it->first) = tmp.sublist(it->first);

  out_list.sublist("State") = TranslateState_();
  out_list.sublist("Cycle Driver") = TranslateCycleDriver_();
  Teuchos::ParameterList& cd_list = out_list.sublist("Cycle Driver");
  out_list.sublist("PKs") = TranslatePKs_(cd_list);

  out_list.sublist("Solvers") = TranslateSolvers_();
  out_list.sublist("Preconditioners") = TranslatePreconditioners_();

  // analysis list
  out_list.sublist("Analysis") = CreateAnalysis_();
  FilterEmptySublists_(out_list);

  // miscalleneous cross-list information
  if (init_filename_.size() > 0) {
    out_list.sublist("State").set<std::string>("initialization filename", init_filename_);
  }

  return out_list;
}
  

/* ******************************************************************
* Extract information of solute components.
****************************************************************** */
void InputConverterU::ParseSolutes_()
{
  bool flag;
  char* tagname;
  char* text_content;

  MemoryManager mm;

  DOMNode* inode = doc_->getElementsByTagName(mm.transcode("phases"))->item(0);
  DOMNode* node = GetUniqueElementByTagsString_(inode, "liquid_phase, dissolved_components, solutes", flag);

  DOMNodeList* children = node->getChildNodes();
  for (int i = 0; i < children->getLength(); ++i) {
    inode = children->item(i);
    tagname = mm.transcode(inode->getNodeName());
    text_content = mm.transcode(inode->getTextContent());

    if (strcmp(tagname, "solute") == 0) {
      phases_["water"].push_back(TrimString_(text_content));
    }
  }
  
  comp_names_all_ = phases_["water"];
}


/* ******************************************************************
* Extract generic verbosity object for all sublists.
****************************************************************** */
Teuchos::ParameterList InputConverterU::TranslateVerbosity_()
{
  Teuchos::ParameterList vlist;

  DOMNodeList* node_list;
  DOMNode* node_attr;
  DOMNamedNodeMap* attr_map;
  char* text_content;

  MemoryManager mm;
    
  // get execution contorls node
  node_list = doc_->getElementsByTagName(mm.transcode("execution_controls"));
  
  for (int i = 0; i < node_list->getLength(); i++) {
    DOMNode* inode = node_list->item(i);

    if (DOMNode::ELEMENT_NODE == inode->getNodeType()) {
      DOMNodeList* children = inode->getChildNodes();

      for (int j = 0; j < children->getLength(); j++) {
        DOMNode* jnode = children->item(j);

        if (DOMNode::ELEMENT_NODE == jnode->getNodeType()) {
          char* tagname = mm.transcode(jnode->getNodeName());
          if (std::string(tagname) == "verbosity") {
            attr_map = jnode->getAttributes();
            node_attr = attr_map->getNamedItem(mm.transcode("level"));
            if (node_attr) {
              text_content = mm.transcode(node_attr->getNodeValue());
              vlist.sublist("VerboseObject").set<std::string>("Verbosity Level", TrimString_(text_content));
              break;
            } else {
              ThrowErrorIllformed_("verbosity", "value", "level");
            }
          }
        }
      }
    }
  }
  return vlist;
}


/* ******************************************************************
* Miscalleneous commands
****************************************************************** */
Teuchos::ParameterList InputConverterU::TranslateMisc_()
{
  Teuchos::ParameterList out_list;
  if (vo_->getVerbLevel() >= Teuchos::VERB_HIGH)
    *vo_->os() << "Translating miscalleneous commands" << std::endl;

  MemoryManager mm;  
  DOMNode* node;
  DOMElement* element;

  bool flag;
  node = GetUniqueElementByTagsString_("misc, echo_translated_input", flag);
  if (flag) {
    element = static_cast<DOMElement*>(node);
    std::string filename = GetAttributeValueS_(element, "file_name", false, "native_v6.xml");

    out_list.set<std::string>("File Name", filename);
  }

  return out_list;
}


/* ******************************************************************
* Analysis list can used by special tools.
****************************************************************** */
Teuchos::ParameterList InputConverterU::CreateAnalysis_()
{
  Teuchos::ParameterList out_list;
  out_list.set<Teuchos::Array<std::string> >("used boundary condition regions", vv_bc_regions_);
  out_list.set<Teuchos::Array<std::string> >("used source regions", vv_src_regions_);
  out_list.set<Teuchos::Array<std::string> >("used observation regions", vv_obs_regions_);

  out_list.sublist("VerboseObject") = verb_list_.sublist("VerboseObject");

  return out_list;
}


/* ******************************************************************
* Filters out empty sublists starting with node "parent".
****************************************************************** */
void InputConverterU::FilterEmptySublists_(Teuchos::ParameterList& plist)
{
  for (Teuchos::ParameterList::ConstIterator it = plist.begin(); it != plist.end(); ++it) {
    if (plist.isSublist(it->first)) {
      Teuchos::ParameterList& slist = plist.sublist(it->first);
      if (slist.numParams() == 0) {
        plist.remove(it->first);
      } else {
        FilterEmptySublists_(slist);
      }
    }
  }
}



/* ******************************************************************
* Output of XML
****************************************************************** */
void InputConverterU::SaveXMLFile(
    Teuchos::ParameterList& out_list, std::string& xmlfilename)
{
  std::string filename(xmlfilename);
  std::string new_extension("_native_v6.xml");
  size_t pos = filename.find(".xml");
  filename.replace(pos, (size_t)4, new_extension, (size_t)0, (size_t)14);

  int precision = out_list.sublist("Miscalleneous").get<int>("output precision", 0);
  if (vo_->getVerbLevel() >= Teuchos::VERB_LOW) {
    Teuchos::OSTab tab = vo_->getOSTab();
    *vo_->os() << "Writing the translated XML to " << filename.c_str() << std::endl;;
  }

  Teuchos::Amanzi_XMLParameterListWriter XMLWriter;
  if (precision > 0) XMLWriter.set_precision(precision);
  Teuchos::XMLObject XMLobj = XMLWriter.toXML(out_list);

  std::ofstream xmlfile;
  xmlfile.open(filename.c_str());
  xmlfile << XMLobj;
}

}  // namespace AmanziInput
}  // namespace Amanzi
