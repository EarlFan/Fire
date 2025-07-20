//       _____   _               
//      |  ___| (_)  _ __    ___ 
//      | |_    | | | '__|  / _ \
//      |  _|   | | | |    |  __/
//      |_|     |_| |_|    \___|
//
//  This file is part of Fire.
//
//  Fire is the legal property of its developers, whose names 
//  are listed in the copyright file included with this source 
//  distribution.
//
//  Fire is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published 
//  by the Free Software Foundation, either version 3 of the License, 
//  or (at your option) any later version.
//  
//  Fire is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
//  GNU General Public License for more details.
//  
//  You should have received a copy of the GNU General Public License
//  along with Fire (file LICENSE).  
//  If not, see <http://www.gnu.org/licenses/>.

//! \file      Input.cpp
//! \author    E Fan, Lisong Shi, Jiaao Hao, Tianhan Zhang, Chih-yung Wen.
//! \version   1.0
//! \date      Jan 01 2025

#include "../libs/ECOGEN/InputOutput/Input.h"
#include "../libs/ECOGEN/InputOutput/HeaderInputOutput.h"
#include "../libs/ECOGEN/Meshes/stretchZone.h"

using namespace tinyxml2;

//***********************************************************************

void Input::entreeModel(std::string casTest)
{
  try{
    //1) Parsing du file XML par la bibliotheque tinyxml2
    //------------------------------------------------------
    std::stringstream fileName(casTest + m_nameModel);
    XMLDocument xmlModel;
    XMLError error(xmlModel.LoadFile(fileName.str().c_str())); //Le file est parse ici
    if (error != XML_SUCCESS) throw ErrorXML(fileName.str(), __FILE__, __LINE__);

    //2) Recuperation des donnees du Model
    //-------------------------------------
    //Recuperation racine du document XML
    XMLNode *xmlNode = xmlModel.FirstChildElement("model");
    if (xmlNode == NULL) throw ErrorXMLRacine("model", fileName.str(), __FILE__, __LINE__);
    
    m_run->m_numberPhases = NS;
    m_run->m_model = new ModFire(m_run->m_numberTransports, m_run->m_numberPhases);
    

    //Thermodynamique
    std::vector<std::string> nameEOS;
    nameEOS.push_back("SG_waterLiq.xml");
    int EOSTrouvee(1);
    {
      m_run->m_numberEos = m_run->m_numberPhases;
      m_run->m_eos = new Eos*[m_run->m_numberEos];
      int numberEos(0);
      for (int i = 0; i < m_run->m_numberEos; i++){ m_run->m_eos[i] = entreeEOS(i, numberEos); }
    }

    //chemical input reading
    XMLElement *element(xmlNode->FirstChildElement("CanteraInput"));
    if (element != NULL) {
      std::string input_file(element->Attribute("input_file"));
      std::string name(element->Attribute("name"));
      std::string transport_type(element->Attribute("transport_type"));
      if(rankCpu==0) 
      {
        std::cout<<"The Cantera mech file is:  "<<input_file<<std::endl;
      }
      m_run->m_sol = Cantera::newSolution(input_file, name, transport_type);
      solution_mech = Cantera::newSolution(input_file, name, transport_type);
      std::shared_ptr<Cantera::ThermoPhase> gas = solution_mech->thermo();
      readMW();

      if (AMRPara::YH2O == 1)
      {
        AMRPara::YH2O = gas->speciesIndex("H2O");
      }

      if (AMRPara::YCO2 == 1)
      {
        AMRPara::YCO2 = gas->speciesIndex("CO2");
      }
      if (AMRPara::YOH == 1)
      {
        AMRPara::YOH = gas->speciesIndex("OH");
      }
    }

    //Reading of the additional physics
    int numberGPA(0);
    int physiqueAddTrouvee(0);
    element = xmlNode->FirstChildElement("additionalPhysic");
    while (element != NULL) {
      physiqueAddTrouvee++;
      //Lecture physique additionnelle
      std::string typeAddPhys(element->Attribute("type"));
      if (typeAddPhys == "") throw ErrorXMLAttribut("type", fileName.str(), __FILE__, __LINE__);
      Tools::uppercase(typeAddPhys);

      if (typeAddPhys == "VISCOSITY") {

        m_run->m_addPhys.push_back(new APFireDiffusion(numberGPA, fileName.str())); 
      }
      else { throw ErrorXMLDev(fileName.str(), __FILE__, __LINE__); }
      element = element->NextSiblingElement("additionalPhysic");
    }
    m_run->m_numberAddPhys = physiqueAddTrouvee;

    //Symmetry terms reading
    m_run->m_symmetryAddPhys = new Symmetry();
    element = xmlNode->FirstChildElement("symmetryTerm");
    if (element == NULL) {
      m_run->m_symm="NOSYMM";
      m_run->m_symmetry = new Symmetry();
    }
    else {
      //Symmetry reading
      std::string typeSym(element->Attribute("type"));
      if (typeSym == "") throw ErrorXMLAttribut("type", fileName.str(), __FILE__, __LINE__);
      Tools::uppercase(typeSym);
      m_run->m_mesh->setSymType(typeSym);
      m_run->m_model->setSymType(typeSym);
      //Switch on the type of symmetry
      if (typeSym == "CYLINDRICAL") { m_run->m_symm="CYLINDRICAL"; m_run->m_symmetry = new SymCylindrical(element, fileName.str()); }
      else if (typeSym == "SPHERICAL") { m_run->m_symm="SPHERICAL";m_run->m_symmetry = new SymSpherical(element, fileName.str()); }
      else { throw ErrorXMLDev(fileName.str(), __FILE__, __LINE__); }
    }
    
    //Source terms reading
    int sourceFound(0);
    element = xmlNode->FirstChildElement("sourceTerms");
    while (element != NULL) {
      sourceFound++;
      //Source reading
      std::string typeSource(element->Attribute("type"));
      if (typeSource == "") throw ErrorXMLAttribut("type", fileName.str(), __FILE__, __LINE__);
      Tools::uppercase(typeSource);
      // std::string orderSource(element->Attribute("order"));
      std::string orderSource("EULER");
      Tools::uppercase(orderSource);
      //Order
      int order = 1;
      //Switch on the type of source
      if (typeSource == "REACTION") 
      { 
        m_run->m_sources.push_back(new SourceReaction(element, order, fileName.str()));
      }
      else { throw ErrorXMLDev(fileName.str(), __FILE__, __LINE__); }
      element = element->NextSiblingElement("sourceTerms");
    }
    m_run->m_numberSources = sourceFound;

  }
  catch (ErrorXML &){ throw; } // Renvoi au niveau suivant
}

//***********************************************************************

Eos* Input::entreeEOS(int index,int &numberEOS)
{

  Eos *eos;
  eos = new EosTPIG(index, numberEOS);
  return eos;
}

//***********************************************************************

void Input::entreeAMRFire(std::string casTest)
{
  try{
    //Methode AMR: initialization des variables
    m_run->m_lvlMax = 0;

    //1) Parsing du file XML par la bibliotheque tinyxml2
    //------------------------------------------------------
    std::stringstream fileName(casTest + m_nameMesh);
    XMLDocument xmlMesh;
    XMLError error(xmlMesh.LoadFile(fileName.str().c_str())); //Le file est parse ici
    if (error != XML_SUCCESS) throw ErrorXML(fileName.str(), __FILE__, __LINE__);

    //2) Recuperation des donnees principales du compute
    //-------------------------------------------------
    //Recuperation racine du document XML
    XMLNode *mesh = xmlMesh.FirstChildElement("mesh");
    if (mesh == NULL) throw ErrorXMLRacine("mesh", fileName.str(), __FILE__, __LINE__);

    XMLElement *element;
    //Recuperation du type de mesh
    element = mesh->FirstChildElement("type");
    if (element == NULL) throw ErrorXMLElement("type", fileName.str(), __FILE__, __LINE__);
    std::string structureMesh(element->Attribute("structure"));
    if (structureMesh == "") throw ErrorXMLAttribut("structure", fileName.str(), __FILE__, __LINE__);
    Tools::uppercase(structureMesh);

    if (structureMesh == "CARTESIAN")
    {
      //----------------MESH CARTESIAN ---------------------
      XMLElement *cartesianMesh;
      cartesianMesh = mesh->FirstChildElement("cartesianMesh");
      if (cartesianMesh == NULL) throw ErrorXMLElement("cartesianMesh", fileName.str(), __FILE__, __LINE__);

      // more AMR control parameters which are declared in the unname namespace in pharos.h. - fane
      element = cartesianMesh->FirstChildElement("AMR_Fire");
      error = element->QueryDoubleAttribute("xi_diff_lim", &AMRPara::xi_diff_lim);
      if (error != XML_NO_ERROR) throw ErrorXMLAttribut("xi_diff_lim", fileName.str(), __FILE__, __LINE__);
      error = element->QueryDoubleAttribute("f1", &AMRPara::f1);
      if (error != XML_NO_ERROR) { AMRPara::f1 = 1.5; }
      error = element->QueryDoubleAttribute("f2", &AMRPara::f2);
      if (error != XML_NO_ERROR) { AMRPara::f2 = 1.5; }
      error = element->QueryDoubleAttribute("f3", &AMRPara::f3);
      if (error != XML_NO_ERROR) throw ErrorXMLAttribut("f3", fileName.str(), __FILE__, __LINE__);
      error = element->QueryDoubleAttribute("dx_expo", &AMRPara::dx_expo);
      if (error != XML_NO_ERROR) {AMRPara::dx_expo = 1.5;}
      error = element->QueryIntAttribute("div_vor_v_flag", &AMRPara::div_vor_v_flag);
      if (error != XML_NO_ERROR) { AMRPara::div_vor_v_flag = 0; }
      error = element->QueryIntAttribute("amr_iter", &AMRPara::amr_iter);
      if (error != XML_NO_ERROR) { AMRPara::amr_iter = 2;}
      error = element->QueryDoubleAttribute("amr_flux", &AMRPara::amr_flux);
      if (error != XML_NO_ERROR) {AMRPara::amr_flux = 0.1;}
      error = element->QueryIntAttribute("YH2O", &AMRPara::YH2O);
      if (error != XML_NO_ERROR) throw ErrorXMLAttribut("YH2O", fileName.str(), __FILE__, __LINE__);
      error = element->QueryIntAttribute("YCO2", &AMRPara::YCO2);
      if (error != XML_NO_ERROR) throw ErrorXMLAttribut("YCO2", fileName.str(), __FILE__, __LINE__);
      error = element->QueryIntAttribute("YOH", &AMRPara::YOH);
      if (error != XML_NO_ERROR) throw ErrorXMLAttribut("YOH", fileName.str(), __FILE__, __LINE__);
      error = element->QueryIntAttribute("gradRho_flag", &AMRPara::gradRho_flag);
      if (error != XML_NO_ERROR) throw ErrorXMLAttribut("gradRho_flag", fileName.str(), __FILE__, __LINE__);
      error = element->QueryIntAttribute("chem_step_gap", &AMRPara::chem_step_gap);
      if (error != XML_NO_ERROR) { AMRPara::chem_step_gap = 1;}
      // if (error != XML_NO_ERROR) throw ErrorXMLAttribut("chem_step_gap", fileName.str(), __FILE__, __LINE__);
      error = element->QueryIntAttribute("riemann_type", &AMRPara::riemann_type);
      if (error != XML_NO_ERROR) throw ErrorXMLAttribut("riemann_type", fileName.str(), __FILE__, __LINE__);
      error = element->QueryIntAttribute("test_prob", &AMRPara::test_prob);
      if (error != XML_NO_ERROR) { AMRPara::test_prob = -1;}
      error = element->QueryIntAttribute("species_diffusion_type", &AMRPara::species_diffusion_type);
      if (error != XML_NO_ERROR) { AMRPara::species_diffusion_type = 0;}

    }
    else{ throw ErrorXMLDev(fileName.str(), __FILE__, __LINE__); }
  }
  catch (ErrorXML &){ throw; } // Renvoi au niveau suivant
}
