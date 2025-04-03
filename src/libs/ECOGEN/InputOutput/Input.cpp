//  
//       ,---.     ,--,    .---.     ,--,    ,---.    .-. .-. 
//       | .-'   .' .')   / .-. )  .' .'     | .-'    |  \| | 
//       | `-.   |  |(_)  | | |(_) |  |  __  | `-.    |   | | 
//       | .-'   \  \     | | | |  \  \ ( _) | .-'    | |\  | 
//       |  `--.  \  `-.  \ `-' /   \  `-) ) |  `--.  | | |)| 
//       /( __.'   \____\  )---'    )\____/  /( __.'  /(  (_) 
//      (__)              (_)      (__)     (__)     (__)     
//
//  This file is part of ECOGEN.
//
//  ECOGEN is the legal property of its developers, whose names 
//  are listed in the copyright file included with this source 
//  distribution.
//
//  ECOGEN is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published 
//  by the Free Software Foundation, either version 3 of the License, 
//  or (at your option) any later version.
//  
//  ECOGEN is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
//  GNU General Public License for more details.
//  
//  You should have received a copy of the GNU General Public License
//  along with ECOGEN (file LICENSE).  
//  If not, see <http://www.gnu.org/licenses/>.

//! \file      Input.cpp
//! \author    F. Petitpas, K. Schmidmayer
//! \version   1.1
//! \date      June 5 2019

#include "Input.h"
#include "HeaderInputOutput.h"
#include "../Meshes/stretchZone.h"

using namespace tinyxml2;

//***********************************************************************

Input::Input(Run *run) : m_run(run)
{
	//Attribution des numbers de version pour la lecture
	m_vMain = 5;
	m_vMesh = 5;
	m_vCI = 4;
	m_vModel = 4;

	m_nameMain = "mainV" + IO::toString(m_vMain) + ".xml";
	m_nameMesh = "meshV" + IO::toString(m_vMesh) + ".xml";
	m_nameCI = "initialConditionsV" + IO::toString(m_vCI) + ".xml";
	m_nameModel = "modelV" + IO::toString(m_vModel) + ".xml";
}

//***********************************************************************

Input::~Input(){}

//***********************************************************************

void Input::lectureInputXML(std::vector<GeometricalDomain*> &domains, std::vector<BoundCond*> &boundCond)
{
  try{
    //1) Parametres generaux du compute
    entreeMain(m_run->m_simulationName);
    //2) Donnees de Mesh
    entreeAMRFire(m_run->m_simulationName);
    entreeMesh(m_run->m_simulationName);
    //3) Donnees models et fluids
    entreeModel(m_run->m_simulationName);
    //4) Lecture des conditions initiales
    entreeConditionsInitiales(m_run->m_simulationName, domains, boundCond);
  }
  catch (ErrorXML &e){
    if(rankCpu==0) std::cerr << e.infoError() << std::endl;
    throw;
  }
}

//***********************************************************************

void Input::entreeMain(std::string casTest)
{
  try{
    //1) Parsing du file XML par la bibliotheque tinyxml2
    //------------------------------------------------------
	  std::stringstream fileName(casTest + m_nameMain);
    XMLDocument xmlMain;
    XMLError error(xmlMain.LoadFile(fileName.str().c_str())); //Le file est parse ici
    if (error != XML_SUCCESS) throw ErrorXML(fileName.str(),__FILE__, __LINE__);
    
    //2) Recuperation des donnees principales du compute
    //-------------------------------------------------
    //Recuperation racine du document XML
    XMLNode *computationParam = xmlMain.FirstChildElement("computationParam");
    if (computationParam == NULL) throw ErrorXMLRacine("computationParam", fileName.str(), __FILE__, __LINE__);

    XMLElement *element, *sousElement;

    //Recuperation name du run
    element = computationParam->FirstChildElement("run");
    if (element == NULL) throw ErrorXMLElement("run", fileName.str(), __FILE__, __LINE__);
    XMLNode* xmlNode2 = element->FirstChild();
    if (xmlNode2 == NULL) throw ErrorXMLElement("run", fileName.str(), __FILE__, __LINE__);
    XMLText* xmlText = xmlNode2->ToText();
    if (xmlText == NULL) throw ErrorXMLElement("run", fileName.str(), __FILE__, __LINE__);

    //Lecture des informations de sorties/prints
    element = computationParam->FirstChildElement("outputMode");
    if (element == NULL) throw ErrorXMLElement("outputMode", fileName.str(), __FILE__, __LINE__);
    //Lecture format sortie
    std::string format(element->Attribute("format"));
    if (format == "") throw ErrorXMLAttribut("format", fileName.str(), __FILE__, __LINE__);
    Tools::uppercase(format);
    if (format == "XML") { m_run->m_outPut = new OutputXML(casTest, xmlText->Value(), element, fileName.str(), this); }
    else if (format == "GNU") { m_run->m_outPut = new OutputGNU(casTest, xmlText->Value(), element, fileName.str(), this); }
    else { throw ErrorXMLDev(fileName.str(), __FILE__, __LINE__); }

    //Lecture des cuts 1D
    element = computationParam->FirstChildElement("cut1D");
    while (element != NULL)
    {
      m_run->m_cuts.push_back(new OutputCutGNU(casTest, xmlText->Value(), element, fileName.str(), LINE, this));
      element = element->NextSiblingElement("cut1D");
    }
    //Lecture des cuts 2D
    element = computationParam->FirstChildElement("cut2D");
    while (element != NULL)
    {
      m_run->m_cuts.push_back(new OutputCutGNU(casTest, xmlText->Value(), element, fileName.str(), PLAN, this));
      element = element->NextSiblingElement("cut2D");
    }
    //Reading probes
    element = computationParam->FirstChildElement("probe");
    while (element != NULL)
    {
      m_run->m_probes.push_back(new OutputProbeGNU(casTest, xmlText->Value(), element, fileName.str(), this));
      element = element->NextSiblingElement("probe");
    }

	//Record global quantity on whole domain (mass, total energy)
	element = computationParam->FirstChildElement("recordGlobalQuantity");
	if (element != NULL) {
		std::string quantity(element->Attribute("quantity"));
		if (quantity == "") throw ErrorXMLAttribut("quantity", fileName.str(), __FILE__, __LINE__);
		Tools::lowercase(quantity);
		if (quantity == "mass" || quantity == "energy") { m_run->m_globalQuantities.push_back(new OutputGlobalGNU(casTest, xmlText->Value(), fileName.str(), this, quantity)); }
		else { throw ErrorXMLAttribut("quantity", fileName.str(), __FILE__, __LINE__); }
	}
    
    //Recuperation Iteration / temps Physique
    element = computationParam->FirstChildElement("timeControlMode");
    if (element == NULL) throw ErrorXMLElement("timeControlMode", fileName.str(), __FILE__, __LINE__);
    error = element->QueryBoolAttribute("iterations", &m_run->m_controleIterations);
    if (error != XML_NO_ERROR) throw ErrorXMLAttribut("iterations", fileName.str(), __FILE__, __LINE__);
    if (m_run->m_controleIterations)
    {
      //Recuperation Iterations / Frequence
      sousElement = element->FirstChildElement("iterations");
      if (sousElement == NULL) throw ErrorXMLElement("iterations", fileName.str(), __FILE__, __LINE__);
      error = sousElement->QueryIntAttribute("number", &m_run->m_nbIte);
      if (error != XML_NO_ERROR) throw ErrorXMLAttribut("number", fileName.str(), __FILE__, __LINE__);
      error = sousElement->QueryIntAttribute("iterFreq", &m_run->m_freq);
      if (error != XML_NO_ERROR) throw ErrorXMLAttribut("iterFreq", fileName.str(), __FILE__, __LINE__);
    }
    else
    {
      //Recuperation Temps / Frequence
      sousElement = element->FirstChildElement("physicalTime");
      if (sousElement == NULL) throw ErrorXMLElement("physicalTime", fileName.str(), __FILE__, __LINE__);
      error = sousElement->QueryFloatAttribute("totalTime", &m_run->m_finalPhysicalTime);
      if (error != XML_NO_ERROR) throw ErrorXMLAttribut("totalTime", fileName.str(), __FILE__, __LINE__);
      error = sousElement->QueryFloatAttribute("timeFreq", &m_run->m_timeFreq);
      if (error != XML_NO_ERROR) throw ErrorXMLAttribut("timeFreq", fileName.str(), __FILE__, __LINE__);
    }

    //Recuperation CFL
    element = computationParam->FirstChildElement("computationControl");
    if (element == NULL) throw ErrorXMLElement("computationControl", fileName.str(), __FILE__, __LINE__);
    error = element->QueryDoubleAttribute("CFL", &m_run->m_cfl);
    if (error != XML_NO_ERROR) throw ErrorXMLAttribut("CFL", fileName.str(), __FILE__, __LINE__);

    //Lecture Ordre2
    element = computationParam->FirstChildElement("secondOrder");
    if(element == NULL) {
      m_run->m_order = "FIRSTORDER";
      m_run->m_globalLimiter = new Limiter;
      m_run->m_interfaceLimiter = new Limiter;
      m_run->m_globalVolumeFractionLimiter = new Limiter;
      m_run->m_interfaceVolumeFractionLimiter = new Limiter;
    }
    else{
      m_run->m_order = "SECONDORDER";
      XMLNode *contenu;
      //Recuperation global limiter
      sousElement = element->FirstChildElement("globalLimiter");
      if (sousElement == NULL) throw ErrorXMLElement("globalLimiter", fileName.str(), __FILE__, __LINE__);
      contenu = sousElement->FirstChild();
      if (contenu == NULL) throw ErrorXMLElement("globalLimiter", fileName.str(), __FILE__, __LINE__);
      std::string globalLimiter = contenu->ToText()->Value();
      Tools::uppercase(globalLimiter);
      if (globalLimiter == "MINMOD") { m_run->m_globalLimiter = new LimiterMinmod; }
      else if (globalLimiter == "VANLEER") { m_run->m_globalLimiter = new LimiterVanLeer; }
      else if (globalLimiter == "VANALBADA") { m_run->m_globalLimiter = new LimiterVanAlbada; }
      else if (globalLimiter == "SUPERBEE") { m_run->m_globalLimiter = new LimiterSuperBee; }
      else if (globalLimiter == "MC") { m_run->m_globalLimiter = new LimiterMC; }
      else if (globalLimiter == "THINC") { throw ErrorXMLAttribut("THINC can only be a volume-fraction limiter", fileName.str(), __FILE__, __LINE__); }
      else { throw ErrorXMLDev(fileName.str(), __FILE__, __LINE__); }
      //Recuperation interface limiter
      std::string interfaceLimiter = globalLimiter;
      sousElement = element->FirstChildElement("interfaceLimiter");
      if (sousElement != NULL) {
        contenu = sousElement->FirstChild();
        if (contenu == NULL) throw ErrorXMLElement("interfaceLimiter", fileName.str(), __FILE__, __LINE__);
        interfaceLimiter = contenu->ToText()->Value();
      }
      Tools::uppercase(interfaceLimiter);
      if (interfaceLimiter == "MINMOD") { m_run->m_interfaceLimiter = new LimiterMinmod; }
      else if (interfaceLimiter == "VANLEER") { m_run->m_interfaceLimiter = new LimiterVanLeer; }
      else if (interfaceLimiter == "VANALBADA") { m_run->m_interfaceLimiter = new LimiterVanAlbada; }
      else if (interfaceLimiter == "SUPERBEE") { m_run->m_interfaceLimiter = new LimiterSuperBee; }
      else if (interfaceLimiter == "MC") { m_run->m_interfaceLimiter = new LimiterMC; }
      else if (interfaceLimiter == "THINC") { throw ErrorXMLAttribut("THINC can only be a volume-fraction limiter", fileName.str(), __FILE__, __LINE__); }
      else { throw ErrorXMLDev(fileName.str(), __FILE__, __LINE__); }
      //Recuperation global volume-fraction limiter
      std::string globalVolumeFractionLimiter = globalLimiter;
      sousElement = element->FirstChildElement("globalVolumeFractionLimiter");
      if (sousElement != NULL) {
        contenu = sousElement->FirstChild();
        if (contenu == NULL) throw ErrorXMLElement("globalVolumeFractionLimiter", fileName.str(), __FILE__, __LINE__);
        globalVolumeFractionLimiter = contenu->ToText()->Value();
      }
      Tools::uppercase(globalVolumeFractionLimiter);
      if (globalVolumeFractionLimiter == "MINMOD") { m_run->m_globalVolumeFractionLimiter = new LimiterMinmod; }
      else if (globalVolumeFractionLimiter == "VANLEER") { m_run->m_globalVolumeFractionLimiter = new LimiterVanLeer; }
      else if (globalVolumeFractionLimiter == "VANALBADA") { m_run->m_globalVolumeFractionLimiter = new LimiterVanAlbada; }
      else if (globalVolumeFractionLimiter == "SUPERBEE") { m_run->m_globalVolumeFractionLimiter = new LimiterSuperBee; }
      else if (globalVolumeFractionLimiter == "MC") { m_run->m_globalVolumeFractionLimiter = new LimiterMC; }
      else if (globalVolumeFractionLimiter == "THINC") { m_run->m_globalVolumeFractionLimiter = new LimiterTHINC; }
      else { throw ErrorXMLDev(fileName.str(), __FILE__, __LINE__); }
      //Recuperation interface volume-fraction limiter
      std::string interfaceVolumeFractionLimiter = interfaceLimiter;
      sousElement = element->FirstChildElement("interfaceVolumeFractionLimiter");
      if (sousElement != NULL) {
        contenu = sousElement->FirstChild();
        if (contenu == NULL) throw ErrorXMLElement("interfaceVolumeFractionLimiter", fileName.str(), __FILE__, __LINE__);
        interfaceVolumeFractionLimiter = contenu->ToText()->Value();
      }
      Tools::uppercase(interfaceVolumeFractionLimiter);
      if (interfaceVolumeFractionLimiter == "MINMOD") { m_run->m_interfaceVolumeFractionLimiter = new LimiterMinmod; }
      else if (interfaceVolumeFractionLimiter == "VANLEER") { m_run->m_interfaceVolumeFractionLimiter = new LimiterVanLeer; }
      else if (interfaceVolumeFractionLimiter == "VANALBADA") { m_run->m_interfaceVolumeFractionLimiter = new LimiterVanAlbada; }
      else if (interfaceVolumeFractionLimiter == "SUPERBEE") { m_run->m_interfaceVolumeFractionLimiter = new LimiterSuperBee; }
      else if (interfaceVolumeFractionLimiter == "MC") { m_run->m_interfaceVolumeFractionLimiter = new LimiterMC; }
      else if (interfaceVolumeFractionLimiter == "THINC") { m_run->m_interfaceVolumeFractionLimiter = new LimiterTHINC; }
      else { throw ErrorXMLDev(fileName.str(), __FILE__, __LINE__); }

    }

    //Reprise de calcul depuis file resultat
    element = computationParam->FirstChildElement("restartSimulation");
    if (element != NULL) { 
      error = element->QueryIntAttribute("restartFileNumber", &m_run->m_restartSimulation);
      if (error != XML_NO_ERROR) throw ErrorXMLAttribut("restartFileNumber", fileName.str(), __FILE__, __LINE__);
      error = element->QueryIntAttribute("AMRsaveFreq", &m_run->m_restartAMRsaveFreq);
      if (error != XML_NO_ERROR) throw ErrorXMLAttribut("AMRsaveFreq", fileName.str(), __FILE__, __LINE__);
    }

  }
  catch (ErrorXML &){ throw; } // Renvoi au niveau suivant
}

//***********************************************************************

void Input::entreeMesh(std::string casTest)
{
  try{
    //Methode AMR: initialization des variables
    m_run->m_lvlMax = 0;
		double criteriaVar(1.e10);
		bool varRho(false), varP(false), varU(false), varAlpha(false),varT(false),varYH2O(false);
		double xiSplit(1.), xiJoin(1.);
    double T_threshhold;

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
    if (structureMesh == "UNSTRUCTURED")
    {
      //----------------MESH NON STRUCTURE ---------------------
      XMLElement *meshNS;
      meshNS = mesh->FirstChildElement("unstructuredMesh");
      if (meshNS == NULL) throw ErrorXMLElement("unstructuredMesh", fileName.str(), __FILE__, __LINE__);
      //Lecture name du file contenant les informations de mesh
      element = meshNS->FirstChildElement("file");
      if (element == NULL) throw ErrorXMLElement("file", fileName.str(), __FILE__, __LINE__);
      std::string meshFile(element->Attribute("name"));
      if (meshFile == "") throw ErrorXMLAttribut("name", fileName.str(), __FILE__, __LINE__);
      //Read format of the mesh file (.msh, .mesh ...)
      std::string meshExtension(MeshUnStruct::readMeshFileExtension(meshFile));
      Tools::uppercase(meshExtension);
      if (meshExtension == "MSH") // Gmsh format
      {
        std::string version(MUSGmsh::readVersion(meshFile));
        if (version == "2.2") {
			    bool switchTags(false);
			    element = meshNS->FirstChildElement("tag");
			    if (element != NULL) {
				    error = element->QueryBoolAttribute("GMSHSwitchTags", &switchTags);
				    if (error != XML_NO_ERROR) throw ErrorXMLAttribut("GMSHSwitchTags", fileName.str(), __FILE__, __LINE__);
			    }
			    m_run->m_mesh = new MUSGmshV2(meshFile, meshExtension, switchTags);
		    }
        else if (version == "4.1") m_run->m_mesh = new MUSGmshV4(meshFile, meshExtension);
        else throw ErrorXML("mesh version not found for file : " + meshFile, __FILE__, __LINE__);
      }
      else if (meshExtension == "MESH") throw ErrorXML("MESH format is not supported for file : " + meshFile, __FILE__, __LINE__);
      else { throw ErrorXML("mesh extension not supported for file : " + meshFile, __FILE__, __LINE__); }
      
      //Recuperation pretraitement parallele
      element = meshNS->FirstChildElement("parallel");
      if (element != NULL) {
        error = element->QueryBoolAttribute("GMSHPretraitement", &m_run->m_parallelPreTreatment);
        if (error != XML_NO_ERROR) throw ErrorXMLAttribut("GMSHPretraitement", fileName.str(), __FILE__, __LINE__);
      }
      //Methode AMR non possible avec mesh non structure
      element = meshNS->FirstChildElement("AMR");
      if (element != NULL) { throw ErrorXMLAttribut("Methode AMR non possible avec mesh non structure", fileName.str(), __FILE__, __LINE__); }
    }
    else if (structureMesh == "CARTESIAN")
    {
      //----------------MESH CARTESIAN ---------------------
      XMLElement *cartesianMesh;
      cartesianMesh = mesh->FirstChildElement("cartesianMesh");
      if (cartesianMesh == NULL) throw ErrorXMLElement("cartesianMesh", fileName.str(), __FILE__, __LINE__);

      //Recuperation des dimensions
      double lX, lY, lZ;
      element = cartesianMesh->FirstChildElement("dimensions");
      if (element == NULL) throw ErrorXMLElement("dimensions", fileName.str(), __FILE__, __LINE__);
      error = element->QueryDoubleAttribute("x", &lX);
      if (error != XML_NO_ERROR) throw ErrorXMLAttribut("x", fileName.str(), __FILE__, __LINE__);
      error = element->QueryDoubleAttribute("y", &lY);
      if (error != XML_NO_ERROR) throw ErrorXMLAttribut("y", fileName.str(), __FILE__, __LINE__);
      error = element->QueryDoubleAttribute("z", &lZ);
      if (error != XML_NO_ERROR) throw ErrorXMLAttribut("z", fileName.str(), __FILE__, __LINE__);
      //Recuperation des numbers de mailles
      int nbX, nbY, nbZ;
      element = cartesianMesh->FirstChildElement("numberCells");
      if (element == NULL) throw ErrorXMLElement("numberCells", fileName.str(), __FILE__, __LINE__);
      error = element->QueryIntAttribute("x", &nbX);
      if (error != XML_NO_ERROR) throw ErrorXMLAttribut("x", fileName.str(), __FILE__, __LINE__);
      if (nbX <= 1)  throw ErrorXMLAttribut("Number of cells in the x-direction has to be superior to 1", fileName.str(), __FILE__, __LINE__);
      error = element->QueryIntAttribute("y", &nbY);
      if (error != XML_NO_ERROR) throw ErrorXMLAttribut("y", fileName.str(), __FILE__, __LINE__);
      error = element->QueryIntAttribute("z", &nbZ);
      if (error != XML_NO_ERROR) throw ErrorXMLAttribut("z", fileName.str(), __FILE__, __LINE__);
      if (nbZ > 1 && nbY <= 1) throw ErrorXMLAttribut("Number of cells in the z-direction can't be superior to 1 if number of cells in the y-direction is equal to 1", fileName.str(), __FILE__, __LINE__);
      
      // set the dimension
      if (nbX > 1 && nbY > 1 && nbZ > 1) {
        AMRPara::dim = 3;
      }
      else if ( nbX > 1 && nbY > 1 && nbZ == 1) {
        AMRPara::dim = 2;
      }
      else if ( nbX > 1 && nbY == 1 && nbZ == 1 ) {
        AMRPara::dim = 1;
      }
      else {
        printf("nbX = %d, nbY = %d, nbZ = %d.\n", nbX, nbY, nbZ);
        throw(" Invalid dimension defined in Cartesian mesh.\n");
      }

      //Stretching data
      std::vector<stretchZone> stretchX, stretchY, stretchZ;
      element = cartesianMesh->FirstChildElement("meshStretching");
      if (element != NULL) {
        double tempStart, tempEnd, tempFactor;
        int tempNumberCells;
        XMLElement *sousElement;
        //X stretching
        if (nbX > 1) {
          sousElement = element->FirstChildElement("XStretching");
          if (sousElement != NULL) {
            XMLElement *stretchElement(sousElement->FirstChildElement("stretch"));
            while (stretchElement != NULL) {
              error = stretchElement->QueryDoubleAttribute("startAt", &tempStart);
              if (error != XML_NO_ERROR) throw ErrorXMLAttribut("startAt", fileName.str(), __FILE__, __LINE__);
              error = stretchElement->QueryDoubleAttribute("endAt", &tempEnd);
              if (error != XML_NO_ERROR) throw ErrorXMLAttribut("endAt", fileName.str(), __FILE__, __LINE__);
              error = stretchElement->QueryDoubleAttribute("factor", &tempFactor);
              if (error != XML_NO_ERROR) throw ErrorXMLAttribut("factor", fileName.str(), __FILE__, __LINE__);
              error = stretchElement->QueryIntAttribute("numberCells", &tempNumberCells);
              if (error != XML_NO_ERROR) throw ErrorXMLAttribut("numberCells", fileName.str(), __FILE__, __LINE__);
              stretchX.push_back(stretchZone(tempStart, tempEnd, tempFactor, tempNumberCells));
              stretchElement = stretchElement->NextSiblingElement("stretch");
            }
          }
        }
        //Y stretching
        if (nbY > 1) {
          sousElement = element->FirstChildElement("YStretching");
          if (sousElement != NULL) {
            XMLElement *stretchElement(sousElement->FirstChildElement("stretch"));
            while (stretchElement != NULL) {
              error = stretchElement->QueryDoubleAttribute("startAt", &tempStart);
              if (error != XML_NO_ERROR) throw ErrorXMLAttribut("startAt", fileName.str(), __FILE__, __LINE__);
              error = stretchElement->QueryDoubleAttribute("endAt", &tempEnd);
              if (error != XML_NO_ERROR) throw ErrorXMLAttribut("endAt", fileName.str(), __FILE__, __LINE__);
              error = stretchElement->QueryDoubleAttribute("factor", &tempFactor);
              if (error != XML_NO_ERROR) throw ErrorXMLAttribut("factor", fileName.str(), __FILE__, __LINE__);
              error = stretchElement->QueryIntAttribute("numberCells", &tempNumberCells);
              if (error != XML_NO_ERROR) throw ErrorXMLAttribut("numberCells", fileName.str(), __FILE__, __LINE__);
              stretchY.push_back(stretchZone(tempStart, tempEnd, tempFactor, tempNumberCells));
              stretchElement = stretchElement->NextSiblingElement("stretch");
            }
          }
        }
        //Z stretching
        if (nbZ > 1) {
          sousElement = element->FirstChildElement("ZStretching");
          if (sousElement != NULL) {
            XMLElement *stretchElement(sousElement->FirstChildElement("stretch"));
            while (stretchElement != NULL) {
              error = stretchElement->QueryDoubleAttribute("startAt", &tempStart);
              if (error != XML_NO_ERROR) throw ErrorXMLAttribut("startAt", fileName.str(), __FILE__, __LINE__);
              error = stretchElement->QueryDoubleAttribute("endAt", &tempEnd);
              if (error != XML_NO_ERROR) throw ErrorXMLAttribut("endAt", fileName.str(), __FILE__, __LINE__);
              error = stretchElement->QueryDoubleAttribute("factor", &tempFactor);
              if (error != XML_NO_ERROR) throw ErrorXMLAttribut("factor", fileName.str(), __FILE__, __LINE__);
              error = stretchElement->QueryIntAttribute("numberCells", &tempNumberCells);
              if (error != XML_NO_ERROR) throw ErrorXMLAttribut("numberCells", fileName.str(), __FILE__, __LINE__);
              stretchZ.push_back(stretchZone(tempStart, tempEnd, tempFactor, tempNumberCells));
              stretchElement = stretchElement->NextSiblingElement("stretch");
            }
          }
        }

        //Stretching errors verifications
        stretchZone::verifyStretching(stretchX, lX, fileName.str());
        stretchZone::verifyStretching(stretchY, lY, fileName.str());
        stretchZone::verifyStretching(stretchZ, lZ, fileName.str());

      } //End stretching

      //Recuperation des variables pour methode AMR
      element = cartesianMesh->FirstChildElement("AMR");
      if (element != NULL) {
        //Interdiction parallele
        error = element->QueryIntAttribute("lvlMax", &m_run->m_lvlMax);
        if (error != XML_NO_ERROR) throw ErrorXMLAttribut("lvlMax", fileName.str(), __FILE__, __LINE__);
        AMRPara::g_lvlmax = m_run->m_lvlMax;
        error = element->QueryIntAttribute("lvlHydro", &m_run->m_lvlHydro);
        if (error != XML_NO_ERROR) throw ErrorXMLAttribut("lvlHydro", fileName.str(), __FILE__, __LINE__);
        error = element->QueryIntAttribute("lvlChem", &m_run->m_lvlChem);
        if (error != XML_NO_ERROR) throw ErrorXMLAttribut("lvlChem", fileName.str(), __FILE__, __LINE__);
        error = element->QueryIntAttribute("lvlGeo", &m_run->m_lvlGeo);
        if (error != XML_NO_ERROR) {m_run->m_lvlGeo = 0;}
        error = element->QueryDoubleAttribute("criteriaVar", &criteriaVar);
        if (error != XML_NO_ERROR) throw ErrorXMLAttribut("criteriaVar", fileName.str(), __FILE__, __LINE__);
        error = element->QueryBoolAttribute("varRho", &varRho);
        if (error != XML_NO_ERROR) throw ErrorXMLAttribut("varRho", fileName.str(), __FILE__, __LINE__);
        error = element->QueryBoolAttribute("varP", &varP);
        if (error != XML_NO_ERROR) throw ErrorXMLAttribut("varP", fileName.str(), __FILE__, __LINE__);
        error = element->QueryBoolAttribute("varU", &varU);
        if (error != XML_NO_ERROR) throw ErrorXMLAttribut("varU", fileName.str(), __FILE__, __LINE__);
        error = element->QueryBoolAttribute("varT", &varT);
        if (error != XML_NO_ERROR) throw ErrorXMLAttribut("varT", fileName.str(), __FILE__, __LINE__);
        error = element->QueryBoolAttribute("varAlpha", &varAlpha);
        if (error != XML_NO_ERROR) throw ErrorXMLAttribut("varAlpha", fileName.str(), __FILE__, __LINE__);
        error = element->QueryDoubleAttribute("T_threshhold", &T_threshhold);
        if (error != XML_NO_ERROR) throw ErrorXMLAttribut("T_threshhold", fileName.str(), __FILE__, __LINE__);
        error = element->QueryDoubleAttribute("xiSplit", &xiSplit);
        if (error != XML_NO_ERROR) throw ErrorXMLAttribut("xiSplit", fileName.str(), __FILE__, __LINE__);
        error = element->QueryDoubleAttribute("xiJoin", &xiJoin);
        if (error != XML_NO_ERROR) throw ErrorXMLAttribut("xiJoin", fileName.str(), __FILE__, __LINE__);
        m_run->m_mesh = new MeshCartesianAMR(lX, nbX, lY, nbY, lZ, nbZ, stretchX, stretchY, stretchZ, T_threshhold, 
          m_run->m_lvlMax, m_run->m_lvlHydro,m_run->m_lvlChem, m_run->m_lvlGeo,criteriaVar, varT, varYH2O, varRho, varP, varU, varAlpha, xiSplit, xiJoin);
      }
      else {
        m_run->m_mesh = new MeshCartesian(lX, nbX, lY, nbY, lZ, nbZ, stretchX, stretchY, stretchZ);
      }

    }
    else{ throw ErrorXMLDev(fileName.str(), __FILE__, __LINE__); }
  }
  catch (ErrorXML &){ throw; } // Renvoi au niveau suivant
}

//***********************************************************************

Eos* Input::entreeEOS(std::string EOS, int &numberEOS)
{
  try{
    //1) Parsing du file XML par la bibliotheque tinyxml2
    //------------------------------------------------------
    std::stringstream fileName("./libEOS/" + EOS);
    XMLDocument xmlEOS;
    XMLError error(xmlEOS.LoadFile(fileName.str().c_str())); //Le file est parse ici
    if (error != XML_SUCCESS) throw ErrorXML(fileName.str(), __FILE__, __LINE__);

    //2) Recuperation des donnees de l'EOS
    //------------------------------------
    //Recuperation racine du document XML
    XMLNode *xmlNode = xmlEOS.FirstChildElement("parametersEOS");
    if (xmlNode == NULL) throw ErrorXMLRacine("parametersEOS", fileName.str(), __FILE__, __LINE__);
    //Recuperation type d'EOS
    XMLElement *element;
    element = xmlNode->FirstChildElement("EOS");
    if (element == NULL) throw ErrorXMLElement("EOS", fileName.str(), __FILE__, __LINE__);
    std::string typeEOS(element->Attribute("type"));
    if (typeEOS == "") throw ErrorXMLAttribut("type", fileName.str(), __FILE__, __LINE__);
    Tools::uppercase(typeEOS);
    //Switch selon EOS
    Eos *eos;
    std::vector<std::string> NamesParametresEos;
    if      (typeEOS == "IG"){ eos = new EosIG(NamesParametresEos, numberEOS); }
    else if (typeEOS == "SG"){ eos = new EosSG(NamesParametresEos, numberEOS); }
    else if (typeEOS == "NASG"){ eos = new EosNASG(NamesParametresEos, numberEOS); }
    else{ throw ErrorXMLEOSInconnue(typeEOS, fileName.str(), __FILE__, __LINE__); } //Cas ou la loi state est inconnue

    //Recuperation des parametres de l'EOS
    element = xmlNode->FirstChildElement("parameters");
    if (element == NULL) throw ErrorXMLElement("parameters", fileName.str(), __FILE__, __LINE__);
    std::vector<double> parametresEos(NamesParametresEos.size());
    for (unsigned int p = 0; p < NamesParametresEos.size(); p++)
    {
      error = element->QueryDoubleAttribute(NamesParametresEos[p].c_str(), &parametresEos[p]);
      if (error != XML_NO_ERROR) throw ErrorXMLAttribut(NamesParametresEos[p].c_str(), fileName.str(), __FILE__, __LINE__);
    }
    eos->assignParametersEos(EOS.c_str(), parametresEos);
    //Lecture des parametres physiques (viscosite, conductivite, etc.)
    eos->readPhysicalParameter(xmlNode, fileName.str());

    return eos;
  }
  catch (ErrorXML &){ throw; } // Renvoi au niveau suivant
}

//***********************************************************************

void Input::entreeConditionsInitiales(std::string casTest, std::vector<GeometricalDomain*> &domains, std::vector<BoundCond*> &boundCond)
{
  try{
    //1) Parsing du file XML par la bibliotheque tinyxml2
    //------------------------------------------------------
    std::stringstream fileName(casTest + m_nameCI);
    XMLDocument xmlModel;
    XMLError error(xmlModel.LoadFile(fileName.str().c_str())); //Le file est parse ici
    if (error != XML_SUCCESS) throw ErrorXML(fileName.str(), __FILE__, __LINE__);

    //------------------------------ CONDITIONS INITIALES -------------------------------

    //2) Recuperation racine CI du document XML
    //-----------------------------------------
    XMLNode *xmlNode = xmlModel.FirstChildElement("CI");
    if (xmlNode == NULL) throw ErrorXMLRacine("CI", fileName.str(), __FILE__, __LINE__);
    XMLElement *elementDomaine(xmlNode->FirstChildElement("physicalDomains"));
    if(elementDomaine == NULL) throw ErrorXMLElement("physicalDomains", fileName.str(), __FILE__, __LINE__);

    //3) Recuperation des domains definis
    //------------------------------------
    std::string nameDomaine, stateDomaine;
    int domaineTrouve(0);
    XMLElement *element(elementDomaine->FirstChildElement("domain"));
    while (element != NULL)
    {
      //A)Lecture name domain
      //*********************
      nameDomaine = element->Attribute("name");
      if (nameDomaine == "") throw ErrorXMLAttribut("name", fileName.str(), __FILE__, __LINE__);
      Tools::uppercase(nameDomaine);

      //B)Lecture state domain
      //**********************
      stateDomaine = element->Attribute("state");
      if (stateDomaine == "") throw ErrorXMLAttribut("state", fileName.str(), __FILE__, __LINE__);
      Tools::uppercase(stateDomaine);
      //Recherche de l state associe au domain
      XMLElement *state(xmlNode->FirstChildElement("state"));
      bool trouve(false); std::string nameEtat;
      while (state != NULL)
      {
        nameEtat = state->Attribute("name");
        Tools::uppercase(nameEtat);
        if (nameEtat == stateDomaine){ trouve = true; break; }
        state = state->NextSiblingElement("state");
      }
      if (!trouve){ throw ErrorXMLEtat(stateDomaine, fileName.str(), __FILE__, __LINE__); }

      //Reading mixture state first
      Mixture *stateMixture(0);
      
      stateMixture = new MixFire(state, fileName.str());

      //Then Reading phases states
      std::vector<Phase*> statesPhases;
      std::string typeMateriau; std::string nameEOS;
      {
        double initialT,initialRho[NS],initialP;
        XMLElement* mixture(state->FirstChildElement("mixture"));
        int nbMateriauxEtat(0);
        while (mixture != NULL)
        {
          // read the dataMix
          XMLElement *dataMix(mixture->FirstChildElement("dataMix"));
          XMLError error;

          error = dataMix->QueryDoubleAttribute("temperature", &initialT);
          if (error != XML_NO_ERROR) throw ErrorXMLAttribut("temperature", fileName.str(), __FILE__, __LINE__);
        
          XMLElement *densityArray(mixture->FirstChildElement("densityArray"));

          struct SpeciesPha *ps;
          for(int ll=0;ll<NS;ll++)
          {
            ps = &Sptr[ll];
            error = densityArray->QueryDoubleAttribute(ps->Name.c_str(), &initialRho[ll]);
            if (error != XML_NO_ERROR) initialRho[ll] = 0;
            // if (error != XML_NO_ERROR) throw ErrorXMLAttribut(ps->Name.c_str(), fileName.str(), __FILE__, __LINE__);
          }

          for(int ll=0;ll<m_run->m_numberPhases;ll++)
          {
            statesPhases.push_back(new PhaseFire(initialT, initialRho[ll], m_run->m_eos[ll]));
          }

          mixture = mixture->NextSiblingElement("mixture");
        }
      }

      //Lecture des variables transportees pour l'state trouve
      std::vector<Transport> statesTransport(m_run->m_numberTransports);
      std::string nameTransport; double valueTransport(0.);
      XMLElement *elementTransport(state->FirstChildElement("transport"));
      int nbTransports(0);
      while (elementTransport != NULL) {
        nameTransport = elementTransport->Attribute("name");
        elementTransport->QueryDoubleAttribute("value", &valueTransport);
        int e(0);
        for (e = 0; e < m_run->m_numberTransports; e++) {
          if (nameTransport == m_run->m_nameGTR[e]) { break; }
        }
        if (e != m_run->m_numberTransports) {
          statesTransport[e].setValue(valueTransport);
          nbTransports++;
        }
        elementTransport = elementTransport->NextSiblingElement("transport");
      }

      //C)Reading optional physical entity
      //**********************************
      int physicalEntity(element->IntAttribute("physicalEntity")); //Default value: 0      

      //C)Reading domain type
      //*********************
      std::string typeDomaine(element->Attribute("type"));
      Tools::uppercase(typeDomaine);
      std::vector<std::string> NamesParametresDomaine;
      if      (typeDomaine == "ENTIREDOMAIN")                    { domains.push_back(new GDEntireDomain(nameDomaine, statesPhases, stateMixture, statesTransport, physicalEntity)); }
      else if (typeDomaine == "ENTIREDOMAINWITHPARTICULARITIES") { domains.push_back(new GDEntireDomainWithParticularities(nameDomaine, statesPhases, stateMixture, statesTransport, physicalEntity)); }
      else if (typeDomaine == "HALFSPACE")                       { domains.push_back(new GDHalfSpace(nameDomaine, statesPhases, stateMixture, statesTransport, element, physicalEntity, fileName.str())); }
      else if (typeDomaine == "DISC")                            { domains.push_back(new GDDisc(nameDomaine, statesPhases, stateMixture, statesTransport, element, physicalEntity, fileName.str())); }
      else if (typeDomaine == "ELLIPSE")                         { domains.push_back(new GDEllipse(nameDomaine, statesPhases, stateMixture, statesTransport, element, physicalEntity, fileName.str())); }
      else if (typeDomaine == "RECTANGLE")                       { domains.push_back(new GDRectangle(nameDomaine, statesPhases, stateMixture, statesTransport, element, physicalEntity, fileName.str())); }
      else if (typeDomaine == "CUBOID")                          { domains.push_back(new GDCuboid(nameDomaine, statesPhases, stateMixture, statesTransport, element, physicalEntity, fileName.str())); }
      else if (typeDomaine == "SPHERE")                          { domains.push_back(new GDSphere(nameDomaine, statesPhases, stateMixture, statesTransport, element, physicalEntity, fileName.str())); }
      else if (typeDomaine == "ELLIPSOID")                       { domains.push_back(new GDEllipsoid(nameDomaine, statesPhases, stateMixture, statesTransport, element, physicalEntity, fileName.str())); }
      else if (typeDomaine == "CYLINDER")                        { domains.push_back(new GDCylinder(nameDomaine, statesPhases, stateMixture, statesTransport, element, physicalEntity, fileName.str())); }
      else{ throw ErrorXMLDomaineInconnu(typeDomaine, fileName.str(), __FILE__, __LINE__); } //Cas ou le domain n a pas ete implemente
      domaineTrouve++;
      //Domaine suivant
      element = element->NextSiblingElement("domain");

      for (int k = 0; k < m_run->m_numberPhases; k++) { delete statesPhases[k]; }
      delete stateMixture;
    }

    //------------------------------ CONDITIONS AUX LIMITES -------------------------------
    //4) Recuperation racine CL du document XML
    //-----------------------------------------
    XMLElement *elementLimites(xmlNode->FirstChildElement("boundaryConditions"));
    if (elementLimites == NULL) throw ErrorXMLElement("boundaryConditions", fileName.str(), __FILE__, __LINE__);

    //5) Recuperation des boundCond definies
    //------------------------------------
    std::string nameBoundCond, stateLimite, typeBoundCond;
    int numBoundCond;
    int boundCondTrouve(0);
    element = elementLimites->FirstChildElement("boundCond");
    while (element != NULL)
    {
      //A)Lecture name boundCond
      //************************
      nameBoundCond = element->Attribute("name");
      if (nameBoundCond == "") throw ErrorXMLAttribut("name", fileName.str(), __FILE__, __LINE__);
      Tools::uppercase(nameBoundCond);

      //B)Lecture type boundCond
      //************************
      typeBoundCond = element->Attribute("type");
      if (typeBoundCond == "") throw ErrorXMLAttribut("type", fileName.str(), __FILE__, __LINE__);
      Tools::uppercase(typeBoundCond);

      //C)Number de la boundCond associe a la geometrie mesh
      //****************************************************
      error = element->QueryIntAttribute("number", &numBoundCond);
      if (error != XML_NO_ERROR) throw ErrorXMLAttribut("number", fileName.str(), __FILE__, __LINE__);

      //D)Reading conLim specific data
      //******************************
      if (typeBoundCond == "NONREFLECTING") { boundCond.push_back(new BoundCondNonReflecting(numBoundCond)); }
      else if (typeBoundCond == "SYMMETRY") { if (m_run->m_order == "FIRSTORDER") { boundCond.push_back(new BoundCondWall(numBoundCond)); } else { boundCond.push_back(new BoundCondSymmetryO2(numBoundCond)); } }
      else if (typeBoundCond == "WALL") { if (m_run->m_order == "FIRSTORDER") { boundCond.push_back(new BoundCondWall(numBoundCond)); } else { boundCond.push_back(new BoundCondWallO2(numBoundCond)); } }
      else if (typeBoundCond == "INJECTION") { boundCond.push_back(new BoundCondInj(numBoundCond, element, m_run->m_numberPhases, m_run->m_numberTransports, m_run->m_nameGTR, m_run->m_eos, fileName.str())); }
      else if (typeBoundCond == "TANK") { boundCond.push_back(new BoundCondTank(numBoundCond, element, m_run->m_numberPhases, m_run->m_numberTransports, m_run->m_nameGTR, m_run->m_eos, fileName.str())); }
      else if (typeBoundCond == "OUTFLOW") { boundCond.push_back(new BoundCondOutflow(numBoundCond, element, m_run->m_numberPhases, m_run->m_numberTransports, m_run->m_nameGTR, fileName.str())); }
      else if (typeBoundCond == "INLET") { boundCond.push_back(new BoundCondInlet(numBoundCond)); }
      else { throw ErrorXMLBoundCondInconnue(typeBoundCond, fileName.str(), __FILE__, __LINE__); } //Cas ou la limite n a pas ete implemente
      boundCondTrouve++;
      //Domaine suivant
      element = element->NextSiblingElement("boundCond");
    }

  }
  catch (ErrorXML &){ throw; } // Renvoi au niveau suivant
}

//***********************************************************************
