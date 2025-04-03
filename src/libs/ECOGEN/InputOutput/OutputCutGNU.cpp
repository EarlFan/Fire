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

//! \file      OutputCutGNU.cpp
//! \author    F. Petitpas, K. Schmidmayer
//! \version   1.1
//! \date      June 5 2019

#include "OutputCutGNU.h"

using namespace tinyxml2;

//***************************************************************

OutputCutGNU::OutputCutGNU(){
  delete m_objet;
}

//***************************************************************

OutputCutGNU::OutputCutGNU(std::string casTest, std::string run, XMLElement *element, std::string fileName, TypeGO type, Input *entree)
{
  try {
    //Modification des attributs
    m_ecritBinaire = false;
    m_simulationName = casTest;
    if (type == LINE) { m_fileNameResults = "cut1D"; }
    else if (type == PLAN) { m_fileNameResults = "cut2D"; }
    else { throw ErrorECOGEN("OutputCutGNU::OutputCutGNU : type de cut inconnu", __FILE__, __LINE__); }
    m_fileNameVisu = "visualization" + m_fileNameResults + ".gnu";
    m_folderOutput = "./results/" + run + "/cuts/";
    m_donneesSeparees = 0;
    // m_numFichier = 0;
    m_input = entree;
    m_run = m_input->getRun();


    XMLElement *sousElement;
    XMLError error;

    double donnee;
    Coord vertex, vecteur;

    //Recuperation donnees de cut 1D
    sousElement = element->FirstChildElement("vertex");
    if (sousElement == NULL) throw ErrorXMLElement("vertex", fileName, __FILE__, __LINE__);
    error = sousElement->QueryDoubleAttribute("x", &donnee); vertex.setX(donnee);
    if (error != XML_NO_ERROR) throw ErrorXMLAttribut("x", fileName, __FILE__, __LINE__);
    error = sousElement->QueryDoubleAttribute("y", &donnee); vertex.setY(donnee);
    if (error != XML_NO_ERROR) throw ErrorXMLAttribut("y", fileName, __FILE__, __LINE__);
    error = sousElement->QueryDoubleAttribute("z", &donnee); vertex.setZ(donnee);
    if (error != XML_NO_ERROR) throw ErrorXMLAttribut("z", fileName, __FILE__, __LINE__);
    if (type == LINE) sousElement = element->FirstChildElement("vecDir");
    else if (type == PLAN) { sousElement = element->FirstChildElement("vecNormal"); }
    else { throw ErrorECOGEN("OutputCutGNU::OutputCutGNU : type de cut inconnu", __FILE__, __LINE__); }
    if (sousElement == NULL) throw ErrorXMLElement("vecDir ou vecNormal", fileName, __FILE__, __LINE__);
    error = sousElement->QueryDoubleAttribute("x", &donnee); vecteur.setX(donnee);
    if (error != XML_NO_ERROR) throw ErrorXMLAttribut("x", fileName, __FILE__, __LINE__);
    error = sousElement->QueryDoubleAttribute("y", &donnee); vecteur.setY(donnee);
    if (error != XML_NO_ERROR) throw ErrorXMLAttribut("y", fileName, __FILE__, __LINE__);
    error = sousElement->QueryDoubleAttribute("z", &donnee); vecteur.setZ(donnee);
    if (error != XML_NO_ERROR) throw ErrorXMLAttribut("z", fileName, __FILE__, __LINE__);

    if (type == LINE) { m_objet = new GOLine(vertex, vecteur); }
    else if (type == PLAN) { m_objet = new GOPlan(vertex, vecteur); }
    else { throw ErrorECOGEN("OutputCutGNU::OutputCutGNU : type de cut inconnu", __FILE__, __LINE__); }

  }
  catch (ErrorECOGEN &) { throw; }
}

//***************************************************************

OutputCutGNU::~OutputCutGNU(){}

//***********************************************************************

void OutputCutGNU::ecritSolution(Mesh *mesh, std::vector<Cell *> *cellsLvl)
{
  std::ofstream fileStream;
  std::string file = m_folderOutput + creationNameFichierGNU(m_fileNameResults.c_str(), -1, rankCpu, m_numFichier);
  fileStream.open(file.c_str());
  mesh->ecritSolutionGnuplot(cellsLvl, fileStream, m_objet);
  fileStream << std::endl;
  fileStream.close();

  try {
    //Creation du file gnuplot pour visualization des resultats
    if (rankCpu == 0) {
      if (m_objet->getType() == LINE) { ecritScriptGnuplot(1); }
      else if (m_objet->getType() == PLAN) { ecritScriptGnuplot(2); }
      else { throw ErrorECOGEN("OutputCutGNU::ecritSolutionSpecifique : type de cut inconnu", __FILE__, __LINE__); }
    }
  }
  catch (ErrorECOGEN &) { throw; }
  m_numFichier++;
}

//***************************************************************

void OutputCutGNU::setNumSortie(int restartTime)
{
  m_numFichier = restartTime;
}