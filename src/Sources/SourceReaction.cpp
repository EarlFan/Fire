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

//! \file      SourceReaction.cpp
//! \author    E Fan, Lisong Shi, Jiaao Hao, Tianhan Zhang, Chih-yung Wen.
//! \version   1.0
//! \date      Jan 01 2025

#include "SourceReaction.h"

using namespace tinyxml2;

const std::string SourceReaction::NAME = "REACTION";

//***********************************************************************

SourceReaction::SourceReaction():
  Source(NAME)
{}

//***********************************************************************

SourceReaction::SourceReaction(XMLElement *element, int order, std::string fileName) : Source(NAME, order)
{
  XMLElement *sousElement(element->FirstChildElement("reaction"));
  if (sousElement == NULL) throw ErrorXMLElement("reaction", fileName, __FILE__, __LINE__);
  //Collecting attributes
  //---------------------
  XMLError error;
  //Specific heat flux value
  error = sousElement->QueryBoolAttribute("reactionFlag", &m_reactionFlag);
  if (error != XML_NO_ERROR) throw ErrorXMLAttribut("reactionFlag", fileName, __FILE__, __LINE__);
}

//***********************************************************************

SourceReaction::~SourceReaction()
{}

//***********************************************************************

void SourceReaction::integrateSourceTerms(Cell *cell, const int &numberPhases, const double &dt, std::shared_ptr<Cantera::Solution> sol_)
{
  //Source terms integration on conservative quantities
  if(m_reactionFlag!=0)
  {
    cell->buildCons(numberPhases);
    cell->getCons()->prepSourceTermsReaction(cell, dt, numberPhases,sol_);
    cell->buildPrim(numberPhases);
  }
}
