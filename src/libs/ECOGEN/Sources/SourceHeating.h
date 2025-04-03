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

#ifndef SOURCEHEATING_H
#define SOURCEHEATING_H

//! \file      SourceHeating.h
//! \author    F. Petitpas, J. Caze
//! \version   1.0
//! \date      October 29 2019

#include "Source.h"

//! \class     SourceHeating
//! \brief     Class for heating source terms
class SourceHeating : public Source
{
public:
  SourceHeating();
  //! \brief     Source constructor from a XML format reading
  //! \details   Reading data from XML file under the following format:
  //!            ex: <dataHeating q="1.d3"/>
  //! \param     element          XML element to read for source term
  //! \param     fileName         string name of readed XML file
  SourceHeating(tinyxml2::XMLElement *element, int order, std::string fileName = "Unknown file");
  virtual ~SourceHeating();

  virtual void prepSourceTerms(Cell *cell, const int &numberPhases, const double &dt, const int i = 0);

private:
  static const std::string NAME;
  double m_q;     //!Specific heat power (W/m3)
};

#endif //SOURCEHEATING_H