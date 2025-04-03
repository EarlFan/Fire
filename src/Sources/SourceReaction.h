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

#ifndef SOURCEREACTION_H
#define SOURCEREACTION_H

//! \file      SourceReaction.h
//! \author    E Fan, Lisong Shi, Jiaao Hao, Tianhan Zhang, Chih-yung Wen.
//! \version   1.0
//! \date      Jan 01 2025

#include "../libs/ECOGEN/Sources/Source.h"

//! \class     SourceReaction
//! \brief     Class for heating source terms
class SourceReaction : public Source
{
public:
  SourceReaction();
  //! \brief     Source constructor from a XML format reading
  //! \details   Reading data from XML file under the following format:
  //!            ex: <dataHeating q="1.d3"/>
  //! \param     element          XML element to read for source term
  //! \param     fileName         string name of readed XML file
  SourceReaction(tinyxml2::XMLElement *element, int order, std::string fileName = "Unknown file");
  virtual ~SourceReaction();

  virtual void integrateSourceTerms(Cell *cell, const int &numberPhases, const double &dt, std::shared_ptr<Cantera::Solution> sol_);

private:
  static const std::string NAME;
  bool m_reactionFlag;     //! reaction flag
};

#endif //SOURCEREACTION_H