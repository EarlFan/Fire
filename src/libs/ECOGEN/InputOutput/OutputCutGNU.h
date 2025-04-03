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

#ifndef OUTPUTCUTGNU_H
#define OUTPUTCUTGNU_H

//! \file      OutputCutGNU.h
//! \author    F. Petitpas, K. Schmidmayer
//! \version   1.0
//! \date      May 03 2018

#include "OutputGNU.h"
#include "../Maths/GOLine.h"
#include "../Maths/GOPlan.h"

class OutputCutGNU : public OutputGNU
{
public:
  OutputCutGNU();
  OutputCutGNU(std::string casTest, std::string run, tinyxml2::XMLElement *element, std::string fileName, TypeGO type, Input *entree);
  virtual ~OutputCutGNU();

  virtual void ecritSolution(Mesh *mesh, std::vector<Cell *> *cellsLvl);

  virtual void prepareOutputInfos() {}; //Aucune infos a ecrire
  virtual void ecritInfos() {};
  virtual void setNumSortie(int restartTime) ;

private:
  GeometricObject *m_objet; //droite ou plan de cut
};

#endif //OUTPUTCUTGNU_H