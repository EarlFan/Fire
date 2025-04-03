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

//! \file      QAPFireDiffusion.cpp
//! \author    E Fan, Lisong Shi, Jiaao Hao, Tianhan Zhang, Chih-yung Wen.
//! \version   1.0
//! \date      Jan 01 2025

#include "QAPFireDiffusion.h"
#include <iostream>

//***********************************************************************

QAPFireDiffusion::QAPFireDiffusion(){}

//***********************************************************************

QAPFireDiffusion::QAPFireDiffusion(AddPhys* addPhys) : QuantitiesAddPhys(addPhys), m_gradsVT(4),m_gradsRhoI(NS)
{
  variableNamesVisc.resize(4);
  numPhasesVisc.resize(4);//-fane
  for (int i = 0; i < 4; ++i) {
    m_gradsVT[i] = 0.;
    numPhasesVisc[i] = -1;//since the mixture variables are used for gradient evaluation, -1 is used for species index.
  }
  for (int i = 0; i < NS; ++i)m_gradsRhoI[i] = 0.;
  variableNamesVisc[0] = velocityU;
  variableNamesVisc[1] = velocityV;
  variableNamesVisc[2] = velocityW;
  variableNamesVisc[3] = temperature;
}

//***********************************************************************

QAPFireDiffusion::~QAPFireDiffusion(){}

//***********************************************************************

void QAPFireDiffusion::computeQuantities(Cell* cell)
{
  cell->computeTransportProperties();
  cell->computeGradient(m_gradsVT, variableNamesVisc, numPhasesVisc);
  // printf("In QAPFireDiffusion::computeQuantities()./\n");
  // for (int i = 0; i < 4; ++i){
  //   printf("m_gradsVT[%d], .x = %le, .y = %le.\n",i,m_gradsVT[i].getX(),m_gradsVT[i].getY());
  // }
  cell->computeRhoIGradient(m_gradsRhoI);
}

//***********************************************************************

void QAPFireDiffusion::setGrad(const Coord &grad, int num)//1:U, 2:V, 3:W, 4:T, 5->NS+4:species
{
  if(num<=4)
  {
    m_gradsVT[num-1] = grad;
  }
  else if(num>4 && num<=NS+4)
  {
    m_gradsRhoI[num-5] = grad;
  }
  else
  {
    printf("Invalid index for QAPFireDiffusion::setGrad(), index = %d.\n",num);
    exit(EXIT_FAILURE);
  }
}

//***********************************************************************

const Coord& QAPFireDiffusion::getGrad(int num) const //1:U, 2:V, 3:W, 4:T, 5->NS+4:species
{
  if(num<=4)
  {
    return m_gradsVT[num-1];
  }
  else if(num>4 && num<=NS+4)
  {
    return m_gradsRhoI[num-5];
  }
  else
  {
    printf("Invalid index for QAPFireDiffusion::getGrad(), index = %d.\n",num);
    exit(EXIT_FAILURE);
  }   
}; 
