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

//! \file      CellO2.cpp
//! \author    E Fan, Lisong Shi, Jiaao Hao, Tianhan Zhang, Chih-yung Wen.
//! \version   1.0
//! \date      Jan 01 2025

#include "../libs/ECOGEN/Order2/CellO2.h"
#include "../libs/ECOGEN/Models/Phase.h"

//***********************************************************************

void CellO2::copyMixture(Mixture *mixture)
{
  m_mixture->copyMixture(*mixture);
  m_mixtureO2->copyMixture(*mixture);
}

//***********************************************************************
// call after timeEvolution, the time evolution should return Q_0 + DQ_1/2 * DT 
void CellO2::addHalfTimeCons(const int &numberPhases)
{
  m_consHalfTime->setBufferFlux2();      //fluxTempXXX receive half time conservative variables
  
  m_cons->addFlux(1., numberPhases);                       // add half time cons to m_cons, now m_cons = Q_0 + DQ_1/2 * DT + Q_1/2

  m_cons->multiply(0.5, numberPhases);                     // Q2 = 0.5*(Q_0 + DQ_1/2*DT + Q_1/2)
}

