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

//! \file      BoundCondInlet.cpp
//! \author    F. Petitpas, K. Schmidmayer
//! \version   1.0
//! \date      February 13 2019

#include "BoundCondInlet.h"

BoundCondInlet BoundCondDefaut2;

//****************************************************************************

BoundCondInlet::BoundCondInlet(){}

//****************************************************************************

BoundCondInlet::BoundCondInlet(const BoundCondInlet& Source, const int lvl) : BoundCond(Source)
{
  m_lvl = lvl;
}

//****************************************************************************

BoundCondInlet::BoundCondInlet(int numPhysique) : BoundCond(numPhysique)
{}

//****************************************************************************

BoundCondInlet::~BoundCondInlet() {}

//****************************************************************************

void BoundCondInlet::creeLimite(TypeMeshContainer<CellInterface *> &cellInterfaces)
{
  cellInterfaces.push_back(new BoundCondInlet(*(this)));
}

//****************************************************************************

void BoundCondInlet::solveRiemannLimite(Cell &cellLeft, const int & numberPhases, const double & dxLeft, double & dtMax)
{
  m_mod->solveRiemannIntern(cellLeft, cellLeft, numberPhases, dxLeft, dxLeft, dtMax);
}

//****************************************************************************

void BoundCondInlet::solveRiemannTransportLimite(Cell &cellLeft, const int & numberTransports) const
{
	m_mod->solveRiemannTransportIntern(cellLeft, cellLeft, numberTransports);
}

//****************************************************************************
//******************************Methode AMR***********************************
//****************************************************************************

void BoundCondInlet::creerCellInterfaceChild()
{
  m_cellInterfacesChildren.push_back(new BoundCondInlet(*this, m_lvl + 1));
}

//***********************************************************************