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

//! \file      PhaseFire.cpp
//! \author    E Fan, Lisong Shi, Jiaao Hao, Tianhan Zhang, Chih-yung Wen.
//! \version   1.0
//! \date      Jan 01 2025

#include "PhaseFire.h"
// #include "../../libs/ECOGEN/Eos/Eos.h"
#include <fstream>

using namespace std;
using namespace tinyxml2;

//***************************************************************************

PhaseFire::PhaseFire() :m_eos(0) {}

//***************************************************************************
// Set temperatrue and density for each species
PhaseFire::PhaseFire(double temperature, double density, Eos *eos) :
  m_eos(eos)
{

}

//***************************************************************************

PhaseFire::~PhaseFire() {}

//***************************************************************************

void PhaseFire::allocateAndCopyPhase(Phase **vecPhase)
{
  *vecPhase = new PhaseFire(*this);
}

//***************************************************************************

void PhaseFire::copyPhase(Phase &phase)
{
  m_eos = phase.getEos();
}
