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

//! \file      readThermo.cpp
//! \author    E Fan, Lisong Shi, Jiaao Hao, Tianhan Zhang, Chih-yung Wen.
//! \version   1.0
//! \date      Jan 01 2025

#include "readThermo.h"
#include "cantera/thermo.h"
#include "cantera/transport.h"
#include "cantera/base/Solution.h"

int LN;

void readMW()
{
    std::shared_ptr<Cantera::ThermoPhase> gas = solution_mech->thermo();
    size_t nsp = gas->nSpecies();
    for(int i=0;i<nsp;i++)
    {
        // struct SpeciesPha *ps;
        // ps = &Sptr[i];
        Sptr[i].Name = gas->speciesName(i);
        Sptr[i].MW = gas->molecularWeight(i);
        Sptr[i].RsPha = Ru / Sptr[i].MW;
    }
}