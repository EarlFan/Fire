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

#ifndef HEADERADDPHYS_H
#define HEADERADDPHYS_H

//! \file      HeaderAddPhys.h
//! \author    K. Schmidmayer
//! \version   1.0
//! \date      December 20 2017
//! \brief     New additional-physics headers should be added here

//Additional physics for Euler model
// #include "APEuler.h"

//Additional physics for Kapila model
// #include "APKapila/APKapila.h"
// #include "APKapila/SurfaceTension/APKSurfaceTension.h"
// #include "APKapila/Viscosity/APKViscosity.h"
// #include "APKapila/Conductivity/APKConductivity.h"

//Additional physics for Fire model
#include "../../../Viscous/APFire.h"
#include "../../../Viscous/DiffusionEffects/APFireDiffusion.h"

//Add here headers of new additional physics

#endif // HEADERADDPHYS_H