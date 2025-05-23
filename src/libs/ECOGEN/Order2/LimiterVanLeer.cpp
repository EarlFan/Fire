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

//! \file      LimiterVanLeer.cpp
//! \author    F. Petitpas, K. Schmidmayer
//! \version   1.0
//! \date      July 19 2018

#include "LimiterVanLeer.h"

//***********************************************************************

LimiterVanLeer::LimiterVanLeer(){}

//***********************************************************************

LimiterVanLeer::~LimiterVanLeer(){}

//***********************************************************************

double LimiterVanLeer::limiteSlope(const double& slope1, const double& slope2) 
{
	double zero(1e-6);
	double slope(0.), produit(slope1*slope2), somme(slope1+slope2);
	if ((std::fabs(somme)>zero) && (produit > zero)) { slope = 2.*produit / somme; }
	return slope;
}

double vanleer(double r) {
	double vanleer;
	if (fabs(1.0 + r) < 1e-20) {
		vanleer = 0;
	}
	else {
		vanleer = (r + fabs(r)) / (1.0 + r);
	}
	return vanleer;
}