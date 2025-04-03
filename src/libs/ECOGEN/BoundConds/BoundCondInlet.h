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

#ifndef BOUNDCONDINLET_H
#define BOUNDCONDINLET_H

//! \file      BoundCondInlet.h
//! \author    F. Petitpas, K. Schmidmayer
//! \version   1.0
//! \date      February 13 2019

#include "BoundCond.h"

class BoundCondInlet : public BoundCond
{
  public:
    BoundCondInlet();
    BoundCondInlet(const BoundCondInlet& Source, const int lvl = 0); //Constructeur de copie (utile pour AMR)
    BoundCondInlet(int numPhysique);
    virtual ~BoundCondInlet();

    virtual void creeLimite(TypeMeshContainer<CellInterface *> &cellInterfaces);
    virtual void solveRiemannLimite(Cell &cellLeft, const int &numberPhases, const double &dxLeft, double &dtMax);
    virtual void solveRiemannTransportLimite(Cell &cellLeft, const int &numberTransports) const;
    
		virtual int whoAmI() const { return 9; };

    //Pour methode AMR
    virtual void creerCellInterfaceChild();  /*!< Creer un child cell interface (non initialize) */

  protected:
  private:
};

extern BoundCondInlet BoundCondDefaut2;

#endif // BOUNDCONDINLET_H
