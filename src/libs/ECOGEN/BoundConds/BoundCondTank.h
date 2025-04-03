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

#ifndef BOUNDCONDTANK_H
#define BOUNDCONDTANK_H

//! \file      BoundCondTank.cpp
//! \author    F. Petitpas, K. Schmidmayer
//! \version   1.0
//! \date      February 13 2019

#include "BoundCond.h"

class BoundCondTank : public BoundCond
{
  public:
    BoundCondTank();
    BoundCondTank(int numPhysique, tinyxml2::XMLElement *element, int &numberPhases, int &numberTransports, std::vector<std::string> nameTransports, Eos **eos, std::string fileName = "Fichier Inconnu");
    BoundCondTank(const BoundCondTank &Source, const int lvl = 0); //Constructeur de copie (utile pour AMR)
    virtual ~BoundCondTank();

    virtual void creeLimite(TypeMeshContainer<CellInterface *> &cellInterfaces);
    virtual void solveRiemannLimite(Cell &cellLeft, const int &numberPhases, const double &dxLeft, double &dtMax);
    virtual void solveRiemannTransportLimite(Cell &cellLeft, const int &numberTransports) const;

    virtual int whoAmI() const { return 5; };
    virtual void printInfo();

    //Pour methode AMR
    virtual void creerCellInterfaceChild();  /*!< Creer un child cell interface (non initialize) */

  protected:
  private:
    int m_numberPhase;
    double *m_ak0;
    double *m_Yk0;
    double *m_rhok0;
    double m_p0;
    double m_T0;
    int m_numberTransports;
    double *m_valueTransport;

};

#endif // BOUNDCONDTANK_H
