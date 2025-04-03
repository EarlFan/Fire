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

#ifndef BOUNDCOND_H
#define BOUNDCOND_H

//! \file      BoundCond.h
//! \author    F. Petitpas, K. Schmidmayer
//! \version   1.1
//! \date      June 5 2019

#include <iostream>

#include "../Order1/CellInterface.h"
#include "../Order2/CellInterfaceO2.h" //Ajouter pour l'AMR, a priori ne pose pas de probleme
#include "../libTierces/tinyxml2.h"
#include "../Errors.h"
#include "../Tools.h"

class BoundCond : public CellInterface
{
  public:
    BoundCond();
    BoundCond(int numPhysique);
    BoundCond(const BoundCond &Source);
    virtual ~BoundCond();

    virtual void creeLimite(TypeMeshContainer<CellInterface *> &cellInterfaces){ Errors::errorMessage("Impossible de creer la limite dans creeLimite"); };
    virtual void creeLimite(TypeMeshContainer<CellInterface *> &cellInterfaces, std::string ordreCalcul) { Errors::errorMessage("Impossible de creer la limite dans creeLimite"); };
    virtual void initialize(Cell *cellLeft, Cell *cellRight);

    virtual void computeFlux(const int &numberPhases, const int &numberTransports, double &dtMax, Limiter &globalLimiter, Limiter &interfaceLimiter, Limiter &globalVolumeFractionLimiter, Limiter &interfaceVolumeFractionLimiter, Prim type = vecPhases);
    virtual void computeFluxAddPhys(const int &numberPhases, AddPhys &addPhys);
    virtual void solveRiemann(const int &numberPhases, const int &numberTransports, double &dtMax, Limiter &globalLimiter, Limiter &interfaceLimiter, Limiter &globalVolumeFractionLimiter, Limiter &interfaceVolumeFractionLimiter, Prim type = vecPhases);
    virtual void addFlux(const int &numberPhases, const int &numberTransports, const double &coefAMR, Prim type = vecPhases) {};  //Ici la fonction ne fait rien car il s agit d une limite a droite et il n y a rien a ajouter a droite.
    virtual void solveRiemannLimite(Cell &cellLeft, const int &numberPhases, const double &dxLeft, double &dtMax) { Errors::errorMessage("Attention solveRiemannLimite non prevu pour limite utilisee"); };
    virtual void solveRiemannTransportLimite(Cell &cellLeft, const int &numberTransports) const { Errors::errorMessage("Attention solveRiemannTransportLimite non prevu pour limite utilisee"); };

    virtual int whoAmI() const { Errors::errorMessage("whoAmI pas prevu pour la limite demandee"); return 0; };
    virtual void printInfo(){};

    virtual const int& getNumPhys() const { return m_numPhysique; };

    //Pour methode AMR
    virtual void computeXi(double lvl, double lvlHydro, double lvlChem, double lvlGeo, double T_threshhold, const double &criteriaVar, const bool &varRho, const bool &varP, const bool &varU, const bool &varT, const bool &varH2O, const bool &varAlpha) {};
    virtual void computeFluxXi();
    virtual void raffineCellInterfaceExterne(const int &nbCellsY, const int &nbCellsZ, const double &dXParent, const double &dYParent, const double &dZParent, Cell *cellRef, const int &dim);
    virtual void deraffineCellInterfaceExterne(Cell *cellRef);

  protected:
    int m_numPhysique;

  private:
};

#endif // BOUNDCOND_H
