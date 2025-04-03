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

#ifndef QUANTITIESADDPHYS_H
#define QUANTITIESADDPHYS_H

//! \file      QuantitiesAddPhys.h
//! \author    K. Schmidmayer
//! \version   1.1
//! \date      June 5 2019

class QuantitiesAddPhys; //pre-declaration of QuantitiesAddPhys class needed for AddPhys.h inclusion

#include "AddPhys.h"

//! \class     QuantitiesAddPhys
//! \brief     General class for additional-physics quantities
//! \details   This is a pure virtual class: can not be instantiated
class QuantitiesAddPhys
{
    public:
    QuantitiesAddPhys();
    //! \brief     Generic model constructor
    //! \param     addPhys              corresponding additional physic
    QuantitiesAddPhys(AddPhys* addPhys);
    virtual ~QuantitiesAddPhys();

    //! \brief     Compute the needed quantities for the additional physic
    //! \param     cell                 corresponding cell
    virtual void computeQuantities(Cell* cell) { Errors::errorMessage("computeQuantities not implemented for used quantities of additional physics"); };
    //! \brief     Compute and send back mass energie linked to the physic (0 if no linked energy)
    double computeEnergyAddPhys();

    //Accessors for the different gradients
    //! \brief     Set the additional-physic gradient with the transmitted values
    //! \param     grad                 transmitted gradient
    //! \param     num                  number to determine the corresponding gradient
    virtual void setGrad(const Coord &grad, int num = -1) { Errors::errorMessage("setGrad not implemented for used quantities of additional physics"); };
    //! \brief     Get the additional-physic gradient
    //! \param     num                  number to determine the corresponding gradient
    virtual const Coord& getGrad(int num = -1) const { Errors::errorMessage("getGrad not implemented for used quantities of additional physics"); return Coord::defaultCoord; };
    //! \brief     Set the gradient of the velocity along the x-direction with the transmitted values
    //! \param     grad                 transmitted gradient
    virtual void setGradU(const Coord &grad) {};
    //! \brief     Set the gradient of the velocity along the y-direction with the transmitted values
    //! \param     grad                 transmitted gradient
    virtual void setGradV(const Coord &grad) {};
    //! \brief     Set the gradient of the velocity along the z-direction with the transmitted values
    //! \param     grad                 transmitted gradient
    virtual void setGradW(const Coord &grad) {};
    //! \brief     Set the gradient of the phase temperature with the transmitted values
    //! \param     phaseNum             number of the corresponding phase
    //! \param     grad                 transmitted gradient
    virtual void setGradTk(int &phaseNum, const Coord &grad) {};
    //! \brief     Return the gradient of the velocity along the x-direction
    virtual Coord getGradU() const { return 0; };
    //! \brief     Return the gradient of the velocity along the y-direction
    virtual Coord getGradV() const { return 0; };
    //! \brief     Return the gradient of the velocity along the z-direction
    virtual Coord getGradW() const { return 0; };
    //! \brief     Return the gradient of the phase temperature
    //! \param     phaseNum             number of the corresponding phase
    virtual const Coord& getGradTk(int &phaseNum) const { return Coord::defaultCoord; };

    //! \brief     Return the corresponding additional-physic class of this quantities class
    AddPhys* getAddPhys() { return m_addPhys; };

    // add by Fan E
    virtual const Coord& getGradVT(int num = -1) const { Errors::errorMessage("getGradVT not implemented for used quantities of additional physics"); return Coord::defaultCoord; };

    virtual const Coord& getGradRhoI(int num = -1) const { Errors::errorMessage("getGradVT not implemented for used quantities of additional physics"); return Coord::defaultCoord; };

    protected:
      AddPhys* m_addPhys;           //!< Corresponding additional-physic class of this quantities class

    private:
};

extern std::vector<Variable> variableNameSurfTens;  //!< Variable name of the corresponding gradient
extern std::vector<int> numPhaseSurfTens;           //!< Number of the phase (here transport)
extern std::vector<Variable> variableNamesVisc;     //!< Variable names of the corresponding gradients
extern std::vector<int> numPhasesVisc;              //!< Number of the phase (-1 for viscosity: works on mixture)
extern std::vector<Variable> variableNamesCond;     //!< Variable names of the corresponding gradients
extern std::vector<int> numPhasesCond;              //!< Number of the phase

#endif // QUANTITIESADDPHYS_H
