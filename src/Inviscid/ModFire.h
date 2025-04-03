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

#ifndef MODFIRE_H
#define MODFIRE_H

//! \file      ModFire.h
//! \author    E Fan, Lisong Shi, Jiaao Hao, Tianhan Zhang, Chih-yung Wen.
//! \version   1.0
//! \date      Jan 01 2025

#include "../libs/ECOGEN/Models/Model.h"
#include "../libs/ECOGEN/Order1/Cell.h"
#include "../Mixture/MixFire.h"

class ModFire;

#include "../Mixture/FluxFire.h"

//! \class     ModFire
//! \brief     Model class for mechanical equilibrium multiphase flows
class ModFire : public Model
{
  public:
    //! \brief     Kapila model constructor
    //! \param     numberTransports    number of additional transport equations
    //! \param     numberPhases        number of phases
    ModFire(int &numberTransports, const int &numberPhases);
    virtual ~ModFire();

    virtual void allocateCons(Flux **cons, const int &numberPhases);
    virtual void allocatePhase(Phase **phase);
    virtual void allocateMixture(Mixture **mixture);

    //! \details    Complete multiphase mechanical equilibrium state from volume fractions, pressure, densities, velocity
    virtual void fulfillState(Phase **phases, Mixture *mixture, const int &numberPhases, Prim type = vecPhases);

    //Hydrodynamic Riemann solvers
    //----------------------------
    virtual void solveRiemannIntern(Cell &cellLeft, Cell &cellRight, const int &numberPhases, const double &dxLeft, const double &dxRight, double &dtMax) const; // Riemann between two computed cells
    virtual void solveRiemannWall(Cell &cellLeft, const int &numberPhases, const double &dxLeft, double &dtMax) const; // Riemann between left cell and wall

    void solveRiemannSymmetryInner(Cell &cellLeft, const int &numberPhases, const double &dxLeft, double &dtMax) const;
    
    void solveRiemannMixture(Mixture &mixtureLeft, Mixture &mixtureRight, const int &numberPhases, const double &dxLeft, const double &dxRight, double &dtMax) const; // Riemann between two computed cells

    void solveRiemannInternHLL(Cell &cellLeft, Cell &cellRight, const int &numberPhases, const double &dxLeft, const double &dxRight, double &dtMax) const;
    void solveRiemannInternHLLC(Cell &cellLeft, Cell &cellRight, const int &numberPhases, const double &dxLeft, const double &dxRight, double &dtMax) const;
    void solveRiemannInternHLLC_LM(Cell &cellLeft, Cell &cellRight, const int &numberPhases, const double &dxLeft, const double &dxRight, double &dtMax) const;

    virtual void reverseProjection(const Coord normal, const Coord tangent, const Coord binormal) const;

    //Accessors
    //---------
    virtual const double& getSM();
    virtual const Coord& getVelocity(const Cell *cell) const { return cell->getMixture()->getVelocity(); };
    virtual Coord& getVelocity(Cell *cell) { return cell->getMixture()->getVelocity(); };

    virtual const std::string& whoAmI() const { return m_name; };

  protected:
  
  private:
    static const std::string NAME;

    friend class FluxFire;
    friend class APFireDiffusion;
};

#endif // MODFIRE_H
