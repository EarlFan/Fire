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

#ifndef FLUXFIRE_H
#define FLUXFIRE_H

//! \file      FluxFire.h
//! \author    E Fan, Lisong Shi, Jiaao Hao, Tianhan Zhang, Chih-yung Wen.
//! \version   1.0
//! \date      Jan 01 2025

#include <iostream>
#include "../libs/ECOGEN/Models/Flux.h"

class FluxFire;

#include "../Inviscid/ModFire.h"

//! \class     FluxFire
//! \brief     Model class for for Kapila system of equations (mechanical equilibrium) flux
class FluxFire : public Flux
{
  public:
    FluxFire();
    FluxFire(ModFire *model); 
    virtual ~FluxFire();

    virtual void printFlux() const;
    virtual void addFlux(double coefA, const int &numberPhases);
    virtual void addFlux(Flux* flux, const int& numberPhases);
    virtual void subtractFlux(double coefA, const int &numberPhases);
    virtual void multiply(double scalar, const int &numberPhases);
    virtual void setBufferFlux(Cell &cell, const int &numberPhases);
    virtual void setBufferFlux2();
    virtual void buildCons(Phase **phases, const int &numberPhases, Mixture *mixture);
    virtual void buildPrim(Phase **phases, Mixture *mixture, const int &numberPhases);
    virtual void setToZero(const int &numberPhases);
    virtual void setToZeroBufferFlux(const int &numberPhases);
    // do nothing in Fire - FanE
    virtual void addNonCons(double coefA, const Cell *cell, const int &numberPhases){};
    virtual void subtractNonCons(double coefA, const Cell *cell, const int &numberPhases){};
    virtual void correctionEnergy(Cell *cell, const int &numberPhases, Prim type = vecPhases) const{};

    virtual void addSymmetricTerms(Phase **phases, Mixture *mixture, const int &numberPhases, const double &r, const double &v);

    // Accessors
    //----------
    virtual const Coord& getQdm() const { return m_qdm; };
    virtual const double& getMasseMix() const {return Errors::defaultDouble; };
    virtual const double& getEnergyMix() const { return m_energ; };
    virtual void setCons(const Flux *cons, const int &numberPhases);

    virtual std::array<double,NS> getMassArray() const { return m_masse;}

    virtual void prepSourceTermsReaction(Cell *cell, const double &dt, const int &numberPhases, std::shared_ptr<Cantera::Solution> sol_) ;

protected:
    std::array<double,NS> m_masse;          //!< mass array - fane
    Coord m_qdm;              //!< momentum array
    double m_energ;          //!< specific internal energy array
    ModFire *m_model;       //!< associated model

  private:

    friend class ModFire;
    // To modify if needed, example: to add a class APKViscosity, add friend class APKViscosity.
    friend class APFireDiffusion;

};

#endif // FLUXFIRE_H
