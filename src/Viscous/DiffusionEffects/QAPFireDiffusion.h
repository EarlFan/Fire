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


#ifndef QAPFIREVISCOSITY_H
#define QAPFIREVISCOSITY_H

//! \file      QAPFireDiffusion.h
//! \author    E Fan, Lisong Shi, Jiaao Hao, Tianhan Zhang, Chih-yung Wen.
//! \version   1.0
//! \date      Jan 01 2025

#include "../../libs/ECOGEN/AdditionalPhysics/QuantitiesAddPhys.h"

//! \class     QAPFireDiffusion
//! \brief     General class for viscous quantities
class QAPFireDiffusion : public QuantitiesAddPhys
{
    public:
    QAPFireDiffusion();
    QAPFireDiffusion(AddPhys* addPhys);
    virtual ~QAPFireDiffusion();

    virtual void computeQuantities(Cell* cell);

    //Accessors
    virtual void setGrad(const Coord &grad, int num = -1);                       //1:U, 2:V, 3:W
    virtual const Coord& getGrad(int num = -1) const; //1:U, 2:V, 3:W, 4:T, 5-14:species
    virtual const Coord& getGradVT(int num = -1) const { return m_gradsVT[num-1]; }; //1:U, 2:V, 3:W, 4:T
    virtual const Coord& getGradRhoI(int num = -1) const { return m_gradsRhoI[num-1]; }; //1-10 species

    protected:
    
    std::vector<double> m_transportProperties; // !<array of transport properties    
    std::vector<Coord> m_gradsVT;                   //!< Gradient vectors of the velocities of the cell in x-, y- and z-directions and T   
    std::vector<Coord> m_gradsRhoI;                   //!< Gradient vectors of the density array of the cell 

    private:
};

#endif // QAPFIREVISCOSITY_H
