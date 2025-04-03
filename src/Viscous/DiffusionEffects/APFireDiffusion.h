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


#ifndef APKFIREVISCOSITY_H
#define APKFIREVISCOSITY_H

//! \file      APFireDiffusion.h
//! \author    E Fan, Lisong Shi, Jiaao Hao, Tianhan Zhang, Chih-yung Wen.
//! \version   1.0
//! \date      Jan 01 2025

#include "../APFire.h"
#include "QAPFireDiffusion.h"
#include "../../libs/ECOGEN/Eos/Eos.h"

//! \class     APFireDiffusion
//! \brief     General class for thermal viscosity for the Kapila model
class APFireDiffusion : public APFire
{
  public:
    APFireDiffusion();
    APFireDiffusion(int& numberQPA, std::string nameFile = "Unknown file");
    virtual ~APFireDiffusion();

    virtual void addQuantityAddPhys(Cell *cell);

    virtual void solveFluxAddPhys(CellInterface *cellInterface, const int &numberPhases);
    virtual void solveFluxAddPhysBoundary(CellInterface *cellInterface, const int &numberPhases);
    //! \brief     Solve the viscosity flux between two cells
    //! \param     velocityLeft         velocity of the left cell
    //! \param     velocityRight        velocity of the right cell
    //! \param     gradULeft            gradient of the velocity in the x-direction of the left cell
    //! \param     gradURight           gradient of the velocity in the x-direction of the right cell
    //! \param     gradVLeft            gradient of the velocity in the y-direction of the left cell
    //! \param     gradVRight           gradient of the velocity in the y-direction of the right cell
    //! \param     gradWLeft            gradient of the velocity in the z-direction of the left cell
    //! \param     gradWRight           gradient of the velocity in the z-direction of the right cell
    //! \param     muMixLeft            dynamic viscosity of the mixture of the left cell
    //! \param     muMixRight           dynamic viscosity of the mixture of the right cell
    //! \param     numberPhases         number of phases
    void solveFluxViscosityInner(Coord &velocityLeft, Coord &velocityRight, 
      std::array<double,NS> &rhoArrayL, std::array<double,NS> &rhoArrayR, 
      Coord &gradULeft, Coord &gradURight, 
      Coord &gradVLeft, Coord &gradVRight, Coord &gradWLeft, Coord &gradWRight, 
      Coord &gradTLeft, Coord &gradTRight, std::array<Coord,NS> &gradRhoLeft, 
      std::array<Coord,NS> &gradRhoRight, std::array<double,NS+2> &m_trans, 
      double m_T_face, int numberPhases, CellInterface *cellInterface) const;
    //! \brief     Solve the viscosity flux at a boundary with an non-reflecting type
    //! \param     velocityLeft         velocity of the left cell
    //! \param     gradULeft            gradient of the velocity in the x-direction of the left cell
    //! \param     gradVLeft            gradient of the velocity in the y-direction of the left cell
    //! \param     gradWLeft            gradient of the velocity in the z-direction of the left cell
    //! \param     muMixLeft            dynamic viscosity of the mixture of the left cell
    //! \param     numberPhases         number of phases
    void solveFluxViscosityNonReflecting(Coord &velocityLeft, 
      std::array<double,NS> &rhoArrayL,
      Coord &gradULeft, 
      Coord &gradVLeft, Coord &gradWLeft, 
      Coord &gradTLeft, std::array<Coord,NS> &gradRhoLeft,
      std::array<double,NS+2> &m_trans, double m_T_face, int numberPhases, CellInterface *cellInterface) const;
    //! \brief     Solve the viscosity flux at a boundary with an wall type
    //! \param     velocityLeft         velocity of the left cell
    //! \param     muMixLeft            dynamic viscosity of the mixture of the left cell
    //! \param     distLeft             distance between the center of the left cell and its corresponding edge
    //! \param     numberPhases         number of phases
    void solveFluxViscosityWall(Coord &velocityLeft, double &muMixLeft, double &distLeft, int numberPhases) const;
    //! \brief     Solve the viscosity flux at a boundary with non-defined type yet
    //! \param     velocityLeft         velocity of the left cell
    //! \param     gradULeft            gradient of the velocity in the x-direction of the left cell
    //! \param     gradVLeft            gradient of the velocity in the y-direction of the left cell
    //! \param     gradWLeft            gradient of the velocity in the z-direction of the left cell
    //! \param     muMixLeft            dynamic viscosity of the mixture of the left cell
    //! \param     numberPhases         number of phases
    void solveFluxViscosityOther(Coord &velocityLeft, Coord &gradULeft, Coord &gradVLeft, Coord &gradWLeft, double &muMixLeft, int numberPhases) const;
    virtual void addNonCons(Cell *cell, const int &numberPhases);
    virtual void addSymmetricTermsRadialAxisOnX(Cell *cell, const int &numberPhases);
    virtual void addSymmetricTermsRadialAxisOnY(Cell *cell, const int &numberPhases);

    virtual void communicationsAddPhys(int numberPhases, const int &dim, const int &lvl);

    void solveFluxViscositySymmetry(Coord &velocityLeft, 
      std::array<double,NS> &rhoArrayL,
      Coord &gradULeft, Coord &gradVLeft, Coord &gradWLeft, 
      Coord &gradTLeft, std::array<Coord,NS> &gradRhoLeft,
      std::array<double,NS+2> &m_trans, double m_T_face, int numberPhases, CellInterface *cellInterface);

  protected:
  
  private:
    // double *m_muk;            //!< Dynamic viscosity (kg/m/s or Pa.s) of each phase (taken from the EOS classes) (buffer)
    std::array<double,NS+2> m_trans;
    // std::array<double,NS+2> m_transd;
    int m_numQPA;             //!< Number of the associated variable for each cell (m_vecGrandeursAddPhys)

    Coord m_velocityLeft;     //!< Left velocity vector for the flux computation (buffer)
    Coord m_gradULeft;        //!< Left gradient of the velocity in the x-direction for the flux computation (buffer)
    Coord m_gradVLeft;        //!< Left gradient of the velocity in the y-direction for the flux computation (buffer)
    Coord m_gradWLeft;        //!< Left gradient of the velocity in the z-direction for the flux computation (buffer)
    Coord m_velocityRight;    //!< Right velocity vector for the flux computation (buffer)
    Coord m_gradURight;       //!< Right gradient of the velocity in the x-direction for the flux computation (buffer)
    Coord m_gradVRight;       //!< Right gradient of the velocity in the y-direction for the flux computation (buffer)
    Coord m_gradWRight;       //!< Right gradient of the velocity in the z-direction for the flux computation (buffer)
    
    std::array<double,NS> m_rhoArrayL;
    double m_TemperatureL;
    Coord m_gradTLeft;
    std::array<Coord,NS> m_gradRhoLeft;
    std::array<double,NS> m_rhoArrayR;
    double m_TemperatureR;
    Coord m_gradTRight;
    std::array<Coord,NS> m_gradRhoRight;

    Coord m_normal;           //!< Normal vector of the corresponding face for the flux computation (buffer)
    Coord m_tangent;          //!< Tangent vector of the corresponding face for the flux computation (buffer)
    Coord m_binormal;         //!< Binormal vector of the corresponding face for the flux computation (buffer)
};

#endif // APKFIREVISCOSITY_H
