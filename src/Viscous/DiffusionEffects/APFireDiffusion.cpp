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

//! \file      APFireDiffusion.cpp
//! \author    E Fan, Lisong Shi, Jiaao Hao, Tianhan Zhang, Chih-yung Wen.
//! \version   1.0
//! \date      Jan 01 2025

#include <iostream>
#include <cmath>
#include <algorithm>
#include "APFireDiffusion.h"

using namespace tinyxml2;

//***********************************************************************

APFireDiffusion::APFireDiffusion(){}

//***********************************************************************

APFireDiffusion::APFireDiffusion(int& numberQPA, std::string fileName)
{
  // m_muk = new double[numberPhases];
  // for (int k = 0; k < numberPhases; k++) {
  //   m_muk[k] = eos[k]->getMu();
  // }
  m_numQPA = numberQPA++;
}

//***********************************************************************

APFireDiffusion::~APFireDiffusion(){  }


//***********************************************************************
// - fane, the size of addQuantitiesAddPhyc
void APFireDiffusion::addQuantityAddPhys(Cell *cell)
{
  cell->getVecQuantitiesAddPhys().push_back(new QAPFireDiffusion(this));
}

//***********************************************************************

void APFireDiffusion::solveFluxAddPhys(CellInterface *cellInterface, const int &numberPhases)
{
  double T_face, p_face;
  // Copy velocities and gradients of left and right cells
  // pointer shall be used for efficiency. - fnae
  m_velocityLeft = cellInterface->getCellGauche()->getMixture()->getVelocity();
  m_velocityRight = cellInterface->getCellDroite()->getMixture()->getVelocity();
  m_rhoArrayL = cellInterface->getCellGauche()->getMixture()->getDensityArray();
  m_rhoArrayR = cellInterface->getCellDroite()->getMixture()->getDensityArray();
  m_XArrayL = cellInterface->getCellGauche()->getMixture()->getMolarFractionArray();
  m_XArrayR = cellInterface->getCellDroite()->getMixture()->getMolarFractionArray();
  m_TemperatureL = cellInterface->getCellGauche()->getMixture()->getTemperature();
  m_TemperatureR = cellInterface->getCellDroite()->getMixture()->getTemperature();
  m_pressureL = cellInterface->getCellGauche()->getMixture()->getPressure();
  m_pressureR = cellInterface->getCellDroite()->getMixture()->getPressure();

  double distanceX = std::fabs(cellInterface->getCellGauche()->distanceX(cellInterface->getCellDroite()));
  double distanceY = std::fabs(cellInterface->getCellGauche()->distanceY(cellInterface->getCellDroite()));
  double distanceZ = std::fabs(cellInterface->getCellGauche()->distanceZ(cellInterface->getCellDroite()));
  double distance = cellInterface->getCellGauche()->distance(cellInterface->getCellDroite());

  m_normal = cellInterface->getFace()->getNormal();
  m_tangent = cellInterface->getFace()->getTangent();
  m_binormal = cellInterface->getFace()->getBinormal();

  T_face = 0.5*(m_TemperatureL+m_TemperatureR); // a better linear interpolation shall be used. - fnae
  p_face = 0.5*(m_pressureL+m_pressureR); // a better linear interpolation shall be used. - fnae

  m_gradULeft = cellInterface->getCellGauche()->getQPA(m_numQPA)->getGradVT(1);
  m_gradURight = cellInterface->getCellDroite()->getQPA(m_numQPA)->getGradVT(1);
  m_gradVLeft = cellInterface->getCellGauche()->getQPA(m_numQPA)->getGradVT(2);
  m_gradVRight = cellInterface->getCellDroite()->getQPA(m_numQPA)->getGradVT(2);
  m_gradWLeft = cellInterface->getCellGauche()->getQPA(m_numQPA)->getGradVT(3);
  m_gradWRight = cellInterface->getCellDroite()->getQPA(m_numQPA)->getGradVT(3);

  m_gradTLeft = cellInterface->getCellGauche()->getQPA(m_numQPA)->getGradVT(4);
  m_gradTRight = cellInterface->getCellDroite()->getQPA(m_numQPA)->getGradVT(4);

  m_gradPLeft = cellInterface->getCellGauche()->getQPA(m_numQPA)->getGradVT(5);
  m_gradPRight = cellInterface->getCellDroite()->getQPA(m_numQPA)->getGradVT(5);

  for(int k = 0 ;k < NS; k++)
  {
    m_gradRhoLeft[k] = cellInterface->getCellGauche()->getQPA(m_numQPA)->getGradRhoI(k+1);
    m_gradRhoRight[k] = cellInterface->getCellDroite()->getQPA(m_numQPA)->getGradRhoI(k+1);
  }

  for(int k = 0 ;k < NS; k++)
  {
    m_gradXLeft[k] = cellInterface->getCellGauche()->getQPA(m_numQPA)->getGradXI(k+1);
    m_gradXRight[k] = cellInterface->getCellDroite()->getQPA(m_numQPA)->getGradXI(k+1);
  }

  Coord S, dTds, correct, gradT_f, gradU_f, gradV_f, gradW_f, gradP_f;
  std::array<Coord,NS> gradRho_f;
  std::array<Coord,NS> gradX_f;

  double term1, term2;

  // S: the unit vector of neighbouring cells centres
  S.setX(distanceX/distance);
  S.setY(distanceY/distance);
  S.setZ(distanceZ/distance);

  gradT_f = 0.5* (m_gradTLeft + m_gradTRight);
  term1 = (m_TemperatureR - m_TemperatureL)/distance;
  term2 = gradT_f.getX()*S.getX() + gradT_f.getY()*S.getY() + gradT_f.getZ()*S.getZ();
  gradT_f += (term1 - term2)*S;

  gradP_f = 0.5* (m_gradPLeft + m_gradPRight);
  term1 = (m_pressureR - m_pressureL)/distance;
  term2 = gradP_f.getX()*S.getX() + gradP_f.getY()*S.getY() + gradP_f.getZ()*S.getZ();
  gradP_f += (term1 - term2)*S;

  gradU_f = 0.5* (m_gradULeft + m_gradURight);
  term1 = (m_velocityRight.getX() - m_velocityLeft.getX())/distance;
  term2 = gradU_f.getX()*S.getX() + gradU_f.getY()*S.getY() + gradU_f.getZ()*S.getZ();

  // printf("gradU, x = %le, y = %le.\n",gradU_f.getX(),gradU_f.getY());
  // printf("inner cell, uL = %le, uR = %le.\n",m_velocityLeft.getX(),m_velocityRight.getX());

  gradU_f += (term1 - term2)*S;  

  gradV_f = 0.5* (m_gradVLeft + m_gradVRight);
  term1 = (m_velocityRight.getY() - m_velocityLeft.getY())/distance;
  term2 = gradV_f.getX()*S.getX() + gradV_f.getY()*S.getY() + gradV_f.getZ()*S.getZ();
  gradV_f += (term1 - term2)*S;  

  gradW_f = 0.5* (m_gradWLeft + m_gradWRight);
  term1 = (m_velocityRight.getZ() - m_velocityLeft.getZ())/distance;
  term2 = gradW_f.getX()*S.getX() + gradW_f.getY()*S.getY() + gradW_f.getZ()*S.getZ();
  gradW_f += (term1 - term2)*S;  

  for(int i=0;i<NS;i++)
  {
    gradRho_f[i] = 0.5* (m_gradRhoLeft[i] + m_gradRhoRight[i]);
    term1 = (m_rhoArrayR[i] - m_rhoArrayL[i])/distance;
    term2 = gradRho_f[i].getX()*S.getX() + gradRho_f[i].getY()*S.getY() + gradRho_f[i].getZ()*S.getZ();
    gradRho_f[i] += (term1 - term2)*S; 
  }

  for(int i=0;i<NS;i++)
  {
    gradX_f[i] = 0.5* (m_gradXLeft[i] + m_gradXRight[i]);
    term1 = (m_XArrayR[i] - m_XArrayL[i])/distance;
    term2 = gradX_f[i].getX()*S.getX() + gradX_f[i].getY()*S.getY() + gradX_f[i].getZ()*S.getZ();
    gradX_f[i] += (term1 - term2)*S; 
  }

  gradT_f.localProjection(m_normal, m_tangent, m_binormal);
  gradU_f.localProjection(m_normal, m_tangent, m_binormal);
  gradV_f.localProjection(m_normal, m_tangent, m_binormal);
  gradW_f.localProjection(m_normal, m_tangent, m_binormal);
  for(int i=0;i<NS;i++){
    gradRho_f[i].localProjection(m_normal, m_tangent, m_binormal);
    gradX_f[i].localProjection(m_normal, m_tangent, m_binormal);
  }

  // // Compute the mixture mu on left and right
  // double muMixLeft(0.), muMixRight(0.);
  // for (int k = 0; k < numberPhases; k++) {
  //   muMixLeft += cellInterface->getCellGauche()->getPhase(k)->getAlpha()*m_muk[k];
  //   muMixRight += cellInterface->getCellDroite()->getPhase(k)->getAlpha()*m_muk[k];
  // }

  // get the transport properties on left and right
  for(int k = 0 ;k < NS+2; k++)
  {
    m_trans[k] = 0.5*(cellInterface->getCellGauche()->getTransportProperties()[k]
      + cellInterface->getCellDroite()->getTransportProperties()[k]);
  }

  double yL,yR;
  yL = cellInterface->getCellGauche()->getPosition().getY();
  yR = cellInterface->getCellDroite()->getPosition().getY();

  // Projection on orientation axes attached to the edge of velocities and gradients - fane
  m_velocityLeft.localProjection(m_normal, m_tangent, m_binormal);
  m_velocityRight.localProjection(m_normal, m_tangent, m_binormal);
  m_gradULeft.localProjection(m_normal, m_tangent, m_binormal);
  m_gradURight.localProjection(m_normal, m_tangent, m_binormal);
  m_gradVLeft.localProjection(m_normal, m_tangent, m_binormal);
  m_gradVRight.localProjection(m_normal, m_tangent, m_binormal);
  m_gradWLeft.localProjection(m_normal, m_tangent, m_binormal);
  m_gradWRight.localProjection(m_normal, m_tangent, m_binormal);
  m_gradTLeft.localProjection(m_normal, m_tangent, m_binormal);
  m_gradTRight.localProjection(m_normal, m_tangent, m_binormal);
  for(int k = 0 ;k < NS; k++)
  {
    m_gradRhoLeft[k].localProjection(m_normal, m_tangent, m_binormal);
    m_gradRhoRight[k].localProjection(m_normal, m_tangent, m_binormal);
  }

  this->solveFluxViscosityInner(m_velocityLeft, m_velocityRight, 
    m_rhoArrayL, m_rhoArrayR, m_XArrayL, m_XArrayR, gradU_f, gradU_f,
    gradV_f, gradV_f, gradW_f, gradW_f, gradT_f, gradT_f, gradP_f, gradP_f,
    gradRho_f, gradRho_f, gradX_f, gradX_f, m_trans,T_face, p_face, numberPhases, cellInterface);

  // Flux projection on the absolute orientation axes - fane
  cellInterface->getMod()->reverseProjection(m_normal, m_tangent, m_binormal);

  // printf("local, drho1dx = %le, dudx = %le, dudy = %le, dvdx = %le, dvdy = %le.\n",
  //   gradRho_f[0].getX(),gradU_f.getX(),gradU_f.getY(),gradV_f.getX(),gradV_f.getY());
  // printf("after reverseProjection, fluxBuff[0] = %le, [1] = %le, [u] = %le, [v] = %le, [E] = %le.\n",
  //   static_cast<FluxFire*> (fluxBuff)->m_masse[0],
  //   static_cast<FluxFire*> (fluxBuff)->m_masse[1],
  //   static_cast<FluxFire*> (fluxBuff)->m_qdm.getX(),
  //   static_cast<FluxFire*> (fluxBuff)->m_qdm.getY(),
  //   static_cast<FluxFire*> (fluxBuff)->m_energ);
}

//***********************************************************************

void APFireDiffusion::solveFluxAddPhysBoundary(CellInterface *cellInterface, const int &numberPhases)
{
  double T_face, p_face;
  double term1, term2;
  ////KS//DEV// BC Injection, Tank, Outflow to do

  // Copy velocities and gradients of left and right cells
  m_velocityLeft = cellInterface->getCellGauche()->getMixture()->getVelocity();
  m_rhoArrayL = cellInterface->getCellGauche()->getMixture()->getDensityArray();
  m_XArrayL = cellInterface->getCellGauche()->getMixture()->getMolarFractionArray();
  m_TemperatureL = cellInterface->getCellGauche()->getMixture()->getTemperature();
  m_pressureL = cellInterface->getCellGauche()->getMixture()->getPressure();

  T_face = m_TemperatureL;
  p_face = m_pressureL;

  m_gradULeft = cellInterface->getCellGauche()->getQPA(m_numQPA)->getGradVT(1);
  m_gradVLeft = cellInterface->getCellGauche()->getQPA(m_numQPA)->getGradVT(2);
  m_gradWLeft = cellInterface->getCellGauche()->getQPA(m_numQPA)->getGradVT(3);
  m_gradTLeft = cellInterface->getCellGauche()->getQPA(m_numQPA)->getGradVT(4);
  m_gradPLeft = cellInterface->getCellGauche()->getQPA(m_numQPA)->getGradVT(5);

  for(int k = 0 ;k < NS; k++)
  {
    m_gradRhoLeft[k] = cellInterface->getCellGauche()->getQPA(m_numQPA)->getGradRhoI(k+1);
  }

  // get the transport properties on left and right
  for(int k = 0 ;k < NS+2; k++)
  {
    m_trans[k] = cellInterface->getCellGauche()->getTransportProperties()[k];
  }

  m_normal = cellInterface->getFace()->getNormal();
  m_tangent = cellInterface->getFace()->getTangent();
  m_binormal = cellInterface->getFace()->getBinormal();

  // Projection on orientation axes attached to the edge of velocities and gradients
  m_velocityLeft.localProjection(m_normal, m_tangent, m_binormal);
  m_gradULeft.localProjection(m_normal, m_tangent, m_binormal);
  m_gradVLeft.localProjection(m_normal, m_tangent, m_binormal);
  m_gradWLeft.localProjection(m_normal, m_tangent, m_binormal);
  m_gradTLeft.localProjection(m_normal, m_tangent, m_binormal);
  for(int k = 0 ;k < NS; k++)
  {
    m_gradRhoLeft[k].localProjection(m_normal, m_tangent, m_binormal);
  }

  // Distances cells/cell interfaces for weighting on the flux
  double distLeft = cellInterface->getCellGauche()->distance(cellInterface);
  double distanceX = std::fabs(cellInterface->getCellGauche()->distanceX(cellInterface));
  double distanceY = std::fabs(cellInterface->getCellGauche()->distanceY(cellInterface));
  double distanceZ = std::fabs(cellInterface->getCellGauche()->distanceZ(cellInterface));

  int typeCellInterface = cellInterface->whoAmI();
  if (typeCellInterface == 3 || typeCellInterface == 1 || typeCellInterface == 4 || typeCellInterface == 5 ) 
  {
    // Cell interface of type NonReflecting, Outflow, Injection, Tank or 
    
    // this->solveFluxViscosityNonReflecting(m_velocityLeft, m_rhoArrayL, m_gradULeft, m_gradVLeft, m_gradWLeft, m_gradTLeft, m_gradRhoLeft, 
    //   m_trans, T_face, numberPhases,cellInterface);

    // initialize ghost cell velocity by reverse normal velocity, not this is in local coordinate
    Coord velocity_ghost = m_velocityLeft;

    Coord gradU_f,gradV_f, gradW_f, gradT_f, S, gradP_f;
    std::array<Coord,NS> gradRho_f;
    std::array<Coord,NS> gradX_f;

  // gradT_f = 0.5* (m_gradTLeft + m_gradTRight);
  // term1 = (m_TemperatureR - m_TemperatureL)/distance;
  // term2 = gradT_f.getX()*S.getX() + gradT_f.getY()*S.getY() + gradT_f.getZ()*S.getZ();
  // gradT_f += (term1 - term2)*S;

    gradU_f = 0;

    // term1 = (velocity_ghost.getX() - m_velocityLeft.getX())/(2*distLeft);
    // term2 = gradU_f.getX()*S.getX() + gradU_f.getY()*S.getY() + gradU_f.getZ()*S.getZ();
    // gradU_f += (term1 - term2)*S; 

    gradV_f = 0; // same value for left/ghost cell

    gradW_f = 0; // same value for left/ghost cell

    gradT_f = 0; // same value for left/ghost cell

    gradP_f = 0; // same value for left/ghost cell

    for(int i=0;i<NS;i++)
    {
      gradRho_f[i] = 0;
    }

    // this->solveFluxViscosityInner(m_velocityLeft, velocity_ghost, 
    //   m_rhoArrayL, m_rhoArrayL, gradU_f, gradU_f,
    //   gradV_f, gradV_f, gradW_f, gradW_f, gradT_f, gradT_f,
    //   gradRho_f, gradRho_f, m_trans,T_face, numberPhases,cellInterface);  

    this->solveFluxViscosityInner(m_velocityLeft, velocity_ghost, 
      m_rhoArrayL, m_rhoArrayL, m_XArrayL, m_XArrayL, gradU_f, gradU_f,
      gradV_f, gradV_f, gradW_f, gradW_f, gradT_f, gradT_f, gradP_f, gradP_f,
      gradRho_f, gradRho_f, gradX_f, gradX_f, m_trans,T_face, p_face, numberPhases, cellInterface);

  }
  else if(typeCellInterface == 9)
  {
    // m_gradULeft = m_gradULeft*2;
    // Coord grad_temp = 0;
    // for(int k = 0 ;k < NS; k++)
    // {
    //   m_gradRhoLeft[k]=0;
    // }
    // this->solveFluxViscosityNonReflecting(m_velocityLeft, m_rhoArrayL, grad_temp, grad_temp, grad_temp, grad_temp, m_gradRhoLeft, 
    //   m_trans, T_face, numberPhases);
    this->solveFluxViscosityNonReflecting(m_velocityLeft, m_rhoArrayL, m_gradULeft, m_gradVLeft, m_gradWLeft, m_gradTLeft, m_gradRhoLeft, 
      m_trans, T_face, numberPhases,cellInterface);
  }
  else if (typeCellInterface == 2)
  {
    // initialize ghost cell velocity by reverse normal velocity, not this is in local coordinate
    Coord velocity_ghost = m_velocityLeft;
    velocity_ghost.setX(-velocity_ghost.getX());

    Coord gradU_f,gradV_f, gradW_f, gradT_f, S, gradP_f;
    std::array<Coord,NS> gradRho_f;
    std::array<Coord,NS> gradX_f;

  // gradT_f = 0.5* (m_gradTLeft + m_gradTRight);
  // term1 = (m_TemperatureR - m_TemperatureL)/distance;
  // term2 = gradT_f.getX()*S.getX() + gradT_f.getY()*S.getY() + gradT_f.getZ()*S.getZ();
  // gradT_f += (term1 - term2)*S;

    gradU_f = m_gradULeft;
    term1 = (velocity_ghost.getX() - m_velocityLeft.getX())/(2*distLeft);
    term2 = gradU_f.getX()*S.getX() + gradU_f.getY()*S.getY() + gradU_f.getZ()*S.getZ();
    gradU_f += (term1 - term2)*S; 

    gradV_f = 0; // same value for left/ghost cell
    term1 = 0; //cancelled
    term2 = gradV_f.getX()*S.getX() + gradV_f.getY()*S.getY() + gradV_f.getZ()*S.getZ();
    gradV_f += (term1 - term2)*S; 

    gradW_f = 0; // same value for left/ghost cell
    term1 = 0; //cancelled
    term2 = gradW_f.getX()*S.getX() + gradW_f.getY()*S.getY() + gradW_f.getZ()*S.getZ();
    gradW_f += (term1 - term2)*S; 

    gradT_f = 0; // same value for left/ghost cell
    term1 = 0; //cancelled
    term2 = gradT_f.getX()*S.getX() + gradT_f.getY()*S.getY() + gradT_f.getZ()*S.getZ();
    gradT_f += (term1 - term2)*S; 

    gradP_f = 0; // same value for left/ghost cell
    term1 = 0; //cancelled
    term2 = gradP_f.getX()*S.getX() + gradP_f.getY()*S.getY() + gradP_f.getZ()*S.getZ();
    gradT_f += (term1 - term2)*S; 

    for(int i=0;i<NS;i++)
    {
      gradRho_f[i] = 0;
      term1 = 0; //cancelled
      term2 = gradRho_f[i].getX()*S.getX() + gradRho_f[i].getY()*S.getY() + gradRho_f[i].getZ()*S.getZ();
      gradRho_f[i] += (term1 - term2)*S; // correction
    }

    for(int i=0;i<NS;i++)
    {
      gradX_f[i] = 0;
      term1 = 0; //cancelled
      term2 = gradX_f[i].getX()*S.getX() + gradX_f[i].getY()*S.getY() + gradX_f[i].getZ()*S.getZ();
      gradX_f[i] += (term1 - term2)*S; // correction
    }

  // printf("gradU_f.x = %le, gradU_f.Y = %le, gradU_f.Z = %le.\n", gradU_f.getX(), gradU_f.getY(), gradU_f.getZ());
  // printf("gradV_f.x = %le, gradV_f.Y = %le, gradV_f.Z = %le.\n", gradV_f.getX(), gradV_f.getY(), gradV_f.getZ());
  // printf("gradUR.x = %le, gradUR.Y = %le, gradUR.Z = %le.\n", m_gradURight.getX(), m_gradURight.getY(), m_gradURight.getZ());
  // printf("avg.x = %le, avg.Y = %le, avg.Z = %le, term1 = %le, term2 = %le.\n", gradU_f.getX(), gradU_f.getY(), gradU_f.getZ(), term1, term2);
  // printf("In term2, UR = %le, UL = %le.\n", m_velocityRight.getX(), m_velocityLeft.getX());

    // this->solveFluxViscosityInner(m_velocityLeft, velocity_ghost, 
    //   m_rhoArrayL, m_rhoArrayL, gradU_f, gradU_f,
    //   gradV_f, gradV_f, gradW_f, gradW_f, gradT_f, gradT_f,
    //   gradRho_f, gradRho_f, m_trans,T_face, numberPhases,cellInterface);
      
    this->solveFluxViscosityInner(m_velocityLeft, velocity_ghost, 
      m_rhoArrayL, m_rhoArrayL, m_XArrayL, m_XArrayL, gradU_f, gradU_f,
      gradV_f, gradV_f, gradW_f, gradW_f, gradT_f, gradT_f, gradP_f, gradP_f,
      gradRho_f, gradRho_f, gradX_f, gradX_f, m_trans,T_face,p_face, numberPhases, cellInterface);


  // printf("dudx = %le, dudy = %le, dvdx = %le, dvdy = %le, dTdx = %le, u = %le, v = %le.\n",dudx,dudy,dvdx,dvdy, dTdx, u,v);
    // printf("fluxBuff[0] = %le, fluxBuff[1] = %le, [2] = %le, [u] = %le, [v] = %le, [E] = %le.\n",
    //   static_cast<FluxFire*> (fluxBuff)->m_masse[0],
    //   static_cast<FluxFire*> (fluxBuff)->m_masse[1],
    //   static_cast<FluxFire*> (fluxBuff)->m_masse[2],
    //   static_cast<FluxFire*> (fluxBuff)->m_qdm.getX(),
    //   static_cast<FluxFire*> (fluxBuff)->m_qdm.getY(),
    //   static_cast<FluxFire*> (fluxBuff)->m_energ);

    // //Writing of viscous terms on each equation of fluxBuffFire
    // for (int i = 0; i<NS; i++)
    // {
    //   static_cast<FluxFire*> (fluxBuff)->m_masse[i] = 0;
    // }
    // static_cast<FluxFire*> (fluxBuff)->m_qdm.setX( 0);
    // static_cast<FluxFire*> (fluxBuff)->m_qdm.setY( 0);
    // static_cast<FluxFire*> (fluxBuff)->m_qdm.setZ( 0);
    // static_cast<FluxFire*> (fluxBuff)->m_energ = 0;
  }
  else if (typeCellInterface == 6) // symmetry
  {
     // initialize ghost cell velocity by reverse normal velocity, not this is in local coordinate
    Coord velocity_ghost = m_velocityLeft;
    velocity_ghost.setX(-velocity_ghost.getX());

    Coord gradU_f,gradV_f, gradW_f, gradT_f, S, gradP_f;
    std::array<Coord,NS> gradRho_f;
    std::array<Coord,NS> gradX_f;

    S.setX(distanceX/distLeft);
    S.setY(distanceY/distLeft);
    S.setZ(distanceZ/distLeft);
    S.localProjection(m_normal, m_tangent, m_binormal);

  // gradT_f = 0.5* (m_gradTLeft + m_gradTRight);
  // term1 = (m_TemperatureR - m_TemperatureL)/distance;
  // term2 = gradT_f.getX()*S.getX() + gradT_f.getY()*S.getY() + gradT_f.getZ()*S.getZ();
  // gradT_f += (term1 - term2)*S;

    gradU_f = m_gradULeft;
    // term1 = (velocity_ghost.getX() - m_velocityLeft.getX())/(2*distLeft);
    // term2 = gradU_f.getX()*S.getX() + gradU_f.getY()*S.getY() + gradU_f.getZ()*S.getZ();
    // gradU_f += (term1 - term2)*S; 

    gradV_f = 0; // same value for left/ghost cell

    gradW_f = 0; // same value for left/ghost cell

    gradT_f = 0; // same value for left/ghost cell

    gradP_f = 0; // same value for left/ghost cell

    for(int i=0;i<NS;i++)
    {
      gradRho_f[i] = 0;
      gradX_f[i] = 0;
    }

  // printf("gradU_f.x = %le, gradU_f.Y = %le, gradU_f.Z = %le.\n", gradU_f.getX(), gradU_f.getY(), gradU_f.getZ());
  // printf("gradV_f.x = %le, gradV_f.Y = %le, gradV_f.Z = %le.\n", gradV_f.getX(), gradV_f.getY(), gradV_f.getZ());
  // printf("gradUR.x = %le, gradUR.Y = %le, gradUR.Z = %le.\n", m_gradURight.getX(), m_gradURight.getY(), m_gradURight.getZ());
  // printf("avg.x = %le, avg.Y = %le, avg.Z = %le, term1 = %le, term2 = %le.\n", gradU_f.getX(), gradU_f.getY(), gradU_f.getZ(), term1, term2);
  // printf("In term2, UR = %le, UL = %le.\n", m_velocityRight.getX(), m_velocityLeft.getX());

    // this->solveFluxViscosityInner(m_velocityLeft, velocity_ghost, 
    //   m_rhoArrayL, m_rhoArrayL, gradU_f, gradU_f,
    //   gradV_f, gradV_f, gradW_f, gradW_f, gradT_f, gradT_f,
    //   gradRho_f, gradRho_f, m_trans,T_face, numberPhases,cellInterface);  

    this->solveFluxViscosityInner(m_velocityLeft, velocity_ghost, 
      m_rhoArrayL, m_rhoArrayL, m_XArrayL, m_XArrayL, gradU_f, gradU_f,
      gradV_f, gradV_f, gradW_f, gradW_f, gradT_f, gradT_f, gradP_f, gradP_f,
      gradRho_f, gradRho_f, gradX_f, gradX_f, m_trans,T_face,p_face, numberPhases, cellInterface);
  }
  else { 
    printf("Viscous B.C. is not managed.\n");
    exit(EXIT_FAILURE);
    // this->solveFluxViscosityOther(m_velocityLeft, m_gradULeft, m_gradVLeft, m_gradWLeft, muMixLeft, numberPhases); 
  }

  // Flux projection on the absolute orientation axes
  cellInterface->getMod()->reverseProjection(m_normal, m_tangent, m_binormal);

  // printf("dudx = %le, dudy = %le, dvdx = %le, dvdy = %le, dTdx = %le, u = %le, v = %le.\n",dudx,dudy,dvdx,dvdy, dTdx, u,v);
  // printf("after reverseProjection, fluxBuff[0] = %le, [1] = %le, [u] = %le, [v] = %le, [E] = %le.\n",
  //   static_cast<FluxFire*> (fluxBuff)->m_masse[0],
  //   static_cast<FluxFire*> (fluxBuff)->m_masse[1],
  //   static_cast<FluxFire*> (fluxBuff)->m_qdm.getX(),
  //   static_cast<FluxFire*> (fluxBuff)->m_qdm.getY(),
  //   static_cast<FluxFire*> (fluxBuff)->m_energ);
}

//***********************************************************************
void APFireDiffusion::solveFluxViscosityInner(Coord &velocityLeft, Coord &velocityRight, 
      std::array<double,NS> &rhoArrayL, std::array<double,NS> &rhoArrayR, 
      std::array<double,NS> &XArrayL, std::array<double,NS> &XArrayR, 
      Coord &gradULeft, Coord &gradURight, 
      Coord &gradVLeft, Coord &gradVRight, 
      Coord &gradWLeft, Coord &gradWRight, 
      Coord &gradTLeft, Coord &gradTRight, 
      Coord &gradPLeft, Coord &gradPRight, 
      std::array<Coord,NS> &gradRhoLeft, std::array<Coord,NS> &gradRhoRight, 
      std::array<Coord,NS> &gradXLeft, std::array<Coord,NS> &gradXRight, 
      std::array<double,NS+2> &trans, double T_face, double P_face, 
      int numberPhases, CellInterface *cellInterface) const
{
  // cout<<"In solveFluxViscosityInner, cellInterface cellInterface->getCellGauche()->x = "<<
  //   cellInterface->getCellGauche()->getPosition().getX()<<endl;

  // string a = cellInterface->getCellGauche()->getSymType();
  // cout<<"a = "<<a<<endl;

  // cout<<"In solveFluxViscosityInner, cellInterface cellInterface->getCellGauche()->getSymType()"<<
  //   cellInterface->getCellGauche()->getSymType()<<endl;
  

  double muMix, kMix, DC[NS];
  muMix = trans[0];// mu_mix
  kMix = trans[1]; // k_mix
  for(int i=0;i<NS;i++) DC[i] = trans[i+2]; // DC_i

  double tau_xx(0.);

	//Extraction of data
	double uL, vL, wL, uR, vR, wR;
	double dudxL, dudyL, dudzL, dudxR, dudyR, dudzR;
  double dTdx, dPdx;
  double dvdxL, dvdyL, dvdxR, dvdyR;
  double dwdxL, dwdzL, dwdxR, dwdzR;
  double drdxL[NS], drdyL[NS], drdzL[NS];
  double drdxR[NS], drdyR[NS], drdzR[NS];
  double drdx[NS],drhoMixdx(0.);

	uL = velocityLeft.getX();
	vL = velocityLeft.getY();
	wL = velocityLeft.getZ();
	uR = velocityRight.getX();
	vR = velocityRight.getY();
	wR = velocityRight.getZ();

	dudxL = gradULeft.getX();
  dudyL = gradULeft.getY();
  dudzL = gradULeft.getZ();
  dudxR = gradURight.getX();
  dudyR = gradURight.getY();
  dudzR = gradURight.getZ();

  dvdxL = gradVLeft.getX();
  dvdyL = gradVLeft.getY();
  dvdxR = gradVRight.getX();
  dvdyR = gradVRight.getY();

  dwdxL = gradWLeft.getX();
  dwdzL = gradWLeft.getZ();
  dwdxR = gradWRight.getX();
  dwdzR = gradWRight.getZ();

  dTdx = 0.5*(gradTLeft.getX()+gradTRight.getX());
  dPdx = 0.5*(gradPLeft.getX()+gradPRight.getX());

  for(int i=0;i<NS;i++)
  {
    drdx[i] = 0.5*(gradRhoLeft[i].getX()+gradRhoRight[i].getX());
    drhoMixdx += drdx[i];
  }

	//Data of the cell interface
  double u, v, w;
  double dudx, dudy, dudz;
  double dvdx, dvdy;
  double dwdx, dwdz;
  // double muMel;

	u = (uL + uR) / 2.;
	v = (vL + vR) / 2.;
	w = (wL + wR) / 2.;

	dudx = (dudxL + dudxR) / 2.;
  dudy = (dudyL + dudyR) / 2.;
  dudz = (dudzL + dudzR) / 2.;

  dvdx = (dvdxL + dvdxR) / 2.;
  dvdy = (dvdyL + dvdyR) / 2.;

  dwdx = (dwdxL + dwdxR) / 2.;
  dwdz = (dwdzL + dwdzR) / 2.;


  // muMel = (muMixLeft + muMixRight) / 2.;

  double rhoMixL(0.), rhoMixR(0.), YsL[NS], YsR[NS],YsM[NS];
  for(int i=0;i<NS;i++){rhoMixL+=rhoArrayL[i];rhoMixR+=rhoArrayR[i];}

  double JAX[NS],JA(0.), thermo[5], JSX[NS];

  double XsM[NS];
  double rhosM[NS];
  double rhoMix_sM(0.);

  // calculate face averaged values
  for(int i=0;i<NS;i++)
  {
    rhosM[i] = 0.5*(rhoArrayL[i] + rhoArrayR[i]);
    rhoMix_sM += rhosM[i];
  }
  for(int i=0;i<NS;i++){YsM[i] = rhosM[i]/rhoMix_sM;}

  std::array<double, NS> molarFractions;
  double totalMoles = 0.0;
  for (int i = 0; i < NS; ++i) {
      totalMoles += rhosM[i] / Sptr[i].MW;
  }
  for (int i = 0; i < NS; ++i) {
      XsM[i] = (rhosM[i] / Sptr[i].MW) / totalMoles;
  }

  for(int i=0;i<NS;i++)
  {
    if(XsM[i] > 1e-15)
    {
      JAX[i] = - DC[i]*(drdx[i] - YsM[i]*drhoMixdx) - rhosM[i]*DC[i]/XsM[i]*(XsM[i]-YsM[i])*dPdx/P_face;
      JA += JAX[i];
    }
    else
    {
      JAX[i] = 0.0; // avoid division by zero
    }
    for(int i=0;i<NS;i++)JSX[i] = JAX[i] - YsM[i]*JA;
  }
  
	//Writing of viscous terms on each equation of fluxBuffFire
	for (int i = 0; i<NS; i++)
	{
    static_cast<FluxFire*> (fluxBuff)->m_masse[i] = JSX[i];
	}
  if(cellInterface->getCellGauche()->getSymType()=="CYLINDRICAL")
  {
    if(cellInterface->getFace()->getPos().getY()>1e-10)
    tau_xx = dudx + dvdy + v/cellInterface->getFace()->getPos().getY() + dwdz;
  }
  else
    tau_xx = dudx + dvdy + dwdz;
  // static_cast<FluxFire*> (fluxBuff)->m_qdm.setX( -muMix / 3. * (4.*dudx - 2.*(dvdy + dwdz)));
  static_cast<FluxFire*> (fluxBuff)->m_qdm.setX( -muMix * ( 2*dudx - 2./3.*tau_xx));
  static_cast<FluxFire*> (fluxBuff)->m_qdm.setY( -muMix * (dvdx + dudy));
  static_cast<FluxFire*> (fluxBuff)->m_qdm.setZ( -muMix * (dwdx + dudz));
  double energyFlux;
  // energyFlux = -muMix * (4./3.*dudx*u + (dvdx + dudy)*v + (dwdx + dudz)*w - 2./3.*(dvdy + dwdz)*u);
  energyFlux = -muMix * (( 2*dudx - 2./3.*tau_xx)*u + (dvdx + dudy)*v + (dwdx + dudz)*w);
  energyFlux += -kMix * dTdx;

  // // old method of calling function cphs()
  // for(int i=0;i<NS;i++) 
  // {
  //   cphs(thermo, i,T_face,"VisF");
  //   energyFlux += (JSX[i]*thermo[2]);
  // }

  // call the method in Cantera
  std::shared_ptr<Cantera::ThermoPhase> gas = solution_mech->thermo();
  double Y[NS] = {0};
  for(int i=0;i<NS;i++)
  {
    Y[i] = 1.0;
    gas->setState_TPY(T_face, 101325, Y);
    energyFlux += (JSX[i]*gas->enthalpy_mass());
    Y[i] = 0.0;
  }

  static_cast<FluxFire*> (fluxBuff)->m_energ = energyFlux;
  
  // printf("x = %le, [u] = %le, [v] = %le, [E] = %le.\n",
  //   cellInterface->getFace()->getPos().getX(),
  //   static_cast<FluxFire*> (fluxBuff)->m_qdm.getX(),
  //   static_cast<FluxFire*> (fluxBuff)->m_qdm.getY(),
  //   static_cast<FluxFire*> (fluxBuff)->m_energ);

  // for(int i=0;i<NS;i++)
  // {
  //   printf("fluxBuff[%d]L = %le.\n", i, static_cast<FluxFire*> (fluxBuff)->m_masse[i]);
  // }

  // printf("Input uL = %le, vL = %le, wL = %le, dudxL = %le, dudyL = %le.\n", uL, vL, wL, dudxL, dudyL);
  // printf("pL = %le, TL = %le.\n", cellInterface->getCellGauche()->getMixture()->getPressure(), cellInterface->getCellGauche()->getMixture()->getTemperature());
  // for(int i=0;i<NS;i++)
  // {
  //   printf("rhoL[%d] = %le.\n", i, cellInterface->getCellGauche()->getMixture()->getDensityArray()[i]);
  // }
  // printf("Input dvdxL = %le, dvdyL = %le, dTdx = %le.\n", dvdxL, dvdyL, dTdx);
  // printf("muMix = %le, kMix = %le.\n", muMix, kMix);
  // for(int i=0;i<NS;i++)
  // {
  //   printf("drho[%d]dxL = %le.\n", i, gradRhoLeft[i].getX());
  // }
}

//***********************************************************************

void APFireDiffusion::solveFluxViscosityNonReflecting(Coord &velocityLeft, 
      std::array<double,NS> &rhoArrayL,
      Coord &gradULeft, 
      Coord &gradVLeft, Coord &gradWLeft, 
      Coord &gradTLeft, std::array<Coord,NS> &gradRhoLeft,
      std::array<double,NS+2> &trans, double T_face, int numberPhases, CellInterface *cellInterface) const
{
  std::array<double,NS> rhoArrayR = rhoArrayL;
  Coord velocityRight = velocityLeft;
  Coord gradURight = gradULeft;
  Coord gradVRight = gradVLeft;
  Coord gradWRight = gradWLeft;
  Coord gradTRight = gradTLeft;
  std::array<Coord,NS> gradRhoRight = gradRhoLeft;

  double muMix, kMix, DC[NS];
  muMix = trans[0];// mu_mix
  kMix = trans[1]; // k_mix
  for(int i=0;i<NS;i++) DC[i] = trans[i+2]; // DC_i

  double tau_xx(0.);

	//Extraction of data
	double uL, vL, wL, uR, vR, wR;
	double dudxL, dudyL, dudzL, dudxR, dudyR, dudzR;
  double dTdx;
  double dvdxL, dvdyL, dvdxR, dvdyR;
  double dwdxL, dwdzL, dwdxR, dwdzR;
  double drdxL[NS], drdyL[NS], drdzL[NS];
  double drdxR[NS], drdyR[NS], drdzR[NS];
  double drdx[NS],drhoMixdx(0.);

	uL = velocityLeft.getX();
	vL = velocityLeft.getY();
	wL = velocityLeft.getZ();
	uR = velocityRight.getX();
	vR = velocityRight.getY();
	wR = velocityRight.getZ();

	dudxL = gradULeft.getX();
  dudyL = gradULeft.getY();
  dudzL = gradULeft.getZ();
  dudxR = gradURight.getX();
  dudyR = gradURight.getY();
  dudzR = gradURight.getZ();

  dvdxL = gradVLeft.getX();
  dvdyL = gradVLeft.getY();
  dvdxR = gradVRight.getX();
  dvdyR = gradVRight.getY();

  dwdxL = gradWLeft.getX();
  dwdzL = gradWLeft.getZ();
  dwdxR = gradWRight.getX();
  dwdzR = gradWRight.getZ();

  dTdx = 0.5*(gradTLeft.getX()+gradTRight.getX());

  for(int i=0;i<NS;i++)
  {
    drdx[i] = 0.5*(gradRhoLeft[i].getX()+gradRhoRight[i].getX());
    drhoMixdx += drdx[i];
  }

	//Data of the cell interface
  double u, v, w;
  double dudx, dudy, dudz;
  double dvdx, dvdy;
  double dwdx, dwdz;
  // double muMel;

	u = (uL + uR) / 2.;
	v = (vL + vR) / 2.;
	w = (wL + wR) / 2.;

	dudx = (dudxL + dudxR) / 2.;
  dudy = (dudyL + dudyR) / 2.;
  dudz = (dudzL + dudzR) / 2.;

  dvdx = (dvdxL + dvdxR) / 2.;
  dvdy = (dvdyL + dvdyR) / 2.;

  dwdx = (dwdxL + dwdxR) / 2.;
  dwdz = (dwdzL + dwdzR) / 2.;


  // muMel = (muMixLeft + muMixRight) / 2.;

  double rhoMixL(0.), rhoMixR(0.), YsL[NS], YsR[NS],YsM[NS];
  for(int i=0;i<NS;i++){rhoMixL+=rhoArrayL[i];rhoMixR+=rhoArrayR[i];}
  for(int i=0;i<NS;i++){YsM[i] = 0.5*(rhoArrayL[i]/rhoMixL + rhoArrayR[i]/rhoMixR);}
  
  double JAX[NS],JA(0.), thermo[5], JSX[NS];
  for(int i=0;i<NS;i++)
  {
    JAX[i] = DC[i]*(drdx[i] - YsM[i]*drhoMixdx); // notes that rho gradient isn't 0
    JA+=JAX[i];
  }

  // for(int i=0;i<NS;i++)
  // {
  //   JAX[i] = DC[i]*drdx[i];
  //   JA+=JAX[i];
  // }

  for(int i=0;i<NS;i++)JSX[i] = -JAX[i]+YsM[i]*JA;
	//Writing of viscous terms on each equation of fluxBuffFire
	for (int i = 0; i<NS; i++)
	{
    static_cast<FluxFire*> (fluxBuff)->m_masse[i] = JSX[i];
	}
  if(cellInterface->getCellGauche()->getSymType()=="CYLINDRICAL")
  {
    if(cellInterface->getFace()->getPos().getY()>1e-10)
    tau_xx = dudx + dvdy + v/cellInterface->getFace()->getPos().getY() + dwdz;
  }
  else
    tau_xx = dudx + dvdy + dwdz;
  // static_cast<FluxFire*> (fluxBuff)->m_qdm.setX( -muMix / 3. * (4.*dudx - 2.*(dvdy + dwdz)));
  static_cast<FluxFire*> (fluxBuff)->m_qdm.setX( -muMix * ( 2*dudx - 2./3.*tau_xx));
  static_cast<FluxFire*> (fluxBuff)->m_qdm.setY( -muMix * (dvdx + dudy));
  static_cast<FluxFire*> (fluxBuff)->m_qdm.setZ( -muMix * (dwdx + dudz));
  double energyFlux;
  // energyFlux = -muMix * (4./3.*dudx*u + (dvdx + dudy)*v + (dwdx + dudz)*w - 2./3.*(dvdy + dwdz)*u);
  energyFlux = -muMix * (( 2*dudx - 2./3.*tau_xx)*u + (dvdx + dudy)*v + (dwdx + dudz)*w);
  energyFlux += -kMix * dTdx;

  // call the method in Cantera
  std::shared_ptr<Cantera::ThermoPhase> gas = solution_mech->thermo();
  double Y[NS] = {0};
  for(int i=0;i<NS;i++)
  {
    Y[i] = 1.0;
    gas->setState_TPY(T_face, 101325, Y);
    energyFlux += (JSX[i]*gas->enthalpy_mass());
    Y[i] = 0.0;
  }

  static_cast<FluxFire*> (fluxBuff)->m_energ = energyFlux;

  // printf("x = %le, [u] = %le, [v] = %le, [E] = %le.\n",
  //   cellInterface->getFace()->getPos().getX(),
  //   static_cast<FluxFire*> (fluxBuff)->m_qdm.getX(),
  //   static_cast<FluxFire*> (fluxBuff)->m_qdm.getY(),
  //   static_cast<FluxFire*> (fluxBuff)->m_energ);

  // for(int i=0;i<NS;i++)
  // {
  //   printf("fluxBuff[%d]L = %le.\n", i, static_cast<FluxFire*> (fluxBuff)->m_masse[i]);
  // }

  // printf("Input uL = %le, vL = %le, wL = %le, dudxL = %le, dudyL = %le.\n", uL, vL, wL, dudxL, dudyL);
  // printf("pL = %le, TL = %le.\n", cellInterface->getCellGauche()->getMixture()->getPressure(), cellInterface->getCellGauche()->getMixture()->getTemperature());
  // for(int i=0;i<NS;i++)
  // {
  //   printf("rhoL[%d] = %le.\n", i, cellInterface->getCellGauche()->getMixture()->getDensityArray()[i]);
  // }
  // printf("Input dvdxL = %le, dvdyL = %le, dTdx = %le.\n", dvdxL, dvdyL, dTdx);
  // printf("muMix = %le, kMix = %le.\n", muMix, kMix);
  // for(int i=0;i<NS;i++)
  // {
  //   printf("drho[%d]dxL = %le.\n", i, gradRhoLeft[i].getX());
  // }
}

//***********************************************************************

void APFireDiffusion::solveFluxViscositySymmetry(Coord &velocityLeft, 
      std::array<double,NS> &rhoArrayL,
      Coord &gradULeft, 
      Coord &gradVLeft, Coord &gradWLeft, 
      Coord &gradTLeft, std::array<Coord,NS> &gradRhoLeft,
      std::array<double,NS+2> &trans, double T_face, int numberPhases, CellInterface *cellInterface)
{

    //Extraction of data
    //Note that dudy, dudz, dvdx and dwdx are nulls
    // double dudx, dvdy, dwdz;
    // dudx = gradULeft.getX();
    // dvdy = gradVLeft.getY();
    // dwdz = gradWLeft.getZ();

    double u = velocityLeft.getX();
    double v = velocityLeft.getY();

    m_normal = cellInterface->getFace()->getNormal();
    m_tangent = cellInterface->getFace()->getTangent();
    m_binormal = cellInterface->getFace()->getBinormal();

    Coord cellLeft_center = cellInterface->getCellGauche()->getElement()->getPosition();
    Coord cellInterface_center = cellInterface->getFace()->getPos();
    Coord cell2face_vec = cellInterface_center - cellLeft_center;
    Coord d = cell2face_vec;
    double distance = d.squaredNorm();

    velocityLeft.localProjection(m_normal, m_tangent, m_binormal);

    double ux = velocityLeft.getX(); // velocity magnitude in normal direction

    double uxnx = ux*m_normal.getX();
    double uxny = ux*m_normal.getY();
    double uxnz = ux*m_normal.getZ();

    double dudx, dudy;


    for (int i = 0; i<NS; i++)
    {
      static_cast<FluxFire*> (fluxBuff)->m_masse[i] = 0;
    }

    double muMix = trans[0];// mu_mix

    //Writing of viscous terms on each equation of fluxBuffEuler
    // static_cast<FluxFire*> (fluxBuff)->m_qdm.setX(-muMix / 3. * (4. * dudx - 2. * (dvdy + dwdz)));
    
    // follow method in Moukalled's book, Eq. 15.149
    // set symmetry velocity to 0, then gradU = 0 - uxnx
    static_cast<FluxFire*> (fluxBuff)->m_qdm.setX(- muMix * 4./3. * (- uxnx )/distance   );
    static_cast<FluxFire*> (fluxBuff)->m_qdm.setX(- muMix * -(uxny) / distance);
    static_cast<FluxFire*> (fluxBuff)->m_qdm.setY(0.);
    static_cast<FluxFire*> (fluxBuff)->m_qdm.setZ(0.);
    static_cast<FluxFire*> (fluxBuff)->m_energ = (- muMix * 4./3. * (- uxnx )/distance   )*u + (- muMix * -(uxny) / distance)*v  ;

    //   double energyFlux;
    // // energyFlux = -muMix * (4./3.*dudx*u + (dvdx + dudy)*v + (dwdx + dudz)*w - 2./3.*(dvdy + dwdz)*u);
    // energyFlux = -muMix * (( 2*dudx - 2./3.*tau_xx)*u + (dvdx + dudy)*v + (dwdx + dudz)*w);

    // incorrect 0 flux
    // static_cast<FluxFire*> (fluxBuff)->m_qdm.setX(0.);
    // static_cast<FluxFire*> (fluxBuff)->m_qdm.setY(0.);
    // static_cast<FluxFire*> (fluxBuff)->m_qdm.setZ(0.);
    // static_cast<FluxFire*> (fluxBuff)->m_energ = 0.;

}

//***********************************************************************

void APFireDiffusion::solveFluxViscosityWall(Coord &velocityLeft, double &muMixLeft, double &distLeft, int numberPhases) const
{
  printf("APFireDiffusion::solveFluxViscosityWall() invalid.\n");
  exit(EXIT_FAILURE);
  double dudx = -velocityLeft.getX() / distLeft;
  double dvdx = -velocityLeft.getY() / distLeft;
  double dwdx = -velocityLeft.getZ() / distLeft;

  // for (int k = 0; k<numberPhases; k++)
  // {
  //   static_cast<FluxFire*> (fluxBuff)->m_alpha[k] = 0.;
  //   static_cast<FluxFire*> (fluxBuff)->m_masse[k] = 0.;
  //   static_cast<FluxFire*> (fluxBuff)->m_energ[k] = 0.;
  // }
  static_cast<FluxFire*> (fluxBuff)->m_qdm.setX(-muMixLeft / 3. * 4. * dudx );
  static_cast<FluxFire*> (fluxBuff)->m_qdm.setY(-muMixLeft * dvdx);
  static_cast<FluxFire*> (fluxBuff)->m_qdm.setZ(-muMixLeft * dwdx);
  static_cast<FluxFire*> (fluxBuff)->m_energ = 0.;
}

//***********************************************************************

void APFireDiffusion::solveFluxViscosityOther(Coord &velocityLeft, Coord &gradULeft, Coord &gradVLeft, Coord &gradWLeft, double &muMixLeft, int numberPhases) const
{
  //Not manage at the moment, just an example
  std::cout << "Viscous boundary is not managed." << std::endl;

  exit(EXIT_FAILURE);

  // To avoid bug when not manage
  // for (int k = 0; k<numberPhases; k++) {
  //   static_cast<FluxFire*> (fluxBuff)->m_alpha[k] = 0.;
  //   static_cast<FluxFire*> (fluxBuff)->m_masse[k] = 0.;
  //   static_cast<FluxFire*> (fluxBuff)->m_energ[k] = 0.;
  // }
  static_cast<FluxFire*> (fluxBuff)->m_qdm = 0.;
  static_cast<FluxFire*> (fluxBuff)->m_energ = 0.;
}

//***********************************************************************

void APFireDiffusion::addNonCons(Cell *cell, const int &numberPhases)
{
  // double dudx = cell->getQPA(m_numQPA)->getGradVT(1).getX();
  // double dudy = cell->getQPA(m_numQPA)->getGradVT(1).getY();
  // double dudz = cell->getQPA(m_numQPA)->getGradVT(1).getZ();
  // double dvdx = cell->getQPA(m_numQPA)->getGradVT(2).getX();
  // double dvdy = cell->getQPA(m_numQPA)->getGradVT(2).getY();
  // double dvdz = cell->getQPA(m_numQPA)->getGradVT(2).getZ();
  // double dwdx = cell->getQPA(m_numQPA)->getGradVT(3).getX();
  // double dwdy = cell->getQPA(m_numQPA)->getGradVT(3).getY();
  // double dwdz = cell->getQPA(m_numQPA)->getGradVT(3).getZ();
  // double termeNonCons = - 2./3.*(dudx + dvdy + dwdz)*(dudx + dvdy + dwdz)
  //                       + 2.*(dudx*dudx + dvdy*dvdy + dwdz*dwdz)
  //                       + (dudy+dvdx)*(dudy+dvdx) + (dudz+dwdx)*(dudz+dwdx) + (dvdz+dwdy)*(dvdz+dwdy);

  // for (int k = 0; k<numberPhases; k++) {
  //   static_cast<FluxFire*> (fluxBuff)->m_alpha[k] = 0.;
  //   static_cast<FluxFire*> (fluxBuff)->m_masse[k] = 0.;
  //   static_cast<FluxFire*> (fluxBuff)->m_energ[k] = cell->getPhase(k)->getAlpha()*m_muk[k] * termeNonCons;
  // }
  // static_cast<FluxFire*> (fluxBuff)->m_qdm = 0.;
  // static_cast<FluxFire*> (fluxBuff)->m_energ = 0.;

  // cell->getCons()->addFlux(1., numberPhases);
}

//***********************************************************************

void APFireDiffusion::addSymmetricTermsRadialAxisOnX(Cell *cell, const int &numberPhases)
{
  printf("APFireDiffusion::addSymmetricTermsRadialAxisOnX() not available.\n");
  exit(EXIT_FAILURE);
  // //Extraction of data
  // double r = cell->getPosition().getX();
  // double dudx = cell->getQPA(m_numQPA)->getGradVT(1).getX();
  // double dudy = cell->getQPA(m_numQPA)->getGradVT(1).getY();
  // double dvdx = cell->getQPA(m_numQPA)->getGradVT(2).getX();
  // double dvdy = cell->getQPA(m_numQPA)->getGradVT(2).getY();
  // double dTdx = cell->getQPA(m_numQPA)->getGradVT(4).getX();
  // double dTdy = cell->getQPA(m_numQPA)->getGradVT(4).getY();
  // double drdx[NS],drdy[NS];
  // for(int i=0;i<NS;i++)
  // {
  //   drdx[i] = cell->getQPA(m_numQPA)->getGradRhoI(i).getX();
  //   drdy[i] = cell->getQPA(m_numQPA)->getGradRhoI(i).getY();
  // }

  // //Compute the mixture mu
  // double muMix(0.);
  // // for (int k = 0; k < numberPhases; k++) {
  // //   muMix += cell->getPhase(k)->getAlpha()*m_muk[k];
  // // }

  // //Writing of symmetrical viscous terms on each equation of fluxBuffKapila
  // // for (int k = 0; k<numberPhases; k++) {
  // //   static_cast<FluxFire*> (fluxBuff)->m_alpha[k] = 0.;
  // //   static_cast<FluxFire*> (fluxBuff)->m_masse[k] = 0.;
  // //   static_cast<FluxFire*> (fluxBuff)->m_energ[k] = 0.;
  // // }
  // static_cast<FluxFire*> (fluxBuff)->m_qdm.setX(muMix * 2. * dudx / r);
  // static_cast<FluxFire*> (fluxBuff)->m_qdm.setY(muMix * (dvdx + dudy) / r);
  // static_cast<FluxFire*> (fluxBuff)->m_energ = muMix * (1./3.*(4.*dudx - 2.*dvdy)*cell->getMixture()->getU() + (dudy+dvdx)*cell->getMixture()->getV()) / r;

  // cell->getCons()->addFlux(1., numberPhases);
}

//***********************************************************************
// the [trans,T] value shall get from cell instead of cellInterface
// and all the gradients are also cell value
void APFireDiffusion::addSymmetricTermsRadialAxisOnY(Cell *cell, const int &numberPhases)
{
  std::array<double,NS+2> trans;
  trans = cell->getTransportProperties();
  double muMix, kMix, DC[NS];
  muMix = trans[0];// mu_mix
  kMix = trans[1]; // k_mix
  for(int i=0;i<NS;i++) DC[i] = m_trans[i+2]; // DC_i

  // extract data
  double r = cell->getPosition().getY();
  double dudx = cell->getQPA(m_numQPA)->getGradVT(1).getX();//u
  double dudy = cell->getQPA(m_numQPA)->getGradVT(1).getY();
  double dvdx = cell->getQPA(m_numQPA)->getGradVT(2).getX();//v
  double dvdy = cell->getQPA(m_numQPA)->getGradVT(2).getY();
  double dTdx = cell->getQPA(m_numQPA)->getGradVT(4).getX();//T
  double dTdy = cell->getQPA(m_numQPA)->getGradVT(4).getY();
  double drdx[NS],drdy[NS];
  for(int i=0;i<NS;i++)
  {
    drdx[i] = cell->getQPA(m_numQPA)->getGradRhoI(i+1).getX();
    drdy[i] = cell->getQPA(m_numQPA)->getGradRhoI(i+1).getY();
  }

  double u, v, T, rhoMix(0.), YsM[NS];
  std::array<double,NS> rhoArray;

  rhoArray = cell->getMixture()->getDensityArray();
  u = cell->getMixture()->getU();
  v = cell->getMixture()->getV();
  T = cell->getMixture()->getTemperature();

  for(int i=0;i<NS;i++)rhoMix+=rhoArray[i];
  for(int i=0;i<NS;i++)YsM[i]=rhoArray[i]/rhoMix;

  double JAY[NS],JA(0.), JSY[NS],Qr;
  for(int i=0;i<NS;i++){JAY[i] = DC[i]*drdy[i]; JA+=JAY[i];}
  for(int i=0;i<NS;i++)JSY[i] = -(JAY[i]-YsM[i]*JA);
  Qr = -kMix*dTdy;

  // for(int i=0;i<NS;i++)
  // {
  //   static_cast<FluxFire*> (fluxBuff)->m_masse[i] = -JSY[i]/r;
  // }

  // static_cast<FluxFire*> (fluxBuff)->m_qdm.setX(1/3*muMix/r*dvdx + muMix/r*dudy);
  // static_cast<FluxFire*> (fluxBuff)->m_qdm.setY(4/3*muMix/r*dvdy - 4/3*muMix/r/r*v);
  // //4/3*muMix/r*dvdy - 4/3*muMix/r/r*v-2/3*muMix/r*dvdx)
  // static_cast<FluxFire*> (fluxBuff)->m_energ = 1/3*dvdx*muMix*u/r + muMix*u/r*dudy 
  //   - 4/3*muMix*v/r*dudx- Qr/r;
  // double energy_Mass(0.), thermo[5];
  // for(int i=0;i<NS;i++)
  // {
  //   cphs(thermo, i,T,"AxiV");
  //   energy_Mass += -JSY[i]*thermo[2]/r;
  // }

  for(int i=0;i<NS;i++)
  {
    static_cast<FluxFire*> (fluxBuff)->m_masse[i] = 0.;
  }

  static_cast<FluxFire*> (fluxBuff)->m_qdm.setX(0.);
  static_cast<FluxFire*> (fluxBuff)->m_qdm.setY((-2*muMix*v/r+2/3*muMix*(dudx+dvdy+v/r))/r);
  //4/3*muMix/r*dvdy - 4/3*muMix/r/r*v-2/3*muMix/r*dvdx)
  static_cast<FluxFire*> (fluxBuff)->m_energ = 0.;

  // printf("rank = %d, fluxBuff[1] = %le, fluxBuff[2] = %le.\n",rankCpu,-JSY[1]/r,-JSY[2]/r);

  // static_cast<FluxFire*> (fluxBuff)->m_qdm.setX(muMix * (dvdx + dudy) / r);
  // static_cast<FluxFire*> (fluxBuff)->m_qdm.setY(muMix * 2. * dvdy / r);
  // static_cast<FluxFire*> (fluxBuff)->m_energ = muMix * (1./3.*(4.*dvdy - 2.*dudx)*cell->getMixture()->getV() + (dudy+dvdx)*cell->getMixture()->getU()) / r;
  
  cell->getCons()->addFlux(1., numberPhases);
}

//***********************************************************************

void APFireDiffusion::communicationsAddPhys(int numberPhases, const int &dim, const int &lvl)
{
  for(int i=1;i<=NS+5+NS;i++)//1:U, 2:V, 3:W, 4:T, 5:p, 6->NS+6:rho_i, NS+7->2NS+6: X_i(molar fraction)
  {
    parallel.communicationsVector(QPA, dim, lvl, m_numQPA, i); //m_grad_i
  }
  MPI_Barrier(MPI_COMM_WORLD);
  parallel.communicationsTransportArray(transportArray, dim, lvl, m_numQPA, -1); //m_grad_i
  MPI_Barrier(MPI_COMM_WORLD);
	// parallel.communicationsVector(QPA, dim, lvl, m_numQPA, 1); //m_gradU
	// parallel.communicationsVector(QPA, dim, lvl, m_numQPA, 2); //m_gradV
	// parallel.communicationsVector(QPA, dim, lvl, m_numQPA, 3); //m_gradW
}

//***********************************************************************
