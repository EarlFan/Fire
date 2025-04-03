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

//! \file      Cell.cpp
//! \author    E Fan, Lisong Shi, Jiaao Hao, Tianhan Zhang, Chih-yung Wen.
//! \version   1.0
//! \date      Jan 01 2025

#include "../libs/ECOGEN/Order1/CellInterface.h"
#include <iostream>

void CellInterface::computeXi(double lvl, double lvlHydro, double lvlChem, double lvlGeo, double T_threshhold, const double &criteriaVar, const bool &varRho, const bool &varP, const bool &varU, const bool &varT, const bool &varH2O, const bool &varAlpha)
{
  if(AMRPara::test_prob == 4) // for /example_cases/2D_Channel_detonation
  {
    // cells in the left domain are not refined after 140 us
    double interface_X = this->getFace()->getPos().getX();
    double velocity_deton = 1400;
    double coarse_init_x = 0.176;
    if(AMRPara::physical_time < 100e-6)
    {
      if(lvl>=2) return;
    }
    else
    {
      // turn-off left domain refinement after 100 us for efficiency
      if(interface_X<(velocity_deton*(AMRPara::physical_time-100e-6)+coarse_init_x))
      {
        return;
      }
    }
  }

  if(AMRPara::test_prob == 5) // for /example_cases/3D_channel_detonation
  {
    // cells in the left domain are not refined after 140 us
    double interface_X = this->getFace()->getPos().getX();
    double velocity_deton = 1400;
    double coarse_init_x = 0.076;
    // turn-off left domain refinement after 100 us for efficiency
    if(interface_X<(velocity_deton*(AMRPara::physical_time-100e-6)+coarse_init_x))
    {
      return;
    }
  }

  if (varRho) { 
    if (m_cellLeft->getXi() < 0.99 || m_cellRight->getXi() < 0.99) this->computeCritereAMR(lvl, lvlHydro, lvlChem, lvlGeo, T_threshhold, criteriaVar, density, -1); 
  }
  if (varP) {
    if (m_cellLeft->getXi() < 0.99 || m_cellRight->getXi() < 0.99) { this->computeCritereAMR(lvl, lvlHydro, lvlChem, lvlGeo, T_threshhold, criteriaVar, pressure); }
  }
  if (varU) {
    if (m_cellLeft->getXi() < 0.99 || m_cellRight->getXi() < 0.99) { this->computeCritereAMR(lvl, lvlHydro, lvlChem, lvlGeo, T_threshhold, criteriaVar, velocityMag); }
  }
  if (varAlpha) {
    if (m_cellLeft->getXi() < 0.99 || m_cellRight->getXi() < 0.99) { this->computeCritereAMR(lvl, lvlHydro, lvlChem, lvlGeo, T_threshhold, criteriaVar, alpha, 1); }
  }
  if (varT) {
    if (m_cellLeft->getXi() < 0.99 || m_cellRight->getXi() < 0.99) { this->computeCritereAMR(lvl, lvlHydro, lvlChem, lvlGeo, T_threshhold, criteriaVar, temperature); }
  }
  if (varH2O) {
    if (m_cellLeft->getXi() < 0.99 || m_cellRight->getXi() < 0.99) { this->computeCritereAMR(lvl, lvlHydro, lvlChem, lvlGeo, T_threshhold, criteriaVar, Y_OH); }
  }
  if (AMRPara::div_vor_v_flag == 1)
  {
    if (m_cellLeft->getXi() < 0.99 || m_cellRight->getXi() < 0.99) { this->computeCritereAMR(lvl, lvlHydro, lvlChem, lvlGeo, T_threshhold, criteriaVar, div_vor_v); }
  }
  if (AMRPara::gradRho_flag == 1)
  {
    if (m_cellLeft->getXi() < 0.99 || m_cellRight->getXi() < 0.99) { this->computeCritereAMR(lvl, lvlHydro, lvlChem, lvlGeo, T_threshhold, criteriaVar, gradRhoAMR); }
  }
  if (AMRPara::YH2O >= 0)
  {
    if (m_cellLeft->getXi() < 0.99 || m_cellRight->getXi() < 0.99) { 
      this->computeCritereAMR(lvl, lvlHydro, lvlChem, lvlGeo, T_threshhold, criteriaVar, density, AMRPara::YH2O); 
    }
  }
  if (AMRPara::YCO2 >= 1)
  {
    if (m_cellLeft->getXi() < 0.99 || m_cellRight->getXi() < 0.99) { 
      this->computeCritereAMR(lvl, lvlHydro, lvlChem, lvlGeo, T_threshhold, criteriaVar, density, AMRPara::YCO2); 
    }
  }
  if (AMRPara::YOH >= 1)
  {
    if (m_cellLeft->getXi() < 0.99 || m_cellRight->getXi() < 0.99) { 
      this->computeCritereAMR(lvl, lvlHydro, lvlChem, lvlGeo, T_threshhold, criteriaVar, density, AMRPara::YOH); 
    }
  }
}

//***********************************************************************
void CellInterface::computeCritereAMR(double lvl, double lvlHydro, double lvlChem, double lvlGeo, double T_threshhold, const double &criteriaVar, Variable nameVariable, int num)
{
  double valueMin, variation, cd, cg;
  double H2O_g, H2O_d, YH2O_g, YH2O_d;
  double rhog = m_cellLeft->selectScalar(density, -1);
  double rhod = m_cellRight->selectScalar(density, -1);
  double Td = m_cellLeft->selectScalar(temperature, num);
  double Tg = m_cellRight->selectScalar(temperature, num);

  // apply AMR strategy 2 to limit the hydrodynamic level
  if(lvlChem > lvlHydro)
  {
    if(Td < T_threshhold and Tg < T_threshhold and lvl <= lvlHydro)
    {
      return;
    }
  }
  
  if (nameVariable == gradRhoAMR)
  {
    double tau_gradRho_l, tau_gradRho_r;
    tau_gradRho_l = m_cellLeft->calcTauGradRho();
    tau_gradRho_r = m_cellRight->calcTauGradRho();

    if(tau_gradRho_l >= AMRPara::f3*AMRPara::sigma_gradRho || tau_gradRho_r >= AMRPara::f3*AMRPara::sigma_gradRho)
    // if(gradRho_l >= AMRPara::f3 || gradRho_r >= AMRPara::f3)
    {
      m_cellLeft->setXi(2.);
      m_cellRight->setXi(2.);
    }
  }
  else  // local gradient
  {
    // special treatment of total density in FIRE model
    if(m_cellLeft->getModel()->whoAmI()=="FIRE")
    {
      if(nameVariable==density)
      {
        if(num == -1)
        {
          cg = m_cellLeft->selectScalar(nameVariable, -1);
          cd = m_cellRight->selectScalar(nameVariable,-1);
        }
        else // mass fraction of {num}_{th} species
        {
          cg = m_cellLeft->selectScalar(nameVariable, num)/rhog;
          cd = m_cellRight->selectScalar(nameVariable, num)/rhod;
        }
      }
      else
      {
        cg = m_cellLeft->selectScalar(nameVariable, num);
        cd = m_cellRight->selectScalar(nameVariable, num);
      }
    }
    else
    {
      printf("This computeCritereAMR only support FIRE.\n");
      exit(EXIT_FAILURE);
    }

    // Valeur de la variation
    valueMin = std::min(std::fabs(cd), std::fabs(cg));
    if (valueMin < 1.e-2) { //Utile pour alpha (quasi-seulement) ou velocity
      if (nameVariable == velocityMag) { valueMin = 0.1; }
      else {                             valueMin = 1.e-2; }
    }
    variation = std::fabs(cd - cg) / valueMin;

    //Mise a jour de xi si la variation est superieure au criteria
    if (variation >= criteriaVar) {
      m_cellLeft->setXi(1.);
      m_cellRight->setXi(1.);
    }
  } // end local gradient

}
