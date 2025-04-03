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

//! \file      Run.cpp
//! \author    E Fan, Lisong Shi, Jiaao Hao, Tianhan Zhang, Chih-yung Wen.
//! \version   1.0
//! \date      Jan 01 2025

#include "../libs/ECOGEN/Run.h"

using namespace tinyxml2;

//***********************************************************************

void Run::calcGradRhoInCell()
{
  calcTotalLeafCell();

  total_gradRho = 0;
  double tau_rho, cell_l, f;

  for(int lvl = m_lvlMax; lvl>=0; lvl--)
  {
    for (unsigned int i = 0; i < m_cellsLvl[lvl].size(); i++) 
    {
      if (!m_cellsLvl[lvl][i]->getSplit())
      {
        tau_rho = m_cellsLvl[lvl][i]->calcTauGradRho();
        total_gradRho += pow(tau_rho, 2);
      }
    }
  }

  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Allreduce(&total_gradRho, &total_gradRho_para, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  
  AMRPara::sigma_gradRho = sqrt(total_gradRho_para/total_leaf_cells_para);

}

//***********************************************************************

void Run::calcTotalLeafCell()
{
  total_leaf_cells = 0;

  for(int lvl = m_lvlMax; lvl>=0; lvl--)
  {
    for (unsigned int i = 0; i < m_cellsLvl[lvl].size(); i++) 
    {
      if (!m_cellsLvl[lvl][i]->getSplit())
      {        
        total_leaf_cells++;
      }
    }
  }

  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Allreduce(&total_leaf_cells, &total_leaf_cells_para, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
}

//***********************************************************************

// reset the MixFire::m_densityArray_change, m_hf_change to zero
void Run::resetDelteRhoHInCell()
{
  for(int lvl = m_lvlMax; lvl>=0; lvl--)
  {
    for (unsigned int i = 0; i < m_cellsLvl[lvl].size(); i++) {
      if (!m_cellsLvl[lvl][i]->getSplit()){
        // reset the MixFire::m_densityArray_change, m_hf_change to zero
        m_cellsLvl[lvl][i]->getMixture()->setDensityArrayChangeToZero();
        m_cellsLvl[lvl][i]->getMixture()->setHfChangeToZero();
      }
    }
  }
}

// calculate MixFire::m_densityArray_change, m_hf_change in the mixture
void Run::delteRhoHInCell(const double &dt_lvl_0)
{
  for(int lvl = m_lvlMax; lvl>=0; lvl--)
  {
    for (unsigned int i = 0; i < m_cellsLvl[lvl].size(); i++) {
      if (!m_cellsLvl[lvl][i]->getSplit()){
        // calculate the rate of species density and enthalpy due to reaction in cell
        m_cellsLvl[lvl][i]->getMixture()->setDensityArrayChangeRate(dt_lvl_0);
        m_cellsLvl[lvl][i]->getMixture()->setHfChangeRate(dt_lvl_0);
      }
    }
  }
}

//***********************************************************************
// calculate the maximum pressure and temperature
void Run::getPeakValues()
{
  Tmax = -1e30;
  Pmax = -1e30;
  Tmin = 1e30;
  Pmin = 1e30;

  double T_temp, p_temp;
  for(int lvl = m_lvlMax; lvl>=0; lvl--)
  {
    for (unsigned int i = 0; i < m_cellsLvl[lvl].size(); i++) 
    {
      T_temp = m_cellsLvl[lvl][i]->getMixture()->getTemperature();
      p_temp = m_cellsLvl[lvl][i]->getMixture()->getPressure();
      if(Tmax < T_temp) Tmax = T_temp;
      if(Pmax < p_temp) Pmax = p_temp;
      if(Tmin > T_temp) Tmin = T_temp;
      if(Pmin > p_temp) Pmin = p_temp;
    }
  }

  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Reduce(&Tmax, &TmaxGlobal, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
  MPI_Reduce(&Pmax, &PmaxGlobal, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
  MPI_Reduce(&Tmin, &TminGlobal, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
  MPI_Reduce(&Pmin, &PminGlobal, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
}

//***********************************************************************
void Run::calcHRR()
{
  for(int lvl = m_lvlMax; lvl>=0; lvl--)
  {
    for (unsigned int i = 0; i < m_cellsLvl[lvl].size(); i++) 
    {
      if (!m_cellsLvl[lvl][i]->getSplit()) 
      {
        m_cellsLvl[lvl][i]->getMixture()->calcHRR();
      }
    }
  }
}