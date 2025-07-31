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

//! \file      Fire.cpp
//! \author    E Fan, Lisong Shi, Jiaao Hao, Tianhan Zhang, Chih-yung Wen.
//! \version   1.0
//! \date      Jan 01 2025

#include "Fire.h"
#include "cantera/thermo.h"
#include "cantera/transport.h"

double AMRPara::sigma_c = 1e30;
double AMRPara::sigma_r = 1e30;
double AMRPara::xi_diff_lim = 1e30;
double AMRPara::f1 = 1e30;
double AMRPara::f2 = 1e30;
double AMRPara::f3 = 1e30;
double AMRPara::dx_expo = 1e10;
int AMRPara::better_unrefine = 100;
int AMRPara::div_vor_v_flag = 100;
int AMRPara::dim = 100;
int AMRPara::amr_iter = 100;
double AMRPara::amr_flux = 1e30;
int AMRPara::YH2O = 100;
int AMRPara::YCO2 = 100;
int AMRPara::YOH = 100;
double AMRPara::sigma_gradRho = 1e30;
int AMRPara::gradRho_flag = 100;
double AMRPara::physical_time = 0;
int AMRPara::chem_step_gap = 1;
int AMRPara::g_lvlmax = -1;
int AMRPara::riemann_type = -1;
int AMRPara::test_prob = -1;
int AMRPara::species_diffusion_type = 0; 
int AMRPara::strang_splitting_flag = 1; 

// return summation of NS elements
double sumNS(double *arr)
{
	double sum = 0;
	for (int i = 0; i < NS; i++)sum += arr[i];
	return sum;
}

double sumNS(std::array<double, NS> arr)
{
	double sum = 0;
	for (int i = 0; i < NS; i++)sum += arr[i];
	return sum;
}

void transport_model(double Ys[NS],double T,double P,double rhoMix,double trans[NS+2])
{
	std::shared_ptr<Cantera::Solution> sol = solution_mech;
	std::shared_ptr<Cantera::ThermoPhase> gas = sol->thermo();
	gas->setState_TPY(T,P,Ys);
	size_t nsp = gas->nSpecies();
	trans[0] = sol->transport()->viscosity();
	trans[1] = sol->transport()->thermalConductivity();
	double diff_coef_temp[NS];
	sol->transport()->getMixDiffCoeffs(diff_coef_temp);
	for(int i=0;i<NS;i++)
	{
		if(AMRPara::species_diffusion_type == 1) // consider the pressure effect
		{
			diff_coef_temp[i] *= 101325/P;
		}
		trans[i+2] = diff_coef_temp[i];
	}     
}

void displayHeaderFire()
{
  std::cout << "************************************************************" << std::endl;
  std::cout << "             _____   _  " << std::endl;
  std::cout << "            |  ___| (_)  _ __    ___  " << std::endl;
  std::cout << "            | |_    | | | '__|  / _ \\" << std::endl;
  std::cout << "            |  _|   | | | |    |  __/" << std::endl;
  std::cout << "            |_|     |_| |_|     \\___| " << std::endl;
  std::cout << std::endl;
  std::cout << "The Fire solver for supersonic combustion and detonation simulations."  << std::endl;
  std::cout << "************************************************************" << std::endl;
}