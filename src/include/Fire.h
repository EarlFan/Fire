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

#ifndef FIRE_H
#define FIRE_H

//! \file      Fire.h
//! \author    E Fan, Lisong Shi, Jiaao Hao, Tianhan Zhang, Chih-yung Wen.
//! \version   1.0
//! \date      Jan 01 2025

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>
#include <array>
#include "cantera/base/Solution.h"

#define FORMULAR_LENGTH 5
#define Ru 8314.4621 //J/(kmol*K)
#define TOLERANCE 1.e-9

class AMRPara
{
public:
	static double sigma_c;
	static double sigma_r;
	static double xi_diff_lim;
	static double f1;
	static double f2;
	static double f3;
	static double dx_expo;
	static int better_unrefine;
	static int div_vor_v_flag;
	static int dim;
	static int amr_iter;
	static double amr_flux;
	static int YH2O;
	static int YCO2;
	static int YOH;
	static double sigma_gradRho;
	static int gradRho_flag;
	static double physical_time;
	static int chem_step_gap;
	static int g_lvlmax;
	static int riemann_type;
	static int test_prob;
	static int species_diffusion_type;
	static int strang_splitting_flag;
};

using namespace std;
//
struct SpeciesPha {
	string Name;					// Species formula, length = 5
	double MW;						// Molecular Weight,
	double RsPha;						// Specific R, J/(kg*K)
};

extern char const * dataBase[NS];

extern struct SpeciesPha Sptr[NS];					//Array of Species Coefs
extern struct SpeciesPha excitedOH;                 //only for OH*
extern int iReaction;
extern char formular[NS][FORMULAR_LENGTH];
extern double YSINF[NS];

extern std::shared_ptr<Cantera::Solution> solution_mech;
extern std::shared_ptr<Cantera::ThermoPhase> gas_global;

// read molecular weight for each species
void readMW(void);

// return summation of NS elements
double sumNS(double *);

// return summation of NS elements
double sumNS(std::array<double, NS> arr);

// calculate transport properties
void transport_model(double Ys[NS],double T,double P,double rhoMix,double trans[NS+2]);

void displayHeaderFire();

#endif //FIRE_H
