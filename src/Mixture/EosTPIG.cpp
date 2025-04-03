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

//! \file      EosTPIG.cpp
//! \author    E Fan, Lisong Shi, Jiaao Hao, Tianhan Zhang, Chih-yung Wen.
//! \version   1.0
//! \date      Jan 01 2025

#include <cmath>
#include <algorithm>
#include <cstring>
#include "EosTPIG.h"

using namespace std;

//***********************************************************************

EosTPIG::EosTPIG(){}

//***********************************************************************

EosTPIG::EosTPIG(vector<string> &nameParameterEos, int &number) :
    Eos(number)
{
  nameParameterEos.push_back("gamma");
  nameParameterEos.push_back("cv");
  nameParameterEos.push_back("energyRef");
  nameParameterEos.push_back("entropyRef");
}

EosTPIG::EosTPIG (int index, int &number):Eos(number)
{

}

//***********************************************************************

EosTPIG::~EosTPIG(){}

//***********************************************************************

//Constant methods
//****************
// p_i = rho_i * R_i * T_i
double EosTPIG::computeTemperature(const double &density, const double &pressure) const
{
  double T = pressure/max(density, epsilonAlphaNull)/m_Rs;
  if(std::isnan(T))
  {
    printf("T isnan in EosTPIG::computeTemperature().\n");
    printf("Index = %d, p = %lf, rho = %lf,epsilonAlphaNull = %le.\n",m_index,pressure,density,epsilonAlphaNull);
    exit(EXIT_FAILURE);
  }
  if(density<1e-30)
  {
    return 500;
  }
  else
  {
    return pressure/max(density, epsilonAlphaNull)/m_Rs;
  }
}

//***********************************************************************

double EosTPIG::computeEnergy(double temperature)const
{
  std::cout<<"Call EosTPIG::computeEnergy(double temperature)" << std::endl;
  exit(EXIT_FAILURE);
}

//***********************************************************************

double EosTPIG::computeEnergy(const double &density,const double &pressure)const
{
  std::cout<<"Call EosTPIG::computeEnergy(const double &density,const double &pressure)" << std::endl;
  exit(EXIT_FAILURE);
}

//***********************************************************************
double EosTPIG::computePressure(const double &density, const double &energy) const
{
  printf("Invalid call of EosTPIG::computePressure().\n");
  exit(EXIT_FAILURE);
}

//***********************************************************************

double EosTPIG::computeDensity(const double &pressure, const double &temperature) const
{
  return pressure/max((m_Rs*temperature), epsilonAlphaNull);
}

//***********************************************************************
double EosTPIG::computeSoundSpeed(const double &density, const double &pressure) const
{
  std::cout<<"Call EosTPIG::computeSoundSpeed(const double &density, const double &pressure)" << std::endl;
  exit(EXIT_FAILURE);
}

//***********************************************************************
double EosTPIG::computeEntropy(const double &temperature, const double &pressure) const
{
  printf("Invalid call of EosTPIG::computeEntropy().\n");
  exit(EXIT_FAILURE);
}

//***********************************************************************
double EosTPIG::computePressureIsentropic(const double &initialPressure, const double &initialDensity, const double &finalDensity) const
{
  printf("Invalid call of EosTPIG::computePressureIsentropic().\n");
  exit(EXIT_FAILURE);
}

//***********************************************************************
double EosTPIG::computePressureHugoniot(const double &initialPressure, const double &initialDensity, const double &finalDensity) const
{
  printf("Invalid call of EosTPIG::computePressureHugoniot().\n");
  exit(EXIT_FAILURE);
}

//***********************************************************************
double EosTPIG::computeDensityIsentropic(const double &initialPressure, const double &initialDensity, const double &finalPressure, double *drhodp) const
{
  printf("Invalid call of EosTPIG::computeDensityIsentropic().\n");
  exit(EXIT_FAILURE);
}

//***********************************************************************
double EosTPIG::computeDensityHugoniot(const double &initialPressure, const double &initialDensity, const double &finalPressure, double *drhodp) const
{
  printf("Invalid call of EosTPIG::computeDensityHugoniot().\n");
  exit(EXIT_FAILURE);
}

//***********************************************************************
double EosTPIG::computeDensityPfinal(const double &initialPressure, const double &initialDensity, const double &finalPressure, double *drhodp) const
{
  printf("Invalid call of EosTPIG::computeDensityPfinal().\n");
  exit(EXIT_FAILURE);  
}

//***********************************************************************
double EosTPIG::computeEnthalpyIsentropic(const double &initialPressure, const double &initialDensity, const double &finalPressure, double *dhdp) const
{
  printf("Invalid call of EosTPIG::computeEnthalpyIsentropic().\n");
  exit(EXIT_FAILURE); 
}

//***********************************************************************
double EosTPIG::computeDensitySaturation(const double &pressure, const double &Tsat, const double &dTsatdP, double *drhodp) const
{
  printf("Invalid call of EosTPIG::computeDensitySaturation().\n");
  exit(EXIT_FAILURE); 
}

//***********************************************************************
double EosTPIG::computeDensityEnergySaturation(const double &pressure, const double &rho, const double &drhodp, double *drhoedp) const
{
  printf("Invalid call of EosTPIG::computeDensityEnergySaturation().\n");
  exit(EXIT_FAILURE); 
}

//***********************************************************************
void EosTPIG::sendSpecialMixtureEos(double &gamPinfOverGamMinusOne, double &eRef, double &oneOverGamMinusOne, double &covolume) const
{
  printf("Invalid call of EosTPIG::sendSpecialMixtureEos().\n");
  exit(EXIT_FAILURE); 
}

//***********************************************************************
double EosTPIG::vfpfh(const double &pressure, const double &enthalpy) const
{
  printf("Invalid call of EosTPIG::vfpfh().\n");
  exit(EXIT_FAILURE); 
}

//***********************************************************************
double EosTPIG::dvdpch(const double &pressure, const double &enthalpy) const
{
  printf("Invalid call of EosTPIG::dvdpch().\n");
  exit(EXIT_FAILURE); 
}

//***********************************************************************
double EosTPIG::dvdhcp(const double &pressure, const double &enthalpy) const
{
  printf("Invalid call of EosTPIG::dvdhcp().\n");
  exit(EXIT_FAILURE); 
}

//***********************************************************************
void EosTPIG::verifyPressure(const double &pressure, const std::string &message) const
{
  if (pressure < 0.) errors.push_back(Errors(message + " : negative pressure in EosTPIG"));
}

void EosTPIG::verifyTemperature(const double &temperature, const std::string &message) const
{
  if (temperature < 1.e-15) 
  {
    printf("EosTPIG::verifyTemperature, T = %lf.\n", temperature);
    errors.push_back(Errors(message + " : too low temperture in EosTPIG"));
  }

}


//***********************************************************************
void EosTPIG::verifyAndModifyPressure(double &pressure) const
{
  if (pressure < 1.e-15) pressure = 1.e-15;
}

//***********************************************************************
const double& EosTPIG::getGamma() const 
{ 

  printf("Invalid call of EosTPIG::getGamma().\n");
  exit(EXIT_FAILURE); 
}

//***********************************************************************
const double& EosTPIG::getCv() const 
{ 
  printf("Invalid call of EosTPIG::getCv().\n");
  exit(EXIT_FAILURE); 
}

//***********************************************************************
const double& EosTPIG::getERef() const 
{  
  printf("Invalid call of EosTPIG::getERef().\n");
  exit(EXIT_FAILURE); 
}

//***********************************************************************
const double& EosTPIG::getSRef() const 
{ 
  printf("Invalid call of EosTPIG::getSRef().\n");
  exit(EXIT_FAILURE); 
}

double EosTPIG::getMW() const
{
  return m_MW;
}

//***********************************************************************