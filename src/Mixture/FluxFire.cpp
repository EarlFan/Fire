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

//! \file      FluxFire.cpp
//! \author    E Fan, Lisong Shi, Jiaao Hao, Tianhan Zhang, Chih-yung Wen.
//! \version   1.0
//! \date      Jan 01 2025

#include <cmath>
#include "FluxFire.h"
#include <algorithm> 
#include <chrono> 

using namespace std;
using namespace std::chrono; 

//***********************************************************************

FluxFire::FluxFire() {}

//***********************************************************************

FluxFire::FluxFire(ModFire *model) : m_model(model)
{

}

//***********************************************************************

FluxFire::~FluxFire(){}

//***********************************************************************

void FluxFire::printFlux() const
{
  // for(int ll=0;ll<NS;ll++) cout << m_masse[ll]<<", "<<endl;
  // cout << m_qdm.getX() << " " << m_energ << endl;
  printf("m_cons, 0 = %le, 1 = %le, 2 = %le.\n",m_masse[0],m_masse[1],m_masse[2]);
  printf("m_cons, u = %30.20le, v = %30.20le, e = %30.20le.\n",m_qdm.getX(),
    m_qdm.getY(),m_energ);
}

//***********************************************************************

void FluxFire::addFlux(double coefA, const int &numberPhases)
{
  for(int ll=0;ll<NS;ll++) m_masse[ll] += coefA*static_cast<FluxFire*> (fluxBuff)->m_masse[ll];
  m_qdm   += coefA*static_cast<FluxFire*> (fluxBuff)->m_qdm;
  m_energ += coefA*static_cast<FluxFire*> (fluxBuff)->m_energ;
}

//***********************************************************************

void FluxFire::addFlux(Flux* flux, const int &numberPhases)
{
  for(int ll=0;ll<NS;ll++) m_masse[ll] += static_cast<FluxFire*> (flux)->m_masse[ll];
  m_qdm   += static_cast<FluxFire*> (flux)->m_qdm;
  m_energ += static_cast<FluxFire*> (flux)->m_energ;
}

//***********************************************************************

void FluxFire::subtractFlux(double coefA, const int &numberPhases)
{
    for(int ll=0;ll<NS;ll++) m_masse[ll] -= coefA*static_cast<FluxFire*> (fluxBuff)->m_masse[ll];
    m_qdm   -= coefA*static_cast<FluxFire*> (fluxBuff)->m_qdm;
    m_energ -= coefA*static_cast<FluxFire*> (fluxBuff)->m_energ;
}

//***********************************************************************

void FluxFire::multiply(double scalar, const int &numberPhases)
{
    for(int ll=0;ll<NS;ll++) m_masse[ll] *= scalar;
    m_qdm   *= scalar;
    m_energ *= scalar;
}

//***********************************************************************

void FluxFire::setBufferFlux(Cell &cell, const int &numberPhases)
{
  static_cast<FluxFire*> (fluxBuff)->buildCons(cell.getPhases(), numberPhases, cell.getMixture());
}

//***********************************************************************

void FluxFire::setBufferFlux2()
{
  static_cast<FluxFire*> (fluxBuff)->m_masse = this->getMassArray();
  static_cast<FluxFire*> (fluxBuff)->m_qdm = this->getQdm();
  static_cast<FluxFire*> (fluxBuff)->m_energ = this->getEnergyMix();
}


//***********************************************************************
// flux: get conservative variables by known mixture
void FluxFire::buildCons(Phase **phases, const int &numberPhases, Mixture *mixture)
{
  
  m_masse = mixture->getDensityArray();

  double totalRho;
  totalRho = sumNS(m_masse);
  if(totalRho<=0)
  {
    printf("In FluxFire::buildCons, negative rho = %le, rho0 = %lf, rho1 = %lf.\n",totalRho,m_masse[0],m_masse[1]);
    exit(EXIT_FAILURE);
  }
  m_qdm = totalRho*mixture->getVelocity();
  m_energ = totalRho*mixture->getTotalEnergy();
}

//***********************************************************************
// from cons(updated in flux) to pris(in phase and mixture)
// Update: rho, rhoArray, velocity(u,v,w)
// Calls con2pri and updates T, p, c

void FluxFire::buildPrim(Phase **phases, Mixture *mixture, const int &numberPhases)
{
  // printf("In FluxFire::buildPrim, rho_N2 = %le, rho_u = %le, m_energy = %le.\n",m_masse[1],m_qdm.getX(),m_energ);
  double rhoMix, T, rhoArray[NS], totalE_mass;
  
  for(int ll=0;ll<numberPhases;ll++)
  {
    if(m_masse[ll]<1e-30)
    {
      //printf("Negative species density in FluxFire::buildPrim, rho0 = %le, rho1 = %le.\n",m_masse[0],m_masse[1]);
      //exit(EXIT_FAILURE);
      m_masse[ll]=0.;
    }
  }

  rhoMix = 0;

  for(int ll=0;ll<numberPhases;ll++) rhoMix+=m_masse[ll];

  if(rhoMix<=0.)
  {
    printf("Invalid rhoMix in FluxFire::buildPrim.\n");
    printf("rho0 = %le, rho1 = %le.\n",m_masse[0],m_masse[1]);
    exit(EXIT_FAILURE);
  }

  for(int ll=0;ll<numberPhases;ll++) rhoArray[ll]=m_masse[ll];

  //printf("rho0 = %lf, rho1 = %lf, rhoE=%lf.\n",m_masse[0],m_masse[1],m_energ);
  //printf("rhoU = %lf, rhoV = %lf.\n",m_qdm.getX(),m_qdm.getY());

  mixture->setDensity(rhoMix);
  mixture->setDensityArray(rhoArray);
  mixture->setVelocity(m_qdm.getX() / rhoMix, m_qdm.getY() / rhoMix, m_qdm.getZ() / rhoMix);
  totalE_mass = m_energ/rhoMix;
  mixture->setTotalEnergy(totalE_mass);

  mixture->con2pri();
}

//***********************************************************************

void FluxFire::setToZero(const int &numberPhases)
{
  for(int ll=0;ll<NS;ll++)m_masse[ll]=0.;
  m_qdm   = 0.;
  m_energ = 0.;
}

//***********************************************************************

void FluxFire::setToZeroBufferFlux(const int &numberPhases)
{
  for (int k = 0; k<numberPhases; k++) {
    static_cast<FluxFire*> (fluxBuff)->m_masse[k] = 0.;
  }
  static_cast<FluxFire*> (fluxBuff)->m_qdm = 0.;
  static_cast<FluxFire*> (fluxBuff)->m_energ = 0.;
}

//***********************************************************************

void FluxFire::addSymmetricTerms(Phase **phases, Mixture *mixture, const int &numberPhases, const double &r, const double &v)
{
  // for axisymmetric assumption in conservative form
  for (int k = 0; k<numberPhases; k++)
  {
    m_masse[k] += 0;
  }
  double m_qdm_y = m_qdm.getY();
  double a = m_qdm_y + mixture->getPressure()/r;
  m_qdm.setY(a);
  // printf("before addSym, m_qdm.Y = %30.20lf,  after addSym, m_qdm.Y = %30.20le, \np/r = %30.20lf, p = %30.20le, r = %30.20le.\n",
  //  m_qdm_y,a,mixture->getPressure()/r,mixture->getPressure(),r);
  
  double totalEnergy = 0;
  m_energ += 0;
}

//***********************************************************************

void FluxFire::setCons(const Flux *cons, const int &numberPhases)
{
  m_masse = cons->getMassArray();
  m_qdm = cons->getQdm();
  m_energ = cons->getEnergyMix();
}

//***********************************************************************

// Old chemical source term
void FluxFire::prepSourceTermsReaction(Cell *cell, const double &dt, const int &numberPhases, std::shared_ptr<Cantera::Solution> sol_)
{
  double Ys[NS],rhoMix(0.),Unity(0),Ys_old[NS],delta_rho[NS];
  std::array<double, NS> rhoArray;
  double rho = cell->getMixture()->getDensity();
  rhoArray = cell->getMixture()->getDensityArray();
  double T = cell->getMixture()->getTemperature();
  double P = cell->getMixture()->getPressure();

  for(int k=0;k<NS;k++) Ys[k] = rhoArray[k]/rho;
  
  for(int k=0;k<NS;k++) Ys_old[k] = Ys[k];

  std::shared_ptr<Cantera::ThermoPhase> gas = sol_->thermo();
	Cantera::IdealGasReactor reactor;
	Cantera::ReactorNet net;

  reactor.insert(sol_);
  net.addReactor(reactor);

  gas->setState_TPY(T,P,Ys);

  reactor.syncState();

	net.setInitialTime(0.0);
	net.advance(dt);// advance the reaction

  // the reaction only change the density array in cons variables
  for(int k=0;k<NS;k++) 
  {
    Ys[k] = gas->massFraction(k);
    if(Ys[k]<1e-30)Ys[k] = 0.;// set inproper mass fraction to 0.
    Unity += Ys[k];
  }

  for(int k=0;k<NS;k++) 
  {
    Ys[k] = Ys[k]/Unity;
    m_masse[k] = Ys[k]*rho;
  }

  T = gas->temperature();

  // set temperature to mixture for better convergence when calculate T from rho[] and e_int
  cell->getMixture()->setTemperature(T);
  cell->getMixture()->setPressure(gas->pressure());
  cell->getMixture()->updateMolarFractionArray();

  // add the chemical species change into mixture
  for(int k=0;k<NS;k++) 
  {
    delta_rho[k]= (Ys[k]-Ys_old[k])*rho;
    cell->getMixture()->addToDensityArrayChange(delta_rho[k],k);
  }

}
  