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

//! \file      ModFire.cpp
//! \author    E Fan, Lisong Shi, Jiaao Hao, Tianhan Zhang, Chih-yung Wen.
//! \version   1.0
//! \date      Jan 01 2025

#include <cmath>
#include <algorithm>
#include "ModFire.h"
#include "../Mixture/PhaseFire.h"
// #include "../../libs/ECOGEN/Relaxations/Relaxation.h"

const std::string ModFire::NAME = "FIRE";

//***********************************************************************

ModFire::ModFire(int &numberTransports, const int &numberPhases) :
  Model(NAME,numberTransports)
{
  fluxBuff = new FluxFire(this);
  // m_relaxations.push_back(new RelaxationP); //Pressure relaxation imposed in this model
  for (int i = 0; i < 4; i++) {
    sourceCons.push_back(new FluxFire(this));
  }
}

//***********************************************************************

ModFire::~ModFire()
{
  delete fluxBuff;
  for (int i = 0; i < 4; i++) {
    delete sourceCons[i];
  }
  sourceCons.clear();
}

//***********************************************************************

void ModFire::allocateCons(Flux **cons, const int &numberPhases)
{
  *cons = new FluxFire(this);
}

//***********************************************************************

void ModFire::allocatePhase(Phase **phase)
{
  *phase = new PhaseFire;
}

//***********************************************************************

void ModFire::allocateMixture(Mixture **mixture)
{
  *mixture = new MixFire;
}

//***********************************************************************
// the rho[], p, V is updated, try to update other variables except totalEnergy
void ModFire::fulfillState(Phase **phases, Mixture *mixture, const int &numberPhases, Prim type)
{
  double rhoMix = 0;
  for(int i=0;i<NS;i++)
  {
    rhoMix+=mixture->getDensityArray()[i];
  }
  if(rhoMix<=0)
  {
    printf("In ModFire::fulfillState, rhoMix = %lf, rankid = %d.\n",rhoMix,rankCpu);
    for(int i=0;i<NS;i++) printf("rho[%d] = %le.\n",i,mixture->getDensityArray()[0]);
    printf("p = %le, T = %le, u = %le.\n",mixture->getPressure(),mixture->getTemperature(),mixture->getU());
    exit(EXIT_FAILURE);
  }
  for(int i=0;i<NS;i++)
  {
    if(mixture->getDensityArray()[i]<1e-30)
    {
      mixture->setDensityI(0.,i);
    }
  }

  //mixture->setTemperature(mixture->computeTem());//update the T - fnae
  mixture->computeMixtureVariables(phases, numberPhases);
  //Specific to restart simulation
  if (type == restart) {
    for (int k = 0; k < numberPhases; k++) { phases[k]->setTemperature(mixture->getTemperature()); }
  }

  mixture->updateMolarFractionArray();
}

void ModFire::solveRiemannIntern(Cell &cellLeft, Cell &cellRight, const int &numberPhases, const double &dxLeft, const double &dxRight, double &dtMax) const
{
  if(AMRPara::riemann_type == 1)
  {
    solveRiemannInternHLL(cellLeft, cellRight, numberPhases, dxLeft, dxRight, dtMax);
  }
  else if(AMRPara::riemann_type == 2)
  {
    solveRiemannInternHLLC(cellLeft, cellRight, numberPhases, dxLeft, dxRight, dtMax);
  }
  else if(AMRPara::riemann_type == 3)
  {
    solveRiemannInternHLLC_LM(cellLeft, cellRight, numberPhases, dxLeft, dxRight, dtMax);
  }
  else
  {
    printf("Invalid riemann_type is given.");
    exit(EXIT_FAILURE);
  }

}

//****************************************************************************
void ModFire::solveRiemannInternHLL(Cell &cellLeft, Cell &cellRight, const int &numberPhases, const double &dxLeft, const double &dxRight, double &dtMax) const
{
  double sL, sR;
  
  double uL = cellLeft.getMixture()->getVelocity().getX(), vL = cellLeft.getMixture()->getVelocity().getY(),
         wL = cellLeft.getMixture()->getVelocity().getZ(), cL = cellLeft.getMixture()->getMixSoundSpeed(), 
         pL = cellLeft.getMixture()->getPressure(), rhoTotalL = cellLeft.getMixture()->getDensity(),
         TL = cellLeft.getMixture()->getTemperature(), intEL=cellLeft.getMixture()->getEnergy();
  std::array<double,NS> rhoArrayL = cellLeft.getMixture()->getDensityArray();
  double uR = cellRight.getMixture()->getVelocity().getX(), vR = cellRight.getMixture()->getVelocity().getY(), 
         wR = cellRight.getMixture()->getVelocity().getZ(), cR = cellRight.getMixture()->getMixSoundSpeed(), 
         pR = cellRight.getMixture()->getPressure(), rhoTotalR = cellRight.getMixture()->getDensity(),
         TR = cellRight.getMixture()->getTemperature(), intER=cellRight.getMixture()->getEnergy();
  std::array<double,NS> rhoArrayR = cellRight.getMixture()->getDensityArray();
  double EL = cellLeft.getMixture()->getEnergy() + 0.5*cellLeft.getMixture()->getVelocity().squaredNorm(),
         ER = cellRight.getMixture()->getEnergy() + 0.5*cellRight.getMixture()->getVelocity().squaredNorm();

  for(int ll=0;ll<NS;ll++)
  {
    if(rhoArrayL[ll]<0. || rhoArrayR[ll]<0.||rhoTotalL<=0.||rhoTotalR<=0.)
    {
      printf("Negative species density in ModEulerTPIG::solveRiemannIntern.\n");
      printf("rhoTotalL = %lf, rhoTotalR = %lf.\n",rhoTotalL,rhoTotalR);
      printf("rhoL0 = %lf, rhoL1 = %lf, rhoR0 = %lf, rhoR1 = %lf.\n",rhoArrayL[0],rhoArrayL[1],rhoArrayR[0],rhoArrayR[1]);
      exit(EXIT_FAILURE);
    }
  }

  //Davies
  sL = min(uL - cL, uR - cR);
  sR = max(uR + cR, uL + cL);

  if (abs(sL)>1.e-3) dtMax = min(dtMax, dxLeft / abs(sL));
  if (abs(sR)>1.e-3) dtMax = min(dtMax, dxRight / abs(sR));

  //compute left and right mass flow rates and sM
  double mL(rhoTotalL*(sL - uL)), mR(rhoTotalR*(sR - uR));
  double sM((pR - pL + mL*uL - mR*uR) / (mL - mR));// = S_STAR

  if (abs(sM)<1.e-8) sM = 0.;

  if (sL >= 0.){
    for(int ll=0;ll<NS;ll++)static_cast<FluxFire*> (fluxBuff)->m_masse[ll] = rhoArrayL[ll]*uL;
    static_cast<FluxFire*> (fluxBuff)->m_qdm.setX(rhoTotalL*uL*uL + pL);
    static_cast<FluxFire*> (fluxBuff)->m_qdm.setY(rhoTotalL*vL*uL);
    static_cast<FluxFire*> (fluxBuff)->m_qdm.setZ(rhoTotalL*wL*uL);
    static_cast<FluxFire*> (fluxBuff)->m_energ = (rhoTotalL*EL + pL)*uL;
  }
  else if (sR < 0.){
    for(int ll=0;ll<NS;ll++)static_cast<FluxFire*> (fluxBuff)->m_masse[ll] = rhoArrayR[ll]*uR;
    static_cast<FluxFire*> (fluxBuff)->m_qdm.setX(rhoTotalR*uR*uR + pR);
    static_cast<FluxFire*> (fluxBuff)->m_qdm.setY(rhoTotalR*vR*uR);
    static_cast<FluxFire*> (fluxBuff)->m_qdm.setZ(rhoTotalR*wR*uR);
    static_cast<FluxFire*> (fluxBuff)->m_energ = (rhoTotalR*ER + pR)*uR;
  }

  //1) Option HLL
  //(S_R*F_L-S_L*F_R+S_L*S_R*(U_R-U_L))/(S_R-S_L) = (S_L*F_R-S_R*F_L+S_R*S_L*(U_L-U_R))/(S_L-S_R)
  else if (abs(sR - sL)>1.e-3)
  {
   for(int ll=0;ll<NS;ll++)static_cast<FluxFire*> (fluxBuff)->m_masse[ll] = 
     (rhoArrayR[ll]*uR*sL - rhoArrayL[ll]*uL*sR + sL*sR*(rhoArrayL[ll] - rhoArrayR[ll])) / (sL - sR);
  //  static_cast<FluxFire*> (fluxBuff)->m_masse = (rhoR*uR*sL - rhoL*uL*sR + sL*sR*(rhoL - rhoR)) / (sL - sR);

   static_cast<FluxFire*> (fluxBuff)->m_qdm.setX(((rhoTotalR*uR*uR + pR)*sL - (rhoTotalL*uL*uL + pL)*sR +
     sL*sR*(rhoTotalL*uL - rhoTotalR*uR)) / (sL - sR));
   
   static_cast<FluxFire*> (fluxBuff)->m_qdm.setY((rhoTotalR*vR*uR*sL - rhoTotalL*vL*uL*sR + sL*sR*(rhoTotalL*vL - rhoTotalR*vR)) / (sL - sR));

   static_cast<FluxFire*> (fluxBuff)->m_qdm.setZ((rhoTotalR*wR*uR*sL - rhoTotalL*wL*uL*sR + sL*sR*(rhoTotalL*wL - rhoTotalR*wR)) / (sL - sR));

   static_cast<FluxFire*> (fluxBuff)->m_energ = ((rhoTotalR*ER + pR)*uR*sL - (rhoTotalL*EL + pL)*uL*sR + sL*sR*(rhoTotalL*EL - rhoTotalR*ER)) / (sL - sR);
  }

  //Contact discontinuity velocity
  static_cast<FluxFire*> (fluxBuff)->m_sM = sM;

}

//****************************************************************************
void ModFire::solveRiemannInternHLLC(Cell &cellLeft, Cell &cellRight, const int &numberPhases, const double &dxLeft, const double &dxRight, double &dtMax) const
{
  double sL, sR;
  
  double uL = cellLeft.getMixture()->getVelocity().getX(), vL = cellLeft.getMixture()->getVelocity().getY(),
         wL = cellLeft.getMixture()->getVelocity().getZ(), cL = cellLeft.getMixture()->getMixSoundSpeed(), 
         pL = cellLeft.getMixture()->getPressure(), rhoTotalL = cellLeft.getMixture()->getDensity(),
         TL = cellLeft.getMixture()->getTemperature(), intEL=cellLeft.getMixture()->getEnergy();
  std::array<double,NS> rhoArrayL = cellLeft.getMixture()->getDensityArray();
  double uR = cellRight.getMixture()->getVelocity().getX(), vR = cellRight.getMixture()->getVelocity().getY(), 
         wR = cellRight.getMixture()->getVelocity().getZ(), cR = cellRight.getMixture()->getMixSoundSpeed(), 
         pR = cellRight.getMixture()->getPressure(), rhoTotalR = cellRight.getMixture()->getDensity(),
         TR = cellRight.getMixture()->getTemperature(), intER=cellRight.getMixture()->getEnergy();
  std::array<double,NS> rhoArrayR = cellRight.getMixture()->getDensityArray();
  double EL = cellLeft.getMixture()->getEnergy() + 0.5*cellLeft.getMixture()->getVelocity().squaredNorm(),
         ER = cellRight.getMixture()->getEnergy() + 0.5*cellRight.getMixture()->getVelocity().squaredNorm();

  for(int ll=0;ll<NS;ll++)
  {
    if(rhoArrayL[ll]<0. || rhoArrayR[ll]<0.||rhoTotalL<=0.||rhoTotalR<=0.)
    {
      printf("Negative species density in ModEulerTPIG::solveRiemannIntern.\n");
      printf("rhoTotalL = %lf, rhoTotalR = %lf.\n",rhoTotalL,rhoTotalR);
      printf("rhoL0 = %lf, rhoL1 = %lf, rhoR0 = %lf, rhoR1 = %lf.\n",rhoArrayL[0],rhoArrayL[1],rhoArrayR[0],rhoArrayR[1]);
      exit(EXIT_FAILURE);
    }
  }

  //Davies
  sL = min(uL - cL, uR - cR);
  sR = max(uR + cR, uL + cL);

  if (abs(sL)>1.e-3) dtMax = min(dtMax, dxLeft / abs(sL));
  if (abs(sR)>1.e-3) dtMax = min(dtMax, dxRight / abs(sR));

  //compute left and right mass flow rates and sM
  double mL(rhoTotalL*(sL - uL)), mR(rhoTotalR*(sR - uR));
  double sM((pR - pL + mL*uL - mR*uR) / (mL - mR));// = S_STAR

  if (abs(sM)<1.e-8) sM = 0.;

  if (sL >= 0.){
    for(int ll=0;ll<NS;ll++)static_cast<FluxFire*> (fluxBuff)->m_masse[ll] = rhoArrayL[ll]*uL;
    static_cast<FluxFire*> (fluxBuff)->m_qdm.setX(rhoTotalL*uL*uL + pL);
    static_cast<FluxFire*> (fluxBuff)->m_qdm.setY(rhoTotalL*vL*uL);
    static_cast<FluxFire*> (fluxBuff)->m_qdm.setZ(rhoTotalL*wL*uL);
    static_cast<FluxFire*> (fluxBuff)->m_energ = (rhoTotalL*EL + pL)*uL;
  }
  else if (sR < 0.){
    for(int ll=0;ll<NS;ll++)static_cast<FluxFire*> (fluxBuff)->m_masse[ll] = rhoArrayR[ll]*uR;
    static_cast<FluxFire*> (fluxBuff)->m_qdm.setX(rhoTotalR*uR*uR + pR);
    static_cast<FluxFire*> (fluxBuff)->m_qdm.setY(rhoTotalR*vR*uR);
    static_cast<FluxFire*> (fluxBuff)->m_qdm.setZ(rhoTotalR*wR*uR);
    static_cast<FluxFire*> (fluxBuff)->m_energ = (rhoTotalR*ER + pR)*uR;
  }
  // 2) Option HLLC
  else if (sM >= 0.) {
    double pStar = mL*(sM - uL) + pL;
    double rhoStar = mL / (sL - sM);
    double Estar = EL + (sM - uL)*(sM + pL / mL);
    for(int ll=0; ll<NS; ll++) static_cast<FluxFire*> (fluxBuff)->m_masse[ll] = rhoStar*rhoArrayL[ll]*sM/rhoTotalL;
    static_cast<FluxFire*> (fluxBuff)->m_qdm.setX(rhoStar*sM*sM+pStar);
    static_cast<FluxFire*> (fluxBuff)->m_qdm.setY(rhoStar*sM*vL);
    static_cast<FluxFire*> (fluxBuff)->m_qdm.setZ(rhoStar*sM*wL);
    static_cast<FluxFire*> (fluxBuff)->m_energ = (rhoStar*Estar + pStar)*sM;
  }
  else {
    double pStar = mR*(sM - uR) + pR;
    double rhoStar = mR / (sR - sM);
    double Estar = ER + (sM - uR)*(sM + pR / mR);
    for(int ll=0; ll<NS; ll++) static_cast<FluxFire*> (fluxBuff)->m_masse[ll] = rhoStar*rhoArrayR[ll]*sM/rhoTotalR;
    static_cast<FluxFire*> (fluxBuff)->m_qdm.setX(rhoStar*sM*sM + pStar);
    static_cast<FluxFire*> (fluxBuff)->m_qdm.setY(rhoStar*sM*vR);
    static_cast<FluxFire*> (fluxBuff)->m_qdm.setZ(rhoStar*sM*wR);
    static_cast<FluxFire*> (fluxBuff)->m_energ = (rhoStar*Estar + pStar)*sM;
  }

  //Contact discontinuity velocity
  static_cast<FluxFire*> (fluxBuff)->m_sM = sM;
}

//****************************************************************************
void ModFire::solveRiemannInternHLLC_LM(Cell &cellLeft, Cell &cellRight, const int &numberPhases, const double &dxLeft, const double &dxRight, double &dtMax) const
{
  double sL, sR;

  double Ma_limit = 0.1;
  
  double uL = cellLeft.getMixture()->getVelocity().getX(), vL = cellLeft.getMixture()->getVelocity().getY(),
         wL = cellLeft.getMixture()->getVelocity().getZ(), cL = cellLeft.getMixture()->getMixSoundSpeed(), 
         pL = cellLeft.getMixture()->getPressure(), rhoTotalL = cellLeft.getMixture()->getDensity(),
         TL = cellLeft.getMixture()->getTemperature(), intEL=cellLeft.getMixture()->getEnergy();
  std::array<double,NS> rhoArrayL = cellLeft.getMixture()->getDensityArray();
  double uR = cellRight.getMixture()->getVelocity().getX(), vR = cellRight.getMixture()->getVelocity().getY(), 
         wR = cellRight.getMixture()->getVelocity().getZ(), cR = cellRight.getMixture()->getMixSoundSpeed(), 
         pR = cellRight.getMixture()->getPressure(), rhoTotalR = cellRight.getMixture()->getDensity(),
         TR = cellRight.getMixture()->getTemperature(), intER=cellRight.getMixture()->getEnergy();
  std::array<double,NS> rhoArrayR = cellRight.getMixture()->getDensityArray();
  double EL = cellLeft.getMixture()->getEnergy() + 0.5*cellLeft.getMixture()->getVelocity().squaredNorm(),
         ER = cellRight.getMixture()->getEnergy() + 0.5*cellRight.getMixture()->getVelocity().squaredNorm();

  for(int ll=0;ll<NS;ll++)
  {
    if(rhoArrayL[ll]<0. || rhoArrayR[ll]<0.||rhoTotalL<=0.||rhoTotalR<=0.)
    {
      printf("Negative species density in ModEulerTPIG::solveRiemannIntern.\n");
      printf("rhoTotalL = %lf, rhoTotalR = %lf.\n",rhoTotalL,rhoTotalR);
      printf("rhoL0 = %lf, rhoL1 = %lf, rhoR0 = %lf, rhoR1 = %lf.\n",rhoArrayL[0],rhoArrayL[1],rhoArrayR[0],rhoArrayR[1]);
      exit(EXIT_FAILURE);
    }
  }

  //Davies
  sL = min(uL - cL, uR - cR);
  sR = max(uR + cR, uL + cL);

  if (abs(sL)>1.e-3) dtMax = min(dtMax, dxLeft / abs(sL));
  if (abs(sR)>1.e-3) dtMax = min(dtMax, dxRight / abs(sR));

  //compute left and right mass flow rates and sM
  double mL(rhoTotalL*(sL - uL)), mR(rhoTotalR*(sR - uR));
  double sM((pR - pL + mL*uL - mR*uR) / (mL - mR));// = S_STAR

  if (abs(sM)<1.e-8) sM = 0.;
  double res_flux[NS+4];
  double flux_central[NS+4];
  double UL[NS+4], UR[NS+4];
  double UL_star[NS+4], UR_star[NS+4];

  if (sL >= 0.){
    for(int ll=0;ll<NS;ll++)static_cast<FluxFire*> (fluxBuff)->m_masse[ll] = rhoArrayL[ll]*uL;
    static_cast<FluxFire*> (fluxBuff)->m_qdm.setX(rhoTotalL*uL*uL + pL);
    static_cast<FluxFire*> (fluxBuff)->m_qdm.setY(rhoTotalL*vL*uL);
    static_cast<FluxFire*> (fluxBuff)->m_qdm.setZ(rhoTotalL*wL*uL);
    static_cast<FluxFire*> (fluxBuff)->m_energ = (rhoTotalL*EL + pL)*uL;
  }
  else if (sR < 0.){
    for(int ll=0;ll<NS;ll++)static_cast<FluxFire*> (fluxBuff)->m_masse[ll] = rhoArrayR[ll]*uR;
    static_cast<FluxFire*> (fluxBuff)->m_qdm.setX(rhoTotalR*uR*uR + pR);
    static_cast<FluxFire*> (fluxBuff)->m_qdm.setY(rhoTotalR*vR*uR);
    static_cast<FluxFire*> (fluxBuff)->m_qdm.setZ(rhoTotalR*wR*uR);
    static_cast<FluxFire*> (fluxBuff)->m_energ = (rhoTotalR*ER + pR)*uR;
  }
  else{
    double Ma_local, sL_LM, sR_LM, phi;
    Ma_local = max(abs(uL/cL), abs(uR/cR));
    phi = sin(min(1., Ma_local/Ma_limit)*M_PI/2);
    sL_LM = phi*sL;
    sR_LM = phi*sR;

    
    for(int ll=0;ll<NS;ll++) flux_central[ll] = 0.5*(rhoArrayL[ll]*uL+rhoArrayR[ll]*uR);
    flux_central[NS] = 0.5*(rhoTotalL*uL*uL + pL + rhoTotalR*uR*uR + pR);
    flux_central[NS+1] = 0.5*(rhoTotalL*vL*uL+rhoTotalR*vR*uR);
    flux_central[NS+2] = 0.5*(rhoTotalL*wL*uL+rhoTotalR*wR*uR);
    flux_central[NS+3] = 0.5*((rhoTotalL*EL + pL)*uL + (rhoTotalR*ER + pR)*uR);

    
    for(int ll=0;ll<NS;ll++)
    {
      UL[ll] = rhoArrayL[ll];
    }
    UL[NS] = rhoTotalL*uL;
    UL[NS+1] = rhoTotalL*vL;
    UL[NS+2] = rhoTotalL*wL;
    UL[NS+3] = rhoTotalL*EL;

    for(int ll=0;ll<NS;ll++)
    {
      UR[ll] = rhoArrayR[ll];
    }
    UR[NS] = rhoTotalR*uR;
    UR[NS+1] = rhoTotalR*vR;
    UR[NS+2] = rhoTotalR*wR;
    UR[NS+3] = rhoTotalR*ER;

    
    
    double pStarL = mL*(sM - uL) + pL;
    double rhoStarL = mL / (sL - sM);
    double EstarL = EL + (sM - uL)*(sM + pL / mL);
    for(int ll=0; ll<NS; ll++) UL_star[ll] = rhoStarL*rhoArrayL[ll]/rhoTotalL;
    UL_star[NS] = rhoStarL*sM;
    UL_star[NS+1] = rhoStarL*vL;
    UL_star[NS+2] = rhoStarL*wL;
    UL_star[NS+3] = rhoStarL*(EL + (sM-uL)*(sM+pL/(rhoTotalL*(sL - uL))));

    double pStarR = mR*(sM - uR) + pR;
    double rhoStarR = mR / (sR - sM);
    double EstarR = ER + (sM - uR)*(sM + pR / mR);
    for(int ll=0; ll<NS; ll++) UR_star[ll] = rhoStarR*rhoArrayR[ll]/rhoTotalR;
    UR_star[NS] = rhoStarR*sM;
    UR_star[NS+1] = rhoStarR*vR;
    UR_star[NS+2] = rhoStarR*wR;
    UR_star[NS+3] = rhoStarR*(ER + (sM-uR)*(sM+pR/(rhoTotalR*(sR - uR))));

    
    for(int ll=0; ll<NS+4; ll++) 
    {
       res_flux[ll] = flux_central[ll] 
        + 0.5*(sL_LM*(UL_star[ll] - UL[ll]) + abs(sM)*(UL_star[ll] - UR_star[ll]) + sR_LM*(UR_star[ll] - UR[ll]) );
    }

    for(int ll=0; ll<NS; ll++) 
    {
      static_cast<FluxFire*> (fluxBuff)->m_masse[ll] = res_flux[ll] ;
    }
    static_cast<FluxFire*> (fluxBuff)->m_qdm.setX(res_flux[NS]);
    static_cast<FluxFire*> (fluxBuff)->m_qdm.setY(res_flux[NS+1]);
    static_cast<FluxFire*> (fluxBuff)->m_qdm.setZ(res_flux[NS+2]);
    static_cast<FluxFire*> (fluxBuff)->m_energ = res_flux[NS+3];

  }

  //Contact discontinuity velocity
  static_cast<FluxFire*> (fluxBuff)->m_sM = sM;
}

//****************************************************************************
void ModFire::solveRiemannMixture(Mixture &mixtureLeft, Mixture &mixtureRight, const int &numberPhases, const double &dxLeft, const double &dxRight, double &dtMax) const
{
  printf("Unwanted use of solveRiemannMixture().\n");
  exit(EXIT_FAILURE);  

}

void ModFire::solveRiemannWall(Cell &cellLeft, const int &numberPhases, const double &dxLeft, double &dtMax) const
{

  // analytical solution

  if(0)
  {
    double sL;
    
    double uL = cellLeft.getMixture()->getVelocity().getX(), vL = cellLeft.getMixture()->getVelocity().getY(),
          wL = cellLeft.getMixture()->getVelocity().getZ(), cL = cellLeft.getMixture()->getMixSoundSpeed(), 
          pL = cellLeft.getMixture()->getPressure(), rhoTotalL = cellLeft.getMixture()->getDensity(),
          TL = cellLeft.getMixture()->getTemperature(), intEL=cellLeft.getMixture()->getEnergy();
    std::array<double,NS> rhoArrayL = cellLeft.getMixture()->getDensityArray();
    double EL = cellLeft.getMixture()->getEnergy() + 0.5*cellLeft.getMixture()->getVelocity().squaredNorm();

    for(int ll=0;ll<NS;ll++)
    {
      if(rhoArrayL[ll]<0. ||rhoTotalL<=0.)
      {
        printf("Negative species density in ModEulerTPIG::solveRiemannIntern.\n");
        printf("rhoTotalL = %lf.\n",rhoTotalL);
        printf("rhoL0 = %lf, rhoL1 = %lf.\n",rhoArrayL[0],rhoArrayL[1]);
        exit(EXIT_FAILURE);
      }
    }

    //Davies
    sL = min(uL - cL, -uL - cL);

    if (abs(sL)>1.e-3) dtMax = min(dtMax, dxLeft / abs(sL));

    double gammaL = rhoTotalL*cL*cL/pL;
    double Cm = uL;
    double pStar;

    if(uL <= 0) // from Toro's book, section 6.3
    {
      pStar = pL*pow((1+0.5*(gammaL-1)*Cm/cL), 2*gammaL/(gammaL-1));
    }
    else
    {
      double Am = 2/(gammaL+1)/rhoTotalL;
      double Bm = (gammaL-1)/(gammaL+1)*pL;
      pStar = pL+Cm/2/Am*(Cm+pow(Cm*Cm+4*Am*(Bm+pL),0.5));
    }

    // double pStar = rhoTotalL*(uL - sL)*uL + pL;

    for(int ll=0;ll<NS;ll++)static_cast<FluxFire*> (fluxBuff)->m_masse[ll] = 0;
    static_cast<FluxFire*> (fluxBuff)->m_qdm.setX(pStar);
    static_cast<FluxFire*> (fluxBuff)->m_qdm.setY(0.);
    static_cast<FluxFire*> (fluxBuff)->m_qdm.setZ(0.);
    static_cast<FluxFire*> (fluxBuff)->m_energ = 0.;

    //Contact discontinuity velocity
    static_cast<FluxFire*> (fluxBuff)->m_sM = 0;
  }


  // same Riemann solver with inner cell
  if(1) 
  {
    solveRiemannSymmetryInner(cellLeft, numberPhases, dxLeft, dtMax);
  }

}


//****************************************************************************

void ModFire::solveRiemannSymmetryInner(Cell &cellLeft, const int &numberPhases, const double &dxLeft, double &dtMax) const
{
   double sL, sR;
  
  double uL = cellLeft.getMixture()->getVelocity().getX(), vL = cellLeft.getMixture()->getVelocity().getY(),
         wL = cellLeft.getMixture()->getVelocity().getZ(), cL = cellLeft.getMixture()->getMixSoundSpeed(), 
         pL = cellLeft.getMixture()->getPressure(), rhoTotalL = cellLeft.getMixture()->getDensity(),
         TL = cellLeft.getMixture()->getTemperature(), intEL=cellLeft.getMixture()->getEnergy();
  std::array<double,NS> rhoArrayL = cellLeft.getMixture()->getDensityArray();

  double uR = -uL, vR = vL, wR = wL, cR = cL, pR = pL, rhoTotalR = rhoTotalL, TR = TL, intER = intEL;

  std::array<double,NS> rhoArrayR = rhoArrayL;

  double EL = cellLeft.getMixture()->getEnergy() + 0.5*cellLeft.getMixture()->getVelocity().squaredNorm(),
         ER = EL;

  for(int ll=0;ll<NS;ll++)
  {
    if(rhoArrayL[ll]<0. || rhoArrayR[ll]<0.||rhoTotalL<=0.||rhoTotalR<=0.)
    {
      printf("Negative species density in ModEulerTPIG::solveRiemannIntern.\n");
      printf("rhoTotalL = %lf, rhoTotalR = %lf.\n",rhoTotalL,rhoTotalR);
      printf("rhoL0 = %lf, rhoL1 = %lf, rhoR0 = %lf, rhoR1 = %lf.\n",rhoArrayL[0],rhoArrayL[1],rhoArrayR[0],rhoArrayR[1]);
      exit(EXIT_FAILURE);
    }
  }

  //Davies
  sL = min(uL - cL, uR - cR);
  sR = max(uR + cR, uL + cL);

  if (abs(sL)>1.e-3) dtMax = min(dtMax, dxLeft / abs(sL));

  //compute left and right mass flow rates and sM
  double mL(rhoTotalL*(sL - uL)), mR(rhoTotalR*(sR - uR));
  double sM((pR - pL + mL*uL - mR*uR) / (mL - mR));// = S_STAR

  if (abs(sM)<1.e-8) sM = 0.;

  if (sL >= 0.){
    for(int ll=0;ll<NS;ll++)static_cast<FluxFire*> (fluxBuff)->m_masse[ll] = rhoArrayL[ll]*uL;
    static_cast<FluxFire*> (fluxBuff)->m_qdm.setX(rhoTotalL*uL*uL + pL);
    static_cast<FluxFire*> (fluxBuff)->m_qdm.setY(rhoTotalL*vL*uL);
    static_cast<FluxFire*> (fluxBuff)->m_qdm.setZ(rhoTotalL*wL*uL);
    static_cast<FluxFire*> (fluxBuff)->m_energ = (rhoTotalL*EL + pL)*uL;
  }
  else if (sR < 0.){
    for(int ll=0;ll<NS;ll++)static_cast<FluxFire*> (fluxBuff)->m_masse[ll] = rhoArrayR[ll]*uR;
    static_cast<FluxFire*> (fluxBuff)->m_qdm.setX(rhoTotalR*uR*uR + pR);
    static_cast<FluxFire*> (fluxBuff)->m_qdm.setY(rhoTotalR*vR*uR);
    static_cast<FluxFire*> (fluxBuff)->m_qdm.setZ(rhoTotalR*wR*uR);
    static_cast<FluxFire*> (fluxBuff)->m_energ = (rhoTotalR*ER + pR)*uR;
  }

  //2) Option HLLC
  else if (sM >= 0.) {
    double pStar = mL*(sM - uL) + pL;
    double rhoStar = mL / (sL - sM);
    double Estar = EL + (sM - uL)*(sM + pL / mL);
    for(int ll=0; ll<NS; ll++) static_cast<FluxFire*> (fluxBuff)->m_masse[ll] = rhoStar*rhoArrayL[ll]*sM/rhoTotalL;
    static_cast<FluxFire*> (fluxBuff)->m_qdm.setX(rhoStar*sM*sM+pStar);
    static_cast<FluxFire*> (fluxBuff)->m_qdm.setY(rhoStar*sM*vL);
    static_cast<FluxFire*> (fluxBuff)->m_qdm.setZ(rhoStar*sM*wL);
    static_cast<FluxFire*> (fluxBuff)->m_energ = (rhoStar*Estar + pStar)*sM;
  }
  else {
    double pStar = mR*(sM - uR) + pR;
    double rhoStar = mR / (sR - sM);
    double Estar = ER + (sM - uR)*(sM + pR / mR);
    for(int ll=0; ll<NS; ll++) static_cast<FluxFire*> (fluxBuff)->m_masse[ll] = rhoStar*rhoArrayR[ll]*sM/rhoTotalR;
    static_cast<FluxFire*> (fluxBuff)->m_qdm.setX(rhoStar*sM*sM + pStar);
    static_cast<FluxFire*> (fluxBuff)->m_qdm.setY(rhoStar*sM*vR);
    static_cast<FluxFire*> (fluxBuff)->m_qdm.setZ(rhoStar*sM*wR);
    static_cast<FluxFire*> (fluxBuff)->m_energ = (rhoStar*Estar + pStar)*sM;
  }

  //Contact discontinuity velocity
  static_cast<FluxFire*> (fluxBuff)->m_sM = sM;
  
}



//****************************************************************************

const double& ModFire::getSM()
{
  return static_cast<FluxFire*> (fluxBuff)->m_sM;
}

//****************************************************************************
//***************************** others methods *******************************
//****************************************************************************

void ModFire::reverseProjection(const Coord normal, const Coord tangent, const Coord binormal) const
{
  Coord fluxProjete;
  fluxProjete.setX(normal.getX()*static_cast<FluxFire*> (fluxBuff)->m_qdm.getX() + tangent.getX()*static_cast<FluxFire*> (fluxBuff)->m_qdm.getY() + binormal.getX()*static_cast<FluxFire*> (fluxBuff)->m_qdm.getZ());
  fluxProjete.setY(normal.getY()*static_cast<FluxFire*> (fluxBuff)->m_qdm.getX() + tangent.getY()*static_cast<FluxFire*> (fluxBuff)->m_qdm.getY() + binormal.getY()*static_cast<FluxFire*> (fluxBuff)->m_qdm.getZ());
  fluxProjete.setZ(normal.getZ()*static_cast<FluxFire*> (fluxBuff)->m_qdm.getX() + tangent.getZ()*static_cast<FluxFire*> (fluxBuff)->m_qdm.getY() + binormal.getZ()*static_cast<FluxFire*> (fluxBuff)->m_qdm.getZ());
  static_cast<FluxFire*> (fluxBuff)->m_qdm.setXYZ(fluxProjete.getX(), fluxProjete.getY(), fluxProjete.getZ());
}

//****************************************************************************