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

#include "../libs/ECOGEN/Order1/Cell.h"

//***********************************************************************

void Cell::fill(std::vector<GeometricalDomain*> &domains, const int &lvlMax)
{
  // normal initial setups
  if(AMRPara::test_prob == -1 || AMRPara::test_prob == 4)
  {
    Coord coordinates;
    coordinates = m_element->getPosition();
    for (unsigned int geom = 0; geom < domains.size(); geom++) {
        domains[geom]->fillIn(this, m_numberPhases, m_numberTransports);
    }
  }
  else // special initial setups
  {
    if(m_model->whoAmI() != "FIRE")
    {
      Coord coordinates;
      coordinates = m_element->getPosition();
      for (unsigned int geom = 0; geom < domains.size(); geom++) {
          domains[geom]->fillIn(this, m_numberPhases, m_numberTransports);
      }
    }
    else // for the fire solver
    {
      double C(2.e4);
      double C1(8.e3);
      Coord coordinates, center;
      Mixture *mixLP, *mixShock, *mixBubble;
      double halfPosition, radius, interface_x;
      coordinates = m_element->getPosition();
      // get the bubble center
      for (unsigned int geom = 0; geom < domains.size(); geom++) 
      {
        // this fillin will copyPhase and setTransport
        domains[geom]->fillIn(this, m_numberPhases, m_numberTransports);
        if(domains[geom]->getName() == "BASE")
        {
          // printf("geom = %d, domain name = chamberLP.\n",geom);
          mixLP = domains[geom]->getMixture();
        }
        else if(domains[geom]->getName() == "HP")
        {
          // printf("geom = %d, domain name = shockedState.\n",geom);
          mixShock = domains[geom]->getMixture();
          halfPosition = domains[geom]->getHalfPosition();
        }
        else if(domains[geom]->getName() == "BUBBLE")
        {
          // printf("geom = %d, domain name = bubble.\n",geom);
          mixBubble = domains[geom]->getMixture();
          center = domains[geom]->centerDisc();
          radius = domains[geom]->getRadius();
        }
        else if(domains[geom]->getName() == "HOTSPOT")
        {
          // printf("geom = %d, domain name = bubble.\n",geom);
          mixBubble = domains[geom]->getMixture();
        }
      }
      
      double T0, density0[NS], alpha, alpha1;

      double celly = coordinates.getY();
      double lambda = 0.0004;
      double A = 0.0005;
      interface_x = A*sin(celly*2*M_PI/lambda)+halfPosition;

      if(coordinates.getX()<halfPosition)
      {
        m_mixture->copyMixture(*mixShock);
        
      }
      else
      {
        m_mixture->copyMixture(*mixLP);// copy from LP state

        // for 3d channel detonation, perturbation in pre-shock region
        if(AMRPara::test_prob == 5) // for example_cases/3D_channel_detonation
        {
          double cell_x = coordinates.getX();
          double cell_y = coordinates.getY();
          double cell_z = coordinates.getZ();

          double p_x(0.03);
          double p_y(0.003);
          double p_z(0.003);

          // pre-shock perturbation
          if(cell_x>0.006 and cell_x<0.006+p_x and cell_y>p_y and cell_y<0.015 - p_y and cell_z>p_z and cell_z<0.015 - p_z)
          {
            std::array<double,NS> rhoInB;
            double rhoMixInB(0.), rhoMixOutB(0.), YsInB[NS], YsOutB[NS], YsCell[NS]; 
            double pInB,TInB,p,T;

            TInB = mixBubble->getTemperature();
            pInB = mixBubble->getPressure();

            rhoInB = mixBubble->getDensityArray();

            for(int i=0;i<NS;i++) {rhoMixInB+=rhoInB[i];}
            for(int i=0;i<NS;i++){YsInB[i]=rhoInB[i]/rhoMixInB;}

            double RMix(0.);
            T = TInB;
            p = pInB;
            for(int i=0;i<NS;i++)RMix+=YsOutB[i]*Ru/Sptr[i].MW;
            // p = rhoMixOutB*RMix*T;

            for(int i=0;i<NS;i++) YsCell[i] = YsInB[i];

            RMix = 0.;
            for(int i=0;i<NS;i++)RMix+=YsCell[i]*Ru/Sptr[i].MW;

            double rhoMixCell = p/RMix/T;
            double rhoArrayCell[NS];
            for(int i=0;i<NS;i++)rhoArrayCell[i]=rhoMixCell*YsCell[i];

            m_mixture->setDensityArray(rhoArrayCell);// reset the density array 
            m_mixture->setDensity(rhoMixCell);
            m_mixture->setTemperature(T);
            m_mixture->setPressure(p);
            m_mixture->updateMolarFractionArray();
            m_mixture->setU(0.);
            m_mixture->setV(0.);
            m_mixture->setW(0.);
          }
        }

        // 1D H2 detonation flame, Paolucci ------ 
        if(AMRPara::test_prob == 1) // for example_cases/1D_detonation
        {
          std::array<double,NS> rhoHigh, rhoLow;
          double pHigh, pLow, uHigh, uLow;
          double YsCell[NS]; 
          double rho, p, u;

          double x_loc = coordinates.getX();
          double x0(0.15);

          alpha = std::tanh(abs(x_loc-x0)/1.2e-4);

          rhoHigh = mixShock->getDensityArray();
          rhoLow = mixLP->getDensityArray();
          pHigh = mixShock->getPressure();
          pLow = mixLP->getPressure();
          uHigh = mixShock->getVelocity().getX();
          uLow = mixLP->getVelocity().getX();

          double rhoArrayCell[NS];
          for(int i=0;i<NS;i++)
          {
            rhoArrayCell[i]=0.5*(rhoHigh[i]+rhoLow[i]-(rhoHigh[i]-rhoLow[i])*alpha);
          }

          // get rho, p, u
          double rhoMixCell(0);
          for(int i=0;i<NS;i++)
          {
            rhoMixCell += rhoArrayCell[i];
          }

          p = 0.5*(pHigh+pLow-(pHigh-pLow)*alpha);
          u = 0.5*(uHigh+uLow-(uHigh-uLow)*alpha);

          for(int i=0;i<NS;i++) 
          {
            YsCell[i] = rhoArrayCell[i]/rhoMixCell;
          }

          double RMix = 0.;
          for(int i=0;i<NS;i++)RMix+=YsCell[i]*Ru/Sptr[i].MW;

          double T = p/rhoMixCell/RMix;

          m_mixture->setDensityArray(rhoArrayCell);// reset the density array 
          m_mixture->setDensity(rhoMixCell);
          m_mixture->setTemperature(T);
          m_mixture->setPressure(p);
          m_mixture->updateMolarFractionArray();
          m_mixture->setU(u);
        }

        // for canonical shock-bubble problems
        if(AMRPara::test_prob == 2 || AMRPara::test_prob == 3) 
        {
          alpha = 1;

          Coord dist = coordinates - center;
          double distNorm = sqrt(dist.getX()*dist.getX()+dist.getY()*dist.getY()+dist.getZ()*dist.getZ());

          //for RSBI
          // if(interface_flag==0) alpha = (tanh(( - radius + distNorm)*C) + 1.0) / 2.0;
          alpha = (tanh(( - radius + distNorm)*C) + 1.0) / 2.0;

          std::array<double,NS> rhoInB, rhoOutB;
          double rhoMixInB(0.), rhoMixOutB(0.), YsInB[NS], YsOutB[NS], YsCell[NS]; 
          double pInB,pOutB,TInB,TOutB,p,T;
          Mixture tempMix; 

          TInB = mixBubble->getTemperature();
          TOutB = mixLP->getTemperature();
          pInB = mixBubble->getPressure();
          pOutB = mixLP->getPressure();

          rhoOutB = mixLP->getDensityArray();
          rhoInB = mixBubble->getDensityArray();

          for(int i=0;i<NS;i++) {rhoMixInB+=rhoInB[i];rhoMixOutB+=rhoOutB[i];}
          for(int i=0;i<NS;i++){YsInB[i]=rhoInB[i]/rhoMixInB; YsOutB[i]=rhoOutB[i]/rhoMixOutB;}

          double RMix(0.);
          T = alpha*TOutB+(1-alpha)*TInB;
          p = alpha*pOutB+(1-alpha)*pInB;
          for(int i=0;i<NS;i++)RMix+=YsOutB[i]*Ru/Sptr[i].MW;
          // p = rhoMixOutB*RMix*T;

          for(int i=0;i<NS;i++) YsCell[i] = alpha*YsOutB[i]+(1-alpha)*YsInB[i];

          RMix = 0.;
          for(int i=0;i<NS;i++)RMix+=YsCell[i]*Ru/Sptr[i].MW;

          double rhoMixCell = p/RMix/T;
          double rhoArrayCell[NS];
          for(int i=0;i<NS;i++)rhoArrayCell[i]=rhoMixCell*YsCell[i];

          m_mixture->setDensityArray(rhoArrayCell);// reset the density array 
          m_mixture->setDensity(rhoMixCell);
          m_mixture->setTemperature(T);
          m_mixture->setPressure(p);
          m_mixture->updateMolarFractionArray();
        }

      }
    }
  }
}


//***********************************************************************

double Cell::selectMassi(int i) const
{
  if(i<NS && i>=0)
  {
    return m_mixture->getDensityArray()[i]; 
  }
  else
  {
    Errors::errorMessage("Invalid species index in Cell::selectMassi().");
  }
}

//***********************************************************************

double Cell::calcTauC()
{
  double res;
  Coord div_v;
  Coord gradU, gradV, gradW;

  gradU = m_vecQuantitiesAddPhys[0]->getGrad(1);
  gradV = m_vecQuantitiesAddPhys[0]->getGrad(2);
  gradW = m_vecQuantitiesAddPhys[0]->getGrad(3);

  div_v.setXYZ(gradU.getX(), gradV.getY(), gradW.getZ());

  double  cell_l, f;
  double expo = AMRPara::dx_expo;
  if (AMRPara::dim == 2) // for 2d square root
  {
    cell_l = sqrt(getElement()->getSizeX()*getElement()->getSizeY());
    f = pow(cell_l, expo);
    res = sqrt(gradU.getX()*gradU.getX() + gradV.getY()*gradV.getY())*f;
  }
  else if (AMRPara::dim == 3)
  {
    cell_l = pow(getElement()->getSizeX()*getElement()->getSizeY()*getElement()->getSizeZ(), 1.0/3);
    f = pow(cell_l, expo);
    res = div_v.norm()*f;
  }

  return res;
}

//***********************************************************************

double Cell::calcTauR()
{
  double res;
  Coord vor_v;
  Coord gradU, gradV, gradW;

  gradU = m_vecQuantitiesAddPhys[0]->getGrad(1);
  gradV = m_vecQuantitiesAddPhys[0]->getGrad(2);
  gradW = m_vecQuantitiesAddPhys[0]->getGrad(3);

  vor_v.setXYZ(gradW.getY()-gradV.getZ(), gradU.getZ()-gradW.getX(), gradV.getX()-gradU.getY());

  double  cell_l, f; 
  double expo = AMRPara::dx_expo;
  if (AMRPara::dim == 2) // for 2d square root
  {
    cell_l = sqrt(getElement()->getSizeX()*getElement()->getSizeY());
    f = pow(cell_l, expo);
    res = abs(gradU.getY() - gradV.getX())*f;
  }
  else if (AMRPara::dim == 3)
  {
    cell_l = pow(getElement()->getSizeX()*getElement()->getSizeY()*getElement()->getSizeZ(), 1.0/3);
    f = pow(cell_l, expo);
    res = vor_v.norm()*f;
  }
  return res;
}

//***********************************************************************
double Cell::calcTauGradRho()
{
  double res;
  double f1 = AMRPara::f1;
  double f2 = AMRPara::f2;

  double gradRhoNorm = this->getDensityGradient();

  double  cell_l;
  double expo = AMRPara::dx_expo;
  if (AMRPara::dim == 2) // for 2d square root
  {
    cell_l = sqrt(getElement()->getSizeX()*getElement()->getSizeY());
    res = gradRhoNorm*pow(cell_l, expo);
  }
  else if (AMRPara::dim == 3)
  {
    cell_l = pow(getElement()->getSizeX()*getElement()->getSizeY()*getElement()->getSizeZ(), 1.0/3);
    res = gradRhoNorm*pow(cell_l, expo);
  }

  return res;
}

//***********************************************************************
// compute gradients of velocity and temperature for viscous flux evaluation.
// add gradients of pressure and molar fraction for pressure-dependent diffusion
void Cell::computeGradient(std::vector<Coord> &grads, std::vector<Variable> &nameVariables, std::vector<int> &numPhases)
{
  int typeCellInterface(0);
  double cg(0.), cd(0.), gradCellInterface(0.);
  double distance(0.), distanceX(0.), distanceY(0.), distanceZ(0.);
  double sumDistanceX(0.), sumDistanceY(0.), sumDistanceZ(0.);
  for (unsigned int g = 0; g < grads.size(); g++) { grads[g] = 0.; }

  for (unsigned int b = 0; b < m_cellInterfaces.size(); b++) {
    if (!m_cellInterfaces[b]->getSplit()) {
      typeCellInterface = m_cellInterfaces[b]->whoAmI();
      if (typeCellInterface == 0) //Cell interface type CellInterface/O2
      {
        // Sum for each cell interface with ponderation using distance in each direction
        // then the cell gradient is normalized by sum of distances.
        distanceX = std::fabs(m_cellInterfaces[b]->getCellGauche()->distanceX(m_cellInterfaces[b]->getCellDroite()));
        distanceY = std::fabs(m_cellInterfaces[b]->getCellGauche()->distanceY(m_cellInterfaces[b]->getCellDroite()));
        distanceZ = std::fabs(m_cellInterfaces[b]->getCellGauche()->distanceZ(m_cellInterfaces[b]->getCellDroite()));
        sumDistanceX += distanceX;
        sumDistanceY += distanceY;
        sumDistanceZ += distanceZ;
        distance = m_cellInterfaces[b]->getCellGauche()->distance(m_cellInterfaces[b]->getCellDroite());

        // Extracting left and right variables values for each cell interface
        // and calculus of the gradients normal to the face
        for (unsigned int g = 0; g < grads.size(); g++) {
          cg = m_cellInterfaces[b]->getCellGauche()->selectScalar(nameVariables[g], numPhases[g]);
          cd = m_cellInterfaces[b]->getCellDroite()->selectScalar(nameVariables[g], numPhases[g]);
          gradCellInterface = (cd - cg) / distance;

          // Projection in the absolute system of coordinate
          grads[g].setXYZ(grads[g].getX() + m_cellInterfaces[b]->getFace()->getNormal().getX()*gradCellInterface*distanceX,
                          grads[g].getY() + m_cellInterfaces[b]->getFace()->getNormal().getY()*gradCellInterface*distanceY,
                          grads[g].getZ() + m_cellInterfaces[b]->getFace()->getNormal().getZ()*gradCellInterface*distanceZ);
        }
      } // end inner cell interface if
      else {
        distanceX = std::fabs(this->distanceX(m_cellInterfaces[b])) * 2.;
        distanceY = std::fabs(this->distanceY(m_cellInterfaces[b])) * 2.;
        distanceZ = std::fabs(this->distanceZ(m_cellInterfaces[b])) * 2.;
        sumDistanceX += distanceX;
        sumDistanceY += distanceY;
        sumDistanceZ += distanceZ;
        distance = this->distance(m_cellInterfaces[b]);
          
        for (unsigned int g = 0; g < grads.size(); g++) {
          if (nameVariables[g] == velocityU || nameVariables[g] == velocityV || nameVariables[g] == velocityW) {
            Coord face_normal, velocity;
            face_normal =  m_cellInterfaces[b]->getFace()->getNormal();
            velocity = this->getMixture()->getVelocity();
            double v_normal = velocity.getX()*face_normal.getX()+velocity.getY()*face_normal.getY()+velocity.getZ()*face_normal.getZ();
            gradCellInterface = -1*v_normal / distance;
            

            if (typeCellInterface == 6 || typeCellInterface == 2) { //Cell Interface of type Symmetry
              // Multiplication of the gradient by the normal direction to guarantee symmetry
              if (nameVariables[g] == velocityU) { gradCellInterface = gradCellInterface * m_cellInterfaces[b]->getFace()->getNormal().getX(); }
              if (nameVariables[g] == velocityV) { gradCellInterface = gradCellInterface * m_cellInterfaces[b]->getFace()->getNormal().getY(); }
              if (nameVariables[g] == velocityW) { gradCellInterface = gradCellInterface * m_cellInterfaces[b]->getFace()->getNormal().getZ(); }

              // Projection in the absolute system of coordinate +
              grads[g].setXYZ(grads[g].getX() + m_cellInterfaces[b]->getFace()->getNormal().getX()*gradCellInterface*distanceX,
                              grads[g].getY() + m_cellInterfaces[b]->getFace()->getNormal().getY()*gradCellInterface*distanceY,
                              grads[g].getZ() + m_cellInterfaces[b]->getFace()->getNormal().getZ()*gradCellInterface*distanceZ);
            }
            else if (typeCellInterface == 3) { //Cell interface of type Wall
              // Projection in the absolute system of coordinate
              grads[g].setXYZ(grads[g].getX() + m_cellInterfaces[b]->getFace()->getNormal().getX()*gradCellInterface*distanceX,
                              grads[g].getY() + m_cellInterfaces[b]->getFace()->getNormal().getY()*gradCellInterface*distanceY,
                              grads[g].getZ() + m_cellInterfaces[b]->getFace()->getNormal().getZ()*gradCellInterface*distanceZ);
            }
            else if (typeCellInterface == 1){ // non-reflecting, no flux at such face
              // double cell_x = this->getElement()->getPosition().getX();
              // grads[g].setXYZ(0,0,0);
              // int ti = 0;
            }
          }
          else // for temperature field, gradient on this interface is 0
          {
            // for temperture, no flux at the B.C. face
            // grads[g].setXYZ(0,0,0);
          }
        } // end grad(g), g for u,v,w,T loop
      } // end boundary interface if
    } // end interface.split() if
  }  // end m_interface loop

  // Verifications in multiD
  if (sumDistanceX <= 1.e-12) { sumDistanceX = 1.; }
  if (sumDistanceY <= 1.e-12) { sumDistanceY = 1.; }
  if (sumDistanceZ <= 1.e-12) { sumDistanceZ = 1.; }

  // Final normalized gradient on the cell
  for (unsigned int g = 0; g < grads.size(); g++) {
    grads[g].setXYZ(grads[g].getX() / sumDistanceX, grads[g].getY() / sumDistanceY, grads[g].getZ() / sumDistanceZ);
  }
}

//***********************************************************************

void Cell::computeRhoIGradient(std::vector<Coord> &grads)
{
  int typeCellInterface(0);
  double cg(0.), cd(0.), gradCellInterface(0.);
  double distance(0.), distanceX(0.), distanceY(0.), distanceZ(0.);
  double sumDistanceX(0.), sumDistanceY(0.), sumDistanceZ(0.);
  Variable flag = density;
  for (unsigned int g = 0; g < grads.size(); g++) { grads[g] = 0.; }

  for (unsigned int b = 0; b < m_cellInterfaces.size(); b++) {
    if (!m_cellInterfaces[b]->getSplit()) {
      typeCellInterface = m_cellInterfaces[b]->whoAmI();
      if (typeCellInterface == 0) //Cell interface type CellInterface/O2, for boundary cellInterface, gradRhoi = 0
      {
        // Sum for each cell interface with ponderation using distance in each direction
        // then the cell gradient is normalized by sum of distances.
        distanceX = std::fabs(m_cellInterfaces[b]->getCellGauche()->distanceX(m_cellInterfaces[b]->getCellDroite()));
        distanceY = std::fabs(m_cellInterfaces[b]->getCellGauche()->distanceY(m_cellInterfaces[b]->getCellDroite()));
        distanceZ = std::fabs(m_cellInterfaces[b]->getCellGauche()->distanceZ(m_cellInterfaces[b]->getCellDroite()));
        sumDistanceX += distanceX;
        sumDistanceY += distanceY;
        sumDistanceZ += distanceZ;
        distance = m_cellInterfaces[b]->getCellGauche()->distance(m_cellInterfaces[b]->getCellDroite());

        // Extracting left and right variables values for each cell interface
        // and calculus of the gradients normal to the face
        for (unsigned int g = 0; g < grads.size(); g++) {
          cg = m_cellInterfaces[b]->getCellGauche()->selectScalar(flag, g);
          cd = m_cellInterfaces[b]->getCellDroite()->selectScalar(flag, g);
          gradCellInterface = (cd - cg) / distance;

          // Projection in the absolute system of coordinate
          grads[g].setXYZ(grads[g].getX() + m_cellInterfaces[b]->getFace()->getNormal().getX()*gradCellInterface*distanceX,
                          grads[g].getY() + m_cellInterfaces[b]->getFace()->getNormal().getY()*gradCellInterface*distanceY,
                          grads[g].getZ() + m_cellInterfaces[b]->getFace()->getNormal().getZ()*gradCellInterface*distanceZ);
        }
      }
      else {
        distanceX = std::fabs(this->distanceX(m_cellInterfaces[b])) * 2.;
        distanceY = std::fabs(this->distanceY(m_cellInterfaces[b])) * 2.;
        distanceZ = std::fabs(this->distanceZ(m_cellInterfaces[b])) * 2.;
        sumDistanceX += distanceX;
        sumDistanceY += distanceY;
        sumDistanceZ += distanceZ;
        distance = this->distance(m_cellInterfaces[b]);

        // for nonreflecting, symmetry, wall, outflow B.C., no flux at such interfaces.
        // gradient on this interface is 0
        if(typeCellInterface == 1 || typeCellInterface == 6 || typeCellInterface == 2 || typeCellInterface == 3)
        {
          for (unsigned int g = 0; g < grads.size(); g++) {
            // grads[g].setXYZ(0, 0, 0);
          }
        }
        else
        {
          // printf("B.C. type = %d, this B.C. is not supported in Cell::computeRhoIGradient().\n", typeCellInterface);
          // exit(EXIT_FAILURE);
        }

      }
    }
  }

  // Verifications in multiD
  if (sumDistanceX <= 1.e-12) { sumDistanceX = 1.; }
  if (sumDistanceY <= 1.e-12) { sumDistanceY = 1.; }
  if (sumDistanceZ <= 1.e-12) { sumDistanceZ = 1.; }

  // Final normalized gradient on the cell
  for (unsigned int g = 0; g < grads.size(); g++) {
    grads[g].setXYZ(grads[g].getX() / sumDistanceX, grads[g].getY() / sumDistanceY, grads[g].getZ() / sumDistanceZ);
  }
}

//***********************************************************************
// compute gradients of molar fraction for pressure-dependent diffusion

void Cell::computeXGradient(std::vector<Coord> &grads)
{
  int typeCellInterface(0);
  double cg(0.), cd(0.), gradCellInterface(0.);
  double distance(0.), distanceX(0.), distanceY(0.), distanceZ(0.);
  double sumDistanceX(0.), sumDistanceY(0.), sumDistanceZ(0.);
  Variable flag = molarFraction;
  for (unsigned int g = 0; g < grads.size(); g++) { grads[g] = 0.; }

  for (unsigned int b = 0; b < m_cellInterfaces.size(); b++) {
    if (!m_cellInterfaces[b]->getSplit()) {
      typeCellInterface = m_cellInterfaces[b]->whoAmI();
      if (typeCellInterface == 0) //Cell interface type CellInterface/O2, for boundary cellInterface, gradRhoi = 0
      {
        // Sum for each cell interface with ponderation using distance in each direction
        // then the cell gradient is normalized by sum of distances.
        distanceX = std::fabs(m_cellInterfaces[b]->getCellGauche()->distanceX(m_cellInterfaces[b]->getCellDroite()));
        distanceY = std::fabs(m_cellInterfaces[b]->getCellGauche()->distanceY(m_cellInterfaces[b]->getCellDroite()));
        distanceZ = std::fabs(m_cellInterfaces[b]->getCellGauche()->distanceZ(m_cellInterfaces[b]->getCellDroite()));
        sumDistanceX += distanceX;
        sumDistanceY += distanceY;
        sumDistanceZ += distanceZ;
        distance = m_cellInterfaces[b]->getCellGauche()->distance(m_cellInterfaces[b]->getCellDroite());

        // Extracting left and right variables values for each cell interface
        // and calculus of the gradients normal to the face
        for (unsigned int g = 0; g < grads.size(); g++) {
          cg = m_cellInterfaces[b]->getCellGauche()->selectScalar(flag, g);
          cd = m_cellInterfaces[b]->getCellDroite()->selectScalar(flag, g);
          gradCellInterface = (cd - cg) / distance;

          // Projection in the absolute system of coordinate
          grads[g].setXYZ(grads[g].getX() + m_cellInterfaces[b]->getFace()->getNormal().getX()*gradCellInterface*distanceX,
                          grads[g].getY() + m_cellInterfaces[b]->getFace()->getNormal().getY()*gradCellInterface*distanceY,
                          grads[g].getZ() + m_cellInterfaces[b]->getFace()->getNormal().getZ()*gradCellInterface*distanceZ);
        }
      }
      else {
        distanceX = std::fabs(this->distanceX(m_cellInterfaces[b])) * 2.;
        distanceY = std::fabs(this->distanceY(m_cellInterfaces[b])) * 2.;
        distanceZ = std::fabs(this->distanceZ(m_cellInterfaces[b])) * 2.;
        sumDistanceX += distanceX;
        sumDistanceY += distanceY;
        sumDistanceZ += distanceZ;
        distance = this->distance(m_cellInterfaces[b]);

        // for nonreflecting, symmetry, wall, outflow B.C., no flux at such interfaces.
        // gradient on this interface is 0
        if(typeCellInterface == 1 || typeCellInterface == 6 || typeCellInterface == 2 || typeCellInterface == 3)
        {
          for (unsigned int g = 0; g < grads.size(); g++) {
            // grads[g].setXYZ(0, 0, 0);
          }
        }
        else
        {
          // printf("B.C. type = %d, this B.C. is not supported in Cell::computeRhoIGradient().\n", typeCellInterface);
          // exit(EXIT_FAILURE);
        }

      }
    }
  }

  // Verifications in multiD
  if (sumDistanceX <= 1.e-12) { sumDistanceX = 1.; }
  if (sumDistanceY <= 1.e-12) { sumDistanceY = 1.; }
  if (sumDistanceZ <= 1.e-12) { sumDistanceZ = 1.; }

  // Final normalized gradient on the cell
  for (unsigned int g = 0; g < grads.size(); g++) {
    grads[g].setXYZ(grads[g].getX() / sumDistanceX, grads[g].getY() / sumDistanceY, grads[g].getZ() / sumDistanceZ);
  }
}

//***********************************************************************
void Cell::computeTransportProperties()
{
  //m_split
  // double Tg(0.), Td(0.), Pg(0.), Pd(0.);
  double T(0.),P(0.),Ys[NS],trans[NS+2],rhoMix(0.);
  std::array<double,NS> mass;
  double distance(0.), distanceX(0.), distanceY(0.), distanceZ(0.);
  if(!m_split)
  {
    T = m_mixture->getTemperature();
    P = m_mixture->getPressure();
    mass = m_mixture->getDensityArray();

    // prepare total density and mass fraction
    for(int i=0;i<NS;i++){rhoMix+=mass[i];}
    for(int i=0;i<NS;i++){Ys[i]=mass[i]/rhoMix;}
    transport_model(Ys,T,P,rhoMix,trans);

    // the transport property at the interface is avaraged from left/right cell values
    for(int i=0;i<NS+2;i++)m_transportProperties[i] = trans[i];
  }
}

//***********************************************************************

std::array<double,NS+2> const & Cell::getTransportProperties() const
{
  return m_transportProperties;
}

//***********************************************************************

void Cell::setTransportProperties(const std::array<double,NS+2> transArray)
{
  for(int i=0;i<NS+2;i++)m_transportProperties[i] = transArray[i];
}

//***********************************************************************
// return gradient of massFraction i
Coord Cell::getGradientYI(int i)
{
  if(i>=NS)
  {
    printf("In Cell::getGradientYI(), Invalid species index = %d.\n",i);
    exit(EXIT_FAILURE);
  }
  
  if(m_model->whoAmI()!="FIRE")
  {
    printf("The Cell::getGradientYI() only supports FIRE model.\n");
    exit(EXIT_FAILURE);
  }
  else
  {
    int var = i; //the ith species
    return this->computeGradient(massFraction, var);
  }
}

//***********************************************************************
// return gradient of total density
double Cell::getDensityGradient()
{
  {
    int var = -1; //for fire model
    return this->computeGradient(density, var).norm();
  }
}

//***********************************************************************
// return gradient of U
Coord Cell::getGradientU()
{  
  if(m_model->whoAmI()!="FIRE")
  {
    printf("The Cell::getGradientU() only supports FIRE model.\n");
    exit(EXIT_FAILURE);
  }
  else
  {
    return this->computeGradient(velocityU, 2);
  }
}

//***********************************************************************
// return gradient of V
Coord Cell::getGradientV()
{  
  if(m_model->whoAmI()!="FIRE")
  {
    printf("The Cell::getGradientV() only supports FIRE model.\n");
    exit(EXIT_FAILURE);
  }
  else
  {
    return this->computeGradient(velocityV, 2);
  }
}

//***********************************************************************
// return gradient of W
Coord Cell::getGradientW()
{  
  if(m_model->whoAmI()!="FIRE")
  {
    printf("The Cell::getGradientW() only supports FIRE model.\n");
    exit(EXIT_FAILURE);
  }
  else
  {
    return this->computeGradient(velocityW, 2);
  }
}

//***********************************************************************
// return gradient of induction time
Coord Cell::getGradientInductionTime()
{  
  if(m_model->whoAmI()!="FIRE")
  {
    printf("The Cell::getGradientInductionTime() only supports FIRE model.\n");
    exit(EXIT_FAILURE);
  }
  else
  {
    return this->computeGradient(inductionTime, 2);
  }
}

//***********************************************************************
// return gradient of P
Coord Cell::getGradientP()
{  
  if(m_model->whoAmI()!="FIRE")
  {
    printf("The Cell::getGradientP() only supports FIRE model.\n");
    exit(EXIT_FAILURE);
  }
  else
  {
    return this->computeGradient(pressure, 2);
  }
}

//***********************************************************************

void Cell::fillBufferVector(double *buffer, int &counter, const int &lvl, const int &neighbour, const int &dim, Variable nameVector, int num, int index) const
{
  if(nameVector != transportArray)
  {
    if (m_lvl == lvl) {
      buffer[++counter] = this->selectVector(nameVector, num, index).getX();
      if (dim > 1) buffer[++counter] = this->selectVector(nameVector, num, index).getY();
      if (dim > 2) buffer[++counter] = this->selectVector(nameVector, num, index).getZ();
    }
    else {
      for (unsigned int i = 0; i < m_childrenCells.size(); i++) {
        if (m_childrenCells[i]->hasNeighboringGhostCellOfCPUneighbour(neighbour)) {
          m_childrenCells[i]->fillBufferVector(buffer, counter, lvl, neighbour, dim, nameVector, num, index);
        }
      }
    }
  }
  else
  {
    if (m_lvl == lvl) {
      for(int i=0;i<NS+2;i++)
      {
        buffer[++counter] = m_transportProperties[i];
      }
    }
    else {
      for (unsigned int i = 0; i < m_childrenCells.size(); i++) {
        if (m_childrenCells[i]->hasNeighboringGhostCellOfCPUneighbour(neighbour)) {
          m_childrenCells[i]->fillBufferVector(buffer, counter, lvl, neighbour, dim, nameVector, num, index);
        }
      }
    }
  }
}

//***********************************************************************

void Cell::getBufferVector(double *buffer, int &counter, const int &lvl, const int &dim, Variable nameVector, int num, int index)
{
  if(nameVector != transportArray)
  {
    if (m_lvl == lvl) {
      Coord temp;

      temp.setX(buffer[++counter]);
      if (dim > 1) temp.setY(buffer[++counter]);
      if (dim > 2) temp.setZ(buffer[++counter]);
      this->setVector(nameVector, temp, num, index);
    }
    else {
      for (unsigned int i = 0; i < m_childrenCells.size(); i++) {
        m_childrenCells[i]->getBufferVector(buffer, counter, lvl, dim, nameVector, num, index);
      }
    }
  }
  else
  {
    if (m_lvl == lvl) {
      for(int i=0;i<NS+2;i++)
      {
        m_transportProperties[i] = buffer[++counter];
      }
    }
    else {
      for (unsigned int i = 0; i < m_childrenCells.size(); i++) {
        m_childrenCells[i]->getBufferVector(buffer, counter, lvl, dim, nameVector, num, index);
      }
    }
  }
}