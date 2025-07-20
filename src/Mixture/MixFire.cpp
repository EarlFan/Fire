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

//! \file      MixFire.cpp
//! \author    E Fan, Lisong Shi, Jiaao Hao, Tianhan Zhang, Chih-yung Wen.
//! \version   1.0
//! \date      Jan 01 2025

#include <cmath>
#include "MixFire.h"

using namespace std;
using namespace tinyxml2;

//***************************************************************************

MixFire::MixFire():m_totalDensity(0.),
    m_pressure(0.),m_temperature(0.),m_energie(0.),m_totalEnergy(0.),m_hf_change(0.),m_rhoExitedOH(0.),m_pMax(0.)
{
  for(int ll=0;ll<NS;ll++)m_densityArray[ll]=0.;
  for(int ll=0;ll<NS;ll++)m_densityArray_change[ll]=0;
}

//***************************************************************************

MixFire::MixFire(XMLElement *state, std::string fileName) :
  m_energie(0.), m_totalEnergy(0.),m_hf_change(0.),m_rhoExitedOH(0.)
{
  XMLElement *sousElement(state->FirstChildElement("mixture"));
  if (sousElement == NULL) throw ErrorXMLElement("mixture", fileName, __FILE__, __LINE__);
  //Attributes reading
  //------------------
  XMLError error;

  XMLElement *dataMix(sousElement->FirstChildElement("dataMix"));
  if (dataMix == NULL) throw ErrorXMLElement("dataMix", fileName, __FILE__, __LINE__);
  //Attributes reading
  //------------------
  //temperature
  error = dataMix->QueryDoubleAttribute("temperature", &m_temperature);
  if (error != XML_NO_ERROR) throw ErrorXMLAttribut("temperature", fileName, __FILE__, __LINE__);
  ////temperature
  //error = dataMix->QueryDoubleAttribute("temperature", &m_temperature);
  //if (error != XML_NO_ERROR) throw ErrorXMLAttribut("temperature", fileName, __FILE__, __LINE__);

  //set velocity
  XMLElement *velocity(sousElement->FirstChildElement("velocity"));
  if (velocity == NULL) throw ErrorXMLElement("velocity", fileName, __FILE__, __LINE__);
  double velocityX(0.), velocityY(0.), velocityZ(0.);
  error = velocity->QueryDoubleAttribute("x", &velocityX);
  if (error != XML_NO_ERROR) throw ErrorXMLAttribut("x", fileName, __FILE__, __LINE__);
  error = velocity->QueryDoubleAttribute("y", &velocityY);
  if (error != XML_NO_ERROR) throw ErrorXMLAttribut("y", fileName, __FILE__, __LINE__);
  error = velocity->QueryDoubleAttribute("z", &velocityZ);
  if (error != XML_NO_ERROR) throw ErrorXMLAttribut("z", fileName, __FILE__, __LINE__);
  m_velocity.setXYZ(velocityX, velocityY, velocityZ);

  //density array
  XMLElement *densityArray(sousElement->FirstChildElement("densityArray"));
  if (densityArray == NULL) throw ErrorXMLElement("densityArray", fileName, __FILE__, __LINE__);
  double temp[NS];
  for(int ll=0;ll<NS;ll++)
  {
    error = densityArray->QueryDoubleAttribute(Sptr[ll].Name.c_str(), &temp[ll]);
    if (error != XML_NO_ERROR) temp[ll] = 0;
    // if (error != XML_NO_ERROR) throw ErrorXMLAttribut(Sptr[ll].Name.c_str(), fileName, __FILE__, __LINE__);
  }
  for(int ll=0;ll<NS;ll++)m_densityArray[ll]=temp[ll];

  //set total density
  m_totalDensity = sumNS(temp);
  
  // set pressure
  m_pressure = this->computePressure();

  updateMolarFractionArray();

  for(int ll=0;ll<NS;ll++)m_densityArray_change[ll]=0;
  m_pMax = m_pressure;
}

MixFire::~MixFire(){}

//***************************************************************************

void MixFire::allocateAndCopyMixture(Mixture **mixture)
{
  *mixture = new MixFire(*this);
}

//***************************************************************************

void MixFire::copyMixture(Mixture &mixture)
{
  m_totalDensity = mixture.getDensity();
  m_densityArray = mixture.getDensityArray();
  m_velocity = mixture.getVelocity();
  m_pressure = mixture.getPressure();
  m_temperature = mixture.getTemperature();           
  m_energie = mixture.getEnergy();
  m_totalEnergy = mixture.getTotalEnergy();
  m_frozenSoundSpeed = mixture.getMixSoundSpeed();
  m_rhoExitedOH = mixture.getRhoExcitedOH();
  m_pMax = mixture.getMaxPressure();
}

//***************************************************************************

double MixFire::computeDensity(const double *alphak, 
    const double *rhok, const int &numberPhases) 
{
  printf("Invalid call of computeDensity.\n");
  exit(EXIT_FAILURE);
};

//***************************************************************************
// p_sum = sum{p_i}
double MixFire::computePressure(const double *alphak, const double *pk, const int &numberPhases)
{
  printf("invalid call of MixFire::computePressure(alpha).\n");
  exit(EXIT_FAILURE);
}


//***************************************************************************

double MixFire::computeInternalEnergy(const double *Yk, const double *ek, const int &numberPhases)
{
  double e(0.);
  for (int k = 0; k<numberPhases; k++)
  {
    e += Yk[k] * ek[k];
  }
  return e;
}

//***************************************************************************
// delete - fane
double MixFire::computeFrozenSoundSpeed(const double *Yk, const double *ck, const int &numberPhases)
{
  printf("invalid call of MixFire::computeFrozenSoundSpeed(alphak).\n");
  exit(EXIT_FAILURE);
}

//***************************************************************************
// temperature is known
void MixFire::computeMixtureVariables(Phase **vecPhase, const int &numberPhases)
{
  m_totalDensity = sumNS(m_densityArray);

  if(m_totalDensity<=0)
  {
    printf("Invalid density at MixFire::computeMixtureVariables(), density = %lf.\n",m_totalDensity);
    exit(EXIT_FAILURE);
  }

  //Mass fraction
  for (int k = 0; k < numberPhases; k++) {
    vecPhase[k]->setDensity(m_densityArray[k]);
    vecPhase[k]->computeMassFraction(m_totalDensity);
  }
  

  m_temperature=this->computeTem();

  m_energie=this->computeInternalEnergy();
  m_frozenSoundSpeed = this->computeFrozenSoundSpeed();
}

//***************************************************************************

void MixFire::internalEnergyToTotalEnergy(vector<QuantitiesAddPhys*> &vecGPA)
{
  m_totalEnergy = m_energie + 0.5*m_velocity.squaredNorm();
  for (unsigned int pa = 0; pa < vecGPA.size(); pa++) {
    m_totalEnergy += vecGPA[pa]->computeEnergyAddPhys() / m_totalDensity; //Caution /m_totalDensity important
  }
}

//***************************************************************************

void MixFire::totalEnergyToInternalEnergy(vector<QuantitiesAddPhys*> &vecGPA)
{

  double rhoMix = sumNS(m_densityArray);

  if(rhoMix<=0)
  {
    printf("Invalid density at MixFire::totalEnergyToInternalEnergy(), density = %lf.\n",rhoMix);
    printf("rho0 = %le, rho1 = %le.\n",m_densityArray[0],m_densityArray[1]);
    exit(EXIT_FAILURE);
  }

  m_energie = m_totalEnergy - 0.5*m_velocity.squaredNorm();
  for (unsigned int pa = 0; pa < vecGPA.size(); pa++) {
    m_energie -= vecGPA[pa]->computeEnergyAddPhys() / m_totalDensity; //Caution /m_totalDensity important
  }
}

//***************************************************************************

void MixFire::localProjection(const Coord &normal, const Coord &tangent, const Coord &binormal)
{
  m_velocity.localProjection(normal, tangent, binormal);
}

//***************************************************************************

void MixFire::reverseProjection(const Coord &normal, const Coord &tangent, const Coord &binormal)
{
  m_velocity.reverseProjection(normal, tangent, binormal);
}

//***************************************************************************

double MixFire::computeTotalDensity()
{
  return sumNS(m_densityArray);
}

//***************************************************************************
// cons to pris
// Updating: T, p, c
void MixFire::con2pri()
{
  // rhoMix, rhoArray, velocity, m_totalEnergy already set in FluxEulerTPIG::buildPrim()
  m_energie = m_totalEnergy - 0.5*m_velocity.squaredNorm();

  double Ys[NS];
  for(int ll=0;ll<NS;ll++)Ys[ll] = m_densityArray[ll]/m_totalDensity;


  // call the method in Cantera
  std::shared_ptr<Cantera::ThermoPhase> gas = solution_mech->thermo();
  gas->setMassFractions(Ys);
  gas->setState_UV(m_energie, 1/m_totalDensity);
  m_temperature = gas->temperature();

  //printf("In MixFire::After con2pri(), T = %le.\n",m_temperature);
  if(m_temperature<=0 || m_totalDensity<=0)
  {
    printf("MixFire::con2pri, T = %lf,rhoMix = %lf.\n",m_temperature,m_totalDensity);
    exit(EXIT_FAILURE);
  }

  m_pressure = this->computePressure();
  m_frozenSoundSpeed = this->computeFrozenSoundSpeed();

  if(m_pressure > m_pMax) m_pMax = m_pressure;

  updateMolarFractionArray();
}

//***************************************************************************
// update m_molarFractionArray
// Updating: T, p, c
void MixFire::updateMolarFractionArray()
{
  std::array<double, NS> molarFractions;
  double totalMoles = 0.0;
  for (int i = 0; i < NS; ++i) {
      totalMoles += m_densityArray[i] / Sptr[i].MW;
  }
  for (int i = 0; i < NS; ++i) {
      molarFractions[i] = (m_densityArray[i] / Sptr[i].MW) / totalMoles;
  }

  m_molarFractionArray = molarFractions;
}

//****************************************************************************
//******************************* THERMOS **********************************
//****************************************************************************
// compute RMix using primitive variables
double MixFire::computeRMix()
{
  for(int i=0;i<NS;i++)
  {
    if(m_densityArray[i]<1e-30)
    {
      m_densityArray[i]=0.;
    }
  }
  double Ys[NS];
  for(int i=0;i<NS;i++)Ys[i] = m_densityArray[i]/m_totalDensity;
  double RMix = 0;
  for(int i=0;i<NS;i++)RMix = RMix + Ys[i]*Ru/Sptr[i].MW;
  if(std::isnan(RMix))
  {
    for(int i=0;i<NS;i++)printf("Ys[%d] = %lf, MW=%lf.\n",i,Ys[i],Sptr[i].MW);
    exit(EXIT_FAILURE);
  }
  return RMix;
}

//***************************************************************************

// compute pressure
double MixFire::computePressure()
{
  return m_totalDensity*this->computeRMix()*m_temperature;
}

//***************************************************************************

// compute temperature using primitive variables
double MixFire::computeTem()
{
  double temp;
  m_totalDensity = sumNS(m_densityArray);
  if(m_totalDensity<=0)
  {
    printf("Invalid total density in computeTem, density= %le, p = %le, rho1 = %le, rho2= %le .\n",
     m_totalDensity,m_pressure,m_densityArray[0],m_densityArray[1]);
    printf("T = %lf, e = %lf, E = %lf.\n",m_temperature,m_energie,m_totalEnergy);
    exit(EXIT_FAILURE);
  }
  temp = m_pressure/(m_totalDensity*this->computeRMix());
  if(std::isnan(temp))
  {
    printf("In MixFire::computeTem():\n");
    printf("T=%lf, p = %lf, totalRho = %lf, RMix = %lf.\n",temp, m_pressure,m_totalDensity,this->computeRMix());
    exit(EXIT_FAILURE);
  }
  return temp;
}

//***************************************************************************

// compute internal energy per mass
double MixFire::computeInternalEnergy()
{
  double intEnergy = 0;
  if(m_totalDensity<=0)
  {
    printf("Invalid density at MixFire::computeInternalEnergy(), density = %lf.",m_totalDensity);
    exit(EXIT_FAILURE);
  }
  double T = this->computeTem();
  double p = this->getPressure();

  // call the method in Cantera
  std::shared_ptr<Cantera::ThermoPhase> gas = solution_mech->thermo();
  double Y[NS];
  for(int i=0;i<NS;i++) Y[i] = m_densityArray[i]/m_totalDensity;
  gas->setState_TPY(T, p, Y);
  intEnergy = gas->intEnergy_mass();

  return intEnergy;
}

//***************************************************************************

// compute sound speed
double MixFire::computeFrozenSoundSpeed()
{
  double gamma = this->computeRMix()/this->computeCv()+1;
  return sqrt(gamma*m_pressure/m_totalDensity);
}

//***************************************************************************

// compute Cv per mass
double MixFire::computeCv()
{
  double Cv = 0;
  double T = this->computeTem();
  double p = this->getPressure();

  // call the method in Cantera
  std::shared_ptr<Cantera::ThermoPhase> gas = solution_mech->thermo();
  double Y[NS];
  for(int i=0;i<NS;i++) Y[i] = m_densityArray[i]/m_totalDensity;
  gas->setState_TPY(T, p, Y);
  Cv = gas->cv_mass();
  return Cv;
}

//***************************************************************************

// compute total energy per mass
double MixFire::computeTotEne()
{
  return (this->computeInternalEnergy()+0.5*m_velocity.squaredNorm());
}

//****************************************************************************
//**************************** DATA PRINTING *********************************
//****************************************************************************

double MixFire::returnScalar(const int &numVar) const
{
  if(numVar == 1) return m_totalDensity;
  else if(numVar == 2) return m_pressure;
  else if(numVar == 3) return m_temperature;
  else if((numVar>3)&&(numVar<=NS+3)) return m_densityArray[numVar-4];
  else if(numVar == (NS+4)) return m_pMax;
  else return 1e30;
}

//***************************************************************************

Coord MixFire::returnVector(const int &numVar) const
{
  switch (numVar)
  {
  case 1:
    return m_velocity; break;
  default:
    return 0; break;
  }
}

//***************************************************************************

string MixFire::returnNameScalar(const int &numVar) const
{
  if(numVar == 1) return "rhoMix";
  else if(numVar == 2) return "Pressure_Mixture";
  else if(numVar == 3) return "Temperature_Mixture";
  else if((numVar>3)&&(numVar<=NS+3))
  {
    string speciesName;
    string rho = "rho_";
    speciesName = Sptr[numVar-4].Name;
    return rho+speciesName;
  } 
  else if(numVar == (NS+4)) return "pMax";
  else return "NoName";
}

//***************************************************************************

string MixFire::returnNameVector(const int &numVar) const
{
  switch (numVar)
  {
  case 1:
    return "Velocity_Mixture"; break;
  default:
    return "NoName"; break;
  }
}

//****************************************************************************
//**************************** DATA READING **********************************
//****************************************************************************

void MixFire::setScalar(const int &numVar, const double &value)
{
  if(numVar==1) m_totalDensity = value;
  else if(numVar==2) m_pressure = value; 
  else if(numVar==3) m_temperature = value;
  else if(numVar>3 && numVar<=NS+3) m_densityArray[numVar-4] = value;
  else if(numVar == (NS+4)) m_pMax = value;
  else Errors::errorMessage("numVar not found in MixFire::setScalar");
}

//***************************************************************************

void MixFire::setVector(const int &numVar, const Coord &value)
{
  switch (numVar)
  {
  case 1:
    m_velocity = value; break;
  default:
    Errors::errorMessage("numVar not found in MixKapila::setVector"); break;
  }
}

//****************************************************************************
//****************************** PARALLEL ************************************
//****************************************************************************

int MixFire::numberOfTransmittedVariables() const
{
  return NS + 9;
}

//***************************************************************************

void MixFire::fillBuffer(double *buffer, int &counter) const
{
  buffer[++counter] = m_velocity.getX();
  buffer[++counter] = m_velocity.getY();
  buffer[++counter] = m_velocity.getZ();
  buffer[++counter] = m_pressure;
  for(int ll=0;ll<NS;ll++)buffer[++counter] = m_densityArray[ll];
  buffer[++counter] = m_totalDensity;
  buffer[++counter] = m_temperature;
  buffer[++counter] = m_energie;
  buffer[++counter] = m_totalEnergy;
  buffer[++counter] = m_pMax;
}

//***************************************************************************

void MixFire::fillBuffer(std::vector<double> &dataToSend) const
{
  dataToSend.push_back(m_velocity.getX());
  dataToSend.push_back(m_velocity.getY());
  dataToSend.push_back(m_velocity.getZ());
  dataToSend.push_back(m_pressure);
  for(int ll=0;ll<NS;ll++)dataToSend.push_back(m_densityArray[ll]);
  dataToSend.push_back(m_totalDensity);
  dataToSend.push_back(m_temperature);
  dataToSend.push_back(m_energie);
  dataToSend.push_back(m_totalEnergy);
  dataToSend.push_back(m_pMax);
}

//***************************************************************************

void MixFire::getBuffer(double *buffer, int &counter)
{
  m_velocity.setX(buffer[++counter]);
  m_velocity.setY(buffer[++counter]);
  m_velocity.setZ(buffer[++counter]);
  m_pressure = buffer[++counter];
  for(int ll=0;ll<NS;ll++)m_densityArray[ll] = buffer[++counter];
  m_totalDensity = buffer[++counter];
  m_temperature = buffer[++counter];
  m_energie = buffer[++counter];
  m_totalEnergy = buffer[++counter];
  m_pMax = buffer[++counter];
}

//***************************************************************************

void MixFire::getBuffer(std::vector<double> &dataToReceive, int &counter)
{
  m_velocity.setX(dataToReceive[counter++]);
  m_velocity.setY(dataToReceive[counter++]);
  m_velocity.setZ(dataToReceive[counter++]);
  m_pressure = dataToReceive[counter++];
  for(int ll=0;ll<NS;ll++)m_densityArray[ll] = dataToReceive[counter++];
  m_totalDensity = dataToReceive[counter++];
  m_temperature = dataToReceive[counter++];
  m_energie = dataToReceive[counter++];
  m_totalEnergy = dataToReceive[counter++];
  m_pMax = dataToReceive[counter++];
}

//****************************************************************************
//******************************* ORDER 2 ************************************
//****************************************************************************

void MixFire::computeSlopesMixture(const Mixture &sLeft, const Mixture &sRight, const double &distance)
{
  m_pressure = (sRight.getPressure() - sLeft.getPressure()) / distance;
  m_velocity.setX((sRight.getVelocity().getX() - sLeft.getVelocity().getX()) / distance);
  m_velocity.setY((sRight.getVelocity().getY() - sLeft.getVelocity().getY()) / distance);
  m_velocity.setZ((sRight.getVelocity().getZ() - sLeft.getVelocity().getZ()) / distance);
  for(int ll=0;ll<NS;ll++)
    m_densityArray[ll] = (sRight.getDensityArray()[ll] - sLeft.getDensityArray()[ll]) / distance;
}

//***************************************************************************

void MixFire::setToZero()
{
  m_pressure = 0.;
  m_velocity.setX(0.); m_velocity.setY(0.); m_velocity.setZ(0.);
  for(int ll=0;ll<NS;ll++) m_densityArray[ll]= 0.;
  //m_totalDensity = 0.;
}

//***************************************************************************

void MixFire::extrapolate(const Mixture &slope, const double &distance)
{
  m_pressure += slope.getPressure()*distance;
  m_velocity.setX(m_velocity.getX() + slope.getVelocity().getX() * distance);
  m_velocity.setY(m_velocity.getY() + slope.getVelocity().getY() * distance);
  m_velocity.setZ(m_velocity.getZ() + slope.getVelocity().getZ() * distance);
  for(int ll=0;ll<NS;ll++)
    m_densityArray[ll] += slope.getDensityArray()[ll]*distance;
}

//***************************************************************************

void MixFire::limitSlopes(const Mixture &slopeGauche, const Mixture &slopeDroite, Limiter &globalLimiter)
{
  m_pressure = globalLimiter.limiteSlope(slopeGauche.getPressure(), slopeDroite.getPressure());
  m_velocity.setX(globalLimiter.limiteSlope(slopeGauche.getVelocity().getX(), slopeDroite.getVelocity().getX()));
  m_velocity.setY(globalLimiter.limiteSlope(slopeGauche.getVelocity().getY(), slopeDroite.getVelocity().getY()));
  m_velocity.setZ(globalLimiter.limiteSlope(slopeGauche.getVelocity().getZ(), slopeDroite.getVelocity().getZ()));
  for(int ll=0;ll<NS;ll++)
    m_densityArray[ll] = globalLimiter.limiteSlope(slopeGauche.getDensityArray()[ll], slopeDroite.getDensityArray()[ll]);
}

//****************************************************************************
//************************** ORDER 2 PARALLEL ********************************
//****************************************************************************

int MixFire::numberOfTransmittedSlopes() const
{
	return NS+4;
}

//***************************************************************************

void MixFire::fillBufferSlopes(double *buffer, int &counter) const
{
	buffer[++counter] = m_velocity.getX();
	buffer[++counter] = m_velocity.getY();
	buffer[++counter] = m_velocity.getZ();
  buffer[++counter] = m_pressure;
  for(int ll=0;ll<NS;ll++)buffer[++counter] = m_densityArray[ll];
}

//***************************************************************************

void MixFire::getBufferSlopes(double *buffer, int &counter)
{
	m_velocity.setX(buffer[++counter]);
	m_velocity.setY(buffer[++counter]);
	m_velocity.setZ(buffer[++counter]);
  m_pressure = buffer[++counter];
  for(int ll=0;ll<NS;ll++)m_densityArray[ll] = buffer[++counter];
}


//****************************************************************************
//******************************* ACCESSORS **********************************
//****************************************************************************

std::array<double,NS> const & MixFire:: getDensityArray() const
{
  return m_densityArray;
}

std::array<double,NS> const & MixFire:: getMolarFractionArray() const
{
  return m_molarFractionArray;
}

//***************************************************************************


//***************************************************************************

const double& MixFire::getTemperature() const
{
  return m_temperature;
}

//***************************************************************************

void MixFire::setDensityArray(const double densityArray[NS])
{
  double rhoMix = 0;
  for(int ll=0;ll<NS;ll++)
  {
    m_densityArray[ll]=densityArray[ll];
    rhoMix += densityArray[ll];
  }
  if(rhoMix<=0.)
  {
    printf("Invalid rhoMix in MixFire::setDensityArray.\n");
    printf("rho0 = %le, rho1 = %le.\n",densityArray[0],densityArray[1]);
    exit(EXIT_FAILURE);
  }
}

//***************************************************************************

void MixFire::setDensity(const double &density){m_totalDensity = density;}

//***************************************************************************

void MixFire::setPressure(const double &p) { m_pressure = p; }

//***************************************************************************

void MixFire::setVelocity(const double &u, const double &v, const double &w) 
{
  m_velocity.setXYZ(u, v, w); 
}

//***************************************************************************

void MixFire::setVelocity(const Coord &vit) 
{
  m_velocity = vit;
}

//***************************************************************************

void MixFire::setU(const double &u) { m_velocity.setX(u); }

//***************************************************************************

void MixFire::setV(const double &v) { m_velocity.setY(v); }

//***************************************************************************

void MixFire::setW(const double &w) { m_velocity.setZ(w); }

//***************************************************************************

void MixFire::setTotalEnergy(double &totalEnergy)
{
  m_totalEnergy = totalEnergy;
}

//****************************************************************************
//***************************** OPERATORS ************************************
//****************************************************************************

void MixFire::changeSign()
{
  m_totalDensity = -m_totalDensity;
  for(int i=0;i<NS;i++) m_densityArray[i]=-m_densityArray[i];
  m_pressure = -m_pressure;
  m_velocity = m_velocity*-1.;
}

//***************************************************************************

void MixFire::multiplyAndAdd(const Mixture &slopesMixtureTemp, const double &coeff)
{
  m_totalDensity += slopesMixtureTemp.getDensity()*coeff;
  for(int i=0;i<NS;i++) 
    m_densityArray[i] += slopesMixtureTemp.getDensityArray()[i]*coeff;
  m_pressure += slopesMixtureTemp.getPressure()*coeff;
  m_velocity += slopesMixtureTemp.getVelocity()*coeff;
}

//***************************************************************************

void MixFire::divide(const double &coeff)
{
  m_totalDensity /= coeff;
  for(int i=0;i<NS;i++) m_densityArray[i] /= coeff;
  m_pressure /= coeff;
  m_velocity /= coeff;
}

//***************************************************************************

void MixFire::setTemperature(const double &T)
{
  m_temperature = T;
}

//***************************************************************************

void MixFire::setDensityI(const double densityI, int i)
{
  m_densityArray[i] = densityI;
}


void MixFire::setInductionTime(double inductionTime_)
{
  m_tau_ind = inductionTime_;
}

double MixFire::getInductionTime()
{
  return m_tau_ind;
}

double MixFire::getMach()
{
  double a = this->computeFrozenSoundSpeed();
  return m_velocity.norm()/a;
}

//***************************************************************************

void MixFire::set_deltaEnthalpy(const double &delta_h_mass)
{
  m_delta_h_mass = delta_h_mass;
}

//***************************************************************************

void MixFire::storeOldStep()
{
  printf("Banish MixFire::storeOldStep().\n");
  exit(EXIT_FAILURE);
  for(int ll=0;ll<NS;ll++) m_densityArray_old[ll]=m_densityArray[ll];
  m_velocity_old = m_velocity;
  m_pressure_old = m_pressure;
  m_temperature_old = m_temperature;
}

//***************************************************************************

void MixFire::chemicalChange()
{
  printf("Banish MixFire::chemicalChange().\n");
  exit(EXIT_FAILURE);
}

//***************************************************************************

// fane
std::array<double,NS> const & MixFire:: getDensityChangeArray() const
{
  return m_densityArray_change;
}

//***************************************************************************

void MixFire::setHfChangeRate(const double &dt_lvl_0)
{
  m_hf_change_rate = m_hf_change/dt_lvl_0;
}

//***************************************************************************

void MixFire::setDensityArrayChangeRate(const double &dt_lvl_0)
{
  for(int ll=0;ll<NS;ll++) m_densityArray_change_rate[ll]=m_densityArray_change[ll]/dt_lvl_0;
}

//***************************************************************************

std::array<double,NS> const & MixFire::getDensityChangeRateArray() const
{
  return m_densityArray_change_rate;
}

//***************************************************************************

void MixFire::addToHfChange(double delta_enthalpy)
{
  m_hf_change += delta_enthalpy;
}

//***************************************************************************

void MixFire::setDensityArrayChangeToZero()
{
  for(int ll=0;ll<NS;ll++) m_densityArray_change[ll]=0;
}

//***************************************************************************

void MixFire::addToDensityArrayChange(double delta_rho_i, int i)
{
  m_densityArray_change[i] += delta_rho_i;
}

//***************************************************************************

void MixFire::calcHRR()
{
  double Ys[NS], Ys_old[NS];
  double dt_hrr(1e-7);

  double hrr_t(0), deltaRho[NS];

  double rho = m_totalDensity;
  std::array<double, NS> rhoArray = m_densityArray;
  double T = m_temperature;
  double P = m_pressure;

  for(int k=0;k<NS;k++) Ys[k] = rhoArray[k]/rho;
  for(int k=0;k<NS;k++) Ys_old[k] = Ys[k];

  std::shared_ptr<Cantera::Solution> sol_ = solution_mech;
  std::shared_ptr<Cantera::ThermoPhase> gas = sol_->thermo();
  Cantera::IdealGasReactor reactor;
	Cantera::ReactorNet net;

  reactor.insert(sol_);
  net.addReactor(reactor);

  gas->setState_TPY(T,P,Ys);
  reactor.syncState();
	net.setInitialTime(0.0);
	net.advance(dt_hrr);

  for(int k=0;k<NS;k++) 
  {
    Ys[k] = gas->massFraction(k);
  }

  std::shared_ptr<Cantera::ThermoPhase> gas2 = solution_mech->thermo();
  double Y[NS] = {0};
  for(int i=0;i<NS;i++)
  {
    Y[i] = 1.0;
    gas2->setState_TPY(298.15, 101325, Y); // formation enthalpy of each species
    hrr_t += (Ys_old[i]-Ys[i])*rho*gas2->enthalpy_mass();
    Y[i] = 0.0;
  }

  m_hrr = hrr_t/dt_hrr;
}