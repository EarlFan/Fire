//  
//       ,---.     ,--,    .---.     ,--,    ,---.    .-. .-. 
//       | .-'   .' .')   / .-. )  .' .'     | .-'    |  \| | 
//       | `-.   |  |(_)  | | |(_) |  |  __  | `-.    |   | | 
//       | .-'   \  \     | | | |  \  \ ( _) | .-'    | |\  | 
//       |  `--.  \  `-.  \ `-' /   \  `-) ) |  `--.  | | |)| 
//       /( __.'   \____\  )---'    )\____/  /( __.'  /(  (_) 
//      (__)              (_)      (__)     (__)     (__)     
//
//  This file is part of ECOGEN.
//
//  ECOGEN is the legal property of its developers, whose names 
//  are listed in the copyright file included with this source 
//  distribution.
//
//  ECOGEN is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published 
//  by the Free Software Foundation, either version 3 of the License, 
//  or (at your option) any later version.
//  
//  ECOGEN is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
//  GNU General Public License for more details.
//  
//  You should have received a copy of the GNU General Public License
//  along with ECOGEN (file LICENSE).  
//  If not, see <http://www.gnu.org/licenses/>.

#ifndef EOS_H
#define EOS_H

//! \file      Eos.h
//! \author    F. Petitpas, K. Schmidmayer, E. Daniel, J. Caze
//! \version   1.1
//! \date      October 11 2019

#include <string>
#include <vector>
#include <cassert>
#include "../Errors.h"
#include "../libTierces/tinyxml2.h"

//! \class     Eos
//! \brief     General class for Equation of State (EOS).
//! \details   This is a pure virtual class: can not be instantiated.
//! \details This base class  and its derived classes contain mainly computational functions to obtain thermodynamical values for each phase: pressure, temperature, internal energy...
//! \details   In this ECOGEN release only three EOS are available : ideal Gas (sub class EosIG), Stiffened (sub class EoSG) and Noble-Abel Stiffened Gas (sub class EosNASG). Appropriate data for both EOS lead to a large choice of gas, liquid or compressible solid.
class Eos
{
    public:
      Eos();
      Eos(int &number);
      virtual ~Eos();


      //Constant methods
      //! \brief display the name of the fluid
      void display() const;

      const std::string& getName() const { return m_name; };
      //! \brief See derived classes 
      virtual std::string getType() const { return "NA"; };
      //! \brief  Return the number associated to the EOS
      //! \return  m_number
      const int& getNumber() const { return m_number; };

      //! \brief Read physical parameters (viscosity, thermal conductivity....)
      void readPhysicalParameter(tinyxml2::XMLNode *element, std::string fileName = "Unknown File");

      //! \brief    compute the total enthalpy of the phase
      //! \param     density           density (\f$\rho \f$)
      //! \param     pressure          pressure (p)
      //! \param     velocity          velocity (m/s)
      //! \return    this -> hTotal
      //! \details  hTotal is computed as :\f$H_{total} =  \epsilon (\rho, p)   + \frac{p}{\rho} +\frac{1}{2}u^2 \f$. 
      double computeTotalEnthalpy(const double &density, const double &pressure, const double &velocity) const;
      //! \brief  Return the dynamic viscosity of the fluid  
      //!return    \f$ \mu \f$ (Unit: Pa.s).
      const double& getMu() const { return m_mu; };
      //! \brief  get the thermal conductivity of the fluid
      //!return    \f$ \lambda \f$ (Unit: W/(m.K)).
      const double& getLambda() const { return m_lambda; };
      //! \brief  Assign the epsilonAlphaNull value for alphaNull option (alpha = 0 => epsilonAlphaNull != 0 ; alpha != 0 => epsilonAlphaNull = 0)
      //! \param     alphaNull         \f$\alpha = 0 \f$ (boolean value)
      //! \param     fileName          Name of the open file
      void assignEpsilonForAlphaNull(bool alphaNull, std::string fileName) const;


      //Gereral methods for all EOS


      //Virtual methods for child classes
      //! \brief See derived classes  
      virtual void assignParametersEos(std::string name, std::vector<double> parametersEos) = 0;
      //! \brief See derived classes 
      virtual double computeTemperature(const double &density, const double &pressure) const=0; //virtual methods 
      //! \brief See derived classes 
      virtual double computeEnergy(const double &density, const double &pressure) const=0;
      //! \brief See derived classes 
      virtual double computePressure(const double &density, const double &energy) const=0;
      //! \brief See derived classes 
      virtual double computeDensity(const double &pressure, const double &temperature) const = 0;
      //! \brief See derived classes 
      virtual double computeSoundSpeed(const double &density, const double &pressure) const=0;
      //! \brief See derived classes 
      virtual double computeEntropy(const double &temperature, const double &pressure) const = 0;
      //! \brief See derived classes 
      virtual double computePressureIsentropic(const double &initialPressure, const double &initialDensity, const double &finalDensity) const=0;
      //! \brief See derived classes 
      virtual double computePressureHugoniot(const double &initialPressure, const double &initialDensity, const double &finalDensity) const=0;
      //! \brief See derived classes 
      virtual double computeDensityIsentropic(const double &initialPressure, const double &initialDensity, const double &finalPressure, double *drhodp=0) const=0;
      //! \brief See derived classes 
      virtual double computeDensityHugoniot(const double &initialPressure, const double &initialDensity, const double &finalPressure, double *drhodp = 0) const = 0;

      virtual double computeDensityPfinal(const double &initialPressure, const double &initialDensity, const double &finalPressure, double *drhodp = 0) const = 0;

      //! \brief See derived classes 
      virtual double computeEnthalpyIsentropic(const double &initialPressure, const double &initialDensity, const double &finalPressure, double *dhdp = 0) const = 0;
      //! \brief See derived classes 
      virtual double computeDensitySaturation(const double &pressure, const double &Tsat, const double &dTsatdP, double *drhodp = 0) const { Errors::errorMessage("computeDensitySaturation not yet programmed for EOS : " + m_name); return 0; };
      //! \brief See derived classes 
      virtual double computeDensityEnergySaturation(const double &pressure, const double &rho, const double &drhodp, double *drhoedp = 0) const { Errors::errorMessage("computeDensityEnergySaturation non prevu pour EOS : " + m_name); return 0; };
      //! \brief See derived classes 
      virtual void sendSpecialMixtureEos(double &gamPinfOverGamMinusOne, double &eRef, double &oneOverGamMinusOne, double &covolume) const = 0;
      //! \brief See derived classes 
      virtual double vfpfh(const double &pressure, const double &enthalpy) const { Errors::errorMessage("vfpfh not yet programmed for EOS : " + m_name); return 0; };

      //Partial derivatives 
      //! \brief See derived classes 
      virtual double dvdpch(const double &pressure, const double &enthalpy) const { Errors::errorMessage("dvdpch not yet programmed for EOS : " + m_name); return 0; };
      //! \brief See derived classes 
      virtual double dvdhcp(const double &pressure, const double &enthalpy) const { Errors::errorMessage("dvdhcp not yet programmed for EOS : " + m_name); return 0; };

      //Checking..
      virtual void verifyPressure(const double &pressure, const std::string &message = " ") const { Errors::errorMessage("verifyPressure not yet programmed Eos : " + m_name); };
      //! \brief See derived classes 
      virtual void verifyAndModifyPressure(double &pressure) const { Errors::errorMessage("verifyAndModifyPressure not yet programmed Eos : "+ m_name); };

      //Mod
      //! \brief See derived classes 
      virtual void sendInfo(double *&data) const = 0; //FP//DEV// Usefull for some exact solvers (not used)

      //Get
      //! \brief See derived classes 
      virtual const double& getGamma() const { return Errors::defaultDouble; };
      //! \brief See derived classes 
      virtual const double& getPInf() const { return Errors::defaultDouble; };
      //! \brief See derived classes 
      virtual const double& getCv() const { return Errors::defaultDouble; };
      //! \brief See derived classes 
      virtual const double& getERef() const { return Errors::defaultDouble; };
      //! \brief See derived classes 
      virtual const double& getSRef() const { return Errors::defaultDouble; };
	
      //Specific method for EosTPIG
      //! \brief  get the the molecular weight
      //!return    \f$ \molarWeight \f$ (Unit: g/mol).
      virtual double getMW() const { Errors::errorMessage("getMW() not yet programmed Eos : " + m_name); return 0;};

      virtual double computeEnergy(double Temperature)const{ Errors::errorMessage("computeEnergy not yet programmed Eos : " + m_name); };

      virtual void verifyTemperature(const double &temperature, const std::string &message = "") const { Errors::errorMessage("verifyTemperature not yet programmed Eos : " + m_name); };

      virtual int getIndex() const { Errors::errorMessage("getIndex not yet programmed Eos : " + m_name); };

    protected:
      int m_number;       //!< Corresponding number of the equation of state
      std::string m_name; //!< Name of the equation of state

      double m_mu;        //!< Dynamic viscosity (kg/m/s or Pa.s)
      double m_lambda;    //!< Thermal conductivity (W/(m.K))
};

extern double epsilonAlphaNull;    //!< Epsilon value to avoid division by 0 when alpha = 0 is activated. If alpha = 0 desactivated, i.e. alpha != 0, then epsilonAlphaNull = 0.

#endif // EOS_H
