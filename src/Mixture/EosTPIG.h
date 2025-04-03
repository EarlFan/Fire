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

#ifndef EOSTPIG_H
#define EOSTPIG_H

//! \file      EosTPIG.h
//! \author    E Fan, Lisong Shi, Jiaao Hao, Tianhan Zhang, Chih-yung Wen.
//! \version   1.0
//! \date      Jan 01 2025

#include "../libs/ECOGEN/Eos/Eos.h"
#include "../include/Fire.h"

//! \class     EosTPIG
//! \brief     Class describing an ideal gas equation of state
class EosTPIG : public Eos
{
    public:
        EosTPIG();
        EosTPIG(std::vector<std::string> &nameParameterEos, int& number);
        virtual ~EosTPIG();


        virtual void assignParametersEos(std::string name, std::vector<double> parametersEos) {};

        virtual double computeTemperature(const double &density,const double &pressure) const;

        virtual double computeEnergy(const double &density,const double &pressure) const;

        virtual double computePressure(const double &density,const double &energy) const;

        virtual double computeDensity(const double &pressure, const double &temperature) const;

        virtual double computeSoundSpeed(const double &density,const double &pressure) const;

        virtual double computeEntropy(const double &temperature, const double &pressure) const;

        virtual double computePressureIsentropic(const double &initialPressure, const double &initialDensity, const double &finalDensity) const; 

        virtual double computePressureHugoniot(const double &initialPressure, const double &initialDensity, const double &finalDensity) const;

        virtual double computeDensityIsentropic(const double &initialPressure, const double &initialDensity, const double &finalPressure, double *drhodp=0) const;

        virtual double computeDensityHugoniot(const double &initialPressure, const double &initialDensity, const double &finalPressure, double *drhodp = 0) const;

        virtual double computeDensityPfinal(const double &initialPressure, const double &initialDensity, const double &finalPressure, double *drhodp = 0) const;

        virtual double computeEnthalpyIsentropic(const double &initialPressure, const double &initialDensity, const double &finalPressure, double *dhdp = 0) const;

        virtual double computeDensitySaturation(const double &pressure, const double &Tsat, const double &dTsatdP, double *drhodp = 0) const;

        virtual double computeDensityEnergySaturation(const double &pressure, const double &rho, const double &drhodp, double *drhoedp = 0) const;

        virtual void sendSpecialMixtureEos(double &gamPinfOverGamMinusOne, double &eRef, double &oneOverGamMinusOne, double &covolume) const;

		virtual double vfpfh(const double &pressure, const double &enthalpy) const;

        virtual double dvdpch(const double &pressure, const double &enthalpy) const;

        virtual double dvdhcp(const double &pressure, const double &enthalpy) const;

        virtual void verifyPressure(const double &pressure, const std::string &message = "") const;

		virtual void verifyAndModifyPressure(double &pressure) const;

		   virtual void sendInfo(double *&data) const {};

        virtual const double& getGamma() const;

        virtual const double& getCv() const;

        virtual const double& getERef() const;

        virtual const double& getSRef() const;

        virtual std::string getType() const { return "TPIG"; };

		EosTPIG (int index, int &number);

	    virtual double getMW() const;

		virtual double computeEnergy(double Temperature)const;

        virtual void verifyTemperature(const double &temperature, const std::string &message) const;

        virtual int getIndex() const { return m_index; };



    protected:
    private:
	    int m_index;
        //std::string m_name;
        double m_MW;
		double m_Rs;
};

#endif // EosTPIG_H
