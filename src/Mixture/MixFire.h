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

#ifndef MIXFIRE_H
#define MIXFIRE_H

//! \file      MixFire.h
//! \author    E Fan, Lisong Shi, Jiaao Hao, Tianhan Zhang, Chih-yung Wen.
//! \version   1.0
//! \date      Jan 01 2025

#include <vector>
#include "../libs/ECOGEN/Models/Mixture.h"
#include "../include/Fire.h"

//! \class     MixFire
//! \brief     Mixture variables for Kapila system of equations (mechanical equilibrium)
class MixFire : public Mixture
{
    public:
      MixFire();
      //! \brief     Mixture constructor from a XML format reading
      //! \details   Reading data from XML file under the following format:
      //!           ex: <mixture>
      //!                 <velocity x = "0." y = "0." z = "0." />
      //!               </mixture>
      //! \param     state           XML element to read for mixture data
      //! \param     fileName       string name of readed XML file
      MixFire(tinyxml2::XMLElement *state, std::string fileName);
      virtual ~MixFire();

      virtual void allocateAndCopyMixture(Mixture **mixture);
      virtual void copyMixture(Mixture &mixture);
      virtual double computeDensity(const double *alphak, const double *rhok, const int &numberPhases);
      virtual double computePressure(const double *alphak, const double *pk, const int &numberPhases);
      virtual double computeInternalEnergy(const double *Yk, const double *ek, const int &numberPhases);
      virtual double computeFrozenSoundSpeed(const double *Yk, const double *ck, const int &numberPhases);
      
      virtual void computeMixtureVariables(Phase **vecPhase, const int &numberPhases);
      virtual void internalEnergyToTotalEnergy(std::vector<QuantitiesAddPhys*> &vecGPA);
      virtual void totalEnergyToInternalEnergy(std::vector<QuantitiesAddPhys*> &vecGPA);

      virtual void localProjection(const Coord &normal, const Coord &tangent, const Coord &binormal);
      virtual void reverseProjection(const Coord &normal, const Coord &tangent, const Coord &binormal);

      //Data printing
      virtual int getNumberScalars() const { return NS+4; };   //For complete output
      //virtual int getNumberScalars() const { return 3; }; //For complete output
      virtual int getNumberVectors() const { return 1; };
      virtual double returnScalar(const int &numVar) const;
      virtual Coord returnVector(const int &numVar) const;
      virtual std::string returnNameScalar(const int &numVar) const;
      virtual std::string returnNameVector(const int &numVar) const;

      //Data reading
      virtual void setScalar(const int &numVar, const double &value);
      virtual void setVector(const int &numVar, const Coord &value);

      //Parallel
      virtual int numberOfTransmittedVariables() const;
      virtual void fillBuffer(double *buffer, int &counter) const;
      virtual void fillBuffer(std::vector<double> &dataToSend) const;
      virtual void getBuffer(double *buffer, int &counter);
      virtual void getBuffer(std::vector<double> &dataToReceive, int &counter);

      //Second order
      virtual void computeSlopesMixture(const Mixture &sLeft, const Mixture &sRight, const double &distance);
      virtual void setToZero();
      virtual void extrapolate(const Mixture &slope, const double &distance);
      virtual void limitSlopes(const Mixture &slopeGauche, const Mixture &slopeDroite, Limiter &globalLimiter);

      //Parallel second order
      virtual int numberOfTransmittedSlopes() const;
      virtual void fillBufferSlopes(double *buffer, int &counter) const;
      virtual void getBufferSlopes(double *buffer, int &counter);

      //Accessors
      virtual const double& getDensity() const { return m_totalDensity; };
      virtual const double& getPressure() const { return m_pressure; };
      virtual const double& getU() const { return m_velocity.getX(); };
      virtual const double& getV() const { return m_velocity.getY(); };
      virtual const double& getW() const { return m_velocity.getZ(); };
      virtual const Coord& getVelocity() const { return m_velocity; };
      virtual Coord& getVelocity() { return m_velocity; };
      virtual const double& getEnergy() const { return m_energie; };
      virtual const double& getTotalEnergy() const { return m_totalEnergy; };
      virtual const double& getMixSoundSpeed() const {return m_frozenSoundSpeed;};
      virtual const double& getFrozenSoundSpeed() const { return m_frozenSoundSpeed; };
      // virtual const double& getWoodSoundSpeed() const { return m_woodSoundSpeed; };

      virtual void setPressure(const double &p);
      virtual void setVelocity(const double &u, const double &v, const double &w);
      virtual void setVelocity(const Coord &vit);
      virtual void setU(const double &u);
      virtual void setV(const double &v);
      virtual void setW(const double &w);
      virtual void setTotalEnergy(double &totalEnergy);

      //Operators
      virtual void changeSign();
      virtual void multiplyAndAdd(const Mixture &slopesMixtureTemp, const double &coeff);
      virtual void divide(const double &coeff);

      //Specific method for Fire
      virtual std::array<double,NS> const & getDensityArray() const;
      
      virtual void updateMolarFractionArray();

      virtual std::array<double,NS> const & getMolarFractionArray() const;

      virtual void setDensityArray(const double densityArray[NS]);

      virtual void setDensity(const double &density);

      virtual double computeTotalDensity();

      virtual double computeRMix();

      virtual double computePressure();

      virtual double computeTem();

      virtual double computeInternalEnergy();

      virtual double computeCv() ;

      virtual double computeTotEne();

      virtual double computeFrozenSoundSpeed() ;

      virtual void pri2con(){};

      virtual void con2pri();

      virtual const double& getTemperature() const ;

      virtual void setTemperature(const double &T) ;

      virtual void setDensityI(const double densityI, int i) ;

      virtual void setInductionTime(double inductionTime_);

      virtual double getInductionTime();

      virtual double getMach();

      virtual void set_deltaEnthalpy(const double &delta_h_mass);

      virtual void storeOldStep();

      virtual void chemicalChange();

      virtual std::array<double,NS> const & getDensityChangeArray() const;

      virtual const double& getHfChange() const { return m_hf_change; };

      virtual void setHfChangeRate(const double &dt_lvl_0);

      virtual const double& getHfChangeRate()const { return m_hf_change_rate; };

      virtual void setDensityArrayChangeRate(const double &dt_lvl_0);

      virtual std::array<double,NS> const & getDensityChangeRateArray() const;

      virtual void setHfChangeToZero(){m_hf_change=0;};

      virtual void addToHfChange(double delta_enthalpy);

      virtual void setDensityArrayChangeToZero();

      virtual void addToDensityArrayChange(double delta_rho_i, int i);

      virtual void setRhoExcitedOH(double rhoExcitedOH){ m_rhoExitedOH =rhoExcitedOH; };

      virtual double getRhoExcitedOH(){ return m_rhoExitedOH;};

      virtual const double& getMaxPressure() const { return m_pMax; };

      virtual void setPMax(double p_larger) { m_pMax = p_larger; };

      virtual void calcHRR();

      virtual double getHRR() { return m_hrr;}

    protected:
    private:
      double m_totalDensity;         //!< total density
      std::array<double, NS> m_densityArray;
      std::array<double, NS> m_molarFractionArray; //!< molar fraction array
      Coord m_velocity;              //!< mixture velocity
      double m_pressure;             //!< mixture pressure
      double m_temperature;           //!< temperature
      double m_energie;              //!< mixture internal specific energy
      double m_totalEnergy;          //!< mixture total specific energy
      double m_frozenSoundSpeed;     //!< frozen sound speed
      double m_tau_ind;              //!< induction time
      double m_delta_h_mass;         //!< delta enthalpy per mass due to chemical reaction

      // storage the old flowfield for later analysis
      std::array<double, NS> m_densityArray_old;
      Coord m_velocity_old;
      double m_pressure_old, m_temperature_old;             //!< mixture pressure
      std::array<double, NS> m_densityArray_change;
      std::array<double, NS> m_densityArray_change_rate;
      double m_hf_change;  //  J/(m^3*s)
      double m_hf_change_rate;  //  J/(m^3*s)
      double m_rhoExitedOH; //
      double m_pMax;
      double m_hrr;  //  J/(m^3*s)
};

#endif // MIXFIRE_H
