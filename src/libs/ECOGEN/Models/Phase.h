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

#ifndef PHASE_H
#define PHASE_H

//! \file      Phase.h
//! \author    F. Petitpas, K. Schmidmayer
//! \version   1.1
//! \date      June 5 2019

#include <fstream>
#include "../Errors.h"
#include "../Eos/Eos.h"
#include "../Maths/Coord.h"
#include "../libTierces/tinyxml2.h"
#include "../Order2/HeaderLimiter.h"

enum Prim { vecPhases, vecPhasesO2, vecSlopes, restart };

//! \class     Phase
//! \brief     Abstract class for a phase
//! \details   Can not be instanciated, variables depend on the model
class Phase
{
  public:
    Phase();
    virtual ~Phase();
    //! \brief     Print phase variables in file stream
    //! \param     fileStream      file stream to write in
    void printPhase(std::ofstream &fileStream) const;
    //! \brief     Copy phase attributes in phase
    //! \param     vecPhase      destination phase variable 
    virtual void allocateAndCopyPhase(Phase **vecPhase) { Errors::errorMessage("allocateAndCopyPhase not available for requested phase type"); };
    //! \brief     Copy phase in phase attributes
    //! \param     vecPhase      source phase to copy
    virtual void copyPhase(Phase &phase) { Errors::errorMessage("copyPhase not available for requested phase type"); };
    //! \brief     Compute extra thermodynammical variables
    //! \details   Computes from velocity, pressure and density
    //! \param     velocity      phase velocity
    virtual void extendedCalculusPhase(const Coord &velocity) { Errors::errorMessage("extendedCalculusPhase not available for requested phase type"); };
    
    virtual void computeMassFraction(const double &density) { Errors::errorMessage("computeMassFraction not available for requested phase type"); };

    virtual void localProjection(const Coord &normal, const Coord &tangent, const Coord &binormal) { Errors::errorMessage("projection not available for requested phase type"); };
    virtual void reverseProjection(const Coord &normal, const Coord &tangent, const Coord &binormal) { Errors::errorMessage("reverseProjection not available for requested phase type"); };

    //Specific methods for data printing
    //----------------------------------
    virtual int getNumberScalars() const { Errors::errorMessage("getNumberScalars not available for requested phase type"); return 0; };
    virtual int getNumberVectors() const { Errors::errorMessage("getNumberVectors not available for requested phase type"); return 0; };
    virtual double returnScalar(const int &numVar) const { Errors::errorMessage("returnScalar not available for requested phase type"); return 0.; };
    virtual Coord returnVector(const int &numVar) const { Errors::errorMessage("returnVector not available for requested phase type"); return 0; };
    virtual std::string returnNameScalar(const int &numVar) const { Errors::errorMessage("returnNameScalar not available for requested phase type"); return 0; };
    virtual std::string returnNameVector(const int &numVar) const { Errors::errorMessage("returnNameVector not available for requested phase type"); return 0; };

    //Specific method for reading from file
    //-------------------------------------
    virtual void setScalar(const int &numVar, const double &value) { Errors::errorMessage("setScalar not available for requested phase type"); };
    virtual void setVector(const int &numVar, const Coord &value) { Errors::errorMessage("setVector not available for requested phase type"); };

    //Specific methods for parallel computing
    //---------------------------------------
    virtual int numberOfTransmittedVariables() const { Errors::errorMessage("numberOfTransmittedVariables not available for requested phase type"); return 0; };
    virtual void fillBuffer(double *buffer, int &counter) const { Errors::errorMessage("fillBuffer not available for requested phase type"); };
    virtual void fillBuffer(std::vector<double> &dataToSend) const { Errors::errorMessage("fillBuffer not available for requested phase type"); };
    virtual void getBuffer(double *buffer, int &counter, Eos **eos) { Errors::errorMessage("getBuffer not available for requested phase type"); };
    virtual void getBuffer(std::vector<double> &dataToReceive, int &counter, Eos **eos) { Errors::errorMessage("getBuffer not available for requested phase type"); };

    //Specific methods for second order
    //---------------------------------
    virtual void computeSlopesPhase(const Phase &sLeft, const Phase &sRight, const double &distance) { Errors::errorMessage("computeSlopesPhase not available for requested phase type"); };
    virtual void setToZero() { Errors::errorMessage("setToZero not available for requested phase type"); };
    virtual void extrapolate(const Phase &slope, const double &distance) { Errors::errorMessage("extrapolate not available for requested phase type"); };
    virtual void limitSlopes(const Phase &slopeGauche, const Phase &slopeDroite, Limiter &globalLimiter, Limiter &volumeFractionLimiter) { Errors::errorMessage("limitSlopes not available for requested phase type"); };

    //Specific methods for parallele computing at second order
    //--------------------------------------------------------
    virtual int numberOfTransmittedSlopes() const { Errors::errorMessage("numberOfTransmittedSlopes not available for requested phase type"); return 0; };
    virtual void fillBufferSlopes(double *buffer, int &counter) const { Errors::errorMessage("fillBufferSlopes not available for requested phase type"); };
    virtual void getBufferSlopes(double *buffer, int &counter) { Errors::errorMessage("getBufferSlopes not available for requested phase type"); };

    //Verifications
    //-------------
    virtual void verifyPhase(const std::string &message = "") const { Errors::errorMessage("verifyPhase non prevu pour type de phase"); };
    virtual void verifyAndCorrectPhase() { Errors::errorMessage("verifyAndCorrectPhase non prevu pour type de phase"); };

    //Accessors
    //---------
    virtual const double& getAlpha() const { Errors::errorMessage("getAlpha not available for requested phase type"); return Errors::defaultDouble; };
    virtual const double& getDensity() const { Errors::errorMessage("getDensity not available for requested phase type"); return Errors::defaultDouble; };
    virtual const double& getPressure() const { Errors::errorMessage("getPressure not available for requested phase type"); return Errors::defaultDouble; };
    virtual const double& getY() const { Errors::errorMessage("getY not available for requested phase type"); return Errors::defaultDouble; };
    virtual const double& getU() const { Errors::errorMessage("getU not available for requested phase type"); return Errors::defaultDouble; };
    virtual const double& getV() const { Errors::errorMessage("getV not available for requested phase type"); return Errors::defaultDouble; };
    virtual const double& getW() const { Errors::errorMessage("getW not available for requested phase type"); return Errors::defaultDouble; };
    virtual Coord& getVelocity() { Errors::errorMessage("getVelocity not available for requested phase type"); return Coord::defaultCoordNonConst; };
    virtual const Coord& getVelocity() const { Errors::errorMessage("getVelocity not available for requested phase type"); return Coord::defaultCoord; };
    virtual Eos* getEos() const { Errors::errorMessage("EOS not available for requested phase type"); return 0; };
    virtual const double& getEnergy() const { Errors::errorMessage("getEnergy impossible avec type de phase demande"); return Errors::defaultDouble; };
    virtual const double& getSoundSpeed() const { Errors::errorMessage("getSoundSpeed impossible avec type de phase demande"); return Errors::defaultDouble; };
    virtual const double& getTotalEnergy() const { Errors::errorMessage("getTotalEnergy impossible avec type de phase demande"); return Errors::defaultDouble; };
    virtual double getTemperature() const { Errors::errorMessage("getT impossible avec type de phase demande"); return 0.; };

    virtual void setAlpha(double alpha) { Errors::errorMessage("setAlpha not available for requested phase type"); };
    virtual void setDensity(double density) { Errors::errorMessage("setDensity not available for requested phase type"); };
    virtual void setPressure(double pressure) { Errors::errorMessage("setPressure not available for requested phase type"); };
    virtual void setVelocity(const double &u, const double &v, const double &w) { Errors::errorMessage("setVelocity not available for requested phase type"); };
    virtual void setVelocity(const Coord &vit) { Errors::errorMessage("setVelocity not available for requested phase type"); };
    virtual void setU(const double &u) { Errors::errorMessage("setU not available for requested phase type"); };
    virtual void setV(const double &v) { Errors::errorMessage("setV not available for requested phase type"); };
    virtual void setW(const double &w) { Errors::errorMessage("setW not available for requested phase type"); };
    virtual void setEos(Eos *eos) { Errors::errorMessage("impossible to associate EOS to the requested phase type"); };
    virtual void setEnergy(double energie) { Errors::errorMessage("setEnergy not available for requested phase type"); };
    virtual void setSoundSpeed(double soundSpeed) { Errors::errorMessage("setSoundSpeed not available for requested phase type"); };
    virtual void setTotalEnergy(double totalEnergy) { Errors::errorMessage("setTotalEnergy not available for requested phase type"); };
    virtual void setTemperature(double temperature) { Errors::errorMessage("setTemperature not available for requested phase type"); };

    //Operators
    //---------
    virtual void changeSign() { Errors::errorMessage("changeSign not available for requested phase type"); };
    virtual void multiplyAndAdd(const Phase &slopesPhasesTemp, const double &coeff) { Errors::errorMessage("multiplyAndAdd not available for requested phase type"); };
    virtual void divide(const double &coeff) { Errors::errorMessage("divide not available for requested phase type"); };

    //specific method for Fire
    virtual void setMassFraction(double){ Errors::errorMessage("setMassFraction not available for requested phase type"); };
  
    virtual void extendedCalculusPhase(const Coord &velocity, double Temperature) { Errors::errorMessage("extendedCalculusPhase not available for requested phase type"); };

  protected:

};

#endif // PHASE_H
