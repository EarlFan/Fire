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

#ifndef RUN_H
#define RUN_H

//! \file      Run.h
//! \author    F. Petitpas, K. Schmidmayer
//! \version   1.1
//! \date      June 5 2019

class Run;

#include <iostream>
#include <fstream>
#include <vector>
#include <list>
#include <cmath>
#include <ctime>
#include <algorithm>
#include <sstream>
#include "Tools.h"
#include "Order1/Cell.h"
#include "Models/HeaderPhase.h"
#include "Order1/CellInterface.h"
#include "Parallel/Parallel.h"
#include "Meshes/HeaderMesh.h"
#include "BoundConds/HeaderBoundCond.h"
#include "Eos/HeaderEquationOfState.h"
#include "Models/HeaderModel.h"
#include "Geometries/HeaderGeometricalDomain.h"
#include "Order2/HeaderLimiter.h"
#include "AdditionalPhysics/HeaderQuantitiesAddPhys.h"
#include "AdditionalPhysics/HeaderAddPhys.h"

#include "InputOutput/Input.h"
#include "InputOutput/Output.h"
#include "timeStats.h"

#include "Relaxations/HeaderRelaxations.h"

//! \class     Run
//! \brief     Class regrouping all information for a simulation
class Run
{
  public:
    Run(std::string nameCasTest, const int &number);
    ~Run();

    //! \brief    Initialization of the simulation
    void initialize(int argc, char* argv[]);
    //! \brief    Hyperbolic resolution + Relaxations + Source terms integration
    //! \details  Hyperbolic part is solved using Finite Volume method : \f[ \frac{U^{n+1}_i-U^{n}_i}{\Delta t} =  -\sum_{faces} \vec{F}^*_f \cdot \vec{n}_f \f]
    void solver();
    //! \brief    Cleaning simulation
    //! \details  Memory desallocations
    void finalize();
    
    void restartSimulation();

    //Accessors
    const int& getNumberPhases() const { return m_numberPhases; };

    void set_numFichier_run(int numFichier) {m_numFichier_run=numFichier;}

    int set_numFichier_run() {return m_numFichier_run;} 

  private:   

    //Specific solvers
    void integrationProcedure(double &dt, int lvl, double &dtMax, int &nbCellsTotalAMR);
    void advancingProcedure(double &dt, int &lvl, double &dtMax);
    void solveHyperbolic(double &dt, int &lvl, double &dtMax);
    void solveHyperbolicO2(double &dt, int &lvl, double &dtMax);
    void solveAdditionalPhysics(double &dt, int &lvl);
    void solveSourceTerms(double &dt, int &lvl);
    void solveRelaxations(int &lvl);
    void verifyErrors() const;
    void prepareAdditionalPhysics(double &dt, int &lvl);
    void resetDelteRhoHInCell();
    void delteRhoHInCell(const double &dt_lvl_0);
    void calcDivVorVInCell();
    void calcGradRhoInCell();
    void calcTotalLeafCell();
    void getPeakValues();
    void calcHRR();

    int m_numTest;                             //!<Number of the simulation

    //Input attributes
    std::string m_simulationName;              //!<Name of the simulation
    bool m_controleIterations;                 //!Choice for time control mode (iteration or physical time)
    int m_nbIte, m_freq;                       //!<Requested number of final time iteration and frequency
    float m_finalPhysicalTime, m_timeFreq;     //!<Requested final physical time of the simulation and time frequency for output printing
    double m_cfl;                              //!<CFL criteria (between 0 and 1)
    int m_numberPhases;                        //!<Number of phases
    int m_numberEos;                           //!<Number of equations of states
    int m_numberTransports;                    //!<Number of additional transport variables
    int m_numberAddPhys;                       //!<Number of additional physical effects
    int m_numberSources;                       //!<Number of additional source terms
    int m_dimension;                           //!<dimension 1, 2 ou 3
    int m_MRF;                                 //!<source term for Moving Reference Frame computation index(in the list of source term)
    std::string m_order;                       //!<Precision scheme order (firstorder or secondOrder)
    std::string m_symm;                        //!<symm flag

    //Specific to AMR method
    int m_lvlMax;                              //!<Maximum AMR level (if 0, then no AMR)
    int m_lvlHydro;                            //!<Maximum AMR level for Hydrodynamic process
    int m_lvlChem;                            //!<Maximum AMR level for Chemical process
    int m_lvlGeo;                            //!<Maximum AMR level for Chemical process
    int m_nbCellsTotalAMR;                     //!<Number de mailles total maximum durant la simulation

    //Geometrical attributes
    bool m_parallelPreTreatment;               //!<Choice for mesh parallel pre-treatment  (needed for first simulation on a new parallel unstructured geometry)
    
    //Calcul attributes
    Mesh *m_mesh;                              //!<Mesh type object: contains all geometrical properties of the simulation
    Model *m_model;                            //!<Model type object: contains the flow model methods
    TypeMeshContainer<Cell *> *m_cellsLvl;                   //!<Array of vectors (one per level) of computational cell objects: Contains physical fluid states.
    TypeMeshContainer<Cell *> *m_cellsLvlGhost;              //!<Array of vectors (one per level) of ghost cell objects.
    TypeMeshContainer<CellInterface *> *m_cellInterfacesLvl; //!<Array of vectors (one per level) of interface objects between cells (or between a cell and a physical domain boundary)
    Eos **m_eos;                               //!<Array of Equations of states: Contains fluid EOS parameters
    std::vector<AddPhys*> m_addPhys;           //!<Vector of Additional physics
    Symmetry *m_symmetry;                      //!<Specific object for symmetry (cylindrical or spherical) if active
    Symmetry *m_symmetryAddPhys;               //!<Object containing the parent class of symmetry to trick the corresponding additional physics argument (avoid taking into account symmetry terms multipled times)
    std::vector<Source*> m_sources;            //!<Vector of source terms
    Limiter *m_globalLimiter;                  //!<Slope limiter type object for second order in space
    Limiter *m_interfaceLimiter;               //!<Slope limiter type object for second order in space specific to interface location
    Limiter *m_globalVolumeFractionLimiter;    //!<Slope limiter type object for second order in space specific to interface advected variables (alpha, transports)
    Limiter *m_interfaceVolumeFractionLimiter; //!<Slope limiter type object for second order in space specific to interface advected variables (alpha, transports) and to interface location
    std::vector<std::string> m_nameGTR;        //!<Vector of transport variable names
    std::vector<std::string> m_nameQPA;        //!<Vector of names of the quantities of additional physics
    std::vector<std::string> m_nameGPH;        //!<Vector of phasic variables name
    double m_dt;                               //!<Explicit time step
    double m_dtNext;                           //!<Next time step
    double m_physicalTime;                     //!<Physical time
    int m_iteration;                           //!<time iteration number
    int m_restartSimulation;                   //!<File number for restarting a simulation
    int m_restartAMRsaveFreq;                  //!<Frequency at which a save to restart a simulation is done (usefull only for AMR)

    //Input/Output attributes
	Input* m_input;						       //!<Input object
    Output* m_outPut;                          //!<Main output object
    std::vector<Output *> m_cuts;              //!<Vector of output objects for cuts
    std::vector<Output *> m_probes;            //!<Vector of output objects for probes
	std::vector<Output*> m_globalQuantities;     //!<Vector of output objects for global quantities (mass, total energy)
    timeStats m_stat;                          //!<Object linked to computational time statistics
    double *m_pMax, *m_pMaxWall;               //!<Maximal pressure found between each written output and its corresponding coordinate (only for few test cases)
    double m_massWanted, m_alphaWanted;        //!<Mass and corresponding volume fraction for special output (only for few test cases)
    double totalMass[NS];
    double totalMassPara[NS];
    double posVorticity, posVorticityPara, negVorticity, negVorticityPara;
    double totalVorticity, totalVorticityPara,enstropy, enstropyPara;
    double Tmax,TmaxGlobal,Pmax,PmaxGlobal,rhoH2O_Max,rhoH2O_MaxGlobal;
    double TmaxA,TmaxAGlobal,PmaxA,PmaxAGlobal,rhoH2OA_Max,rhoH2OA_MaxGlobal;
    double Tmin, TminGlobal, Pmin, PminGlobal;
    double XL, XLGlobal, XR, XRGlobal, YR, YRGlobal;
    double flameArea,flameAreaGlobal;
    double bubbleVolume, bubbleVolumeGlobal;
    double totalFormationEnthalpy, totalFormationEnthalpyPara;
    double totalChemicalRate[NS],totalChemicalRatePara[NS];
    double heatReleaseRate, heatReleaseRatePara;
    double rhoH2Orate_max, rhoH2Orate_maxPara, YH2O_max, YH2O_maxPara;
    double magGradOH_max, magGradOH_maxPara;
    double flameVolume, flameVolumeGlobal;
    int m_numFichier_run;
    std::shared_ptr<Cantera::Solution> m_sol;

    double total_tau_c;
    double total_tau_r;
    double total_tau_c_para, total_tau_r_para;
    int total_leaf_cells, total_leaf_cells_para;
    Mixture* m_mixShock;
    //double sigma_c, sigma_r;

    double total_gradRho;
    double total_gradRho_para;

    double chem_dt;
    bool flag_calc_chem_gap;

    friend class Input;
    friend class Output;
    friend class OutputXML;
    friend class OutputGNU;
    friend class OutputProbeGNU;
	  friend class OutputGlobalGNU;
    friend class Mesh;
};

#endif // RUN_H
