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

//! \file      Run.cpp
//! \author    F. Petitpas, K. Schmidmayer
//! \version   1.1
//! \date      June 5 2019

#include "Run.h"

using namespace tinyxml2;

//***********************************************************************

Run::Run(std::string nameCasTest, const int &number) : m_numberTransports(0), m_restartSimulation(0), m_restartAMRsaveFreq(0),
  m_dt(1e-15), m_physicalTime(0.), m_iteration(0), m_simulationName(nameCasTest), m_numTest(number), m_MRF(-1)
{
  m_stat.initialize();
  AMRPara::sigma_c = 1e30; AMRPara::sigma_r = 1e30;
  chem_dt = 0;
  flag_calc_chem_gap = 0;
}

//***********************************************************************

Run::~Run(){}

//***********************************************************************

void Run::initialize(int argc, char* argv[])
{

  //1) Reading input file (XML format)
  //----------------------------------
  std::vector<GeometricalDomain*> domains;
  std::vector<BoundCond*> boundCond;
  try {
    m_input = new Input(this);
    m_input->lectureInputXML(domains, boundCond);
  }
  catch (ErrorXML &) { throw; }
  TB = new Tools(m_numberPhases);

  //2) Initialization of parallel computing (also needed for 1 CPU)
  //---------------------------------------------------------------
  parallel.initialization(argc, argv);
  if (Ncpu > 1){
    MPI_Barrier(MPI_COMM_WORLD);
    if (rankCpu == 0) std::cout << "T" << m_numTest << " | Number of CPU: " << Ncpu << std::endl;
  }
  //3) Mesh data initialization
  //---------------------------
  m_mesh->attributLimites(boundCond);
  m_cellsLvl = new TypeMeshContainer<Cell *>[m_lvlMax + 1];
  m_cellsLvlGhost = new TypeMeshContainer<Cell *>[m_lvlMax + 1];
  m_cellInterfacesLvl = new TypeMeshContainer<CellInterface *>[m_lvlMax + 1];
  try {
    if (m_restartSimulation > 0) {
      AMRPara::sigma_c = 1e30; AMRPara::sigma_r = 1e30;
      if (rankCpu == 0) std::cout << "Restarting simulation from result file number: " << m_restartSimulation << "..." << std::endl;
      m_outPut->readInfos();
      if (m_mesh->getType() == AMR) {
        if (m_restartSimulation % m_restartAMRsaveFreq == 0) {
          m_outPut->readDomainDecompostion(m_mesh, m_restartSimulation);
        }
        else {
          Errors::errorMessage("Run::restartSimulation: Restart files not available");
        }
      }
    }
    // the CellO2/CellO2Ghost is determined by push_back() method in this function. - fane
    m_dimension = m_mesh->initializeGeometrie(m_cellsLvl[0], m_cellsLvlGhost[0], m_cellInterfacesLvl[0], m_restartSimulation, m_parallelPreTreatment, m_order);
  }
  catch (ErrorECOGEN &) { throw; }

  //4) Main array initialization using model and phase number
  //---------------------------------------------------------
  for (int i = 0; i < m_cellsLvl[0].size(); i++) { m_cellsLvl[0][i]->allocate(m_numberPhases, m_numberTransports, m_addPhys, m_model); }
  for (int i = 0; i < m_cellsLvlGhost[0].size(); i++) { m_cellsLvlGhost[0][i]->allocate(m_numberPhases, m_numberTransports, m_addPhys, m_model); }
  //Attribution model and slopes to faces
  for (int i = 0; i < m_mesh->getNumberFaces(); i++) { m_cellInterfacesLvl[0][i]->associeModel(m_model); }

  //5) Physical data initialization: filling fluid states
  //-----------------------------------------------------
  {
    Mixture *mixLP, *mixShock, *mixBubble;
    double shockPosition, radius;
    Coord center;
    for (unsigned int geom = 0; geom < domains.size(); geom++) 
    {
      // this fillin will copyPhase and setTransport
      if(domains[geom]->getName() == "BASE")
      {
        mixLP = domains[geom]->getMixture();
      }
      else if(domains[geom]->getName() == "HP")
      {
        mixShock = domains[geom]->getMixture();
        shockPosition = domains[geom]->getHalfPosition();
      }
      else if(domains[geom]->getName() == "BUBBLE")
      {
        mixBubble = domains[geom]->getMixture();
        // if(rankCpu==0) printf("In input, bubble T = %le.\n",mixBubble->getTemperature());
        center = domains[geom]->centerDisc();
        radius = domains[geom]->getRadius();
      }
    }
    for (int i = 0; i < m_cellsLvl[0].size(); i++) 
    // { m_cellsLvl[0][i]->fillFire(domains, m_lvlMax, mixLP, mixShock, mixBubble, shockPosition, center, radius); }
    { m_cellsLvl[0][i]->fill(domains, m_lvlMax); }
    for (int i = 0; i < m_cellsLvlGhost[0].size(); i++) 
    // { m_cellsLvlGhost[0][i]->fillFire(domains, m_lvlMax, mixLP, mixShock, mixBubble, shockPosition, center, radius);}
    { m_cellsLvlGhost[0][i]->fill(domains, m_lvlMax);}
  }

  //EOS filling
  m_cellsLvl[0][0]->allocateEos(m_numberPhases, m_model);
  //Complete fluid state with additional calculations (sound speed, energies, mixture variables, etc.)
  for (int i = 0; i < m_cellsLvl[0].size(); i++) { 
    m_cellsLvl[0][i]->completeFulfillState(); 
    }

  //6) Allocate Sloped and buffer Cells for Riemann problems
  //--------------------------------------------------------
  int allocateSlopeLocal = 0;
  for (int i = 0; i < m_mesh->getNumberFaces(); i++) { m_cellInterfacesLvl[0][i]->allocateSlopes(m_numberPhases, m_numberTransports, allocateSlopeLocal); }
  cellLeft = new CellO2; cellRight = new CellO2;
  cellLeft->allocate(m_numberPhases, m_numberTransports, m_addPhys, m_model);
  cellRight->allocate(m_numberPhases, m_numberTransports, m_addPhys, m_model);
  domains[0]->fillIn2(cellLeft, m_numberPhases, m_numberTransports);
  domains[0]->fillIn2(cellRight, m_numberPhases, m_numberTransports);

  //7) Initialization of persistant communications for parallel computing
  //--------------------------------------------------------------------
  MPI_Barrier(MPI_COMM_WORLD);
  m_mesh->initializePersistentCommunications(m_numberPhases, m_numberTransports, m_cellsLvl[0], m_order);
  if (Ncpu > 1) { parallel.communicationsPrimitives(m_eos, 0); }
  
  //8) AMR initialization
  //---------------------
  m_mesh->procedureRaffinementInitialization(m_cellsLvl, m_cellsLvlGhost, m_cellInterfacesLvl, m_addPhys, m_model,
    m_nbCellsTotalAMR, domains, m_eos, m_restartSimulation, m_order, m_numberPhases, m_numberTransports);

  // // for restart the simulation from DME spherical flame. - fane
  // for (unsigned int geom = 0; geom < domains.size(); geom++) 
  // {
  //   if(domains[geom]->getName() == "HP")
  //   {
  //     m_mixShock = domains[geom]->getMixture();
  //   }
  // }

  //9) Output file preparation
  //--------------------------
  m_outPut->prepareOutput(*cellLeft);
  for (unsigned int c = 0; c < m_cuts.size(); c++) m_cuts[c]->prepareOutput(*cellLeft);
  for (unsigned int p = 0; p < m_probes.size(); p++) m_probes[p]->prepareOutput(*cellLeft);
  for (unsigned int g = 0; g < m_globalQuantities.size(); g++) m_globalQuantities[g]->prepareOutput(*cellLeft);

  //10) Restart simulation
  //----------------------
  if (m_restartSimulation > 0) {
    try { this->restartSimulation(); }
    catch (ErrorECOGEN &) { throw; }

  }

  // delete the domains after restarting
  for (unsigned int d = 0; d < domains.size(); d++) { delete domains[d]; }

  //11) Printing t0 solution
  //------------------------
  if (m_restartSimulation == 0) {
    try {
      //Only for few test cases
      //-----
      //m_pMax = new double[4];
      //m_pMaxWall = new double[4];
      //m_pMax[0] = 0.;
      //m_pMaxWall[0] = 0.;
      //for (unsigned int c = 0; c < m_cellsLvl[0].size(); c++) { m_cellsLvl[0][c]->lookForPmax(m_pMax, m_pMaxWall); }
      //if (Ncpu > 1) { parallel.computePMax(m_pMax[0], m_pMaxWall[0]); }
      // m_massWanted = 0.; //Initial mass
      // m_alphaWanted = -1.;
      // for (unsigned int c = 0; c < m_cellsLvl[0].size(); c++) { m_cellsLvl[0][c]->computeMass(m_massWanted, m_alphaWanted); }
      // if (Ncpu > 1) { parallel.computeMassTotal(m_massWanted); }
      // m_massWanted = 0.9*m_massWanted; //Percentage of the initial mass we want
      // m_alphaWanted = 0.;
      //-----
      m_outPut->prepareOutputInfos();
      // print the output info in the console. - fane
      if (rankCpu == 0) m_outPut->ecritInfos();
      m_outPut->saveInfosMailles();
      if (m_mesh->getType() == AMR) m_outPut->printTree(m_mesh, m_cellsLvl, m_restartAMRsaveFreq);
      for (unsigned int c = 0; c < m_cuts.size(); c++) 
      { 
        m_cuts[c]->setNumSortie(m_outPut->getNumSortie());
        m_cuts[c]->ecritSolution(m_mesh, m_cellsLvl); 
      }
      for (unsigned int p = 0; p < m_probes.size(); p++) { if (m_probes[p]->possesses()) m_probes[p]->ecritSolution(m_mesh, m_cellsLvl); }
	  for (unsigned int g = 0; g < m_globalQuantities.size(); g++) { m_globalQuantities[g]->ecritSolution(m_mesh, m_cellsLvl); }
	  m_outPut->ecritSolution(m_mesh, m_cellsLvl);

	  //double totalMass(0.), massBuff(0.);
	  //for (int c = 0; c < m_cellsLvl[0].size(); c++) {
		 // massBuff += m_cellsLvl[0][c]->getPhase(0)->getDensity() * m_cellsLvl[0][c]->getElement()->getVolume();
	  //}

	  //MPI_Allreduce(&massBuff, &totalMass, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

	  //if (rankCpu == 0) {
		 // std::cout << "Total mass : " << totalMass << std::endl;
		 // totalMass = 0.;
	  //}
	  
    }
    catch (ErrorXML &) { throw; }
  }

  MPI_Barrier(MPI_COMM_WORLD);
  if (rankCpu == 0) std::cout << "Run::initialize OK" << std::endl;
}

//***********************************************************************

void Run::restartSimulation()
{
  std::ifstream fileStream;

  //Reconstruct the mesh and get physical data from restart point
  try {
    if (m_mesh->getType() == AMR) {
      if (m_restartSimulation % m_restartAMRsaveFreq == 0) {
        m_outPut->readTree(m_mesh, m_cellsLvl, m_cellsLvlGhost, m_cellInterfacesLvl, m_addPhys, m_model, m_eos, m_nbCellsTotalAMR);
      }
    }
    m_outPut->readResults(m_mesh, m_cellsLvl);

    // for DME shock-cool flame interaction study
    // if(m_restartSimulation == 1)
    // {
    //   // revise the post-shock states
    //   for(int lvl = m_lvlMax; lvl>=0; lvl--)
    //   {
    //     double x_c,y_c;
    //     for (unsigned int i = 0; i < m_cellsLvl[lvl].size(); i++) 
    //     {
    //       x_c = m_cellsLvl[lvl][i]->getPosition().getX();
    //       if(x_c<0.11)
    //       {
    //         m_cellsLvl[lvl][i]->copyMixture(m_mixShock);
    //       }
    //     }
    //   }
    // }
  }
  catch (ErrorECOGEN &) { fileStream.close(); throw; }
  fileStream.close();

  //Communicate physical data between processors and complete fluid state with additional calculations (sound speed, energies, mixture variables, etc.)
  if (Ncpu > 1) {
    for (int lvl = 0; lvl <= m_lvlMax; lvl++) {
      parallel.communicationsPrimitives(m_eos, lvl);
      parallel.communicationsTransports(lvl);
    }
  }
  for (int lvl = 0; lvl <= m_lvlMax; lvl++) { //With reduced output //KS//FP// To implement dynamically
    for (unsigned int i = 0; i < m_cellsLvl[lvl].size(); i++) { m_cellsLvl[lvl][i]->completeFulfillState(restart); }
  }
  // for (int lvl = 0; lvl <= m_lvlMax; lvl++) { //With complete output //KS//FP// To implement dynamically
  //   for (unsigned int i = 0; i < m_cellsLvl[lvl].size(); i++) { m_cellsLvl[lvl][i]->fulfillState(restart); }
  //   if (m_numberAddPhys) {
  //       for (unsigned int i = 0; i < m_cellsLvl[lvl].size(); i++) { if (!m_cellsLvl[lvl][i]->getSplit()) { m_cellsLvl[lvl][i]->prepareAddPhys(); } }
  //       if (Ncpu > 1) {
  //         m_stat.startCommunicationTime();
  //         for (unsigned int pa = 0; pa < m_addPhys.size(); pa++) { m_addPhys[pa]->communicationsAddPhys(m_numberPhases, m_dimension, lvl); }
  //         m_stat.endCommunicationTime();
  //       }
  //   }
  // }
  if (m_mesh->getType() == AMR) {
    for (int lvl = 0; lvl < m_lvlMax; lvl++) {
      for (unsigned int i = 0; i < m_cellsLvl[lvl].size(); i++) { m_cellsLvl[lvl][i]->averageChildrenInParent(); }
    }
  }
  if (Ncpu > 1) {
    for (int lvl = 0; lvl <= m_lvlMax; lvl++) { parallel.communicationsPrimitives(m_eos, lvl); }
  }
  //KS//FP//DEV// Apparemment fulfillState avec Prim::restart n'est pas a jour dans tous les modeles (seulement pour Kapila)

  MPI_Barrier(MPI_COMM_WORLD);
  if (rankCpu == 0) std::cout << "restartSimulation OK" << std::endl;
}

//***********************************************************************

void Run::solver()
{
  int nbCellsTotalAMRMax = m_nbCellsTotalAMR;
  double dtMax;

  //-------------------
  //Time iterative loop
  //-------------------
  bool computeFini(false); bool print(false);
  double printSuivante(m_physicalTime+m_timeFreq);
  while (!computeFini) {
    //Errors checking
    try {
      this->verifyErrors();
    }
    catch (ErrorECOGEN &) { throw; }
		
    //------------------- INTEGRATION PROCEDURE -------------------
    
    //Setting cons variable to zero for spatial scheme on dU/dt: no need for time step at this point
    for (unsigned int i = 0; i < m_cellsLvl[0].size(); i++) { m_cellsLvl[0][i]->setToZeroConsGlobal(m_numberPhases, m_numberTransports); }
    dtMax = 1.e10;
    int lvlDep = 0;
    resetDelteRhoHInCell();

    // the main integration function.
    this->integrationProcedure(m_dt, lvlDep, dtMax, m_nbCellsTotalAMR);

    if (AMRPara::div_vor_v_flag == 1) {
      if(m_iteration == 0) {
        AMRPara::sigma_c = 1e30; AMRPara::sigma_r = 1e30;
      }
      else {
        calcDivVorVInCell();
      }
    }

    if (AMRPara::gradRho_flag == 1) {
      if(m_iteration == 0) {
        calcTotalLeafCell();
        AMRPara::sigma_gradRho = 1e30;
      }
      else {
        // check whether AMR is turned on
        if (m_lvlMax>0) calcGradRhoInCell();
      }
    }
    else
    {
      calcTotalLeafCell();
    }

    AMRPara::physical_time = m_physicalTime + m_dt;

    getPeakValues();
    if(rankCpu==0)
    {
      printf("**********************************************\n");
      printf("After step %d, t = %le, dt = %le, lvlMax = %d, lvlChem = %d, lvlGeo = %d.\n",m_iteration,m_physicalTime+m_dt,m_dt,m_lvlMax,m_lvlChem,m_lvlGeo);
      if (AMRPara::gradRho_flag == 1) 
      {
        if(m_iteration == 0) {
          if(rankCpu == 0) printf("total_leaf_cells = %d.\n",total_leaf_cells_para);
        }
        else {
          if(rankCpu == 0)
          {
            printf("total_leaf_cells = %d.\n",  total_leaf_cells_para);
          }
        }
      }
      else
      {
        if(rankCpu == 0) printf("total_leaf_cells_para = %d.\n",total_leaf_cells_para);
      }
      // printf("T_max = %le K, T_min = %le K, p_max = %le Pa, p_min = %le Pa.\n", TmaxGlobal, TminGlobal, PmaxGlobal, PminGlobal);
      printf("T_range = [%.2f, %.2f] K, p_range = [%.4E, %.4E] Pa\n", 
        TminGlobal, TmaxGlobal, PminGlobal, PmaxGlobal); 
      printf("**********************************************\n");
    }
    //-------------------- CONTROL ITERATIONS/TIME ---------------------

    //Still alive...
    if (m_iteration !=0 && m_iteration % 1000 == 0 && rankCpu == 0) { 
      std::cout << "T" << m_numTest << " | Iteration " << m_iteration << " / Timestep " << m_dt << " / Progress " << m_physicalTime / m_finalPhysicalTime * 100. << "%" << std::endl;
    }

    m_physicalTime += m_dt;
    m_iteration++;
    //Managing output files printing / End of time iterative loop
    if (m_controleIterations) {
      if (m_iteration%m_freq == 0) { print = true; }
      if (m_iteration >= m_nbIte) { computeFini = true; }
    }
    else {
      if (m_physicalTime >= printSuivante) { print = true; printSuivante += m_timeFreq; }
      if (m_physicalTime >= m_finalPhysicalTime) 
      { 
        print = true; 
        computeFini = true; 
      }
    }
    //Managing Sources evolutions
    for (unsigned int s = 0; s < m_sources.size(); s++) { m_sources[s]->sourceEvolution(m_physicalTime); }

    #ifdef DEBUG
      if (rankCpu == 0) std::cout << m_iteration << ", dt = " << m_dt << std::endl;
    #endif // DEBUG


    //------------------------ OUTPUT FILES PRINTING -------------------------
    nbCellsTotalAMRMax = std::max(nbCellsTotalAMRMax, m_nbCellsTotalAMR);
    m_dtNext = m_cfl * dtMax;
    if (Ncpu > 1) { parallel.computeDt(m_dtNext); }
    if (print) {
      calcHRR();
      m_stat.updateComputationTime();
      //General printings
      //Only for few test case
      //-----
      //if (Ncpu > 1) { parallel.computePMax(m_pMax[0], m_pMaxWall[0]); }
      //m_pMax[0] = 0.;
      //m_pMaxWall[0] = 0.;
      // m_alphaWanted = 1.; //Volume fraction corresponding to the wanted mass
      // double mass(0.);
      // do {
      //   m_alphaWanted -= 0.01;
      //   mass = 0.;
      //   for (unsigned int c = 0; c < m_cellsLvl[0].size(); c++) { m_cellsLvl[0][c]->computeMass(mass, m_alphaWanted); }
      //   if (Ncpu > 1) { parallel.computeMassTotal(mass); }
      // } while (mass < 0.999*m_massWanted && m_alphaWanted > 0.001);
      // if (m_alphaWanted < 1.e-10) m_alphaWanted = 0.;
      //-----
      if (rankCpu == 0) m_outPut->ecritInfos();
      m_outPut->saveInfosMailles();
      if (m_mesh->getType() == AMR) m_outPut->printTree(m_mesh, m_cellsLvl, m_restartAMRsaveFreq);
      // add setNumSortie to correctly set m_numFichier_run - fane 20210427
      for (unsigned int c = 0; c < m_cuts.size(); c++) 
      { 
        m_cuts[c]->setNumSortie(m_outPut->getNumSortie());
        m_cuts[c]->ecritSolution(m_mesh, m_cellsLvl); 
      }
      for (unsigned int g = 0; g < m_globalQuantities.size(); g++) { m_globalQuantities[g]->ecritSolution(m_mesh, m_cellsLvl); }
      m_outPut->ecritSolution(m_mesh, m_cellsLvl);

      // this->set_numFichier_run(m_outPut->getNumSortie());

      if (rankCpu == 0) std::cout << "OK" << std::endl;
      print = false;
	}
    //Printing probes data
    for (unsigned int p = 0; p < m_probes.size(); p++) { 
      if((m_probes[p]->possesses()) && m_probes[p]->getNextTime()<=m_physicalTime) m_probes[p]->ecritSolution(m_mesh, m_cellsLvl);
    }

    //-------------------------- TIME STEP UPDATING --------------------------

    m_dt = m_dtNext;
    // m_dt = 1e-4;

  } //time iterative loop end
  if (rankCpu == 0) std::cout << "T" << m_numTest << " | -------------------------------------------" << std::endl;
  MPI_Barrier(MPI_COMM_WORLD);
  if (m_mesh->getType() == AMR) {
    double localLoad(0.);
    for (int lvl = m_lvlMax; lvl >= 0; lvl--) {
      for (unsigned int i = 0; i < m_cellsLvl[0].size(); i++) {
      m_cellsLvl[0][i]->computeLoad(localLoad, lvl);
      }
    }
    std::cout << "T" << m_numTest << " | Final local load on CPU " << rankCpu << " : " << localLoad << std::endl;
  }
}

//***********************************************************************
void Run::calcDivVorVInCell()
{
  total_tau_c = 0;
  total_tau_r = 0;
  total_leaf_cells = 0;
  double tau_c_t, tau_r_t, cell_l, f;

  for(int lvl = m_lvlMax; lvl>=0; lvl--)
  {
    for (unsigned int i = 0; i < m_cellsLvl[lvl].size(); i++) 
    {
      if (!m_cellsLvl[lvl][i]->getSplit())
      {
        tau_c_t = m_cellsLvl[lvl][i]->calcTauC();
        tau_r_t = m_cellsLvl[lvl][i]->calcTauR();
        
        total_leaf_cells++;
        total_tau_c += pow(tau_c_t, 2);
        total_tau_r += pow(tau_r_t, 2);
      }
    }
  }

  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Allreduce(&total_tau_c, &total_tau_c_para, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&total_tau_r, &total_tau_r_para, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&total_leaf_cells, &total_leaf_cells_para, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  AMRPara::sigma_c = sqrt(total_tau_c_para/total_leaf_cells_para);
  AMRPara::sigma_r = sqrt(total_tau_r_para/total_leaf_cells_para);

  if(rankCpu == 0)
  {
    printf("On cpu = %d, sigma_c = %le, sigma_r = %le.\n", rankCpu, AMRPara::sigma_c, AMRPara::sigma_r);
  }
}


//***********************************************************************

void Run::integrationProcedure(double &dt, int lvl, double &dtMax, int &nbCellsTotalAMR)
{
  //1) AMR Level time step determination
  double dtLvl = dt * std::pow(2., -(double)lvl); 
  
  //2) Refinement procedure
  if (m_lvlMax > 0) {
    // communicate gradients before AMR refinement
    // shall be revised in the future as only gradRho, gradUV is required for AMR process.
    for (unsigned int i = 0; i < m_cellsLvl[lvl].size(); i++) { if (!m_cellsLvl[lvl][i]->getSplit()) { m_cellsLvl[lvl][i]->prepareAddPhys(); } }
    if (Ncpu > 1) {
      m_stat.startCommunicationTime();
      for (unsigned int pa = 0; pa < m_addPhys.size(); pa++) { m_addPhys[pa]->communicationsAddPhys(m_numberPhases, m_dimension, lvl); }
      m_stat.endCommunicationTime();
    }

    m_stat.startAMRTime();

    m_mesh->procedureRaffinement(m_cellsLvl, m_cellsLvlGhost, m_cellInterfacesLvl, lvl, m_addPhys, m_model, nbCellsTotalAMR, m_eos);

    if (Ncpu > 1) { if (lvl == 0) { if (m_iteration % (static_cast<int>(1./m_cfl/0.6) + 1) == 0) {
      m_mesh->parallelLoadBalancingAMR(m_cellsLvl, m_cellsLvlGhost, m_cellInterfacesLvl, m_order, m_numberPhases, m_numberTransports, m_addPhys, m_model, m_eos, nbCellsTotalAMR);
      for (unsigned int p = 0; p < m_probes.size(); p++) {
        m_probes[p]->locateProbeInMesh(m_cellsLvl[0], m_mesh->getNumberCells()); //Locate new probes CPU after Load Balancing
      }
    } } }

    m_stat.endAMRTime();
  }

  //3) Slopes determination for second order and gradients for additional physics
  //Fait ici pour avoir une mise a jour d'effectuer lors de l'execution de la procedure de niveau lvl+1 (donc pour les slopes plus besoin de les faire au debut de resolHyperboliqueO2)
  if (m_order == "SECONDORDER") {
    for (unsigned int i = 0; i < m_cellInterfacesLvl[lvl].size(); i++) { if (!m_cellInterfacesLvl[lvl][i]->getSplit()) { m_cellInterfacesLvl[lvl][i]->computeSlopes(m_numberPhases, m_numberTransports); } }
    if (Ncpu > 1) {
      m_stat.startCommunicationTime();
      parallel.communicationsSlopes(lvl);
      if (lvl > 0) { parallel.communicationsSlopes(lvl - 1); }
      m_stat.endCommunicationTime();
    }
  }
  if (lvl < m_lvlMax) {
    if (m_numberAddPhys) {
      for (unsigned int i = 0; i < m_cellsLvl[lvl].size(); i++) { if (!m_cellsLvl[lvl][i]->getSplit()) { m_cellsLvl[lvl][i]->prepareAddPhys(); } }
      if (Ncpu > 1) {
        m_stat.startCommunicationTime();
        for (unsigned int pa = 0; pa < m_addPhys.size(); pa++) { m_addPhys[pa]->communicationsAddPhys(m_numberPhases, m_dimension, lvl); }
        m_stat.endCommunicationTime();
      }
    }
    //4) Recursive call for level up integration procedure
    this->integrationProcedure(dt, lvl + 1, dtMax, nbCellsTotalAMR);
  }

  //5) Advancement procedure
  this->advancingProcedure(dtLvl, lvl, dtMax);

  //6) Additional calculations for AMR levels > 0
  if (lvl > 0) {
    if (m_order == "SECONDORDER") {
      for (unsigned int i = 0; i < m_cellInterfacesLvl[lvl].size(); i++) { 
        if (!m_cellInterfacesLvl[lvl][i]->getSplit()) { 
          m_cellInterfacesLvl[lvl][i]->computeSlopes(m_numberPhases, m_numberTransports); 
        }
      }
      if (Ncpu > 1) {
        m_stat.startCommunicationTime();
        parallel.communicationsSlopes(lvl);
        if (lvl > 0) { parallel.communicationsSlopes(lvl - 1); }
        m_stat.endCommunicationTime();
      }
    }
    if (lvl < m_lvlMax) { this->integrationProcedure(dt, lvl + 1, dtMax, nbCellsTotalAMR); }
    this->advancingProcedure(dtLvl, lvl, dtMax);
  }
}

//***********************************************************************

void Run::advancingProcedure(double &dt, int &lvl, double &dtMax)
{
  //1) Finite volume scheme for hyperbolic systems (Godunov or MUSCL)
  if (m_order == "FIRSTORDER") { this->solveHyperbolic(dt, lvl, dtMax); }
  else { this->solveHyperbolicO2(dt, lvl, dtMax); }


  if (m_numberSources) this->solveSourceTerms(dt, lvl);

  //4) Relaxations to equilibria
  if (m_numberPhases > 1) this->solveRelaxations(lvl);

  //5) Averaging childs cells in mother cell (if AMR)
  if (lvl < m_lvlMax) { for (unsigned int i = 0; i < m_cellsLvl[lvl].size(); i++) { m_cellsLvl[lvl][i]->averageChildrenInParent(); } }

  //6) Final communications
  if (Ncpu > 1) {
    m_stat.startCommunicationTime();
    parallel.communicationsPrimitives(m_eos, lvl);
    m_stat.endCommunicationTime();
  }
  //Only for few test case
  //-----
  //for (unsigned int i = 0; i < m_cellsLvl[lvl].size(); i++) { m_cellsLvl[lvl][i]->lookForPmax(m_pMax, m_pMaxWall); }
  //-----
}

//***********************************************************************

void Run::solveHyperbolicO2(double &dt, int &lvl, double &dtMax)
{
  //1) m_cons saves for AMR/second order combination
  //------------------------------------------------
  for (unsigned int i = 0; i < m_cellsLvl[lvl].size(); i++) {
    if (!m_cellsLvl[lvl][i]->getSplit()) { 
      m_cellsLvl[lvl][i]->saveCons(m_numberPhases, m_numberTransports); 
    } 
  }

  //2) Spatial second order scheme
  //------------------------------
  //Fluxes are determined at each cells interfaces and stored in the m_cons variable of corresponding cells. Hyperbolic maximum time step determination
  // compute and add inviscid flux to fluxBuffer, no dt
  for (unsigned int i = 0; i < m_cellInterfacesLvl[lvl].size(); i++) 
  { 
    if (!m_cellInterfacesLvl[lvl][i]->getSplit()) 
    { 
      m_cellInterfacesLvl[lvl][i]->computeFlux(m_numberPhases, m_numberTransports, 
        dtMax, *m_globalLimiter, *m_interfaceLimiter, *m_globalVolumeFractionLimiter, 
        *m_interfaceVolumeFractionLimiter); 
    } 
  }

  //2.1) if viscous is on
  if (m_numberAddPhys) 
  {
    //1) Preparation of variables for additional (gradients computations, etc) and communications
    //-------------------------------------------------------------------------------------------
    if (Ncpu > 1) {
      m_stat.startCommunicationTime();
      parallel.communicationsPrimitives(m_eos, lvl);
      m_stat.endCommunicationTime();
    }
    for (unsigned int i = 0; i < m_cellsLvl[lvl].size(); i++) { 
      if (!m_cellsLvl[lvl][i]->getSplit()) 
      {
        m_cellsLvl[lvl][i]->prepareAddPhys(); 
      } 
    }
    // exit(EXIT_SUCCESS);
    if (Ncpu > 1) {
      m_stat.startCommunicationTime();
      for (unsigned int pa = 0; pa < m_addPhys.size(); pa++) { 
        m_addPhys[pa]->communicationsAddPhys(m_numberPhases, m_dimension, lvl); 
      }
      m_stat.endCommunicationTime();
    }

    //2) Additional physics fluxes determination (Surface tensions, viscosity, conductivity, ...)
    //-------------------------------------------------------------------------------------------
    //Calcul de la somme des flux des physiques additionnelles que l on stock dans m_cons de chaque cell
    for (unsigned int pa = 0; pa < m_addPhys.size(); pa++) {
      for (unsigned int i = 0; i < m_cellInterfacesLvl[lvl].size(); i++) { 
        if (!m_cellInterfacesLvl[lvl][i]->getSplit()) {
        // printf("cellInterface i = %d.\n",i);
        m_cellInterfacesLvl[lvl][i]->computeFluxAddPhys(m_numberPhases, *m_addPhys[pa]); 
      } 
    }
      for (unsigned int i = 0; i < m_cellsLvl[lvl].size(); i++) { 
        if (!m_cellsLvl[lvl][i]->getSplit()) { 
          m_cellsLvl[lvl][i]->addNonConsAddPhys(m_numberPhases, *m_addPhys[pa], m_symmetry); 
        } 
      }
    }

    //3) Time evolution for additional physics
    // this operation are moved to in predictionOrdre2()
    //----------------------------------------
    for (unsigned int i = 0; i < m_cellsLvl[lvl].size(); i++) {
      if (!m_cellsLvl[lvl][i]->getSplit()) {
        // m_cellsLvl[lvl][i]->timeEvolutionAddPhys(dt, m_numberPhases, m_numberTransports);   //Obtention des cons pour shema sur (Un+1-Un)/dt
        // m_cellsLvl[lvl][i]->buildPrim(m_numberPhases);                                      //On peut reconstruire Prim a partir de m_cons
        // m_cellsLvl[lvl][i]->setToZeroCons(m_numberPhases, m_numberTransports);              //Mise a zero des cons pour shema spatial sur dU/dt : permet de s affranchir du pas de temps
      }
    }
  }

  //3)Prediction step using slopes
  //------------------------------
  for (unsigned int i = 0; i < m_cellsLvl[lvl].size(); i++) { 
    if (!m_cellsLvl[lvl][i]->getSplit()) { 
      m_cellsLvl[lvl][i]->predictionOrdre2(dt, m_numberPhases, m_numberTransports, m_symmetry); 
    } 
  }
  //3b) Option: Activate relaxation during prediction //KS//FP// To implement dynamically
  //3c) Option: Activate additional physics during prediction //KS//FP// To implement dynamically
  //3d) Option: Activate source terms during prediction //KS//FP// To implement dynamically

  //4) m_cons recovery for AMR/second order combination (substotute to setToZeroCons)
  //---------------------------------------------------------------------------------
  for (unsigned int i = 0; i < m_cellsLvl[lvl].size(); i++) { 
    if (!m_cellsLvl[lvl][i]->getSplit()) { 
      m_cellsLvl[lvl][i]->recuperationCons(m_numberPhases, m_numberTransports); 
    } 
  }

  //5) vecPhasesO2 communications
  //-----------------------------
  if (Ncpu > 1) {
    m_stat.startCommunicationTime();
    parallel.communicationsPrimitives(m_eos, lvl, vecPhasesO2);
    m_stat.endCommunicationTime();
  }

  //6) Optional new slopes determination (improves code stability)
  //--------------------------------------------------------------
  for (unsigned int i = 0; i < m_cellInterfacesLvl[lvl].size(); i++) { 
    if (!m_cellInterfacesLvl[lvl][i]->getSplit()) { 
      m_cellInterfacesLvl[lvl][i]->computeSlopes(m_numberPhases, m_numberTransports, vecPhasesO2); 
    } 
  }

  if (Ncpu > 1) {
    m_stat.startCommunicationTime();
    parallel.communicationsSlopes(lvl);
    if (lvl > 0) { parallel.communicationsSlopes(lvl - 1); }
    m_stat.endCommunicationTime();
  }

  //7) Spatial scheme on predicted variables
  //----------------------------------------
  //Fluxes are determined at each cells interfaces and stored in the m_cons variableof corresponding cells. Hyperbolic maximum time step determination
  for (unsigned int i = 0; i < m_cellInterfacesLvl[lvl].size(); i++) { 
    if (!m_cellInterfacesLvl[lvl][i]->getSplit()) { 
      m_cellInterfacesLvl[lvl][i]->computeFlux(m_numberPhases, m_numberTransports, dtMax, *m_globalLimiter, *m_interfaceLimiter, *m_globalVolumeFractionLimiter, *m_interfaceVolumeFractionLimiter, vecPhasesO2); 
    } 
  }

//7.1) if viscous is on
  if (m_numberAddPhys) 
  {
    for (unsigned int i = 0; i < m_cellsLvl[lvl].size(); i++) { 
      if (!m_cellsLvl[lvl][i]->getSplit()) 
      {
        m_cellsLvl[lvl][i]->prepareAddPhys(); 
      } 
    }
    // exit(EXIT_SUCCESS);
    if (Ncpu > 1) {
      m_stat.startCommunicationTime();
      for (unsigned int pa = 0; pa < m_addPhys.size(); pa++) { m_addPhys[pa]->communicationsAddPhys(m_numberPhases, m_dimension, lvl); }
      m_stat.endCommunicationTime();
    }

    //2) Additional physics fluxes determination (Surface tensions, viscosity, conductivity, ...)
    //-------------------------------------------------------------------------------------------
    //Calcul de la somme des flux des physiques additionnelles que l on stock dans m_cons de chaque cell
    for (unsigned int pa = 0; pa < m_addPhys.size(); pa++) {
      for (unsigned int i = 0; i < m_cellInterfacesLvl[lvl].size(); i++) { if (!m_cellInterfacesLvl[lvl][i]->getSplit()) {
        // printf("cellInterface i = %d.\n",i);
        m_cellInterfacesLvl[lvl][i]->computeFluxAddPhys(m_numberPhases, *m_addPhys[pa]); } }
      for (unsigned int i = 0; i < m_cellsLvl[lvl].size(); i++) { if (!m_cellsLvl[lvl][i]->getSplit()) { m_cellsLvl[lvl][i]->addNonConsAddPhys(m_numberPhases, *m_addPhys[pa], m_symmetry); } }
    }

    //3) Time evolution for additional physics
    // this operation are performed in predictionOrdre2()
    //----------------------------------------
    for (unsigned int i = 0; i < m_cellsLvl[lvl].size(); i++) {
      if (!m_cellsLvl[lvl][i]->getSplit()) {
        // m_cellsLvl[lvl][i]->timeEvolutionAddPhys(dt, m_numberPhases, m_numberTransports);   //Obtention des cons pour shema sur (Un+1-Un)/dt
        // m_cellsLvl[lvl][i]->buildPrim(m_numberPhases);                                      //On peut reconstruire Prim a partir de m_cons
        // m_cellsLvl[lvl][i]->setToZeroCons(m_numberPhases, m_numberTransports);              //Mise a zero des cons pour shema spatial sur dU/dt : permet de s affranchir du pas de temps
      }
    }
  }

  //8) Time evolution
  //-----------------
  for (unsigned int i = 0; i < m_cellsLvl[lvl].size(); i++) {
    if (!m_cellsLvl[lvl][i]->getSplit()) {
      m_cellsLvl[lvl][i]->timeEvolution(dt, m_numberPhases, m_numberTransports, m_symmetry, vecPhasesO2);   //Obtention des cons pour shema sur (Un+1-Un)/dt
      m_cellsLvl[lvl][i]->addHalfTimeCons(m_numberPhases);
      m_cellsLvl[lvl][i]->buildPrim(m_numberPhases);                                                        //On peut reconstruire Prim a partir de m_cons
      m_cellsLvl[lvl][i]->setToZeroCons(m_numberPhases, m_numberTransports);                                //Mise a zero des cons pour shema spatial sur dU/dt : permet de s affranchir du pas de temps
    }
  }

}

//***********************************************************************

void Run::solveHyperbolic(double &dt, int &lvl, double &dtMax)
{
  //1) Spatial scheme
  //-----------------
  //Fluxes are determined at each cells interfaces and stored in the m_cons variableof corresponding cells. Hyperbolic maximum time step determination
  for (unsigned int i = 0; i < m_cellInterfacesLvl[lvl].size(); i++) 
  { if (!m_cellInterfacesLvl[lvl][i]->getSplit()) 
    { 
      m_cellInterfacesLvl[lvl][i]->computeFlux(m_numberPhases, m_numberTransports, dtMax, *m_globalLimiter, *m_interfaceLimiter, *m_globalVolumeFractionLimiter, *m_interfaceVolumeFractionLimiter); 
    } 
  }

  if (m_numberAddPhys) 
  {
    //1) Preparation of variables for additional (gradients computations, etc) and communications
    //-------------------------------------------------------------------------------------------
    if (Ncpu > 1) {
      m_stat.startCommunicationTime();
      parallel.communicationsPrimitives(m_eos, lvl);
      m_stat.endCommunicationTime();
    }
    for (unsigned int i = 0; i < m_cellsLvl[lvl].size(); i++) { 
      if (!m_cellsLvl[lvl][i]->getSplit()) 
      {
        m_cellsLvl[lvl][i]->prepareAddPhys(); 
      } 
    }
    // exit(EXIT_SUCCESS);
    if (Ncpu > 1) {
      m_stat.startCommunicationTime();
      for (unsigned int pa = 0; pa < m_addPhys.size(); pa++) { 
        m_addPhys[pa]->communicationsAddPhys(m_numberPhases, m_dimension, lvl); 
      }
      m_stat.endCommunicationTime();
    }

    //2) Additional physics fluxes determination (Surface tensions, viscosity, conductivity, ...)
    //-------------------------------------------------------------------------------------------
    //Calcul de la somme des flux des physiques additionnelles que l on stock dans m_cons de chaque cell
    for (unsigned int pa = 0; pa < m_addPhys.size(); pa++) {
      for (unsigned int i = 0; i < m_cellInterfacesLvl[lvl].size(); i++) { 
        if (!m_cellInterfacesLvl[lvl][i]->getSplit()) {
        m_cellInterfacesLvl[lvl][i]->computeFluxAddPhys(m_numberPhases, *m_addPhys[pa]); 
        } 
      }
      for (unsigned int i = 0; i < m_cellsLvl[lvl].size(); i++) { 
        if (!m_cellsLvl[lvl][i]->getSplit()) { 
          m_cellsLvl[lvl][i]->addNonConsAddPhys(m_numberPhases, *m_addPhys[pa], m_symmetry); 
          
        } 
      }
    }

  }

  //3) build primitive variable and recover m_con to Zero

  for (unsigned int i = 0; i < m_cellsLvl[lvl].size(); i++) {
    if (!m_cellsLvl[lvl][i]->getSplit()) {
      m_cellsLvl[lvl][i]->timeEvolution(dt, m_numberPhases, m_numberTransports, m_symmetry);   //Obtention des cons pour shema sur (Un+1-Un)/dt
      m_cellsLvl[lvl][i]->buildPrim(m_numberPhases);                                           //On peut reconstruire Prim a partir de m_cons
      m_cellsLvl[lvl][i]->setToZeroCons(m_numberPhases, m_numberTransports);                   //Mise a zero des cons pour shema spatial sur dU/dt : permet de s affranchir du pas de temps
    }
  }

}

//***********************************************************************

void Run::prepareAdditionalPhysics(double &dt, int &lvl)
{
  //1) Preparation of variables for additional (gradients computations, etc) and communications
  //-------------------------------------------------------------------------------------------
  if (Ncpu > 1) {
    m_stat.startCommunicationTime();
    parallel.communicationsPrimitives(m_eos, lvl);
    m_stat.endCommunicationTime();
  }
  for (unsigned int i = 0; i < m_cellsLvl[lvl].size(); i++) { 
    if (!m_cellsLvl[lvl][i]->getSplit()) 
    {
      printf("cell i = %d, x = %le, y = %le.\n",i,m_cellsLvl[lvl][i]->getPosition().getX(),m_cellsLvl[lvl][i]->getPosition().getY());
      m_cellsLvl[lvl][i]->prepareAddPhys(); 
    } 
  }
  if (Ncpu > 1) {
    m_stat.startCommunicationTime();
    for (unsigned int pa = 0; pa < m_addPhys.size(); pa++) { m_addPhys[pa]->communicationsAddPhys(m_numberPhases, m_dimension, lvl); }
    m_stat.endCommunicationTime();
  }
}

//***********************************************************************

void Run::solveAdditionalPhysics(double &dt, int &lvl)
{
  //1) Preparation of variables for additional (gradients computations, etc) and communications
  //-------------------------------------------------------------------------------------------
  if (Ncpu > 1) {
    m_stat.startCommunicationTime();
    parallel.communicationsPrimitives(m_eos, lvl);
    m_stat.endCommunicationTime();
  }
  for (unsigned int i = 0; i < m_cellsLvl[lvl].size(); i++) { 
    if (!m_cellsLvl[lvl][i]->getSplit()) 
    {
      m_cellsLvl[lvl][i]->prepareAddPhys(); 
    } 
  }
  // exit(EXIT_SUCCESS);
  if (Ncpu > 1) {
    m_stat.startCommunicationTime();
    for (unsigned int pa = 0; pa < m_addPhys.size(); pa++) { m_addPhys[pa]->communicationsAddPhys(m_numberPhases, m_dimension, lvl); }
    m_stat.endCommunicationTime();
  }

  //2) Additional physics fluxes determination (Surface tensions, viscosity, conductivity, ...)
  //-------------------------------------------------------------------------------------------
  //Calcul de la somme des flux des physiques additionnelles que l on stock dans m_cons de chaque cell
  for (unsigned int pa = 0; pa < m_addPhys.size(); pa++) {
    for (unsigned int i = 0; i < m_cellInterfacesLvl[lvl].size(); i++) { if (!m_cellInterfacesLvl[lvl][i]->getSplit()) {
      m_cellInterfacesLvl[lvl][i]->computeFluxAddPhys(m_numberPhases, *m_addPhys[pa]); } }
    for (unsigned int i = 0; i < m_cellsLvl[lvl].size(); i++) { if (!m_cellsLvl[lvl][i]->getSplit()) { m_cellsLvl[lvl][i]->addNonConsAddPhys(m_numberPhases, *m_addPhys[pa], m_symmetry); } }
  }

  //3) Time evolution for additional physics
  //----------------------------------------
  for (unsigned int i = 0; i < m_cellsLvl[lvl].size(); i++) {
    if (!m_cellsLvl[lvl][i]->getSplit()) {
      m_cellsLvl[lvl][i]->timeEvolutionAddPhys(dt, m_numberPhases, m_numberTransports);   //Obtention des cons pour shema sur (Un+1-Un)/dt
      m_cellsLvl[lvl][i]->buildPrim(m_numberPhases);                                      //On peut reconstruire Prim a partir de m_cons
      m_cellsLvl[lvl][i]->setToZeroCons(m_numberPhases, m_numberTransports);              //Mise a zero des cons pour shema spatial sur dU/dt : permet de s affranchir du pas de temps
    }
  }
}

//***********************************************************************

void Run::solveSourceTerms(double &dt, int &lvl)
{
  
  std::shared_ptr<Cantera::Solution> sol = m_sol;

  for (unsigned int i = 0; i < m_cellsLvl[lvl].size(); i++) {
    if (!m_cellsLvl[lvl][i]->getSplit()) {
      for (unsigned int s = 0; s < m_sources.size(); s++) { m_sources[s]->integrateSourceTerms(m_cellsLvl[lvl][i], m_numberPhases, dt, sol); }
      m_cellsLvl[lvl][i]->setToZeroCons(m_numberPhases, m_numberTransports);
    }
  }
}

//***********************************************************************

void Run::solveRelaxations(int &lvl)
{
  for (unsigned int i = 0; i < m_cellsLvl[lvl].size(); i++) { 
    if (!m_cellsLvl[lvl][i]->getSplit()) { 
		m_model->relaxations(m_cellsLvl[lvl][i], m_numberPhases);
    } 
  }
  //Reset of colour function (transports) using volume fraction
  for (unsigned int pa = 0; pa < m_addPhys.size(); pa++) {
    if (m_addPhys[pa]->reinitializationActivated()) {
      m_addPhys[pa]->reinitializeColorFunction(m_cellsLvl, lvl);
      if (Ncpu > 1) {
        m_stat.startCommunicationTime();
        parallel.communicationsTransports(lvl);
        m_stat.endCommunicationTime();
      }
    }
  }
  for (unsigned int i = 0; i < m_cellsLvl[lvl].size(); i++) { if (!m_cellsLvl[lvl][i]->getSplit()) { m_cellsLvl[lvl][i]->prepareAddPhys(); } }
  if (Ncpu > 1) {
    m_stat.startCommunicationTime();
    for (unsigned int pa = 0; pa < m_addPhys.size(); pa++) { m_addPhys[pa]->communicationsAddPhys(m_numberPhases, m_dimension, lvl); }
    m_stat.endCommunicationTime();
  }
  //Optional energy corrections and other relaxations
  for (unsigned int i = 0; i < m_cellsLvl[lvl].size(); i++) {
    if (!m_cellsLvl[lvl][i]->getSplit()) {
      m_cellsLvl[lvl][i]->correctionEnergy(m_numberPhases);               //Correction des energies
      //if (m_evaporation) m_cellsLvl[lvl][i]->relaxPTMu(m_numberPhases); //Relaxation des pressures, temperatures et potentiels chimiques //KS//FP// To implement dynamically
    }
  }
}

//***********************************************************************

void Run::verifyErrors() const
{
  try {
    if (Ncpu > 1) {
      parallel.verifyStateCPUs();
    }
    else if (errors.size() != 0) {
      for (unsigned int e = 0; e < errors.size(); e++) {
        errors[e].afficheError();
      }
      throw ErrorECOGEN("Stop code after error... not managed");
    }
  }
  catch (ErrorECOGEN &) { throw; }
}

//***********************************************************************

void Run::finalize()
{
  //Global desallocations (some are recursives)
  for (int i = 0; i < m_cellInterfacesLvl[0].size(); i++) { delete m_cellInterfacesLvl[0][i]; }
  for (int i = 0; i < m_cellsLvl[0].size(); i++) { delete m_cellsLvl[0][i]; }
  for (int i = 0; i < m_cellsLvlGhost[0].size(); i++) { delete m_cellsLvlGhost[0][i]; }
  
  for (int i = 0; i < m_numberEos; i++) { delete m_eos[i]; } delete[] m_eos;
  //Additional physics desallocations
  for (unsigned int pa = 0; pa < m_addPhys.size(); pa++) { delete m_addPhys[pa]; }
  for (unsigned int s = 0; s < m_sources.size(); s++) { delete m_sources[s]; }
  //Second order desallocations
  if (m_order == "SECONDORDER") {
    for (int k = 0; k < m_numberPhases; k++) { delete slopesPhasesLocal1[k]; }
    for (int k = 0; k < m_numberPhases; k++) { delete slopesPhasesLocal2[k]; }
    delete[] slopesPhasesLocal1;
    delete[] slopesPhasesLocal2;
    delete slopesMixtureLocal1;
    delete slopesMixtureLocal2;
    delete[] slopesTransportLocal1;
    delete[] slopesTransportLocal2;
  }
  //Parallel desaloccations
	m_mesh->finalizeParallele(m_lvlMax);
  //Desallocations others
  delete TB;
  delete cellLeft; delete cellRight;
  delete m_mesh;
  delete m_model;
  delete m_globalLimiter; delete m_interfaceLimiter; delete m_globalVolumeFractionLimiter; delete m_interfaceVolumeFractionLimiter;
  delete m_input; 
  delete m_outPut;
  for (unsigned int s = 0; s < m_cuts.size(); s++) { delete m_cuts[s]; }
  //Desallocations AMR
  delete[] m_cellsLvl;
  delete[] m_cellInterfacesLvl;
  delete[] m_cellsLvlGhost;
}