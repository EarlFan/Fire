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

//! \file      MeshCartesianAMR.cpp
//! \author    K. Schmidmayer, F. Petitpas, B. Dorschner
//! \version   1.1
//! \date      June 5 2019

#include <algorithm>
#include <unordered_map> //For hash
//#include <map> //For ordered map

#include "MeshCartesianAMR.h"

//***********************************************************************

MeshCartesianAMR::MeshCartesianAMR(double lX, int numberCellsX, double lY, int numberCellsY, double lZ, int numberCellsZ,
  std::vector<stretchZone> stretchX, std::vector<stretchZone> stretchY, std::vector<stretchZone> stretchZ, double T_threshhold,
	int lvlMax, int lvlHydro, int lvlChem, int lvlGeo, double criteriaVar, bool varT, bool varYH2O, bool varRho, bool varP, bool varU, bool varAlpha, double xiSplit, double xiJoin) :
  MeshCartesian(lX, numberCellsX, lY, numberCellsY, lZ, numberCellsZ, stretchX, stretchY, stretchZ),
  m_lvlMax(lvlMax), m_lvlHydro(lvlHydro),m_lvlChem(lvlChem), m_lvlGeo(lvlGeo), m_T_threshhold(T_threshhold),m_criteriaVar(criteriaVar), m_varT(varT), m_varYH2O(varYH2O),m_varRho(varRho), m_varP(varP), m_varU(varU), m_varAlpha(varAlpha), m_xiSplit(xiSplit), m_xiJoin(xiJoin)
{
  m_type = AMR;
}

//***********************************************************************

MeshCartesianAMR::~MeshCartesianAMR(){}

//***********************************************************************

int MeshCartesianAMR::initializeGeometrie(TypeMeshContainer<Cell *> &cells, TypeMeshContainer<Cell *> &cellsGhost, TypeMeshContainer<CellInterface *> &cellInterfaces,
  const int &restartSimulation, bool pretraitementParallele, std::string ordreCalcul)
{
  this->meshStretching();
  this->initializeGeometrieAMR(cells, cellsGhost, cellInterfaces, restartSimulation, ordreCalcul);
  return m_geometrie;
}


//***********************************************************************

void MeshCartesianAMR::initializeGeometrieAMR(TypeMeshContainer<Cell *> &cells, TypeMeshContainer<Cell *> &cellsGhost, TypeMeshContainer<CellInterface *> &cellInterfaces, const int &restartSimulation, std::string ordreCalcul)
{
  int ix, iy, iz;
  m_numberCellsX = m_numberCellsXGlobal;
  m_numberCellsY = m_numberCellsYGlobal;
  m_numberCellsZ = m_numberCellsZGlobal;

  //Domain decomposition
  //--------------------
  std::array<int,3> physicalDomainSizes={{m_numberCellsXGlobal,m_numberCellsYGlobal,m_numberCellsZGlobal}};
  if (restartSimulation == 0) { m_decomp = decomposition::Decomposition(physicalDomainSizes); }
  else { m_decomp.updatePhysicalDomainSizes(physicalDomainSizes); }
  auto keys = m_decomp.initialize(Ncpu, rankCpu, restartSimulation);

  for(unsigned int i = 0; i < keys.size(); ++i)
  {
    if (ordreCalcul == "FIRSTORDER") { cells.push_back(new Cell); }
    else { cells.push_back(new CellO2); }
    m_elements.push_back(new ElementCartesian());
    m_elements[i]->setKey(keys[i]);
    cells[i]->setElement(m_elements[i], i);
  }

  //Create cells and elements
  //-------------------------
  this->assignElementProperties(cells, keys);

  //Create cell interfaces, faces and ghost cells
  //---------------------------------------------
  m_numberCellsCalcul = cells.size();
  createCellInterfacesFacesAndGhostCells(cells, cellsGhost, cellInterfaces, ordreCalcul);
  m_numberCellsTotal = cells.size() + cellsGhost.size();
  m_numberFacesTotal = cellInterfaces.size();

  // printf("cells.size() = %d, cellsTotal = %d, facesTotal = %d.\n",
  //   m_numberCellsCalcul,
  //   m_numberCellsTotal,
  //   m_numberFacesTotal);
  // std::cout<<"cpu "<<rankCpu
  //   << " m_numberCellsCalcul "<<m_numberCellsCalcul
  //   << " m_numberCellsTotal "<<m_numberCellsTotal
  //   << " m_numberFacesTotal "<<m_numberFacesTotal
  //   <<std::endl;
}

//***********************************************************************

void MeshCartesianAMR::assignElementProperties(TypeMeshContainer<Cell *> &cells, std::vector<decomposition::Key<3>> &keys)
{
  int ix, iy, iz;
  double volume(0.);
  for(unsigned int i = 0; i < keys.size(); ++i)
  {
    auto coord = keys[i].coordinate();
    ix = coord.x(); iy = coord.y(); iz = coord.z();
    volume = m_dXi[ix] * m_dYj[iy] * m_dZk[iz];
    if(this->getSymType()=="CYLINDRICAL")volume=volume*m_posYj[iy];
    cells[i]->getElement()->setVolume(volume);

    //CFL length
    double lCFL(1.e10);
    if (m_numberCellsX != 1) { lCFL = std::min(lCFL, m_dXi[ix]); }
    if (m_numberCellsY != 1) { lCFL = std::min(lCFL, m_dYj[iy]); }
    if (m_numberCellsZ != 1) { lCFL = std::min(lCFL, m_dZk[iz]); }
    if (m_geometrie > 1) lCFL *= 0.6;

    cells[i]->getElement()->setLCFL(lCFL);
    cells[i]->getElement()->setPos(m_posXi[ix], m_posYj[iy], m_posZk[iz]);
    cells[i]->getElement()->setSize(m_dXi[ix], m_dYj[iy], m_dZk[iz]);
  }
}

//***********************************************************************

void MeshCartesianAMR::createCellInterfacesFacesAndGhostCells(TypeMeshContainer<Cell *> &cells, TypeMeshContainer<Cell *> &cellsGhost,     
TypeMeshContainer<CellInterface*> &cellInterfaces, std::string ordreCalcul)
{
  double posX=0, posY=0., posZ=0.;
  using key_type=decomposition::Key<3>;
  using coordinate_type = key_type::coordinate_type;
  std::array<coordinate_type,6> offsets;
  std::fill(offsets.begin(), offsets.end(), coordinate_type(0));

  std::unordered_map<key_type, Cell*, key_type::hash_functor> cell_map; //For hash
  //std::map<key_type, Cell*> cell_map; //For ordered map
  for (auto c: cells) {
    cell_map.insert(std::make_pair(c->getElement()->getKey(), c));
  }

  for (int d = 0; d < 3; d++)
  {
    offsets[2*d][d] =-1;
    offsets[2*d+1][d] =+1;
  }

  for (unsigned int i = 0; i < cells.size(); ++i)
  {
    const auto coord = cells[i]->getElement()->getKey().coordinate();
    const auto ix = coord.x(), iy = coord.y(), iz = coord.z();
    for (int idx = 0; idx < 2*m_geometrie; idx++)
    {
      const auto offset=offsets[idx];

      posX = m_posXi[ix] + 0.5*m_dXi[ix]*offset[0];
      posY = m_posYj[iy] + 0.5*m_dYj[iy]*offset[1];
      posZ = m_posZk[iz] + 0.5*m_dZk[iz]*offset[2];

      // printf("m_dXi[ix] = %le.\n",m_dXi[ix]);
      // exit(EXIT_FAILURE);

      Coord normal, tangent,binormal;
      normal.setXYZ(static_cast<double>(offset[0]), 
                    static_cast<double>(offset[1]), 
                    static_cast<double>(offset[2])); 

      //Xdir
      if (offset[0] == 1) 
      {
        tangent.setXYZ( 0.,1.,0.); 
        binormal.setXYZ(0.,0.,1.); 
      }
      if (offset[0] == -1)
      {
        tangent.setXYZ( 0.,-1.,0.); 
        binormal.setXYZ(0.,0.,1.); 
      }

      //Ydir
      if (offset[1] == 1)
      {
        tangent.setXYZ( -1.,0.,0.); 
        binormal.setXYZ(0.,0.,1.); 
      }
      if (offset[1] == -1)
      {
        tangent.setXYZ( 1.,0.,0.); 
        binormal.setXYZ(0.,0.,1.); 
      }

      //Zdir
      if (offset[2] == 1) 
      {
        tangent.setXYZ( 1.,0.,0.); 
        binormal.setXYZ(0.,1.,0.); 
      }
      if (offset[2] == -1)
      {
        tangent.setXYZ(-1.,0.,0.); 
        binormal.setXYZ(0.,1.,0.); 
      }

      auto neighborCell = cells[i]->getElement()->getKey().coordinate() + offset;
      if (!m_decomp.is_inside(neighborCell)) //Offset is at a physical boundary
      {
        //Create boundary cell interface
        if (offset[0] == 1) //xDir=N
          m_limXp->creeLimite(cellInterfaces);
        if (offset[0] == -1) //xDir=0
          m_limXm->creeLimite(cellInterfaces);
        if (offset[1] == 1) //yDir=N
          m_limYp->creeLimite(cellInterfaces);
        if (offset[1] == -1) //yDir=0
          m_limYm->creeLimite(cellInterfaces);
        if (offset[2] == 1) //zDir=N
          m_limZp->creeLimite(cellInterfaces);
        if (offset[2] == -1) //zDir=0
          m_limZm->creeLimite(cellInterfaces);

        cellInterfaces.back()->initialize(cells[i], nullptr);

        cells[i]->addCellInterface(cellInterfaces.back());
        m_faces.push_back(new FaceCartesian());
        cellInterfaces.back()->setFace(m_faces.back());

        if (offset[0])//xDir=N
        {
          m_faces.back()->setSize(0.0, m_dYj[iy], m_dZk[iz]);
          if(this->getSymType()=="CYLINDRICAL")
          {
            // the xlim face area shall be revised - fnae
            // xDir = 0 or N? 
            m_faces.back()->initializeAutres(m_dYj[iy] * m_dZk[iz]*posY, normal, tangent, binormal);
          }
          else
            m_faces.back()->initializeAutres(m_dYj[iy] * m_dZk[iz], normal, tangent, binormal);
        }
        if (offset[1])
        {
          m_faces.back()->setSize(m_dXi[ix], 0.0, m_dZk[iz]);
          if(this->getSymType()=="CYLINDRICAL")
            m_faces.back()->initializeAutres(m_dXi[ix] * m_dZk[iz]*posY, normal, tangent, binormal);
          else
            m_faces.back()->initializeAutres(m_dXi[ix] * m_dZk[iz], normal, tangent, binormal);
        }
        if (offset[2])
        {
          m_faces.back()->setSize(m_dXi[ix], m_dYj[iy], 0.0);
          m_faces.back()->initializeAutres(m_dYj[iy] * m_dXi[ix], normal, tangent, binormal);
        }
        m_faces.back()->setPos(posX, posY, posZ);

      }
      else //Offset is an internal cell (ghost or not)
      {
        //Get neighbor key
        auto nKey = cells[i]->getElement()->getKey().neighbor(offset);
        int neighbour = m_decomp.get_rank(nKey);

        if (offset[0]>0 || offset[1]>0 || offset[2]>0) //Positive offset
        {
          //Create cell interface
          if (ordreCalcul == "FIRSTORDER") { cellInterfaces.push_back(new CellInterface); }
          else { cellInterfaces.push_back(new CellInterfaceO2); }     

          m_faces.push_back(new FaceCartesian());
          cellInterfaces.back()->setFace(m_faces.back());

          if (offset[0])
          {
            m_faces.back()->setSize(0.0, m_dYj[iy], m_dZk[iz]);
            if(this->getSymType()=="CYLINDRICAL")
              m_faces.back()->initializeAutres(m_dYj[iy] * m_dZk[iz]*posY, normal, tangent, binormal);
            else
              m_faces.back()->initializeAutres(m_dYj[iy] * m_dZk[iz], normal, tangent, binormal);
          }
          if (offset[1])
          {
            m_faces.back()->setSize(m_dXi[ix], 0.0, m_dZk[iz]);
            if(this->getSymType()=="CYLINDRICAL")
              m_faces.back()->initializeAutres(m_dXi[ix] * m_dZk[iz]*posY, normal, tangent, binormal);
            else
              m_faces.back()->initializeAutres(m_dXi[ix] * m_dZk[iz], normal, tangent, binormal);
          }
          if (offset[2])
          {
            m_faces.back()->setSize(m_dXi[ix], m_dYj[iy], 0.0);
            m_faces.back()->initializeAutres(m_dYj[iy] * m_dXi[ix], normal, tangent, binormal);
          }
          m_faces.back()->setPos(posX, posY, posZ);

          //Try to find the neighbor cell into the non-ghost cells
          auto it_pair = cell_map.find(nKey);

          if (it_pair != cell_map.end()) //Neighbor cell is a non-ghost cell
          {

            auto it = it_pair->second;
            //Update cell interface
            cellInterfaces.back()->initialize(cells[i], it);
            cells[i]->addCellInterface(cellInterfaces.back());
            it->addCellInterface(cellInterfaces.back());
          }
          else //Neighbor cell is a ghost cell
          {
            //Try to find the neighbor cell in the already created ghost cells
            auto it2 = std::find_if(cellsGhost.begin(), cellsGhost.end(),
                    [&nKey](Cell* _k0){ 
                    return _k0->getElement()->getKey() == nKey;
                    });

            if (it2 == cellsGhost.end()) //Ghost cell does not exist
            {
              //Create ghost cell and update cell interface
              if (ordreCalcul == "FIRSTORDER") { cellsGhost.push_back(new CellGhost); }
              else { cellsGhost.push_back(new CellO2Ghost); }
              m_elements.push_back(new ElementCartesian());
              m_elements.back()->setKey(nKey);
              cellsGhost.back()->setElement(m_elements.back(), cellsGhost.size()-1);
              cellsGhost.back()->pushBackSlope();
              parallel.addSlopesToSend(neighbour);
              parallel.addSlopesToReceive(neighbour);

              //Update parallel communications
              parallel.setNeighbour(neighbour);

              //Try to find the current cell into the already added cells to send
              auto cKey = cells[i]->getElement()->getKey();
              auto it3 = std::find_if(parallel.getElementsToSend(neighbour).begin(), parallel.getElementsToSend(neighbour).end(),
                      [&cKey](Cell* _k0){ 
                      return _k0->getElement()->getKey() == cKey;
                      });

              if (it3 == parallel.getElementsToSend(neighbour).end()) //Current cell not added in send vector
              {
                // add the cell into vector that send to neighbour
                parallel.addElementToSend(neighbour, cells[i]);
              }
              parallel.addElementToReceive(neighbour, cellsGhost.back());
              cellsGhost.back()->setRankOfNeighborCPU(neighbour);// -  stop here

              const auto coord = nKey.coordinate();
              const auto nix = coord.x(), niy = coord.y(), niz = coord.z();

              double volume = m_dXi[nix] * m_dYj[niy] * m_dZk[niz];
              if(this->getSymType()=="CYLINDRICAL")volume=volume*m_posYj[niy];
              cellsGhost.back()->getElement()->setVolume(volume);

              double lCFL(1.e10);
              if (m_numberCellsX != 1) { lCFL = std::min(lCFL, m_dXi[nix]); }
              if (m_numberCellsY != 1) { lCFL = std::min(lCFL, m_dYj[niy]); }
              if (m_numberCellsZ != 1) { lCFL = std::min(lCFL, m_dZk[niz]); }
              if (m_geometrie > 1) lCFL *= 0.6;

              cellsGhost.back()->getElement()->setLCFL(lCFL);
              cellsGhost.back()->getElement()->setPos(m_posXi[nix], m_posYj[niy], m_posZk[niz]);
              cellsGhost.back()->getElement()->setSize(m_dXi[nix], m_dYj[niy], m_dZk[niz]);

              //Update pointers cells <-> cell interfaces
              cellInterfaces.back()->initialize(cells[i], cellsGhost.back());
              cells[i]->addCellInterface(cellInterfaces.back());
              cellsGhost.back()->addCellInterface(cellInterfaces.back());
            }
            else { //Ghost cell exists
              //Update parallel communications
              //Try to find the current cell into the already added cells to send
              auto cKey = cells[i]->getElement()->getKey();
              auto it3 = std::find_if(parallel.getElementsToSend(neighbour).begin(), parallel.getElementsToSend(neighbour).end(),
                      [&cKey](Cell* _k0){ 
                      return _k0->getElement()->getKey() == cKey;
                      });
              if (it3 == parallel.getElementsToSend(neighbour).end()) //Current cell not added in send vector
              {
                parallel.addElementToSend(neighbour, cells[i]);
              }

              //Update pointers cells <-> cell interfaces
              cellInterfaces.back()->initialize(cells[i], *it2);
              cells[i]->addCellInterface(cellInterfaces.back());
              (*it2)->addCellInterface(cellInterfaces.back());
              (*it2)->pushBackSlope();
              parallel.addSlopesToSend(neighbour);
              parallel.addSlopesToReceive(neighbour);
            }
          }
        }
        else //Negative offset
        {
          //Try to find the neighbor cell into the non-ghost cells
          auto it_pair = cell_map.find(nKey);

          if (it_pair == cell_map.end()) //Neighbor cell is a ghost cell
          {
            //Create cell interface related to the ghost cell
            if (ordreCalcul == "FIRSTORDER") { cellInterfaces.push_back(new CellInterface); }
            else { cellInterfaces.push_back(new CellInterfaceO2); }     

            m_faces.push_back(new FaceCartesian());
            cellInterfaces.back()->setFace(m_faces.back());

            if (offset[0])
            {
              normal.setXYZ( 1.,0.,0.); 
              tangent.setXYZ( 0.,1.,0.); 
              binormal.setXYZ(0.,0.,1.); 
              m_faces.back()->setSize(0.0, m_dYj[iy], m_dZk[iz]);
              if(this->getSymType()=="CYLINDRICAL")
                m_faces.back()->initializeAutres(m_dYj[iy] * m_dZk[iz]*posY, normal, tangent, binormal);
              else 
                m_faces.back()->initializeAutres(m_dYj[iy] * m_dZk[iz], normal, tangent, binormal);
            }
            if (offset[1])
            {
              normal.setXYZ( 0.,1.,0.); 
              tangent.setXYZ( -1.,0.,0.); 
              binormal.setXYZ(0.,0.,1.); 
              m_faces.back()->setSize(m_dXi[ix], 0.0, m_dZk[iz]);
              if(this->getSymType()=="CYLINDRICAL")
                m_faces.back()->initializeAutres(m_dXi[ix] * m_dZk[iz]*posY, normal, tangent, binormal);
              else
                m_faces.back()->initializeAutres(m_dXi[ix] * m_dZk[iz], normal, tangent, binormal);
            }
            if (offset[2])
            {
              normal.setXYZ( 0.,0.,1.); 
              tangent.setXYZ( 1.,0.,0.); 
              binormal.setXYZ(0.,1.,0.); 
              m_faces.back()->setSize(m_dXi[ix], m_dYj[iy], 0.0);
              m_faces.back()->initializeAutres(m_dYj[iy] * m_dXi[ix], normal, tangent, binormal);
            }
            m_faces.back()->setPos(posX, posY, posZ);

            //Try to find the neighbor cell into the already created ghost cells
            auto it2 = std::find_if(cellsGhost.begin(), cellsGhost.end(),
                    [&nKey](Cell* _k0){ 
                    return _k0->getElement()->getKey() == nKey;
                     });

            if (it2 == cellsGhost.end()) //Ghost cell does not exist
            {
              //Create ghost cell
              if (ordreCalcul == "FIRSTORDER") { cellsGhost.push_back(new CellGhost); }
              else { cellsGhost.push_back(new CellO2Ghost); }
              m_elements.push_back(new ElementCartesian());
              m_elements.back()->setKey(nKey);
              cellsGhost.back()->setElement(m_elements.back(), cellsGhost.size()-1);
              cellsGhost.back()->pushBackSlope();
              parallel.addSlopesToSend(neighbour);
              parallel.addSlopesToReceive(neighbour);

              //Update parallel communications
              parallel.setNeighbour(neighbour);
              //Try to find the current cell into the already added cells to send
              auto cKey = cells[i]->getElement()->getKey();
              auto it3 = std::find_if(parallel.getElementsToSend(neighbour).begin(), parallel.getElementsToSend(neighbour).end(),
                      [&cKey](Cell* _k0){ 
                      return _k0->getElement()->getKey() == cKey;
                      });
              if (it3 == parallel.getElementsToSend(neighbour).end()) //Current cell not added in send vector
              {
                parallel.addElementToSend(neighbour, cells[i]);
              }
              parallel.addElementToReceive(neighbour, cellsGhost.back());
              cellsGhost.back()->setRankOfNeighborCPU(neighbour);

              const auto coord = nKey.coordinate();
              const auto nix = coord.x(), niy = coord.y(), niz = coord.z();

              double volume = m_dXi[nix] * m_dYj[niy] * m_dZk[niz];
              if(this->getSymType()=="CYLINDRICAL")volume=volume*m_posYj[niy];
              cellsGhost.back()->getElement()->setVolume(volume);

              double lCFL(1.e10);
              if (m_numberCellsX != 1) { lCFL = std::min(lCFL, m_dXi[nix]); }
              if (m_numberCellsY != 1) { lCFL = std::min(lCFL, m_dYj[niy]); }
              if (m_numberCellsZ != 1) { lCFL = std::min(lCFL, m_dZk[niz]); }
              if (m_geometrie > 1) lCFL *= 0.6;

              cellsGhost.back()->getElement()->setLCFL(lCFL);
              cellsGhost.back()->getElement()->setPos(m_posXi[nix], m_posYj[niy], m_posZk[niz]);
              cellsGhost.back()->getElement()->setSize(m_dXi[nix], m_dYj[niy], m_dZk[niz]);

              //Update pointers cells <-> cell interfaces
              cellInterfaces.back()->initialize(cellsGhost.back(), cells[i]);
              cells[i]->addCellInterface(cellInterfaces.back());
              cellsGhost.back()->addCellInterface(cellInterfaces.back());
            }
            else //Ghost cell exists
            {
              //Update parallel communications
              //Try to find the current cell into the already added cells to send
              auto cKey = cells[i]->getElement()->getKey();
              auto it3 = std::find_if(parallel.getElementsToSend(neighbour).begin(), parallel.getElementsToSend(neighbour).end(),
                      [&cKey](Cell* _k0){ 
                      return _k0->getElement()->getKey() == cKey;
                      });
              if (it3 == parallel.getElementsToSend(neighbour).end()) //Current cell not added in send vector
              {
                parallel.addElementToSend(neighbour, cells[i]);
              }

              //Update pointers cells <-> cell interfaces
              cellInterfaces.back()->initialize(*it2, cells[i]);
              cells[i]->addCellInterface(cellInterfaces.back());
              (*it2)->addCellInterface(cellInterfaces.back());
              (*it2)->pushBackSlope();
              parallel.addSlopesToSend(neighbour);
              parallel.addSlopesToReceive(neighbour);
            }
          }
        } //Negative offset
      } //Offset is an internal cell (ghost or not)
    } //Offsets
  } //Internal, non-ghost cells

  if (Ncpu > 1) {
    for (int i=0;i<Ncpu;++i)
    {
      // std::sort(parallel.getElementsToSend(i).begin(),parallel.getElementsToSend(i).end(),[&](Cell* child0, Cell* child1)
      // {
      //   return child0->getElement()->getKey() < child1->getElement()->getKey();
      // });
      std::sort(parallel.getElementsToReceive(i).begin(),parallel.getElementsToReceive(i).end(),[&](Cell* child0, Cell* child1)
      {
        return child0->getElement()->getKey() < child1->getElement()->getKey();
      });
    }
  }
}

//***********************************************************************

void MeshCartesianAMR::procedureRaffinementInitialization(TypeMeshContainer<Cell *> *cellsLvl, TypeMeshContainer<Cell *> *cellsLvlGhost, 
  TypeMeshContainer<CellInterface *> *cellInterfacesLvl,
  const std::vector<AddPhys*> &addPhys, Model *model, int &nbCellsTotalAMR, std::vector<GeometricalDomain*> &domains,
  Eos **eos, const int &restartSimulation, std::string ordreCalcul, const int &numberPhases, const int &numberTransports)
{
  nbCellsTotalAMR = m_numberCellsCalcul;

  if (restartSimulation == 0) { //Only for simulation from input files
    for (int iterInit = 0; iterInit < 2; iterInit++) {
      for (int lvl = 0; lvl < m_lvlMax; lvl++) {
        if (Ncpu > 1) { parallel.communicationsPrimitives(eos, lvl); }
        this->procedureRaffinement(cellsLvl, cellsLvlGhost, cellInterfacesLvl, lvl, addPhys, model, nbCellsTotalAMR, eos);
        for (unsigned int i = 0; i < cellsLvl[lvl + 1].size(); i++) {
          cellsLvl[lvl + 1][i]->fill(domains, m_lvlMax);
        }
        for (unsigned int i = 0; i < cellsLvl[lvl + 1].size(); i++) {
          cellsLvl[lvl + 1][i]->completeFulfillState();
        }
        for (unsigned int i = 0; i < cellsLvl[lvl].size(); i++) {
          cellsLvl[lvl][i]->averageChildrenInParent();
        }
      }
    }
    for (int lvl = 0; lvl <= m_lvlMax; lvl++) {
      
      if (Ncpu > 1) { 
        parallel.communicationsPrimitives(eos, lvl); 
        parallel.communicationsPrimitives(eos, lvl,vecPhasesO2);
      }
      for (unsigned int i = 0; i < cellsLvl[lvl].size(); i++) {
        if (!cellsLvl[lvl][i]->getSplit()) { cellsLvl[lvl][i]->completeFulfillState(); }
      }
    }
    if (Ncpu > 1) { this->parallelLoadBalancingAMR(cellsLvl, cellsLvlGhost, cellInterfacesLvl, 
      ordreCalcul, numberPhases, numberTransports, addPhys, model, eos, nbCellsTotalAMR, true); }
  }
}

//***********************************************************************

void MeshCartesianAMR::procedureRaffinement(TypeMeshContainer<Cell *> *cellsLvl, TypeMeshContainer<Cell *> *cellsLvlGhost, TypeMeshContainer<CellInterface *> *cellInterfacesLvl, const int &lvl,
  const std::vector<AddPhys*> &addPhys, Model *model, int &nbCellsTotalAMR, Eos **eos)
{
  //1) Calcul de Xi dans chaque cell de niveau lvl
  //-------------------------------------------------
  for (unsigned int i = 0; i < cellsLvl[lvl].size(); i++) { cellsLvl[lvl][i]->setToZeroXi(); }
  for (unsigned int i = 0; i < cellInterfacesLvl[lvl].size(); i++) 
  { 
    cellInterfacesLvl[lvl][i]->computeXi(lvl, m_lvlHydro,m_lvlChem, m_lvlGeo, m_T_threshhold, m_criteriaVar, m_varRho, m_varP, m_varU, m_varT, m_varYH2O, m_varAlpha); 
  }
  // bool varP2 = m_varP;
  // if (lvl >= 5) { varP2 = false; }
  // for (unsigned int i = 0; i < cellInterfacesLvl[lvl].size(); i++) { cellInterfacesLvl[lvl][i]->computeXi(m_criteriaVar, m_varRho, varP2, m_varU, m_varAlpha); }
  // bool varU2 = m_varU;
  // if (lvl >= 2) { varU2 = false; }
  // for (unsigned int i = 0; i < cellInterfacesLvl[lvl].size(); i++) { cellInterfacesLvl[lvl][i]->computeXi(m_criteriaVar, m_varRho, m_varP, varU2, m_varAlpha); }
  // for (unsigned int i = 0; i < cellsLvl[lvl].size(); i++) {
  //  double x(0.), y(0.), z(0.);
  //  x = cellsLvl[lvl][i]->getPosition().getX();
  //  y = cellsLvl[lvl][i]->getPosition().getY();
  //  //z = cellsLvl[lvl][i]->getPosition().getZ();
  //  //if (std::pow((x*x + y*y + z*z), 0.5) > 500.e-6) {
  //  //if (std::pow((x*x + y*y), 0.5) > 6.e-4) {
  //  //if ((x > 250e-6) || (y > 200.e-6)) {
  //  //if (x > 15.) {
  //  if (std::pow((x*x + y * y), 0.5) > 5.) {
  //      cellsLvl[lvl][i]->setToZeroXi();
  //  }
  // }
  if (Ncpu > 1) { parallel.communicationsXi( lvl); }
  
  //2) Smoothing de Xi
  //------------------
  for (int iterDiff = 0; iterDiff < AMRPara::amr_iter; iterDiff++) { //Arbitrary number of iterations
		//Mise a zero cons xi
    for (unsigned int i = 0; i < cellsLvl[lvl].size(); i++) { cellsLvl[lvl][i]->setToZeroConsXi(); }

    //Calcul des "flux"
    // for (unsigned int i = 0; i < cellInterfacesLvl[lvl].size(); i++) { 
    //   if (!cellInterfacesLvl[lvl][i]->getSplit()) {
    //     cellInterfacesLvl[lvl][i]->computeFluxXi();
    //   }
    // }

    for (unsigned int i = 0; i < cellInterfacesLvl[lvl].size(); i++) { cellInterfacesLvl[lvl][i]->computeFluxXi(); }

    //Evolution temporelle
    for (unsigned int i = 0; i < cellsLvl[lvl].size(); i++) { cellsLvl[lvl][i]->timeEvolutionXi(); }
		if (Ncpu > 1) { parallel.communicationsXi( lvl); }
  }

	if (lvl < m_lvlMax) {
    int lvlPlus1 = lvl + 1;
    //3) Raffinement des cells et cell interfaces
    //-------------------------------------------
    for (unsigned int i = 0; i < cellsLvl[lvl].size(); i++) { cellsLvl[lvl][i]->chooseRefine(m_xiSplit, m_numberCellsY, m_numberCellsZ, addPhys, model, nbCellsTotalAMR); }

    //4) Deraffinement des cells et cell interfaces
    //---------------------------------------------
    for (unsigned int i = 0; i < cellsLvl[lvl].size(); i++) { cellsLvl[lvl][i]->chooseUnrefine(m_xiJoin, nbCellsTotalAMR, m_lvlMax); }

    if (Ncpu > 1) {
      //5) Raffinement et deraffinement des cells fantomes
      //-----------------------------------------------------
      //Communication split + Raffinement et deraffinement des cells fantomes + Reconstruction du tableau de cells fantomes de niveau lvl + 1
      // Split communication + Refinement and derefinement of phantom cells + Reconstruction of the array of phantom cells at level lvl + 1
      parallel.communicationsSplit(lvl);
      cellsLvlGhost[lvlPlus1].clear();
      // fane
      for (unsigned int i = 0; i < cellsLvlGhost[lvl].size(); i++) { cellsLvlGhost[lvl][i]->chooseRefineDeraffineGhost(m_numberCellsY, m_numberCellsZ, addPhys, model, cellsLvlGhost); }
      //Communications primitives pour mettre a jour les cells deraffinees
      // Primitive communications for updating derafined cells
      parallel.communicationsPrimitives(eos, lvl);

      //6) Mise a jour des communications persistantes au niveau lvl + 1
      // 6) Update of persistent communications at lvl + 1 level
      //----------------------------------------------------------------
      parallel.communicationsNumberGhostCells(lvlPlus1);	
      //Communication des numbers d'elements a envoyer et a recevoir de chaque cote de la limite parallele
      // Communication of the numbers of elements to be sent and received on each side of the parallel limit
      parallel.updatePersistentCommunicationsLvlAMR(lvlPlus1, m_geometrie);
    }

    //7) Reconstruction des tableaux de cells et cell interfaces lvl + 1
    //7) Reconstruction of arrays of cells and cell interfaces lvl + 1
    //------------------------------------------------------------------
    cellsLvl[lvlPlus1].clear();
    cellInterfacesLvl[lvlPlus1].clear();
    for (unsigned int i = 0; i < cellsLvl[lvl].size(); i++) { cellsLvl[lvl][i]->buildLvlCellsAndLvlInternalCellInterfacesArrays(cellsLvl, cellInterfacesLvl); }
    for (unsigned int i = 0; i < cellInterfacesLvl[lvl].size(); i++) { cellInterfacesLvl[lvl][i]->constructionTableauCellInterfacesExternesLvl(cellInterfacesLvl); }
  }
}

//***********************************************************************

std::string MeshCartesianAMR::whoAmI() const
{
  return "CARTESIAN_AMR";
}

//**************************************************************************
//******************************** PRINTING ********************************
//**************************************************************************

void MeshCartesianAMR::ecritHeaderPiece(std::ofstream &fileStream, TypeMeshContainer<Cell *> *cellsLvl) const
{
  int numberCells = 0, numberPointsParMaille = 4;
  for (int lvl = 0; lvl <= m_lvlMax; lvl++) {
    for (unsigned int i = 0; i < cellsLvl[lvl].size(); i++) {
      if (!cellsLvl[lvl][i]->getSplit()) {
        numberCells += 1;
      }
    }
  }
  if (m_numberCellsZ > 1) { numberPointsParMaille = 8; }

  fileStream << "    <Piece NumberOfPoints=\"" << numberPointsParMaille*numberCells << "\" NumberOfCells=\"" << numberCells << "\">" << std::endl;
}

//***********************************************************************
// recupereNoeuds: recover knots, used for output.
void MeshCartesianAMR::recupereNoeuds(std::vector<double> &jeuDonnees, std::vector<Cell *> *cellsLvl) const
{
  int dimZ = 0;
  if (m_numberCellsZ > 1) dimZ = 1;

  double dXsur2(0.), dYsur2(0.), dZsur2(0.);
  for (int lvl = 0; lvl <= m_lvlMax; lvl++) {
    dXsur2 = 0.; dYsur2 = 0.; dZsur2 = 0.;
    for (unsigned int i = 0; i < cellsLvl[lvl].size(); i++) {
      if (!cellsLvl[lvl][i]->getSplit()) {
        dXsur2 = 0.5*cellsLvl[lvl][i]->getSizeX();
        dYsur2 = 0.5*cellsLvl[lvl][i]->getSizeY();
        dZsur2 = 0.5*cellsLvl[lvl][i]->getSizeZ();
        //Point 0
        jeuDonnees.push_back(cellsLvl[lvl][i]->getPosition().getX() - dXsur2);
        jeuDonnees.push_back(cellsLvl[lvl][i]->getPosition().getY() - dYsur2);
        jeuDonnees.push_back(cellsLvl[lvl][i]->getPosition().getZ() - dZsur2*dimZ);
        //Point 1
        jeuDonnees.push_back(cellsLvl[lvl][i]->getPosition().getX() + dXsur2);
        jeuDonnees.push_back(cellsLvl[lvl][i]->getPosition().getY() - dYsur2);
        jeuDonnees.push_back(cellsLvl[lvl][i]->getPosition().getZ() - dZsur2*dimZ);
        //Point 2
        jeuDonnees.push_back(cellsLvl[lvl][i]->getPosition().getX() + dXsur2);
        jeuDonnees.push_back(cellsLvl[lvl][i]->getPosition().getY() + dYsur2);
        jeuDonnees.push_back(cellsLvl[lvl][i]->getPosition().getZ() - dZsur2*dimZ);
        //Point 3
        jeuDonnees.push_back(cellsLvl[lvl][i]->getPosition().getX() - dXsur2);
        jeuDonnees.push_back(cellsLvl[lvl][i]->getPosition().getY() + dYsur2);
        jeuDonnees.push_back(cellsLvl[lvl][i]->getPosition().getZ() - dZsur2*dimZ);

        if (dimZ > 0.99) {
          //Point 4
          jeuDonnees.push_back(cellsLvl[lvl][i]->getPosition().getX() - dXsur2);
          jeuDonnees.push_back(cellsLvl[lvl][i]->getPosition().getY() - dYsur2);
          jeuDonnees.push_back(cellsLvl[lvl][i]->getPosition().getZ() + dZsur2);
          //Point 5
          jeuDonnees.push_back(cellsLvl[lvl][i]->getPosition().getX() + dXsur2);
          jeuDonnees.push_back(cellsLvl[lvl][i]->getPosition().getY() - dYsur2);
          jeuDonnees.push_back(cellsLvl[lvl][i]->getPosition().getZ() + dZsur2);
          //Point 6
          jeuDonnees.push_back(cellsLvl[lvl][i]->getPosition().getX() + dXsur2);
          jeuDonnees.push_back(cellsLvl[lvl][i]->getPosition().getY() + dYsur2);
          jeuDonnees.push_back(cellsLvl[lvl][i]->getPosition().getZ() + dZsur2);
          //Point 7
          jeuDonnees.push_back(cellsLvl[lvl][i]->getPosition().getX() - dXsur2);
          jeuDonnees.push_back(cellsLvl[lvl][i]->getPosition().getY() + dYsur2);
          jeuDonnees.push_back(cellsLvl[lvl][i]->getPosition().getZ() + dZsur2);
        }
      } //Fin cell non split
    } //Fin Cells
  } //Fin Levels
}

//***********************************************************************

void MeshCartesianAMR::recupereConnectivite(std::vector<double> &jeuDonnees, std::vector<Cell *> *cellsLvl) const
{
  int dimZ(0);
  int numberPointsParMaille(4);
  if (m_numberCellsZ > 1) { dimZ = 1; numberPointsParMaille = 8; }

  if (dimZ < 0.99) {
    int numCell(0);
    for (int lvl = 0; lvl <= m_lvlMax; lvl++) {
      for (unsigned int i = 0; i < cellsLvl[lvl].size(); i++) {
        if (!cellsLvl[lvl][i]->getSplit()) {
          jeuDonnees.push_back(numCell*numberPointsParMaille);
          jeuDonnees.push_back(numCell*numberPointsParMaille+1);
          jeuDonnees.push_back(numCell*numberPointsParMaille+2);
          jeuDonnees.push_back(numCell*numberPointsParMaille+3);
          numCell++;
        }
      }
    }
  }
  else {
    int numCell(0);
    for (int lvl = 0; lvl <= m_lvlMax; lvl++) {
      for (unsigned int i = 0; i < cellsLvl[lvl].size(); i++) {
        if (!cellsLvl[lvl][i]->getSplit()) {
          jeuDonnees.push_back(numCell*numberPointsParMaille);
          jeuDonnees.push_back(numCell*numberPointsParMaille + 1);
          jeuDonnees.push_back(numCell*numberPointsParMaille + 2);
          jeuDonnees.push_back(numCell*numberPointsParMaille + 3);
          jeuDonnees.push_back(numCell*numberPointsParMaille + 4);
          jeuDonnees.push_back(numCell*numberPointsParMaille + 5);
          jeuDonnees.push_back(numCell*numberPointsParMaille + 6);
          jeuDonnees.push_back(numCell*numberPointsParMaille + 7);
          numCell++;
        }
      }
    }
  }
}

//***********************************************************************

void MeshCartesianAMR::recupereOffsets(std::vector<double> &jeuDonnees, std::vector<Cell *> *cellsLvl) const
{
  int numberPointsParMaille(4);
  if (m_numberCellsZ > 1) { numberPointsParMaille = 8; }
  int numCell(0);
  for (int lvl = 0; lvl <= m_lvlMax; lvl++) {
    for (unsigned int i = 0; i < cellsLvl[lvl].size(); i++) {
      if (!cellsLvl[lvl][i]->getSplit()) {
        jeuDonnees.push_back((numCell + 1)*numberPointsParMaille);
        numCell++;
      }
    }
  }
}

//****************************************************************************

void MeshCartesianAMR::recupereTypeCell(std::vector<double> &jeuDonnees, std::vector<Cell *> *cellsLvl) const
{
  int type(9);
  if (m_numberCellsZ > 1) { type = 12; }
  int numCell(0);
  for (int lvl = 0; lvl <= m_lvlMax; lvl++) {
    for (unsigned int i = 0; i < cellsLvl[lvl].size(); i++) {
      if (!cellsLvl[lvl][i]->getSplit()) {
        jeuDonnees.push_back(type);
        numCell++;
      }
    }
  }
}

//***********************************************************************

void MeshCartesianAMR::recupereDonnees(TypeMeshContainer<Cell *> *cellsLvl, std::vector<double> &jeuDonnees, const int var, int phase) const
{
  jeuDonnees.clear();
  double transport(0.);
  for (int lvl = 0; lvl <= m_lvlMax; lvl++) {
    for (unsigned int i = 0; i < cellsLvl[lvl].size(); i++) {
      if (!cellsLvl[lvl][i]->getSplit()) {
        if (var > 0) { //On veut recuperer les donnees scalars
          if (phase >= 0) { jeuDonnees.push_back(cellsLvl[lvl][i]->getPhase(phase)->returnScalar(var)); }      //Donnees de phases
          else if (phase == -1) { jeuDonnees.push_back(cellsLvl[lvl][i]->getMixture()->returnScalar(var)); }   //Donnees de mixture
          else if (phase == -2) {
            transport = cellsLvl[lvl][i]->getTransport(var - 1).getValue();
            if (transport < 1.e-20) { transport = 0.; }
            jeuDonnees.push_back(transport);
          }
          else if (phase == -3) { jeuDonnees.push_back(cellsLvl[lvl][i]->getXi()); }
          else if (phase == -4) { jeuDonnees.push_back(cellsLvl[lvl][i]->getDensityGradient()); }
          else if (phase == -5) { jeuDonnees.push_back(cellsLvl[lvl][i]->getMixture()->getInductionTime()); }
          else if (phase == -6) { jeuDonnees.push_back(cellsLvl[lvl][i]->getMixture()->getMach()); }
          else if (phase == -7) { jeuDonnees.push_back(cellsLvl[lvl][i]->getMixture()->getDensityChangeRateArray()[5]); }//density change rate of H2O
          else if (phase == -8) { jeuDonnees.push_back(cellsLvl[lvl][i]->getMixture()->getHRR()); }//hrr
          else if (phase == -9) { jeuDonnees.push_back(cellsLvl[lvl][i]->getMixture()->getDensityChangeRateArray()[1]); }//density change rate of H2
          else if (phase == -10) { jeuDonnees.push_back(cellsLvl[lvl][i]->getMixture()->getDensityChangeRateArray()[2]); }//density change rate of H
          else if (phase == -11) { jeuDonnees.push_back(cellsLvl[lvl][i]->getMixture()->getDensityChangeRateArray()[3]); }//density change rate of OH
          else if (phase == -12) { jeuDonnees.push_back(cellsLvl[lvl][i]->getMixture()->getDensityChangeRateArray()[4]); }//density change rate of O
          else if (phase == -13) { jeuDonnees.push_back(cellsLvl[lvl][i]->getMixture()->getDensityChangeRateArray()[6]); }//density change rate of HO2
          else if (phase == -14) { jeuDonnees.push_back(cellsLvl[lvl][i]->getMixture()->getDensityChangeRateArray()[7]); }//density change rate of H2O2
          else if (phase == -15) { jeuDonnees.push_back(cellsLvl[lvl][i]->getMixture()->getDensityChangeRateArray()[9]); }//density change rate of O2
          else if (phase == -16) { jeuDonnees.push_back(cellsLvl[lvl][i]->getMixture()->getRhoExcitedOH()); }//numerical distribution of OH*

          // else if (phase == -5) { jeuDonnees.push_back(static_cast<double>(rankCpu)); }
          // else if (phase == -6) { jeuDonnees.push_back(static_cast<double>(cellsLvl[lvl][i]->getElement()->getKey().getIndex())); }
          // else if (phase == -7) { jeuDonnees.push_back(cellsLvl[lvl][i]->getGradientYI(0)); }// get the gradient of Y_SF6
          else { Errors::errorMessage("MeshCartesianAMR::recupereDonnees: unknown number of phase: ", phase); }
        }
        else { //On veut recuperer les donnees vectorielles
          if (phase >= 0) { //Donnees de phases
            jeuDonnees.push_back(cellsLvl[lvl][i]->getPhase(phase)->returnVector(-var).getX());
            jeuDonnees.push_back(cellsLvl[lvl][i]->getPhase(phase)->returnVector(-var).getY());
            jeuDonnees.push_back(cellsLvl[lvl][i]->getPhase(phase)->returnVector(-var).getZ());
          }
          else if (phase == -1){  // Velocity
            jeuDonnees.push_back(cellsLvl[lvl][i]->getMixture()->returnVector(-var).getX());
            jeuDonnees.push_back(cellsLvl[lvl][i]->getMixture()->returnVector(-var).getY());
            jeuDonnees.push_back(cellsLvl[lvl][i]->getMixture()->returnVector(-var).getZ());
          }
          else if (phase == -2) { //gradient U, these gradients shall be stored to avoid repeated computation. - fnae, 230404
            Coord temp = cellsLvl[lvl][i]->getGradientU();
            jeuDonnees.push_back(temp.getX());
            jeuDonnees.push_back(temp.getY());
            jeuDonnees.push_back(temp.getZ());
          }
          else if (phase == -3) { //gradient V
            Coord temp = cellsLvl[lvl][i]->getGradientV();
            jeuDonnees.push_back(temp.getX());
            jeuDonnees.push_back(temp.getY());
            jeuDonnees.push_back(temp.getZ());
          }
          else if (phase == -4) { //gradient Y[0th]
            Coord temp = cellsLvl[lvl][i]->getGradientYI(0);
            jeuDonnees.push_back(temp.getX());
            jeuDonnees.push_back(temp.getY());
            jeuDonnees.push_back(temp.getZ());
          }
          else if (phase == -5) { //gradient Y[3rd]
            Coord temp = cellsLvl[lvl][i]->getGradientYI(3);
            jeuDonnees.push_back(temp.getX());
            jeuDonnees.push_back(temp.getY());
            jeuDonnees.push_back(temp.getZ());
          }
          else if (phase == -6) { //gradient of induction time
            Coord temp = cellsLvl[lvl][i]->getGradientInductionTime();
            jeuDonnees.push_back(temp.getX());
            jeuDonnees.push_back(temp.getY());
            jeuDonnees.push_back(temp.getZ());
          }
          else if (phase == -7) { //gradient of Pressure
            Coord temp = cellsLvl[lvl][i]->getGradientP();
            jeuDonnees.push_back(temp.getX());
            jeuDonnees.push_back(temp.getY());
            jeuDonnees.push_back(temp.getZ());
          }
          else { Errors::errorMessage("MeshCartesianAMR::recupereDonnees: unknown number of phase: ", phase); }
        } //Fin vecteur
      } //Fin split
    } //fin lvl
  } //fin levels
}

//****************************************************************************

void MeshCartesianAMR::setDataSet(std::vector<double> &jeuDonnees, TypeMeshContainer<Cell *> *cellsLvl, const int var, int phase) const
{
  int iterDataSet(0);
  Coord vec;
  for (int lvl = 0; lvl <= m_lvlMax; lvl++) {
    for (unsigned int i = 0; i < cellsLvl[lvl].size(); i++) {
      if (!cellsLvl[lvl][i]->getSplit()) {
        if (var > 0) { //Scalars data are first set
          if (phase >= 0) { cellsLvl[lvl][i]->getPhase(phase)->setScalar(var, jeuDonnees[iterDataSet++]); } //phases data
          else if (phase == -1) { cellsLvl[lvl][i]->getMixture()->setScalar(var, jeuDonnees[iterDataSet++]); }  //mixture data
          else if (phase == -2) { cellsLvl[lvl][i]->getTransport(var - 1).setValue(jeuDonnees[iterDataSet++]); } //transport data
          else if (phase == -3) { cellsLvl[lvl][i]->setXi(jeuDonnees[iterDataSet++]); } //xi indicator
          else { Errors::errorMessage("MeshCartesianAMR::setDataSet: unknown phase number: ", phase); }
        }
        else { //On veut recuperer les donnees vectorielles
          if (phase >= 0) { //Phases data
            vec.setXYZ(jeuDonnees[iterDataSet], jeuDonnees[iterDataSet + 1], jeuDonnees[iterDataSet + 2]);
            cellsLvl[lvl][i]->getPhase(phase)->setVector(-var, vec);
            iterDataSet += 3;
          }
          else if (phase == -1) {  //Mixture data
            vec.setXYZ(jeuDonnees[iterDataSet], jeuDonnees[iterDataSet + 1], jeuDonnees[iterDataSet + 2]);
            cellsLvl[lvl][i]->getMixture()->setVector(-var, vec);
            iterDataSet += 3;
          }
          else { Errors::errorMessage("MeshCartesianAMR::setDataSet: unknown phase number: ", phase); }
        } //Fin vecteur
      } // Fin split
    } // Fin lvl
  } // Fin levels
}

//***********************************************************************

void MeshCartesianAMR::refineCellAndCellInterfaces(Cell *cell, const std::vector<AddPhys*> &addPhys, Model *model, int &nbCellsTotalAMR)
{
  bool refineExternalCellInterfaces(true);
  cell->refineCellAndCellInterfaces(m_numberCellsY, m_numberCellsZ, addPhys, model, refineExternalCellInterfaces);
  nbCellsTotalAMR += cell->getNumberCellsChildren() - 1;
}

//***********************************************************************

void MeshCartesianAMR::printDomainDecomposition(std::ofstream &fileStream)
{
  m_decomp.printDomainDecomposition(fileStream);
}

//***********************************************************************

void MeshCartesianAMR::readDomainDecomposition(std::ifstream &fileStream)
{
  std::array<int,3> temp={{m_numberCellsXGlobal,m_numberCellsYGlobal,m_numberCellsZGlobal}};
  m_decomp = decomposition::Decomposition(temp);
  m_decomp.readDomainDecomposition(fileStream);
}

//****************************************************************************
//****************************** Parallele ***********************************
//****************************************************************************

void MeshCartesianAMR::initializePersistentCommunications(const int numberPhases, const int numberTransports, const TypeMeshContainer<Cell *> &cells, std::string ordreCalcul)
{
	m_numberPhases = numberPhases;
	m_numberTransports = numberTransports;
	int numberVariablesPhaseATransmettre = cells[0]->getPhase(0)->numberOfTransmittedVariables();
	numberVariablesPhaseATransmettre *= m_numberPhases;
	int numberVariablesMixtureATransmettre = cells[0]->getMixture()->numberOfTransmittedVariables();
	int m_numberPrimitiveVariables = numberVariablesPhaseATransmettre + numberVariablesMixtureATransmettre + m_numberTransports;
  int m_numberSlopeVariables(0);
  if (ordreCalcul == "SECONDORDER") {
    int numberSlopesPhaseATransmettre = cells[0]->getPhase(0)->numberOfTransmittedSlopes();
    numberSlopesPhaseATransmettre *= m_numberPhases;
    int numberSlopesMixtureATransmettre = cells[0]->getMixture()->numberOfTransmittedSlopes();
    m_numberSlopeVariables = numberSlopesPhaseATransmettre + numberSlopesMixtureATransmettre + m_numberTransports + 1 + 1; //+1 for the interface detection + 1 for slope index
  }
	parallel.initializePersistentCommunicationsAMR(m_numberPrimitiveVariables, m_numberSlopeVariables, m_numberTransports, m_geometrie, m_lvlMax);  
}

//***********************************************************************

void MeshCartesianAMR::finalizeParallele(const int &lvlMax)
{
	parallel.finalizeAMR(lvlMax);
}

//***********************************************************************

void MeshCartesianAMR::parallelLoadBalancingAMR(TypeMeshContainer<Cell *> *cellsLvl, TypeMeshContainer<Cell *> *cellsLvlGhost, 
  TypeMeshContainer<CellInterface *> *cellInterfacesLvl, std::string ordreCalcul,
  const int &numberPhases, const int &numberTransports, const std::vector<AddPhys*> &addPhys, 
  Model *model, Eos **eos, int &nbCellsTotalAMR, bool init)
{
//return; //KS//BD//
  bool balance(false);
  do {
    balance = false;
    std::vector<typename decomposition::Key<3>::value_type> indicesSendStartGlobal;
    std::vector<typename decomposition::Key<3>::value_type> indicesSendEndGlobal;
    std::vector<typename decomposition::Key<3>::value_type> indicesReceiveStartGlobal;
    std::vector<typename decomposition::Key<3>::value_type> indicesReceiveEndGlobal;
    //for (int lvl = 0; lvl <= m_lvlMax; ++lvl) { //For levelwise balancing
    int lvl(0); //For global balancing
      this->computePotentialBalancing(cellsLvl, init, lvl, balance, ordreCalcul,
        indicesSendStartGlobal, indicesSendEndGlobal, 
        indicesReceiveStartGlobal, indicesReceiveEndGlobal);
    //} //For levelwise balancing

    if (balance) {
      this->balance(cellsLvl, cellsLvlGhost, cellInterfacesLvl, ordreCalcul, numberPhases, numberTransports, addPhys, model, eos,
        nbCellsTotalAMR, indicesSendStartGlobal, indicesSendEndGlobal, indicesReceiveStartGlobal, indicesReceiveEndGlobal);
    }
  } while (balance);

  //Update gradients for level max only (others are updated within the recursive time-stepping loop)
  for (unsigned int i = 0; i < cellsLvl[m_lvlMax].size(); i++) { if (!cellsLvl[m_lvlMax][i]->getSplit()) { cellsLvl[m_lvlMax][i]->prepareAddPhys(); } }
  for (unsigned int pa = 0; pa < addPhys.size(); pa++) { addPhys[pa]->communicationsAddPhys(m_numberPhases, m_geometrie, m_lvlMax); }
}

//***********************************************************************

void MeshCartesianAMR::computePotentialBalancing(TypeMeshContainer<Cell *> *cellsLvl, bool init, int lvl, bool &balance, std::string ordreCalcul,
    std::vector<typename decomposition::Key<3>::value_type> &indicesSendStartGlobal, std::vector<typename decomposition::Key<3>::value_type> &indicesSendEndGlobal,
    std::vector<typename decomposition::Key<3>::value_type> &indicesReceiveStartGlobal, std::vector<typename decomposition::Key<3>::value_type> &indicesReceiveEndGlobal)
{
  //Typical example of one situation for CPU 1 (during 1 balance iteration):
  //
  //                CPU 0       CPU 1                                 CPU 2                           
  // Current load |-------|---------------|-----------------------------------------------------------|
  //                      |   localLoad   |
  //                      |      localLoadEndPosition
  //             localLoadStartPosition
  //
  //                          CPU 0                       CPU 1                       CPU 2
  // Ideal load   |---------------------------|---------------------------|---------------------------|
  //                                 idealLoadStartPosition      idealLoadEndPosition
  //
  // Ideal load shifts    |------------------->
  //                                      |------------------------------->
  //                       idealLoadShiftStart    idealLoadShiftEnd
  //
  // Possible load shifts |--------------->
  //                                      |------------------------------->
  //                    possibleLoadShiftStart   possibleLoadShiftEnd
  //
  //                        CPU 0                       CPU 1                         CPU 2
  // Final load   |-----------------------|-------------------------------|---------------------------|
  //                                      |        finalLocalLoad         |
  //                                      |                   finalLocalLoadEndPosition
  //                          finalLocalLoadStartPosition
  //
  // Note that the possible load shift start gives birth of a new start that may
  // constrain the possible load shift end (situation not present in this example).

  int counter(0);
  MPI_Request req_neighborP1;
  MPI_Request req_neighborM1;
  MPI_Status status;

  //1) Compute and communicate loads and positions
  //----------------------------------------------

  //Compute local load
  double localLoad(0.);
  for (unsigned int i = 0; i < cellsLvl[0].size(); i++) {
    cellsLvl[0][i]->computeLoad(localLoad, lvl);
  }

  //Communicate overall loads
  double *loadPerCPU = new double[Ncpu];
  MPI_Allgather(&localLoad, 1, MPI_DOUBLE, loadPerCPU, 1, MPI_DOUBLE, MPI_COMM_WORLD);

  //Compute ideal load end position
  double idealLoadEndPosition(0.);
  for (int i = 0; i < Ncpu; ++i) { idealLoadEndPosition += loadPerCPU[i]; }
  idealLoadEndPosition *= static_cast<double>(rankCpu + 1) / Ncpu;

  //Compute local load end position
  double localLoadEndPosition(0.);
  for (int i = 0; i <= rankCpu; ++i) { localLoadEndPosition += loadPerCPU[i]; }
  delete[] loadPerCPU;

  //2) Compute and communicate ideal shifts
  //---------------------------------------

  //Determine load I should send/receive to/from/ CPU P1, using end positions of local load and ideal load
  //if (diff > 0), I wish to receive loads from CPU P1
  //else,          I wish to send loads to CPU P1
  double idealLoadShiftEnd(0.);
  idealLoadShiftEnd = idealLoadEndPosition - localLoadEndPosition;

  //Communicate what I wish to send/receive to/from neighbours
  double idealLoadShiftStart(0.);
  if (rankCpu != Ncpu - 1) {
    MPI_Isend(&idealLoadShiftEnd, 1, MPI_DOUBLE, rankCpu+1, rankCpu+1, MPI_COMM_WORLD, &req_neighborP1);
    MPI_Wait(&req_neighborP1, &status);
  }
  if (rankCpu != 0) {
    MPI_Irecv(&idealLoadShiftStart, 1, MPI_DOUBLE, rankCpu-1, rankCpu, MPI_COMM_WORLD, &req_neighborM1);
    MPI_Wait(&req_neighborM1, &status);
  }

  //3) Compute and communicate possible shifts (real balance)
  //---------------------------------------------------------

  //Determine and communicate what I can send/receive to/from neighbours (limited by current local cells/load)
  double possibleLoadShiftStart(0.), possibleLoadShiftEnd(0.);
  int lvlMax(0);
  int numberOfCellsToSendStart(0), numberOfCellsToSendEnd(0);
  int numberOfCellsToReceiveStart(0), numberOfCellsToReceiveEnd(0);

  //For load shift start
  if (std::fabs(idealLoadShiftStart) > 0.5) {
    if (idealLoadShiftStart > 0.5) {
      //Determine and send possible load shift start
      for (unsigned int i = 0; i < cellsLvl[0].size() - 1; i++) {
        lvlMax = 0;
        cellsLvl[0][i]->computeLvlMax(lvlMax);
        //if (lvlMax == lvl) { //For levelwise balancing
          cellsLvl[0][i]->computeLoad(possibleLoadShiftStart, lvl);
          ++numberOfCellsToSendStart;
        //} //For levelwise balancing
        if (static_cast<int>(std::round(possibleLoadShiftStart)) >= static_cast<int>(std::round(idealLoadShiftStart))) break;
      }
      if (numberOfCellsToSendStart != 0) --numberOfCellsToSendStart;
      MPI_Isend(&numberOfCellsToSendStart, 1, MPI_INT, rankCpu-1, rankCpu, MPI_COMM_WORLD, &req_neighborM1);
      MPI_Wait(&req_neighborM1, &status);
    }
    else {
      //Receive possible load shift start
      MPI_Irecv(&numberOfCellsToReceiveStart, 1, MPI_INT, rankCpu-1, rankCpu, MPI_COMM_WORLD, &req_neighborM1);
      MPI_Wait(&req_neighborM1, &status);
    }
  }

  //For load shift end
  if (std::fabs(idealLoadShiftEnd) > 0.5) {
    if (idealLoadShiftEnd < -0.5) {
      //Determine and send possible load shift end
      for (int i = cellsLvl[0].size() - 1; i >= 0; --i) {
        lvlMax = 0;
        cellsLvl[0][i]->computeLvlMax(lvlMax);
         //if (lvlMax == lvl) { //For levelwise balancing
          cellsLvl[0][i]->computeLoad(possibleLoadShiftEnd, lvl);
          ++numberOfCellsToSendEnd;
         //} //For levelwise balancing
        if (static_cast<int>(std::round(-possibleLoadShiftEnd)) <= static_cast<int>(std::round(idealLoadShiftEnd)) ||
            (numberOfCellsToSendEnd + numberOfCellsToSendStart + 1) == cellsLvl[0].size()
            ) break;
      }
      if (numberOfCellsToSendEnd != 0) --numberOfCellsToSendEnd;
      MPI_Isend(&numberOfCellsToSendEnd, 1, MPI_INT, rankCpu+1, rankCpu+1, MPI_COMM_WORLD, &req_neighborP1);
      MPI_Wait(&req_neighborP1, &status);
    }
    else {
      //Receive possible load shift end
      MPI_Irecv(&numberOfCellsToReceiveEnd, 1, MPI_INT, rankCpu+1, rankCpu+1, MPI_COMM_WORLD, &req_neighborP1);
      MPI_Wait(&req_neighborP1, &status);
    }
  }
  // double localLoadStartPosition = localLoadEndPosition - localLoad;
  // localLoadStartPosition += possibleLoadShiftStart;
  // localLoadEndPosition += possibleLoadShiftEnd;
  // double finalLocalLoad =localLoadEndPosition-localLoadStartPosition;
  // std::cout<<"cpu "<<rankCpu<<" lvl "<<lvl<<" localLoadStartPosition "<<localLoadStartPosition<<" localLoadEndPosition "<<localLoadEndPosition
  // <<" initialLocalLoad "<<localLoad
  // <<" finalLocalLoad "<<finalLocalLoad<<std::endl;

  //4) Update criterion to balance
  //------------------------------
  double relativePossibleLoadShiftMax(0.), relativePossibleLoadShiftLocal(0.);
  relativePossibleLoadShiftLocal = std::max(std::max(numberOfCellsToSendStart, numberOfCellsToReceiveStart), std::max(numberOfCellsToSendEnd,numberOfCellsToReceiveEnd));
  if (localLoad > 1.e-8) { relativePossibleLoadShiftLocal /= localLoad; }
  MPI_Allreduce(&relativePossibleLoadShiftLocal, &relativePossibleLoadShiftMax, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

  if (init) {
    if (relativePossibleLoadShiftMax > 1.e-10) { balance = true; }
  }
  else {
    if (relativePossibleLoadShiftMax > 0.05) { balance = true; }
  }

  if (!balance ) return;

  //5) Send/Receive indices of base cell keys (lvl = 0) and partially update the domain decomposition with new starts and ends
  //--------------------------------------------------------------------------------------------------------------------------

  //Send/receive indices of the base cell keys
  std::vector<typename decomposition::Key<3>::value_type> indicesSendStart;
  if (numberOfCellsToSendStart > 0) {
    for (unsigned int i = 0; i < cellsLvl[0].size(); ++i) {
      lvlMax = 0;
      cellsLvl[0][i]->computeLvlMax(lvlMax);
      //if (lvlMax == lvl) { //For levelwise balancing
        indicesSendStart.push_back(cellsLvl[0][i]->getElement()->getKey().getIndex());
        if (indicesSendStart.size() == numberOfCellsToSendStart) { break; }
      //} //For levelwise balancing
    }
    MPI_Isend(&indicesSendStart[0], numberOfCellsToSendStart, MPI_UNSIGNED_LONG_LONG, rankCpu-1, rankCpu, MPI_COMM_WORLD, &req_neighborM1);
    MPI_Wait(&req_neighborM1, &status);
  }

  std::vector<typename decomposition::Key<3>::value_type> indicesReceiveEnd(numberOfCellsToReceiveEnd);
  if (numberOfCellsToReceiveEnd > 0) {
    MPI_Irecv(&indicesReceiveEnd[0], numberOfCellsToReceiveEnd, MPI_UNSIGNED_LONG_LONG, rankCpu+1, rankCpu+1, MPI_COMM_WORLD, &req_neighborP1);
    MPI_Wait(&req_neighborP1, &status);
  }

  std::vector<typename decomposition::Key<3>::value_type> indicesSendEnd;
  if (numberOfCellsToSendEnd > 0) {
    for (int i = cellsLvl[0].size() - 1; i >= 0 ; --i) {
      lvlMax = 0;
      cellsLvl[0][i]->computeLvlMax(lvlMax);
      //if (lvlMax == lvl) { //For levelwise balancing
        indicesSendEnd.push_back(cellsLvl[0][i]->getElement()->getKey().getIndex());
        if (indicesSendEnd.size() == numberOfCellsToSendEnd) { break; }
      //} //For levelwise balancing
    }
    MPI_Isend(&indicesSendEnd[0], numberOfCellsToSendEnd, MPI_UNSIGNED_LONG_LONG, rankCpu+1, rankCpu+1, MPI_COMM_WORLD, &req_neighborP1);
    MPI_Wait(&req_neighborP1, &status);
  }

  std::vector<typename decomposition::Key<3>::value_type> indicesReceiveStart(numberOfCellsToReceiveStart);
  if (numberOfCellsToReceiveStart > 0) {
    MPI_Irecv(&indicesReceiveStart[0], numberOfCellsToReceiveStart, MPI_UNSIGNED_LONG_LONG, rankCpu-1, rankCpu, MPI_COMM_WORLD, &req_neighborM1);
    MPI_Wait(&req_neighborM1, &status);
  }

  //Accumulation of send/receive of each level
  indicesSendStartGlobal.insert(indicesSendStartGlobal.end(), indicesSendStart.begin(), indicesSendStart.end());
  indicesSendEndGlobal.insert(indicesSendEndGlobal.end(), indicesSendEnd.begin(), indicesSendEnd.end());
  indicesReceiveStartGlobal.insert(indicesReceiveStartGlobal.end(), indicesReceiveStart.begin(), indicesReceiveStart.end());
  indicesReceiveEndGlobal.insert(indicesReceiveEndGlobal.end(), indicesReceiveEnd.begin(), indicesReceiveEnd.end());
}

//***********************************************************************

void MeshCartesianAMR::balance(TypeMeshContainer<Cell *> *cellsLvl, TypeMeshContainer<Cell *> *cellsLvlGhost, 
  TypeMeshContainer<CellInterface *> *cellInterfacesLvl, std::string ordreCalcul,
  const int &numberPhases, const int &numberTransports, const std::vector<AddPhys*> &addPhys, Model *model, Eos **eos, int &nbCellsTotalAMR,
    std::vector<typename decomposition::Key<3>::value_type> &indicesSendStartGlobal, std::vector<typename decomposition::Key<3>::value_type> &indicesSendEndGlobal,
    std::vector<typename decomposition::Key<3>::value_type> &indicesReceiveStartGlobal, std::vector<typename decomposition::Key<3>::value_type> &indicesReceiveEndGlobal)
{
  int counter(0), counterSplit(0);
  MPI_Request req_neighborP1;
  MPI_Request req_neighborM1;
  MPI_Status status;

  std::sort(indicesSendStartGlobal.begin(), indicesSendStartGlobal.end());
  std::sort(indicesSendEndGlobal.begin(), indicesSendEndGlobal.end());
  std::sort(indicesReceiveStartGlobal.begin(), indicesReceiveStartGlobal.end());
  std::sort(indicesReceiveEndGlobal.begin(),   indicesReceiveEndGlobal.end());

  const auto numberOfCellsToReceiveStartGlobal=indicesReceiveStartGlobal.size();
  const auto numberOfCellsToReceiveEndGlobal=indicesReceiveEndGlobal.size();
  const auto numberOfCellsToSendStartGlobal=indicesSendStartGlobal.size();
  const auto numberOfCellsToSendEndGlobal=indicesSendEndGlobal.size();

  //1) Create the corresponding base cells (received)
  //-------------------------------------------------

  //Create and insert cells and elements
  TypeMeshContainer<Cell *> bufferReceiveCellsStart(numberOfCellsToReceiveStartGlobal);
  TypeMeshContainer<Cell *> bufferReceiveCellsEnd(numberOfCellsToReceiveEndGlobal);
  TypeMeshContainer<Cell *> bufferReceiveCells;
  std::vector<decomposition::Key<3>> keysReceiveStart(numberOfCellsToReceiveStartGlobal);
  std::vector<decomposition::Key<3>> keysReceiveEnd(numberOfCellsToReceiveEndGlobal);
  m_elements.clear();
  for(unsigned int i = 0; i < numberOfCellsToReceiveStartGlobal; ++i) {
    if (ordreCalcul == "FIRSTORDER") { bufferReceiveCellsStart[i] = new Cell; }
    else { bufferReceiveCellsStart[i] = new CellO2; }
    m_elements.push_back(new ElementCartesian());
    keysReceiveStart[i] = decomposition::Key<3>(indicesReceiveStartGlobal[i]);
    m_elements.back()->setKey(keysReceiveStart[i]);
    bufferReceiveCellsStart[i]->setElement(m_elements.back(), i);
  }
  for(unsigned int i = 0; i < numberOfCellsToReceiveEndGlobal; ++i) {
    if (ordreCalcul == "FIRSTORDER") { bufferReceiveCellsEnd[i] = new Cell; }
    else { bufferReceiveCellsEnd[i] = new CellO2; }
    m_elements.push_back(new ElementCartesian());
    keysReceiveEnd[i] = decomposition::Key<3>(indicesReceiveEndGlobal[i]);
    m_elements.back()->setKey(keysReceiveEnd[i]);
    bufferReceiveCellsEnd[i]->setElement(m_elements.back(), i);
  }

  cellsLvl[0].insert(cellsLvl[0].begin(), bufferReceiveCellsStart.begin(), bufferReceiveCellsStart.end());
  cellsLvl[0].insert(cellsLvl[0].end(), bufferReceiveCellsEnd.begin(), bufferReceiveCellsEnd.end());
  std::sort(cellsLvl[0].begin(),cellsLvl[0].end(),[&](Cell* child0, Cell* child1)
      {
        return child0->getElement()->getKey() < child1->getElement()->getKey();
      });

  //Assigning element properties
  this->assignElementProperties(bufferReceiveCellsStart, keysReceiveStart);
  this->assignElementProperties(bufferReceiveCellsEnd, keysReceiveEnd);

  bufferReceiveCells.insert(bufferReceiveCells.end(), bufferReceiveCellsStart.begin(), bufferReceiveCellsStart.end());
  bufferReceiveCells.insert(bufferReceiveCells.end(), bufferReceiveCellsEnd.begin(), bufferReceiveCellsEnd.end());

  //2) Erase the pointers to the corresponding base cells from cellsLvl (cells are not deleted yet)
  //-----------------------------------------------------------------------------------------------

  TypeMeshContainer<Cell *> bufferSendCellsStart(numberOfCellsToSendStartGlobal);
  TypeMeshContainer<Cell *> bufferSendCellsEnd(numberOfCellsToSendEndGlobal);
  TypeMeshContainer<Cell *> bufferSendCells;
  if (indicesSendStartGlobal.size() != 0) {
    counter = 0;
    auto it = cellsLvl[0].begin();
    while (it != cellsLvl[0].end()) {
      if (indicesSendStartGlobal[counter] == (*it)->getElement()->getKey().getIndex()) 
      // match indices in current cellLvl with indicesSend
      {
        bufferSendCellsStart[counter] = *it;
        ++counter;
        it = cellsLvl[0].erase(it);
        if (counter == numberOfCellsToSendStartGlobal) break;
      }
      else {
        ++it;
      }
    }
  }

  if (indicesSendEndGlobal.size() != 0) {
    counter = 0;
    auto rit = cellsLvl[0].rbegin();
    while (rit != cellsLvl[0].rend()) {
      if (indicesSendEndGlobal[numberOfCellsToSendEndGlobal - 1 - counter] == (*rit)->getElement()->getKey().getIndex()) {
        bufferSendCellsEnd[numberOfCellsToSendEndGlobal - 1 - counter] = *rit;
        ++counter;
        auto it = cellsLvl[0].erase(--rit.base());
        rit = std::reverse_iterator<std::vector<Cell*>::iterator>(it);
        if (counter == numberOfCellsToSendEndGlobal)break;
      }
      else {
        ++rit;
      }
    }
  }

  bufferSendCells.insert(bufferSendCells.end(), bufferSendCellsStart.begin(), bufferSendCellsStart.end());
  bufferSendCells.insert(bufferSendCells.end(), bufferSendCellsEnd.begin(), bufferSendCellsEnd.end());

  //3) Update the domain decomposition with new starts and ends
  //-----------------------------------------------------------
  //Communicate between each CPU their respective new starts
  std::vector<typename decomposition::Key<3>::value_type> localKeys;
  std::vector<int> localRanks;
  localKeys.push_back(cellsLvl[0][0]->getElement()->getKey().getIndex());
  localRanks.push_back(rankCpu);
  for (unsigned int i = 1; i < cellsLvl[0].size(); i++)
  {
      if (!m_decomp.areConsecutive(cellsLvl[0][i-1]->getElement()->getKey().getIndex(), 
                  cellsLvl[0][i]->getElement()->getKey().getIndex())
         )
      // if not consecutive
      {
          localKeys.push_back(cellsLvl[0][i]->getElement()->getKey().getIndex());
          localRanks.push_back(rankCpu);
      }
  }
  m_decomp.communicateMaps(Ncpu, localKeys, localRanks, rankCpu);//- fane, not understand

  //Recombine starts of all CPU
  m_decomp.recombineStarts();//- fane, not understand

  //4) Create cell interfaces, faces and ghost cells of level 0
  //-----------------------------------------------------------
  for (int b = 0; b < cellInterfacesLvl[0].size(); b++) { delete cellInterfacesLvl[0][b]; }
  for (int i = 0; i < cellsLvl[0].size(); i++) { cellsLvl[0][i]->clearExternalCellInterfaces(m_numberCellsY, m_numberCellsZ); }
  for (int i = 0; i < cellsLvlGhost[0].size(); i++) { delete cellsLvlGhost[0][i]; }
  cellInterfacesLvl[0].clear();
  m_faces.clear();
  cellsLvlGhost[0].clear();
  parallel.clearElementsAndSlopesToSendAndReceivePLusNeighbour();
  m_numberCellsCalcul = cellsLvl[0].size();
  // fane, create cell, cellGhost, cellInterface and relation, but not set physical values
  createCellInterfacesFacesAndGhostCells(cellsLvl[0], cellsLvlGhost[0], cellInterfacesLvl[0], ordreCalcul);
  m_numberCellsTotal = cellsLvl[0].size() + cellsLvlGhost[0].size();//-fane
  m_numberFacesTotal = cellInterfacesLvl[0].size();
//   std::cout<<"cpu "<<rankCpu
//   << " m_numberCellsCalcul "<<m_numberCellsCalcul
//   << " m_numberCellsTotal "<<m_numberCellsTotal
//   << " m_numberFacesTotal "<<m_numberFacesTotal
//   <<std::endl;

  //5) Allocate physical variables of cells and cell interfaces level 0
  //-------------------------------------------------------------------
  for (int i = 0; i < bufferReceiveCells.size(); i++) { bufferReceiveCells[i]->allocate(numberPhases, numberTransports, addPhys, model); }
  for (int i = 0; i < cellsLvlGhost[0].size(); i++) { cellsLvlGhost[0][i]->allocate(numberPhases, numberTransports, addPhys, model); }
  //Attribution model and slopes to faces
  int allocateSlopeLocal = 1;
  for (int b = 0; b < cellInterfacesLvl[0].size(); b++) {
    cellInterfacesLvl[0][b]->associeModel(model);
    cellInterfacesLvl[0][b]->allocateSlopes(numberPhases, numberTransports, allocateSlopeLocal);
  }

  //6) Send/Receive physical values of cells lvl >= 0 and create new cells and new internal cell interfaces of lvl > 0
  //------------------------------------------------------------------------------------------------------------------
  //Count and communicate number of data (physical values and cell trees) to send/receive + fill corresponding vectors
  int numberSendStart(0), numberSendEnd(0), numberReceiveStart(0), numberReceiveEnd(0);
  int numberSplitSendStart(0), numberSplitSendEnd(0), numberSplitReceiveStart(0), numberSplitReceiveEnd(0);
  std::vector<double> dataToSendStart, dataToSendEnd;
  std::vector<int> dataSplitToSendStart, dataSplitToSendEnd;
  if (numberOfCellsToSendStartGlobal > 0) {
    //Count and fill buffer vector send
    for (int lvl = 0; lvl <= m_lvlMax; lvl++) {
      for (int i = 0; i < bufferSendCellsStart.size(); i++) {
        // fane
        bufferSendCellsStart[i]->fillDataToSend(dataToSendStart, dataSplitToSendStart, lvl);
      }
    }
    numberSendStart = dataToSendStart.size();
    numberSplitSendStart = dataSplitToSendStart.size();
    MPI_Isend(&numberSendStart, 1, MPI_INT, rankCpu-1, rankCpu, MPI_COMM_WORLD, &req_neighborM1);
    MPI_Wait(&req_neighborM1, &status);
    MPI_Isend(&numberSplitSendStart, 1, MPI_INT, rankCpu-1, rankCpu, MPI_COMM_WORLD, &req_neighborM1);
    MPI_Wait(&req_neighborM1, &status);
  }
  if (numberOfCellsToReceiveEndGlobal > 0) {
    MPI_Irecv(&numberReceiveEnd, 1, MPI_INT, rankCpu+1, rankCpu+1, MPI_COMM_WORLD, &req_neighborP1);
    MPI_Wait(&req_neighborP1, &status);
    MPI_Irecv(&numberSplitReceiveEnd, 1, MPI_INT, rankCpu+1, rankCpu+1, MPI_COMM_WORLD, &req_neighborP1);
    MPI_Wait(&req_neighborP1, &status);
  }
  if (numberOfCellsToSendEndGlobal > 0) {
    //Count and fill buffer vector send
    for (int lvl = 0; lvl <= m_lvlMax; lvl++) {
      for (int i = 0; i < bufferSendCellsEnd.size(); i++) {
        bufferSendCellsEnd[i]->fillDataToSend(dataToSendEnd, dataSplitToSendEnd, lvl);
      }
    }
    numberSendEnd = dataToSendEnd.size();
    numberSplitSendEnd = dataSplitToSendEnd.size();
    MPI_Isend(&numberSendEnd, 1, MPI_INT, rankCpu+1, rankCpu+1, MPI_COMM_WORLD, &req_neighborP1);
    MPI_Wait(&req_neighborP1, &status);
    MPI_Isend(&numberSplitSendEnd, 1, MPI_INT, rankCpu+1, rankCpu+1, MPI_COMM_WORLD, &req_neighborP1);
    MPI_Wait(&req_neighborP1, &status);
  }
  if (numberOfCellsToReceiveStartGlobal > 0) {
    MPI_Irecv(&numberReceiveStart, 1, MPI_INT, rankCpu-1, rankCpu, MPI_COMM_WORLD, &req_neighborM1);
    MPI_Wait(&req_neighborM1, &status);
    MPI_Irecv(&numberSplitReceiveStart, 1, MPI_INT, rankCpu-1, rankCpu, MPI_COMM_WORLD, &req_neighborM1);
    MPI_Wait(&req_neighborM1, &status);
  }

  //Send/Receive data
  std::vector<double> dataToReceiveStart(numberReceiveStart), dataToReceiveEnd(numberReceiveEnd);
  std::vector<int> dataSplitToReceiveStart(numberSplitReceiveStart), dataSplitToReceiveEnd(numberSplitReceiveEnd);
  if (numberOfCellsToSendStartGlobal > 0) {
    MPI_Isend(&dataToSendStart[0], numberSendStart, MPI_DOUBLE, rankCpu-1, rankCpu, MPI_COMM_WORLD, &req_neighborM1);
    MPI_Wait(&req_neighborM1, &status);
    MPI_Isend(&dataSplitToSendStart[0], numberSplitSendStart, MPI_INT, rankCpu-1, rankCpu, MPI_COMM_WORLD, &req_neighborM1);
    MPI_Wait(&req_neighborM1, &status);
  }
  if (numberOfCellsToReceiveEndGlobal > 0) {
    MPI_Irecv(&dataToReceiveEnd[0], numberReceiveEnd, MPI_DOUBLE, rankCpu+1, rankCpu+1, MPI_COMM_WORLD, &req_neighborP1);
    MPI_Wait(&req_neighborP1, &status);
    MPI_Irecv(&dataSplitToReceiveEnd[0], numberSplitReceiveEnd, MPI_INT, rankCpu+1, rankCpu+1, MPI_COMM_WORLD, &req_neighborP1);
    MPI_Wait(&req_neighborP1, &status);
    //Get buffer vector receive + Refine cells and internal cell interfaces
    // fane - receive the data and process, fnae, works but not solve the problem
    counter = 0; counterSplit = 0;
    for (int lvl = 0; lvl <= m_lvlMax; lvl++) {
      for (int i = 0; i < bufferReceiveCellsEnd.size(); i++) {
        bufferReceiveCellsEnd[i]->getDataToReceiveAndRefine(dataToReceiveEnd, dataSplitToReceiveEnd, lvl, eos, counter, counterSplit, m_numberCellsY, m_numberCellsZ, addPhys, model);
      }
    }
  }
  if (numberOfCellsToSendEndGlobal > 0) {
    MPI_Isend(&dataToSendEnd[0], numberSendEnd, MPI_DOUBLE, rankCpu+1, rankCpu+1, MPI_COMM_WORLD, &req_neighborP1);
    MPI_Wait(&req_neighborP1, &status);
    MPI_Isend(&dataSplitToSendEnd[0], numberSplitSendEnd, MPI_INT, rankCpu+1, rankCpu+1, MPI_COMM_WORLD, &req_neighborP1);
    MPI_Wait(&req_neighborP1, &status);
  }
  if (numberOfCellsToReceiveStartGlobal > 0) {
    MPI_Irecv(&dataToReceiveStart[0], numberReceiveStart, MPI_DOUBLE, rankCpu-1, rankCpu, MPI_COMM_WORLD, &req_neighborM1);
    MPI_Wait(&req_neighborM1, &status);
    MPI_Irecv(&dataSplitToReceiveStart[0], numberSplitReceiveStart, MPI_INT, rankCpu-1, rankCpu, MPI_COMM_WORLD, &req_neighborM1);
    MPI_Wait(&req_neighborM1, &status);
    //Get buffer vector receive + Refine cells and internal cell interfaces
    counter = 0; counterSplit = 0;
    for (int lvl = 0; lvl <= m_lvlMax; lvl++) {
      for (int i = 0; i < bufferReceiveCellsStart.size(); i++) {
        bufferReceiveCellsStart[i]->getDataToReceiveAndRefine(dataToReceiveStart, dataSplitToReceiveStart, lvl, eos, counter, counterSplit, m_numberCellsY, m_numberCellsZ, addPhys, model);
      }
    }
  }

  //Delete sent cells
  for (int i = 0; i < bufferSendCells.size(); i++) { delete bufferSendCells[i]; }

  //7) Update persistent communications of cells lvl 0
  //--------------------------------------------------
  int dim(1);
  if (m_numberCellsZ != 1) { dim = 3; }
  else if (m_numberCellsY != 1) { dim = 2; }
  parallel.updatePersistentCommunicationsAMR(dim);
  
  //8) Refine external cell interfaces, refine ghost cells and update persistent communications of lvl > 0
  //------------------------------------------------------------------------------------------------------
  for (int lvl = 0; lvl < m_lvlMax; lvl++) {
    //Refine external cell interfaces
    for (int i = 0; i < cellsLvl[lvl].size(); i++) {
      if (cellsLvl[lvl][i]->getSplit()) {
        for (int b = 0; b < cellsLvl[lvl][i]->getCellInterfacesSize(); b++) {
          if (!cellsLvl[lvl][i]->getCellInterface(b)->getSplit()) {
            cellsLvl[lvl][i]->getCellInterface(b)->raffineCellInterfaceExterne(m_numberCellsY, m_numberCellsZ, 
              cellsLvl[lvl][i]->getElement()->getSizeX(), cellsLvl[lvl][i]->getElement()->getSizeY(), cellsLvl[lvl][i]->getElement()->getSizeZ(), cellsLvl[lvl][i], dim);
          }
        }
        cellsLvl[lvl][i]->updatePointersInternalCellInterfaces();
      }
      cellsLvl[lvl][i]->setToZeroCons(numberPhases, numberTransports);
    }

    //Refine ghost cells
    parallel.communicationsSplit(lvl);
    cellsLvlGhost[lvl + 1].clear();
    for (unsigned int i = 0; i < cellsLvlGhost[lvl].size(); i++) { cellsLvlGhost[lvl][i]->chooseRefineDeraffineGhost(m_numberCellsY, m_numberCellsZ, addPhys, model, cellsLvlGhost); }
    parallel.communicationsPrimitives(eos, lvl);

    //Update of persistent communications of cells lvl + 1
    parallel.communicationsNumberGhostCells(lvl + 1);
    parallel.updatePersistentCommunicationsLvlAMR(lvl + 1, m_geometrie);

    //Reconstruction of the arrays of cells and cell interfaces of lvl + 1
    cellsLvl[lvl + 1].clear();
    cellInterfacesLvl[lvl + 1].clear();
    for (unsigned int i = 0; i < cellsLvl[lvl].size(); i++) { cellsLvl[lvl][i]->buildLvlCellsAndLvlInternalCellInterfacesArrays(cellsLvl, cellInterfacesLvl); }
    for (unsigned int i = 0; i < cellInterfacesLvl[lvl].size(); i++) { cellInterfacesLvl[lvl][i]->constructionTableauCellInterfacesExternesLvl(cellInterfacesLvl); }
  }
  parallel.communicationsPrimitives(eos, m_lvlMax);
  parallel.communicationsPrimitives(eos, m_lvlMax,vecPhasesO2);
  nbCellsTotalAMR = 0;
  for (int i = 0; i < cellsLvl[0].size(); i++) { cellsLvl[0][i]->updateNbCellsTotalAMR(nbCellsTotalAMR); }
}

//***********************************************************************
