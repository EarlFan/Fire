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

//! \file      CellInterface.cpp
//! \author    F. Petitpas, K. Schmidmayer, S. Le Martelot, B. Dorschner
//! \version   1.1
//! \date      June 5 2019

#include "CellInterface.h"
#include <iostream>

//Utile pour la resolution des problemes de Riemann
Cell *cellLeft;
Cell *cellRight;

//***********************************************************************

CellInterface::CellInterface() : m_mod(0), m_cellLeft(0), m_cellRight(0), m_face(0), m_cellInterfacesChildren(0)
{
  m_lvl = 0;
}

//***********************************************************************

CellInterface::CellInterface(int lvl) : m_mod(0), m_cellLeft(0), m_cellRight(0), m_face(0), m_cellInterfacesChildren(0)
{
  m_lvl = lvl;
}

//***********************************************************************

CellInterface::~CellInterface()
{
  for (unsigned int i = 0; i < m_cellInterfacesChildren.size(); i++) {
    delete m_cellInterfacesChildren[i];
  }
  m_cellInterfacesChildren.clear();
  delete m_face;
}

//***********************************************************************

void CellInterface::initialize(Cell *cellLeft, Cell *cellRight)
{
  m_cellLeft = cellLeft;
  m_cellRight = cellRight;
}

//***********************************************************************

void CellInterface::initializeGauche(Cell *cellLeft)
{
  m_cellLeft = cellLeft;
}

//***********************************************************************

void CellInterface::initializeDroite(Cell *cellRight)
{
  m_cellRight = cellRight;
}

//***********************************************************************

void CellInterface::setFace(Face *face)
{
  m_face = face;
}

//***********************************************************************

void CellInterface::computeFlux(const int &numberPhases, const int &numberTransports, double &dtMax, Limiter &globalLimiter, Limiter &interfaceLimiter, Limiter &globalVolumeFractionLimiter, Limiter &interfaceVolumeFractionLimiter, Prim type)
{
  this->solveRiemann(numberPhases, numberTransports, dtMax, globalLimiter, interfaceLimiter, globalVolumeFractionLimiter, interfaceVolumeFractionLimiter, type);

  if (m_cellLeft->getLvl() == m_cellRight->getLvl()) {     //CoefAMR = 1 pour les deux
    this->addFlux(numberPhases, numberTransports, 1.);      //Ajout du flux sur maille droite
    this->subtractFlux(numberPhases, numberTransports, 1.);     //Retrait du flux sur maille gauche
  }
  else if (m_cellLeft->getLvl() > m_cellRight->getLvl()) { //CoefAMR = 1 pour la gauche et 0.5 pour la droite
    this->addFlux(numberPhases, numberTransports, 0.5);     //Ajout du flux sur maille droite
    this->subtractFlux(numberPhases, numberTransports, 1.);     //Retrait du flux sur maille gauche
  }
  else {                                                      //CoefAMR = 0.5 pour la gauche et 1 pour la droite
    this->addFlux(numberPhases, numberTransports, 1.);      //Ajout du flux sur maille droite
    this->subtractFlux(numberPhases, numberTransports, 0.5);    //Retrait du flux sur maille gauche
  }
}

//***********************************************************************

void CellInterface::computeFluxAddPhys(const int &numberPhases, AddPhys &addPhys)
{
  addPhys.computeFluxAddPhys(this, numberPhases);
}

//***********************************************************************

void CellInterface::solveRiemann(const int &numberPhases, const int &numberTransports, double &dtMax, Limiter &globalLimiter, Limiter &interfaceLimiter, Limiter &globalVolumeFractionLimiter, Limiter &interfaceVolumeFractionLimiter, Prim type)
{
  //Projection des velocities sur repere attache a la face
  //Projection of velocities on landmark attached to the face
  m_cellLeft->localProjection(m_face->getNormal(), m_face->getTangent(), m_face->getBinormal(), numberPhases);
  m_cellRight->localProjection(m_face->getNormal(), m_face->getTangent(), m_face->getBinormal(), numberPhases);

  //Probleme de Riemann
  double dxLeft(m_cellLeft->getElement()->getLCFL());
  double dxRight(m_cellRight->getElement()->getLCFL());
  dxLeft = dxLeft*std::pow(2., (double)m_lvl);// recover to dx at 0th level
  dxRight = dxRight*std::pow(2., (double)m_lvl);
  m_mod->solveRiemannIntern(*m_cellLeft, *m_cellRight, numberPhases, dxLeft, dxRight, dtMax);
  //Traitement des fonctions de transport (m_Sm connu : doit etre place apres l appel au Solveur de Riemann)
  // Processing of transport functions (known m_Sm: must be placed after calling the Riemann Solver)
  if (numberTransports > 0) { m_mod->solveRiemannTransportIntern(*m_cellLeft, *m_cellRight, numberTransports); }

  //Projection du flux sur le repere absolu
  m_mod->reverseProjection(m_face->getNormal(), m_face->getTangent(), m_face->getBinormal());
  m_cellLeft->reverseProjection(m_face->getNormal(), m_face->getTangent(), m_face->getBinormal(), numberPhases);
  m_cellRight->reverseProjection(m_face->getNormal(), m_face->getTangent(), m_face->getBinormal(), numberPhases);
}

//***********************************************************************

void CellInterface::addFlux(const int &numberPhases, const int &numberTransports, const double &coefAMR)
{
  //No "time step"
  double coefA = m_face->getSurface() / m_cellRight->getElement()->getVolume() * coefAMR;
  m_cellRight->getCons()->addFlux(coefA, numberPhases);
  m_cellRight->getCons()->addNonCons(coefA, m_cellRight, numberPhases);
  for (int k = 0; k < numberTransports; k++) {
    m_cellRight->getConsTransport(k)->addFlux(coefA, k);
    m_cellRight->getConsTransport(k)->addNonCons(coefA, m_cellRight->getTransport(k).getValue(), m_mod->getSM());
  }
}

//***********************************************************************

void CellInterface::subtractFlux(const int &numberPhases, const int &numberTransports, const double &coefAMR)
{
  //No "time step"
  double coefA = m_face->getSurface() / m_cellLeft->getElement()->getVolume() * coefAMR;
  m_cellLeft->getCons()->subtractFlux(coefA, numberPhases);
  m_cellLeft->getCons()->subtractNonCons(coefA, m_cellLeft, numberPhases);
  for (int k = 0; k < numberTransports; k++) {
    m_cellLeft->getConsTransport(k)->subtractFlux(coefA, k);
    m_cellLeft->getConsTransport(k)->subtractNonCons(coefA, m_cellLeft->getTransport(k).getValue(), m_mod->getSM());
  }
}

//***********************************************************************

double CellInterface::distance(Cell *c)
{
  return m_face->distance(c->getElement());
}

//***********************************************************************
void CellInterface::EffetsSurface1D(const int &numberPhases)
{
  if (m_cellRight != NULL){ m_cellRight->getCons()->addTuyere1D(m_face->getNormal(), m_face->getSurface(), m_cellRight, numberPhases); }
  m_cellLeft->getCons()->subtractTuyere1D(m_face->getNormal(), m_face->getSurface(), m_cellLeft, numberPhases);
}

//***********************************************************************

void CellInterface::associeModel(Model *model)
{
  m_mod = model;
}

//***********************************************************************

Face *CellInterface::getFace()
{
  return m_face;
}

//***********************************************************************

Model *CellInterface::getMod() const
{
  return m_mod;
}

//***********************************************************************

Cell *CellInterface::getCellGauche() const
{
  return m_cellLeft;
}

//***********************************************************************

Cell *CellInterface::getCellDroite() const
{
  return m_cellRight;
}

//****************************************************************************
//******************************Method AMR***********************************
//****************************************************************************



//***********************************************************************

void CellInterface::computeFluxXi()
{
  double Xi_flux = AMRPara::amr_flux;
  double Xi_lower = AMRPara::xi_diff_lim;
  if (m_cellLeft->getLvl() > m_cellRight->getLvl()) {
    if ((m_cellLeft->getXi() > Xi_lower)) {
      m_cellLeft->addFluxXi(Xi_flux);
      m_cellRight->addFluxXi(Xi_flux);
    }
  }
  else if (m_cellLeft->getLvl() < m_cellRight->getLvl()) {
    if ((m_cellRight->getXi() > Xi_lower)) {
      m_cellLeft->addFluxXi(Xi_flux);
      m_cellRight->addFluxXi(Xi_flux);
    }
  }
  else {
    if ((m_cellLeft->getXi() > Xi_lower) || (m_cellRight->getXi() > Xi_lower)) {
      m_cellLeft->addFluxXi(Xi_flux);
      m_cellRight->addFluxXi(Xi_flux);
    }
  }
}

//***********************************************************************

void CellInterface::creerCellInterfaceChild()
{
  m_cellInterfacesChildren.push_back(new CellInterface(m_lvl + 1));
}

//***********************************************************************

void CellInterface::creerCellInterfaceChildInterne(const int &lvl, std::vector<CellInterface*> *childrenInternalCellInterfaces)
{
  (*childrenInternalCellInterfaces).push_back(new CellInterface(lvl + 1));
}

//***********************************************************************

void CellInterface::creerFaceChild(CellInterface *cellInterfaceParent)
{
  m_face = cellInterfaceParent->m_face->creerNouvelleFace();
}

//***********************************************************************

void CellInterface::raffineCellInterfaceExterne(const int &nbCellsY, const int &nbCellsZ, const double &dXParent, const double &dYParent, const double &dZParent, Cell *cellRef, const int &dim)
{
  //La creation des children cell interfaces n'est pas systematique, on regarde d'abord si ces children cell interfaces ne sont pas deja crees.
  //Dans tous les cas on re-attribut les liaisons cells/cell interfaces.

  double epsilon(1.e-6);
  int allocateSlopeLocal(1);
  double surfaceChild(std::pow(0.5,dim-1.)*m_face->getSurface());

  if (nbCellsZ == 1) {
    if (nbCellsY == 1) {

      //--------------------------------------------------
      //--------------------- Cas 1D ---------------------
      //--------------------------------------------------

      //Cell interface pas encore split -> Creation des children cell interfaces
      //------------------------------------------------------------------------
      if (m_lvl == cellRef->getLvl()) {

        this->creerCellInterfaceChild();
        m_cellInterfacesChildren[0]->m_face = m_face->creerNouvelleFace();
        m_cellInterfacesChildren[0]->m_face->initializeAutres(surfaceChild, m_face->getNormal(), m_face->getTangent(), m_face->getBinormal());
        m_cellInterfacesChildren[0]->m_face->setPos(m_face->getPos().getX(), m_face->getPos().getY(), m_face->getPos().getZ());
        m_cellInterfacesChildren[0]->m_face->setSize(m_face->getSize());
        if (m_face->getPos().getX() < cellRef->getElement()->getPosition().getX()) {
          //Cell interface number 1 (gauche)
          m_cellInterfacesChildren[0]->initializeGauche(m_cellLeft);
          m_cellInterfacesChildren[0]->initializeDroite(cellRef->getCellChild(0));
          m_cellLeft->addCellInterface(m_cellInterfacesChildren[0]);
          cellRef->getCellChild(0)->addCellInterface(m_cellInterfacesChildren[0]);
        }
        else {
          //Cell interface number 2 (droite)
          m_cellInterfacesChildren[0]->initializeGauche(cellRef->getCellChild(1));
          m_cellInterfacesChildren[0]->initializeDroite(m_cellRight);
          cellRef->getCellChild(1)->addCellInterface(m_cellInterfacesChildren[0]);
          m_cellRight->addCellInterface(m_cellInterfacesChildren[0]);
        }
        m_cellInterfacesChildren[0]->associeModel(m_mod);
        m_cellInterfacesChildren[0]->allocateSlopes(cellRef->getNumberPhases(), cellRef->getNumberTransports(), allocateSlopeLocal);
      }

      //Cell interface deja split -> on met seulement a jour les liaisons cells/cell interfaces
      //---------------------------------------------------------------------------------------
      else {
        if (m_face->getPos().getX() < cellRef->getElement()->getPosition().getX()) {
          //Cell interface number 1 (gauche)
          m_cellRight = cellRef->getCellChild(0);
          cellRef->getCellChild(0)->addCellInterface(this);
        }
        else {
          //Cell interface number 2 (droite)
          m_cellLeft = cellRef->getCellChild(1);
          cellRef->getCellChild(1)->addCellInterface(this);
        }
      }
    }
    else {
      //--------------------------------------------------
      //--------------------- Cas 2D ---------------------
      //--------------------------------------------------

      //Cell interface pas encore split -> Creation des children cell interfaces
      //------------------------------------------------------------------------
      if (m_lvl == cellRef->getLvl()) {

        //Creation des cell interfaces et faces enfants avec premiere initialization
        //--------------------------------------------------------------------------
        for (int i = 0; i < 2; i++) {
          this->creerCellInterfaceChild();
          m_cellInterfacesChildren[i]->m_face = m_face->creerNouvelleFace();
          m_cellInterfacesChildren[i]->m_face->initializeAutres(surfaceChild, m_face->getNormal(), m_face->getTangent(), m_face->getBinormal());
        }

        //Face selon X
        //------------
        if (std::fabs(m_face->getNormal().getX()) > epsilon) {
          if (m_face->getPos().getX() < cellRef->getElement()->getPosition().getX()) {
            for (int i = 0; i < 2; i++) {
              //First face
              if (i == 0) {
                m_cellInterfacesChildren[i]->m_face->setPos(m_face->getPos().getX(), m_face->getPos().getY() - 0.25*dYParent, m_face->getPos().getZ());
                m_cellInterfacesChildren[i]->initializeGauche(m_cellLeft);
                m_cellInterfacesChildren[i]->initializeDroite(cellRef->getCellChild(0));
                m_cellLeft->addCellInterface(m_cellInterfacesChildren[i]);
                cellRef->getCellChild(0)->addCellInterface(m_cellInterfacesChildren[i]);
                if(m_mod->getSymType()=="CYLINDRICAL")
                {
                  m_cellInterfacesChildren[i]->getFace()->setSurface(0.5*m_face->getSizeY()*m_face->getSizeZ()*(m_face->getPos().getY() - 0.25*dYParent));
                }
              }
              //Second face
              else {
                m_cellInterfacesChildren[i]->m_face->setPos(m_face->getPos().getX(), m_face->getPos().getY() + 0.25*dYParent, m_face->getPos().getZ());
                m_cellInterfacesChildren[i]->initializeGauche(m_cellLeft);
                m_cellInterfacesChildren[i]->initializeDroite(cellRef->getCellChild(2));
                m_cellLeft->addCellInterface(m_cellInterfacesChildren[i]);
                cellRef->getCellChild(2)->addCellInterface(m_cellInterfacesChildren[i]);
                if(m_mod->getSymType()=="CYLINDRICAL")
                  m_cellInterfacesChildren[i]->getFace()->setSurface(0.5*m_face->getSizeY()*m_face->getSizeZ()*(m_face->getPos().getY() + 0.25*dYParent));
              }
              m_cellInterfacesChildren[i]->m_face->setSize(0., 0.5*m_face->getSizeY(), m_face->getSizeZ());
            }
          }
          //Cote droite, right side.
          else {
            for (int i = 0; i < 2; i++) {
              //First face
              if (i == 0) {
                m_cellInterfacesChildren[i]->m_face->setPos(m_face->getPos().getX(), m_face->getPos().getY() - 0.25*dYParent, m_face->getPos().getZ());
                m_cellInterfacesChildren[i]->initializeGauche(cellRef->getCellChild(1));
                m_cellInterfacesChildren[i]->initializeDroite(m_cellRight);
                cellRef->getCellChild(1)->addCellInterface(m_cellInterfacesChildren[i]);
                m_cellRight->addCellInterface(m_cellInterfacesChildren[i]);
                if(m_mod->getSymType()=="CYLINDRICAL")
                  m_cellInterfacesChildren[i]->getFace()->setSurface(0.5*m_face->getSizeY()*m_face->getSizeZ()*(m_face->getPos().getY() - 0.25*dYParent));
              }
              //Second face
              else {
                m_cellInterfacesChildren[i]->m_face->setPos(m_face->getPos().getX(), m_face->getPos().getY() + 0.25*dYParent, m_face->getPos().getZ());
                m_cellInterfacesChildren[i]->initializeGauche(cellRef->getCellChild(3));
                m_cellInterfacesChildren[i]->initializeDroite(m_cellRight);
                cellRef->getCellChild(3)->addCellInterface(m_cellInterfacesChildren[i]);
                m_cellRight->addCellInterface(m_cellInterfacesChildren[i]);
                if(m_mod->getSymType()=="CYLINDRICAL")
                  m_cellInterfacesChildren[i]->getFace()->setSurface(0.5*m_face->getSizeY()*m_face->getSizeZ()*(m_face->getPos().getY() + 0.25*dYParent));
              }
              m_cellInterfacesChildren[i]->m_face->setSize(0., 0.5*m_face->getSizeY(), m_face->getSizeZ());
            }
          }
        }

        //Face selon Y
        //------------
        else {
          //Cote bas
          if (m_face->getPos().getY() < cellRef->getElement()->getPosition().getY()) {
            for (int i = 0; i < 2; i++) {
              //First face
              if (i == 0) {
                m_cellInterfacesChildren[i]->m_face->setPos(m_face->getPos().getX() - 0.25*dXParent, m_face->getPos().getY(), m_face->getPos().getZ());
                m_cellInterfacesChildren[i]->initializeGauche(m_cellLeft);
                m_cellInterfacesChildren[i]->initializeDroite(cellRef->getCellChild(0));
                m_cellLeft->addCellInterface(m_cellInterfacesChildren[i]);
                cellRef->getCellChild(0)->addCellInterface(m_cellInterfacesChildren[i]);
                if(m_mod->getSymType()=="CYLINDRICAL")
                  m_cellInterfacesChildren[i]->getFace()->setSurface(0.5*m_face->getSizeX()*m_face->getSizeZ()*m_face->getPos().getY());
              }
              //Second face
              else {
                m_cellInterfacesChildren[i]->m_face->setPos(m_face->getPos().getX() + 0.25*dXParent, m_face->getPos().getY(), m_face->getPos().getZ());
                m_cellInterfacesChildren[i]->initializeGauche(m_cellLeft);
                m_cellInterfacesChildren[i]->initializeDroite(cellRef->getCellChild(1));
                m_cellLeft->addCellInterface(m_cellInterfacesChildren[i]);
                cellRef->getCellChild(1)->addCellInterface(m_cellInterfacesChildren[i]);
                if(m_mod->getSymType()=="CYLINDRICAL")
                  m_cellInterfacesChildren[i]->getFace()->setSurface(0.5*m_face->getSizeX()*m_face->getSizeZ()*m_face->getPos().getY());
              }
              m_cellInterfacesChildren[i]->m_face->setSize(0.5*m_face->getSizeX(), 0., m_face->getSizeZ());
            }
          }
          //Cote haut
          else {
            for (int i = 0; i < 2; i++) {
              //First face
              if (i == 0) {
                m_cellInterfacesChildren[i]->m_face->setPos(m_face->getPos().getX() - 0.25*dXParent, m_face->getPos().getY(), m_face->getPos().getZ());
                m_cellInterfacesChildren[i]->initializeGauche(cellRef->getCellChild(2));
                m_cellInterfacesChildren[i]->initializeDroite(m_cellRight);
                cellRef->getCellChild(2)->addCellInterface(m_cellInterfacesChildren[i]);
                m_cellRight->addCellInterface(m_cellInterfacesChildren[i]);
                if(m_mod->getSymType()=="CYLINDRICAL")
                  m_cellInterfacesChildren[i]->getFace()->setSurface(0.5*m_face->getSizeX()*m_face->getSizeZ()*m_face->getPos().getY());
              }
              //Second face
              else {
                m_cellInterfacesChildren[i]->m_face->setPos(m_face->getPos().getX() + 0.25*dXParent, m_face->getPos().getY(), m_face->getPos().getZ());
                m_cellInterfacesChildren[i]->initializeGauche(cellRef->getCellChild(3));
                m_cellInterfacesChildren[i]->initializeDroite(m_cellRight);
                cellRef->getCellChild(3)->addCellInterface(m_cellInterfacesChildren[i]);
                m_cellRight->addCellInterface(m_cellInterfacesChildren[i]);
                if(m_mod->getSymType()=="CYLINDRICAL")
                  m_cellInterfacesChildren[i]->getFace()->setSurface(0.5*m_face->getSizeX()*m_face->getSizeZ()*m_face->getPos().getY());
              }
              m_cellInterfacesChildren[i]->m_face->setSize(0.5*m_face->getSizeX(), 0., m_face->getSizeZ());
            }
          }
        }

        //Association du model et des slopes
        //----------------------------------
        for (int i = 0; i < 2; i++) {
          m_cellInterfacesChildren[i]->associeModel(m_mod);
          m_cellInterfacesChildren[i]->allocateSlopes(cellRef->getNumberPhases(), cellRef->getNumberTransports(), allocateSlopeLocal);
        }

      }

      //Cell interface deja split -> on met seulement a jour les liaisons cells/cell interfaces
      // we only update the cells / cell interfaces links
      //---------------------------------------------------------------------------------------
      else {

        //Face selon X
        //------------
        if (std::fabs(m_face->getNormal().getX()) > epsilon) {
          //Cote gauche
          if (m_face->getPos().getX() < cellRef->getElement()->getPosition().getX()) {
            //First face
            if (m_face->getPos().getY() < cellRef->getElement()->getPosition().getY()) {
              m_cellRight = cellRef->getCellChild(0);
              cellRef->getCellChild(0)->addCellInterface(this);
            }
            //Second face
            else {
              m_cellRight = cellRef->getCellChild(2);
              cellRef->getCellChild(2)->addCellInterface(this);
            }
          }
          //Cote droite
          else {
            //First face
            if (m_face->getPos().getY() < cellRef->getElement()->getPosition().getY()) {
              m_cellLeft = cellRef->getCellChild(1);
              cellRef->getCellChild(1)->addCellInterface(this);
            }
            //Second face
            else {
              m_cellLeft = cellRef->getCellChild(3);
              cellRef->getCellChild(3)->addCellInterface(this);
            }
          }
        }

        //Face selon Y
        //------------
        else {
          //Cote bas
          if (m_face->getPos().getY() < cellRef->getElement()->getPosition().getY()) {
            //First face
            if (m_face->getPos().getX() < cellRef->getElement()->getPosition().getX()) {
              m_cellRight = cellRef->getCellChild(0);
              cellRef->getCellChild(0)->addCellInterface(this);
            }
            //Second face
            else {
              m_cellRight = cellRef->getCellChild(1);
              cellRef->getCellChild(1)->addCellInterface(this);
            }
          }
          //Cote haut
          else {
            //First face
            if (m_face->getPos().getX() < cellRef->getElement()->getPosition().getX()) {
              m_cellLeft = cellRef->getCellChild(2);
              cellRef->getCellChild(2)->addCellInterface(this);
            }
            //Second face
            else {
              m_cellLeft = cellRef->getCellChild(3);
              cellRef->getCellChild(3)->addCellInterface(this);
            }
          }
        }
      }
    }
  }
  else {

    //--------------------------------------------------
    //--------------------- Cas 3D ---------------------
    //--------------------------------------------------

    //Cell interface pas encore split -> Creation des children cell interfaces
    //------------------------------------------------------------------------
    if (m_lvl == cellRef->getLvl()) {

      //Creation des cell interfaces et faces enfants avec premiere initialization
      //--------------------------------------------------------------------------
      for (int i = 0; i < 4; i++) {
        this->creerCellInterfaceChild();
        m_cellInterfacesChildren[i]->m_face = m_face->creerNouvelleFace();
        m_cellInterfacesChildren[i]->m_face->initializeAutres(surfaceChild, m_face->getNormal(), m_face->getTangent(), m_face->getBinormal());
        m_cellInterfacesChildren[i]->m_face->setSize(0.5*m_face->getSize());
      }

      //Face selon X
      //------------
      if (std::fabs(m_face->getNormal().getX()) > epsilon) {
        //Cote gauche
        if (m_face->getPos().getX() < cellRef->getElement()->getPosition().getX()) {
          for (int i = 0; i < 4; i++) {
            //First face
            if (i == 0) {
              m_cellInterfacesChildren[i]->m_face->setPos(m_face->getPos().getX(), m_face->getPos().getY() - 0.25*dYParent, m_face->getPos().getZ() - 0.25*dZParent);
              m_cellInterfacesChildren[i]->initializeGauche(m_cellLeft);
              m_cellInterfacesChildren[i]->initializeDroite(cellRef->getCellChild(0));
              m_cellLeft->addCellInterface(m_cellInterfacesChildren[i]);
              cellRef->getCellChild(0)->addCellInterface(m_cellInterfacesChildren[i]);
            }
            //Second face
            else if (i == 1) {
              m_cellInterfacesChildren[i]->m_face->setPos(m_face->getPos().getX(), m_face->getPos().getY() - 0.25*dYParent, m_face->getPos().getZ() + 0.25*dZParent);
              m_cellInterfacesChildren[i]->initializeGauche(m_cellLeft);
              m_cellInterfacesChildren[i]->initializeDroite(cellRef->getCellChild(4));
              m_cellLeft->addCellInterface(m_cellInterfacesChildren[i]);
              cellRef->getCellChild(4)->addCellInterface(m_cellInterfacesChildren[i]);
            }
            //Third face
            else if (i == 2) {
              m_cellInterfacesChildren[i]->m_face->setPos(m_face->getPos().getX(), m_face->getPos().getY() + 0.25*dYParent, m_face->getPos().getZ() - 0.25*dZParent);
              m_cellInterfacesChildren[i]->initializeGauche(m_cellLeft);
              m_cellInterfacesChildren[i]->initializeDroite(cellRef->getCellChild(2));
              m_cellLeft->addCellInterface(m_cellInterfacesChildren[i]);
              cellRef->getCellChild(2)->addCellInterface(m_cellInterfacesChildren[i]);
            }
            //Fourth face
            else {
              m_cellInterfacesChildren[i]->m_face->setPos(m_face->getPos().getX(), m_face->getPos().getY() + 0.25*dYParent, m_face->getPos().getZ() + 0.25*dZParent);
              m_cellInterfacesChildren[i]->initializeGauche(m_cellLeft);
              m_cellInterfacesChildren[i]->initializeDroite(cellRef->getCellChild(6));
              m_cellLeft->addCellInterface(m_cellInterfacesChildren[i]);
              cellRef->getCellChild(6)->addCellInterface(m_cellInterfacesChildren[i]);
            }
          }
        }
        //Cote droite
        else {
          for (int i = 0; i < 4; i++) {
            //First face
            if (i == 0) {
              m_cellInterfacesChildren[i]->m_face->setPos(m_face->getPos().getX(), m_face->getPos().getY() - 0.25*dYParent, m_face->getPos().getZ() - 0.25*dZParent);
              m_cellInterfacesChildren[i]->initializeGauche(cellRef->getCellChild(1));
              m_cellInterfacesChildren[i]->initializeDroite(m_cellRight);
              cellRef->getCellChild(1)->addCellInterface(m_cellInterfacesChildren[i]);
              m_cellRight->addCellInterface(m_cellInterfacesChildren[i]);
            }
            //Second face
            else if (i == 1) {
              m_cellInterfacesChildren[i]->m_face->setPos(m_face->getPos().getX(), m_face->getPos().getY() - 0.25*dYParent, m_face->getPos().getZ() + 0.25*dZParent);
              m_cellInterfacesChildren[i]->initializeGauche(cellRef->getCellChild(5));
              m_cellInterfacesChildren[i]->initializeDroite(m_cellRight);
              cellRef->getCellChild(5)->addCellInterface(m_cellInterfacesChildren[i]);
              m_cellRight->addCellInterface(m_cellInterfacesChildren[i]);
            }
            //Third face
            else if (i == 2) {
              m_cellInterfacesChildren[i]->m_face->setPos(m_face->getPos().getX(), m_face->getPos().getY() + 0.25*dYParent, m_face->getPos().getZ() - 0.25*dZParent);
              m_cellInterfacesChildren[i]->initializeGauche(cellRef->getCellChild(3));
              m_cellInterfacesChildren[i]->initializeDroite(m_cellRight);
              cellRef->getCellChild(3)->addCellInterface(m_cellInterfacesChildren[i]);
              m_cellRight->addCellInterface(m_cellInterfacesChildren[i]);
            }
            //Fourth face
            else {
              m_cellInterfacesChildren[i]->m_face->setPos(m_face->getPos().getX(), m_face->getPos().getY() + 0.25*dYParent, m_face->getPos().getZ() + 0.25*dZParent);
              m_cellInterfacesChildren[i]->initializeGauche(cellRef->getCellChild(7));
              m_cellInterfacesChildren[i]->initializeDroite(m_cellRight);
              cellRef->getCellChild(7)->addCellInterface(m_cellInterfacesChildren[i]);
              m_cellRight->addCellInterface(m_cellInterfacesChildren[i]);
            }
          }
        }
      }

      //Face selon Y
      //------------
      else if (std::fabs(m_face->getNormal().getY()) > epsilon) {
        //Cote bas
        if (m_face->getPos().getY() < cellRef->getElement()->getPosition().getY()) {
          for (int i = 0; i < 4; i++) {
            //First face
            if (i == 0) {
              m_cellInterfacesChildren[i]->m_face->setPos(m_face->getPos().getX() - 0.25*dXParent, m_face->getPos().getY(), m_face->getPos().getZ() - 0.25*dZParent);
              m_cellInterfacesChildren[i]->initializeGauche(m_cellLeft);
              m_cellInterfacesChildren[i]->initializeDroite(cellRef->getCellChild(0));
              m_cellLeft->addCellInterface(m_cellInterfacesChildren[i]);
              cellRef->getCellChild(0)->addCellInterface(m_cellInterfacesChildren[i]);
            }
            //Second face
            else if (i == 1) {
              m_cellInterfacesChildren[i]->m_face->setPos(m_face->getPos().getX() + 0.25*dXParent, m_face->getPos().getY(), m_face->getPos().getZ() - 0.25*dZParent);
              m_cellInterfacesChildren[i]->initializeGauche(m_cellLeft);
              m_cellInterfacesChildren[i]->initializeDroite(cellRef->getCellChild(1));
              m_cellLeft->addCellInterface(m_cellInterfacesChildren[i]);
              cellRef->getCellChild(1)->addCellInterface(m_cellInterfacesChildren[i]);
            }
            //Third face
            else if (i == 2) {
              m_cellInterfacesChildren[i]->m_face->setPos(m_face->getPos().getX() - 0.25*dXParent, m_face->getPos().getY(), m_face->getPos().getZ() + 0.25*dZParent);
              m_cellInterfacesChildren[i]->initializeGauche(m_cellLeft);
              m_cellInterfacesChildren[i]->initializeDroite(cellRef->getCellChild(4));
              m_cellLeft->addCellInterface(m_cellInterfacesChildren[i]);
              cellRef->getCellChild(4)->addCellInterface(m_cellInterfacesChildren[i]);
            }
            //Fourth face
            else {
              m_cellInterfacesChildren[i]->m_face->setPos(m_face->getPos().getX() + 0.25*dXParent, m_face->getPos().getY(), m_face->getPos().getZ() + 0.25*dZParent);
              m_cellInterfacesChildren[i]->initializeGauche(m_cellLeft);
              m_cellInterfacesChildren[i]->initializeDroite(cellRef->getCellChild(5));
              m_cellLeft->addCellInterface(m_cellInterfacesChildren[i]);
              cellRef->getCellChild(5)->addCellInterface(m_cellInterfacesChildren[i]);
            }
          }
        }
        //Cote haut
        else {
          for (int i = 0; i < 4; i++) {
            //First face
            if (i == 0) {
              m_cellInterfacesChildren[i]->m_face->setPos(m_face->getPos().getX() - 0.25*dXParent, m_face->getPos().getY(), m_face->getPos().getZ() - 0.25*dZParent);
              m_cellInterfacesChildren[i]->initializeGauche(cellRef->getCellChild(2));
              m_cellInterfacesChildren[i]->initializeDroite(m_cellRight);
              cellRef->getCellChild(2)->addCellInterface(m_cellInterfacesChildren[i]);
              m_cellRight->addCellInterface(m_cellInterfacesChildren[i]);
            }
            //Second face
            else if (i == 1) {
              m_cellInterfacesChildren[i]->m_face->setPos(m_face->getPos().getX() + 0.25*dXParent, m_face->getPos().getY(), m_face->getPos().getZ() - 0.25*dZParent);
              m_cellInterfacesChildren[i]->initializeGauche(cellRef->getCellChild(3));
              m_cellInterfacesChildren[i]->initializeDroite(m_cellRight);
              cellRef->getCellChild(3)->addCellInterface(m_cellInterfacesChildren[i]);
              m_cellRight->addCellInterface(m_cellInterfacesChildren[i]);
            }
            //Third face
            else if (i == 2) {
              m_cellInterfacesChildren[i]->m_face->setPos(m_face->getPos().getX() - 0.25*dXParent, m_face->getPos().getY(), m_face->getPos().getZ() + 0.25*dZParent);
              m_cellInterfacesChildren[i]->initializeGauche(cellRef->getCellChild(6));
              m_cellInterfacesChildren[i]->initializeDroite(m_cellRight);
              cellRef->getCellChild(6)->addCellInterface(m_cellInterfacesChildren[i]);
              m_cellRight->addCellInterface(m_cellInterfacesChildren[i]);
            }
            //Fourth face
            else {
              m_cellInterfacesChildren[i]->m_face->setPos(m_face->getPos().getX() + 0.25*dXParent, m_face->getPos().getY(), m_face->getPos().getZ() + 0.25*dZParent);
              m_cellInterfacesChildren[i]->initializeGauche(cellRef->getCellChild(7));
              m_cellInterfacesChildren[i]->initializeDroite(m_cellRight);
              cellRef->getCellChild(7)->addCellInterface(m_cellInterfacesChildren[i]);
              m_cellRight->addCellInterface(m_cellInterfacesChildren[i]);
            }
          }
        }
      }

      //Face selon Z
      //------------
      else {
        //Cote devant
        if (m_face->getPos().getZ() < cellRef->getElement()->getPosition().getZ()) {
          for (int i = 0; i < 4; i++) {
            //First face
            if (i == 0) {
              m_cellInterfacesChildren[i]->m_face->setPos(m_face->getPos().getX() - 0.25*dXParent, m_face->getPos().getY() - 0.25*dYParent, m_face->getPos().getZ());
              m_cellInterfacesChildren[i]->initializeGauche(m_cellLeft);
              m_cellInterfacesChildren[i]->initializeDroite(cellRef->getCellChild(0));
              m_cellLeft->addCellInterface(m_cellInterfacesChildren[i]);
              cellRef->getCellChild(0)->addCellInterface(m_cellInterfacesChildren[i]);
            }
            //Second face
            else if (i == 1) {
              m_cellInterfacesChildren[i]->m_face->setPos(m_face->getPos().getX() + 0.25*dXParent, m_face->getPos().getY() - 0.25*dYParent, m_face->getPos().getZ());
              m_cellInterfacesChildren[i]->initializeGauche(m_cellLeft);
              m_cellInterfacesChildren[i]->initializeDroite(cellRef->getCellChild(1));
              m_cellLeft->addCellInterface(m_cellInterfacesChildren[i]);
              cellRef->getCellChild(1)->addCellInterface(m_cellInterfacesChildren[i]);
            }
            //Third face
            else if (i == 2) {
              m_cellInterfacesChildren[i]->m_face->setPos(m_face->getPos().getX() - 0.25*dXParent, m_face->getPos().getY() + 0.25*dYParent, m_face->getPos().getZ());
              m_cellInterfacesChildren[i]->initializeGauche(m_cellLeft);
              m_cellInterfacesChildren[i]->initializeDroite(cellRef->getCellChild(2));
              m_cellLeft->addCellInterface(m_cellInterfacesChildren[i]);
              cellRef->getCellChild(2)->addCellInterface(m_cellInterfacesChildren[i]);
            }
            //Fourth face
            else {
              m_cellInterfacesChildren[i]->m_face->setPos(m_face->getPos().getX() + 0.25*dXParent, m_face->getPos().getY() + 0.25*dYParent, m_face->getPos().getZ());
              m_cellInterfacesChildren[i]->initializeGauche(m_cellLeft);
              m_cellInterfacesChildren[i]->initializeDroite(cellRef->getCellChild(3));
              m_cellLeft->addCellInterface(m_cellInterfacesChildren[i]);
              cellRef->getCellChild(3)->addCellInterface(m_cellInterfacesChildren[i]);
            }
          }
        }
        //Cote derriere
        else {
          for (int i = 0; i < 4; i++) {
            //First face
            if (i == 0) {
              m_cellInterfacesChildren[i]->m_face->setPos(m_face->getPos().getX() - 0.25*dXParent, m_face->getPos().getY() - 0.25*dYParent, m_face->getPos().getZ());
              m_cellInterfacesChildren[i]->initializeGauche(cellRef->getCellChild(4));
              m_cellInterfacesChildren[i]->initializeDroite(m_cellRight);
              cellRef->getCellChild(4)->addCellInterface(m_cellInterfacesChildren[i]);
              m_cellRight->addCellInterface(m_cellInterfacesChildren[i]);
            }
            //Second face
            else if (i == 1) {
              m_cellInterfacesChildren[i]->m_face->setPos(m_face->getPos().getX() + 0.25*dXParent, m_face->getPos().getY() - 0.25*dYParent, m_face->getPos().getZ());
              m_cellInterfacesChildren[i]->initializeGauche(cellRef->getCellChild(5));
              m_cellInterfacesChildren[i]->initializeDroite(m_cellRight);
              cellRef->getCellChild(5)->addCellInterface(m_cellInterfacesChildren[i]);
              m_cellRight->addCellInterface(m_cellInterfacesChildren[i]);
            }
            //Third face
            else if (i == 2) {
              m_cellInterfacesChildren[i]->m_face->setPos(m_face->getPos().getX() - 0.25*dXParent, m_face->getPos().getY() + 0.25*dYParent, m_face->getPos().getZ());
              m_cellInterfacesChildren[i]->initializeGauche(cellRef->getCellChild(6));
              m_cellInterfacesChildren[i]->initializeDroite(m_cellRight);
              cellRef->getCellChild(6)->addCellInterface(m_cellInterfacesChildren[i]);
              m_cellRight->addCellInterface(m_cellInterfacesChildren[i]);
            }
            //Fourth face
            else {
              m_cellInterfacesChildren[i]->m_face->setPos(m_face->getPos().getX() + 0.25*dXParent, m_face->getPos().getY() + 0.25*dYParent, m_face->getPos().getZ());
              m_cellInterfacesChildren[i]->initializeGauche(cellRef->getCellChild(7));
              m_cellInterfacesChildren[i]->initializeDroite(m_cellRight);
              cellRef->getCellChild(7)->addCellInterface(m_cellInterfacesChildren[i]);
              m_cellRight->addCellInterface(m_cellInterfacesChildren[i]);
            }
          }
        }
      }

      //Association du model et des slopes
      //----------------------------------
      for (int i = 0; i < 4; i++) {
        m_cellInterfacesChildren[i]->associeModel(m_mod);
        m_cellInterfacesChildren[i]->allocateSlopes(cellRef->getNumberPhases(), cellRef->getNumberTransports(), allocateSlopeLocal);
      }

    }

    //Cell interface deja split -> on met seulement a jour les liaisons cells/cell interfaces
    //---------------------------------------------------------------------------------------
    else {

      //Face selon X
      //------------
      if (std::fabs(std::fabs(m_face->getNormal().getX()) - 1.) < epsilon) {
        //Cote gauche
        if (m_face->getPos().getX() < cellRef->getElement()->getPosition().getX()) {
          if (m_face->getPos().getY() < cellRef->getElement()->getPosition().getY()) {
            if (m_face->getPos().getZ() < cellRef->getElement()->getPosition().getZ()) {
              //First face
              m_cellRight = cellRef->getCellChild(0);
              cellRef->getCellChild(0)->addCellInterface(this);
            }
            else {
              //Second face
              m_cellRight = cellRef->getCellChild(4);
              cellRef->getCellChild(4)->addCellInterface(this);
            }
          }
          else {
            if (m_face->getPos().getZ() < cellRef->getElement()->getPosition().getZ()) {
              //Third face
              m_cellRight = cellRef->getCellChild(2);
              cellRef->getCellChild(2)->addCellInterface(this);
            }
            else {
              //Fourth face
              m_cellRight = cellRef->getCellChild(6);
              cellRef->getCellChild(6)->addCellInterface(this);
            }
          }
        }
        //Cote droite
        else {
          if (m_face->getPos().getY() < cellRef->getElement()->getPosition().getY()) {
            if (m_face->getPos().getZ() < cellRef->getElement()->getPosition().getZ()) {
              //First face
              m_cellLeft = cellRef->getCellChild(1);
              cellRef->getCellChild(1)->addCellInterface(this);
            }
            else {
              //Second face
              m_cellLeft = cellRef->getCellChild(5);
              cellRef->getCellChild(5)->addCellInterface(this);
            }
          }
          else {
            if (m_face->getPos().getZ() < cellRef->getElement()->getPosition().getZ()) {
              //Third face
              m_cellLeft = cellRef->getCellChild(3);
              cellRef->getCellChild(3)->addCellInterface(this);
            }
            else {
              //Fourth face
              m_cellLeft = cellRef->getCellChild(7);
              cellRef->getCellChild(7)->addCellInterface(this);
            }
          }
        }
      }

      //Face selon Y
      //------------
      else if (std::fabs(std::fabs(m_face->getNormal().getY()) - 1.) < epsilon) {
        //Cote bas
        if (m_face->getPos().getY() < cellRef->getElement()->getPosition().getY()) {
          if (m_face->getPos().getZ() < cellRef->getElement()->getPosition().getZ()) {
            if (m_face->getPos().getX() < cellRef->getElement()->getPosition().getX()) {
              //First face
              m_cellRight = cellRef->getCellChild(0);
              cellRef->getCellChild(0)->addCellInterface(this);
            }
            else {
              //Second face
              m_cellRight = cellRef->getCellChild(1);
              cellRef->getCellChild(1)->addCellInterface(this);
            }
          }
          else {
            if (m_face->getPos().getX() < cellRef->getElement()->getPosition().getX()) {
              //Third face
              m_cellRight = cellRef->getCellChild(4);
              cellRef->getCellChild(4)->addCellInterface(this);
            }
            else {
              //Fourth face
              m_cellRight = cellRef->getCellChild(5);
              cellRef->getCellChild(5)->addCellInterface(this);
            }
          }
        }
        //Cote haut
        else {
          if (m_face->getPos().getZ() < cellRef->getElement()->getPosition().getZ()) {
            if (m_face->getPos().getX() < cellRef->getElement()->getPosition().getX()) {
              //First face
              m_cellLeft = cellRef->getCellChild(2);
              cellRef->getCellChild(2)->addCellInterface(this);
            }
            else {
              //Second face
              m_cellLeft = cellRef->getCellChild(3);
              cellRef->getCellChild(3)->addCellInterface(this);
            }
          }
          else {
            if (m_face->getPos().getX() < cellRef->getElement()->getPosition().getX()) {
              //Third face
              m_cellLeft = cellRef->getCellChild(6);
              cellRef->getCellChild(6)->addCellInterface(this);
            }
            else {
              //Fourth face
              m_cellLeft = cellRef->getCellChild(7);
              cellRef->getCellChild(7)->addCellInterface(this);
            }
          }
        }
      }

      //Face selon Z
      //------------
      else {
        //Cote devant
        if (m_face->getPos().getZ() < cellRef->getElement()->getPosition().getZ()) {
          if (m_face->getPos().getY() < cellRef->getElement()->getPosition().getY()) {
            if (m_face->getPos().getX() < cellRef->getElement()->getPosition().getX()) {
              //First face
              m_cellRight = cellRef->getCellChild(0);
              cellRef->getCellChild(0)->addCellInterface(this);
            }
            else {
              //Second face
              m_cellRight = cellRef->getCellChild(1);
              cellRef->getCellChild(1)->addCellInterface(this);
            }
          }
          else {
            if (m_face->getPos().getX() < cellRef->getElement()->getPosition().getX()) {
              //Third face
              m_cellRight = cellRef->getCellChild(2);
              cellRef->getCellChild(2)->addCellInterface(this);
            }
            else {
              //Fourth face
              m_cellRight = cellRef->getCellChild(3);
              cellRef->getCellChild(3)->addCellInterface(this);
            }
          }
        }
        //Cote derriere
        else {
          if (m_face->getPos().getY() < cellRef->getElement()->getPosition().getY()) {
            if (m_face->getPos().getX() < cellRef->getElement()->getPosition().getX()) {
              //First face
              m_cellLeft = cellRef->getCellChild(4);
              cellRef->getCellChild(4)->addCellInterface(this);
            }
            else {
              //Second face
              m_cellLeft = cellRef->getCellChild(5);
              cellRef->getCellChild(5)->addCellInterface(this);
            }
          }
          else {
            if (m_face->getPos().getX() < cellRef->getElement()->getPosition().getX()) {
              //Third face
              m_cellLeft = cellRef->getCellChild(6);
              cellRef->getCellChild(6)->addCellInterface(this);
            }
            else {
              //Fourth face
              m_cellLeft = cellRef->getCellChild(7);
              cellRef->getCellChild(7)->addCellInterface(this);
            }
          }
        }
      }
    }
  }
}

//***********************************************************************

void CellInterface::deraffineCellInterfaceExterne(Cell *cellRef)
{
  //On parcourt seulement les parent cell interfaces pour regarder si la cell voisine de celle de reference a des enfants,
  //si oui (enfants), on ne peut pas deraffiner le cell interface, on reaffecte donc les liaisons cells/cell interfaces des children cell interfaces,
  //si non (pas enfants), on peut deraffiner le cell interface et mettre a jour les neighbouring cells.
  //Plus, si je suis un parent cell interface, qui a donc des enfants, mais que la cell de reference ne les connait pas encore, on les ajoute a ses cell interfaces.

  //Parcours les parent cell interfaces
  if (cellRef->getLvl() == m_lvl) {
    //cellRef est la cell gauche
    if (m_cellLeft == cellRef) {
      //Si la cell voisine (droite) a des enfants, on reaffecte les liaisons cells/cell interfaces des children cell interfaces
      if (m_cellRight->getSplit()) {
        for (unsigned int cellInterfaceChild = 0; cellInterfaceChild < m_cellInterfacesChildren.size(); cellInterfaceChild++) {
          m_cellInterfacesChildren[cellInterfaceChild]->initializeGauche(cellRef);
        }
      }
      //La cell voisine (droite) n'a pas d'enfants, on deraffine le cell interface et met a jour les neighbouring cells
      else {
        //Il faut aussi enlever ses children cell interfaces des cell interfaces de la cell droite et de mes cell interfaces
        for (unsigned int cellInterfaceChild = 0; cellInterfaceChild < m_cellInterfacesChildren.size(); cellInterfaceChild++) {
          m_cellRight->deleteCellInterface(m_cellInterfacesChildren[cellInterfaceChild]);
          cellRef->deleteCellInterface(m_cellInterfacesChildren[cellInterfaceChild]);
        }
        this->deraffineCellInterfacesChildren();
      }
    }
    //cellRef est la cell droite
    else {
      //Si la cell voisine (gauche) a des enfants, on reaffecte les liaisons cells/cell interfaces des children cell interfaces
      if (m_cellLeft->getSplit()) {
        for (unsigned int cellInterfaceChild = 0; cellInterfaceChild < m_cellInterfacesChildren.size(); cellInterfaceChild++) {
          m_cellInterfacesChildren[cellInterfaceChild]->initializeDroite(cellRef);
        }
      }
      //La cell voisine (gauche) n'a pas d'enfants, on deraffine le cell interface et met a jour les neighbouring cells
      else {
        //Il faut aussi enlever ses children cell interfaces des cell interfaces de la cell gauche et de mes cell interfaces
        for (unsigned int cellInterfaceChild = 0; cellInterfaceChild < m_cellInterfacesChildren.size(); cellInterfaceChild++) {
          m_cellLeft->deleteCellInterface(m_cellInterfacesChildren[cellInterfaceChild]);
          cellRef->deleteCellInterface(m_cellInterfacesChildren[cellInterfaceChild]);
        }
        this->deraffineCellInterfacesChildren();
      }
    }

    //Si je suis un cell interface parent, qui a donc des enfants, mais que la cell de reference ne les connait pas encore, on les ajoute a ses cell interfaces.
    for (unsigned int cellInterfaceChild = 0; cellInterfaceChild < m_cellInterfacesChildren.size(); cellInterfaceChild++) {
      bool ajoutChildAuxCellInterfacesCellRef(true);
      for (int i = 0; i < cellRef->getCellInterfacesSize(); i++) {
        if (cellRef->getCellInterface(i) == m_cellInterfacesChildren[cellInterfaceChild]) { ajoutChildAuxCellInterfacesCellRef = false; break; }
      }
      if (ajoutChildAuxCellInterfacesCellRef) { cellRef->addCellInterface(m_cellInterfacesChildren[cellInterfaceChild]); }
    }
  }
}

//***********************************************************************

void CellInterface::deraffineCellInterfacesChildren()
{
  for (unsigned int i = 0; i < m_cellInterfacesChildren.size(); i++) {
    delete m_cellInterfacesChildren[i];
  }
  m_cellInterfacesChildren.clear();
}

//***********************************************************************

void CellInterface::constructionTableauCellInterfacesExternesLvl(std::vector<CellInterface *> *cellInterfacesLvl)
{
  for (unsigned int i = 0; i < m_cellInterfacesChildren.size(); i++) {
    cellInterfacesLvl[m_lvl + 1].push_back(m_cellInterfacesChildren[i]);
  }
}

//***********************************************************************

bool CellInterface::getSplit() const
{
  bool split = false;
  if (m_cellInterfacesChildren.size() > 0) { split = true; }
  return split;
}

//***********************************************************************

int CellInterface::getNumberCellInterfacesChildren() const
{
  return m_cellInterfacesChildren.size();
}

//***********************************************************************

CellInterface* CellInterface::getCellInterfaceChild(const int &numChild)
{
  return m_cellInterfacesChildren[numChild];
}

//***********************************************************************

CellInterface* CellInterface::getCellInterfaceChildBack()
{
  return m_cellInterfacesChildren.back();
}

//***************************************************************************

void CellInterface::updatePointersInternalCellInterfaces()
{
  //Check if the cells pointed by me also points to me
  bool foundCellInterface(false);
  //Left
  for (int b = 0; b < m_cellLeft->getCellInterfacesSize(); b++) {
    if (m_cellLeft->getCellInterface(b) == this) {
      foundCellInterface = true; break;
    }
  }
  if (!foundCellInterface) m_cellLeft->addCellInterface(this);

  //Right
  foundCellInterface = false;
  for (int b = 0; b < m_cellRight->getCellInterfacesSize(); b++) {
    if (m_cellRight->getCellInterface(b) == this) {
      foundCellInterface = true; break;
    }
  }
  if (!foundCellInterface) m_cellRight->addCellInterface(this);

  //Also check my children
  for (int b = 0; b < m_cellInterfacesChildren.size(); b++) {
    m_cellInterfacesChildren[b]->updatePointersInternalCellInterfaces();
  }
}

//***********************************************************************
