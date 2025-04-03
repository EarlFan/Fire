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

//! \file      ViscousParallel.cpp
//! \author    E Fan, Lisong Shi, Jiaao Hao, Tianhan Zhang, Chih-yung Wen.
//! \version   1.0
//! \date      Jan 01 2025

#include "../libs/ECOGEN/Parallel/Parallel.h"
#include "../libs/ECOGEN/Eos/Eos.h"

//****************************************************************************
//*********************** Methods for the transport array ****************************
//****************************************************************************

void Parallel::initializePersistentCommunicationsTransportArray(const int &size)
{
  for (int neighbour = 0; neighbour < Ncpu; neighbour++) {
    if (m_isNeighbour[neighbour]) {
      //Determination of the number of variables to communicate, as much variables as the dimension (1,2 or 3)
      int numberSend = (size)*m_numberElementsToSendToNeighbour[neighbour];
      int numberReceive = (size)*m_numberElementsToReceiveFromNeighbour[neighbour];

      //New sending request and its associated buffer
      m_reqSendTransportArray[0][neighbour] = new MPI_Request;
      m_bufferSendTransportArray[0][neighbour] = new double[numberSend];
      MPI_Send_init(m_bufferSendTransportArray[0][neighbour], numberSend, MPI_DOUBLE, neighbour, neighbour, MPI_COMM_WORLD, m_reqSendTransportArray[0][neighbour]);

      //New receiving request and its associated buffer
      m_reqReceiveTransportArray[0][neighbour] = new MPI_Request;
      m_bufferReceiveTransportArray[0][neighbour] = new double[numberReceive];
      MPI_Recv_init(m_bufferReceiveTransportArray[0][neighbour], numberReceive, MPI_DOUBLE, neighbour, rankCpu, MPI_COMM_WORLD, m_reqReceiveTransportArray[0][neighbour]);
    }
  }
}

//***********************************************************************

void Parallel::finalizePersistentCommunicationsTransportArray(const int &lvlMax)
{
  for (int lvl = 0; lvl <= lvlMax; lvl++) {
    for (int neighbour = 0; neighbour < Ncpu; neighbour++)	{
      if (m_isNeighbour[neighbour]) {
        MPI_Request_free(m_reqSendTransportArray[lvl][neighbour]);
        MPI_Request_free(m_reqReceiveTransportArray[lvl][neighbour]);
        delete m_reqSendTransportArray[lvl][neighbour];
        delete[] m_bufferSendTransportArray[lvl][neighbour];
        delete m_reqReceiveTransportArray[lvl][neighbour];
        delete[] m_bufferReceiveTransportArray[lvl][neighbour];
      }
    }
    delete[] m_reqSendTransportArray[lvl];
    delete[] m_bufferSendTransportArray[lvl];
    delete[] m_reqReceiveTransportArray[lvl];
    delete[] m_bufferReceiveTransportArray[lvl];
  }
  m_reqSendTransportArray.clear();
  m_bufferSendTransportArray.clear();
  m_reqReceiveTransportArray.clear();
  m_bufferReceiveTransportArray.clear();
}

//***********************************************************************

void Parallel::communicationsTransportArray(Variable nameVector, const int &dim, int lvl, int num, int index)
{
  int count(0);
  MPI_Status status;

  for (int neighbour = 0; neighbour < Ncpu; neighbour++) {
    if (m_isNeighbour[neighbour]) {
      //Prepation of sendings
      count = -1;
      for (int i = 0; i < m_numberElementsToSendToNeighbour[neighbour]; i++) {
        //Automatic filing of m_bufferSendVector function of gradient coordinates
        m_elementsToSend[neighbour][i]->fillBufferVector(m_bufferSendTransportArray[lvl][neighbour], count, lvl, neighbour, dim, nameVector, num, index);
      }

      //Sending request
      MPI_Start(m_reqSendTransportArray[lvl][neighbour]);

      //Receiving request
      MPI_Start(m_reqReceiveTransportArray[lvl][neighbour]);
    }
  }
  for (int neighbour = 0; neighbour < Ncpu; neighbour++) {
    if (m_isNeighbour[neighbour]) {
      //Waiting
      MPI_Wait(m_reqSendTransportArray[lvl][neighbour], &status);
      MPI_Wait(m_reqReceiveTransportArray[lvl][neighbour], &status);
      //Receivings
      count = -1;
      for (int i = 0; i < m_numberElementsToReceiveFromNeighbour[neighbour]; i++) {
        //Automatic filing of m_bufferReceiveVector function of gradient coordinates
        m_elementsToReceive[neighbour][i]->getBufferVector(m_bufferReceiveTransportArray[lvl][neighbour], count, lvl, dim, nameVector, num, index);
      }
    }
  }
}