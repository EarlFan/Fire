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

//! \file      main.cpp
//! \author    E Fan, Lisong Shi, Jiaao Hao, Tianhan Zhang, Chih-yung Wen.
//! \version   1.0
//! \date      Jan 01 2025

// Updated include paths for new location
#include "libs/ECOGEN/Run.h"
#include "include/Fire.h"
#include "libs/ECOGEN/Errors.h"
#include "libs/ECOGEN/libTierces/tinyxml2.h"

#include <time.h>

using namespace tinyxml2;

std::shared_ptr<Cantera::Solution> solution_mech;
std::shared_ptr<Cantera::ThermoPhase> gas_global;

//***********************************************************************

int main(int argc, char* argv[])
{
  Cantera::suppress_deprecation_warnings();
  
  Run* run(0);

  //Error redirection
  std::ofstream Out("Error.txt");
  std::streambuf* strBackup(std::cerr.rdbuf());
  std::cerr.rdbuf(Out.rdbuf());

  //Parallel initialization
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rankCpu);
  MPI_Comm_size(MPI_COMM_WORLD, &Ncpu);

  if(rankCpu == 0) {
    displayHeaderFire();
    // Get the directory name of the executable
    namespace fs = std::filesystem;
    fs::path current_path = fs::current_path();

    std::string last_dir;
    if (!current_path.filename().empty()) {
        last_dir = current_path.filename().string();
    } else {
        last_dir = current_path.parent_path().filename().string();
    }
    
    std::cout << "Running Fire in the case folder: " << last_dir << std::endl;
  }
  MPI_Barrier(MPI_COMM_WORLD);

  // run one test case from the input_files folder
  try {
    run = new Run("input_files/", 0);
    MPI_Barrier(MPI_COMM_WORLD);

    MPI_Barrier(MPI_COMM_WORLD);
    //2) Execution of the test case
    run->initialize(argc, argv);
    run->solver();
    run->finalize();
    //3) Removal of the test case
    delete run;
  }

  //Gestion of the exceptions
  //-------------------------
  catch (ErrorXML &e) {
    if (rankCpu == 0) std::cout << e.infoError() << std::endl;
    if (run) {
      run->finalize();
      delete run;
    }
  }
  catch (ErrorECOGEN &e) {
    if(rankCpu==0) std::cerr << e.infoError() << std::endl;
    if (run) {
      run->finalize();
      delete run;
    }
  }
  
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Finalize();
  std::cerr.rdbuf(strBackup);
  Out.close();
  return 0;
}

//***********************************************************************