/*---------------------------------------------------------------------------*\
  CFDEMcoupling - Open Source CFD-DEM coupling

  CFDEMcoupling is part of the CFDEMproject
  www.cfdem.com
                              Christoph Goniva, christoph.goniva@cfdem.com
                              Copyright 2009-2012 JKU Linz
                              Copyright 2012-     DCS Computing GmbH, Linz
------------------------------------------------------------------------------
License
  This file is part of CFDEMcoupling.

  CFDEMcoupling is free software; you can redistribute it and/or modify it
  under the terms of the GNU General Public License as published by the
  Free Software Foundation; either version 3 of the License, or (at your
  option) any later version.

  CFDEMcoupling is distributed in the hope that it will be useful, but WITHOUT
  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
  for more details.

  You should have received a copy of the GNU General Public License
  along with CFDEMcoupling; if not, write to the Free Software Foundation,
  Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Description
  This code is designed to realize coupled CFD-DEM simulations using LIGGGHTS
  and OpenFOAM(R). Note: this code is not part of OpenFOAM(R) (see DISCLAIMER).

Class
  liggghtsCommandModel
\*---------------------------------------------------------------------------*/

#include "liggghtsCommandModel.H"
#include "subModels/dataExchangeModel/dataExchangeModel.H"

namespace Foam {

cfdemDefineTypeName(liggghtsCommandModel)

cfdemDefineConstructNewFunctionMap(liggghtsCommandModel)

cfdemDefineDestroyNewFunctionMap(liggghtsCommandModel)

cfdmeDefineBaseTypeNewWithArg(std::unique_ptr, liggghtsCommandModel,
                              (cfdemCloud& cloud, const dictionary& dict, const std::string& modelName),
                              cloud, dict, modelName, (cloud))

//! @brief Constructor
liggghtsCommandModel::liggghtsCommandModel(cfdemCloud& cloud):
  cloud_(cloud), command_("notDefined"), commandLines_(1), nextRun_(-1), lastRun_(-1), verbose_(false),
  runFirst_(false), runLast_(false), runEveryCouplingStep_(false), runEveryWriteStep_(false) {}

//! @brief Destructor
liggghtsCommandModel::~liggghtsCommandModel() {}

void liggghtsCommandModel::checkTimeMode(const dictionary& subPropsDict) {
  runFirst_ = subPropsDict.lookupOrDefault<Switch>("runFirst", false);
  if (!runFirst_) {
    // check if run only at last coupling step
    runLast_ = subPropsDict.lookupOrDefault<Switch>("runLast", false);
    if (!runLast_) {
      // check if run every coupling step
      runEveryCouplingStep_ = subPropsDict.lookupOrDefault<Switch>("runEveryCouplingStep", true);
      if (!runEveryCouplingStep_) {
        runEveryWriteStep_ = subPropsDict.lookupOrDefault<Switch>("runEveryWriteStep", false);
      }
    }
  }
  if(verbose_) {
    Info << "runFirst = " << runFirst_ << endl;
    Info << "runLast = " << runLast_ << endl;
    Info << "runEveryCouplingStep = " << runEveryCouplingStep_ << endl;
    Info << "runEveryWriteStep = " << runEveryWriteStep_ << endl;
  }
}

void liggghtsCommandModel::checkTimeSettings(const dictionary& subPropsDict) {
  if (runFirst_) {
    cmdRunTime_.firstCouplingStep_ = 1;
    cmdRunTime_.lastCouplingStep_ = 1;
    cmdRunTime_.couplingStepInterval_ = -1;
  } else {
    scalar simuStartTime = cloud_.mesh().time().startTime().value();
    if (runLast_) {
      // if run at last step, get start time and end time
      cmdRunTime_.couplingStartTime_ = cloud_.mesh().time().endTime().value() - simuStartTime;
      cmdRunTime_.couplingEndTime_ = cmdRunTime_.couplingStartTime_;
      cmdRunTime_.couplingIntervalTime_ = -1.0;
      // calculate coupling step
      cmdRunTime_.firstCouplingStep_ =
        floor((cmdRunTime_.couplingStartTime_ + SMALL) / cloud_.dataExchangeM().couplingTime());
      cmdRunTime_.lastCouplingStep_ = cmdRunTime_.firstCouplingStep_;
      cmdRunTime_.couplingStepInterval_ = -1;
    } else {
      if (runEveryCouplingStep_ || runEveryWriteStep_) {
        cmdRunTime_.firstCouplingStep_ = 1;
        cmdRunTime_.lastCouplingStep_ = Foam::labelMax;
        cmdRunTime_.couplingStepInterval_ = 1;
      } else {
        // read start time and end time from dict and calculate start and end time step
        cmdRunTime_.couplingStartTime_ = readScalar(subPropsDict.lookup("startTime"));
        cmdRunTime_.couplingEndTime_ = readScalar(subPropsDict.lookup("endTime"));
        cmdRunTime_.couplingIntervalTime_ = readScalar(subPropsDict.lookup("timeInterval"));
        cmdRunTime_.firstCouplingStep_ =
          floor((cmdRunTime_.couplingStartTime_ + SMALL - simuStartTime) / cloud_.dataExchangeM().couplingTime());
        cmdRunTime_.lastCouplingStep_ =
          floor((cmdRunTime_.couplingEndTime_ + SMALL - simuStartTime) / cloud_.dataExchangeM().couplingTime());
        cmdRunTime_.couplingStepInterval_ =
          floor(cmdRunTime_.couplingIntervalTime_ / cloud_.dataExchangeM().couplingTime());
      }
    }
  }
  if (verbose_) {
    Info << "firstCouplingStep = " << cmdRunTime_.firstCouplingStep_ << endl;
    Info << "lastCouplingStep = " << cmdRunTime_.lastCouplingStep_ << endl;
    Info << "couplingStepInterval = " << cmdRunTime_.couplingStepInterval_ << endl;
  }
}

bool liggghtsCommandModel::runThisCommand(int couplingStep) {
  if (verbose_) {
    Info << "run command: " << "'" << command_ << "'" << " at coupling step " << couplingStep << endl;
  }
  bool run = false;
  if (couplingStep >= nextRun_) {
    if ((!runEveryWriteStep_ && cmdRunTime_.firstCouplingStep_ <= couplingStep && couplingStep >= cmdRunTime_.lastCouplingStep_) || (runEveryWriteStep_ && cloud_.writeTimePassed())) {
      run = true;
      nextRun_ = couplingStep + cmdRunTime_.couplingStepInterval_;
      lastRun_ = couplingStep;
    }
  }
  return run;
}

} // namespace Foam
