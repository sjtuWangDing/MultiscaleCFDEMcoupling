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

#ifndef __LIGGGHTS_COMMAND_MODEL_H__
#define __LIGGGHTS_COMMAND_MODEL_H__

#include "cloud/cfdemCloud.H"
#include "base/runTimeSelectionTables.H"

namespace Foam {

class liggghtsCommandModel {

public:
  /* --------------------------------- base class and functions -------------------------- */

  //! \brief Runtime type information
  cfdemTypeName("liggghtsCommandModel");

  //! \brief Declare runtime constructor selection
  cfdemDeclareRunTimeSelection(std::unique_ptr, liggghtsCommandModel, (cfdemCloud& cloud), (cloud));

  //! \brief Selector
  static std::unique_ptr<liggghtsCommandModel> New(cfdemCloud& cloud,
                                                   const dictionary& dict,
                                                   const std::string& modelName);

  //! \brief Constructor
  liggghtsCommandModel(cfdemCloud& cloud);

  //! \brief Destructor
  virtual ~liggghtsCommandModel();

  /* --------------------------------- sub class ----------------------------------------- */

  struct CmdRunTime {
    scalar couplingStartTime_;
    scalar couplingEndTime_;
    scalar couplingIntervalTime_;
    int firstCouplingStep_;
    int lastCouplingStep_;
    int couplingStepInterval_;
    CmdRunTime():
      couplingStartTime_(-1.0), couplingEndTime_(-1.0), couplingIntervalTime_(0.0),
      firstCouplingStep_(-1), lastCouplingStep_(-1), couplingStepInterval_(0) {}
  };

  /* --------------------------------- virtual functions --------------------------------- */

  //! \brief 根据 index 获取 command 的内容
  virtual std::string getCommand(int index) const = 0;

  //! \param couplingStep 当前耦合步计数
  virtual bool runCommand(int couplingStep) = 0;

  /* --------------------------------- member functions ---------------------------------- */

  void checkTimeMode(const dictionary& subPropsDict);

  void checkTimeSettings(const dictionary& subPropsDict);

  bool runThisCommand(int couplingStep);

  /* --------------------------------- public interfaces --------------------------------- */

  inline int commandLines() const { return commandLines_; }

protected:

  cfdemCloud& cloud_;

  std::string command_;

  int commandLines_;

  Switch verbose_;

  Switch runFirst_;

  Switch runLast_;

  Switch runEveryCouplingStep_;

  Switch runEveryWriteStep_;

  CmdRunTime cmdRunTime_;

  int nextRun_;

  int lastRun_;
};

} // namespace Foam

#endif // __LIGGGHTS_COMMAND_MODEL_H__
