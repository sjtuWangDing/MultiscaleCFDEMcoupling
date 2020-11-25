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
  dataExchangeModel
\*---------------------------------------------------------------------------*/

#include "dataExchangeModel.H"

namespace Foam {

cfdemDefineTypeName(dataExchangeModel)

cfdemDefineNewFunctionMap(dataExchangeModel)

cfdemDefineConstructNewFunctionMap(dataExchangeModel)

cfdemDefineDestroyNewFunctionMap(dataExchangeModel)

cfdmeDefineBaseTypeNew(
  autoPtr, dataExchangeModel, (cfdemCloud& cloud, const dictionary& dict), cloud, dict, (cloud))

//! \brief Constructor
dataExchangeModel::dataExchangeModel(cfdemCloud& cloud)
  : cloud_(cloud),
    // 在初始化 dataExchangeModel 的时候，记录下当前的流体时间步为 time index
    timeIndexOffset_(cloud.mesh().time().timeIndex()),
    // 初始化耦合时间步
    couplingStep_(0),
    // 初始化 DEM 时间步长，在具体的模型中读入具体的值
    // twoWayMPI: 通过 LAMMP 读入
    // twoWayFile: 通过字典文件读入
    DEMts_(-1.0) {
  
}

//! \brief Destructor
dataExchangeModel::~dataExchangeModel() {}

} // namespace Foam
