/*---------------------------------------------------------------------------*\
    CFDEMcoupling - Open Source CFD-DEM coupling

    CFDEMcoupling is part of the CFDEMproject
    www.cfdem.com
                                Christoph Goniva, christoph.goniva@cfdem.com
                                Copyright 2009-2012 JKU Linz
                                Copyright 2012-     DCS Computing GmbH, Linz
                                Copyright (C) 2013-     Graz University of  
                                                        Technology, IPPT
-------------------------------------------------------------------------------
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
    mixConstDiffSmoothing

SourceFiles
    mixConstDiffSmoothing.C

\*---------------------------------------------------------------------------*/

#ifndef mixConstDiffSmoothing_H
#define mixConstDiffSmoothing_H

#include "smoothingModel.H"

namespace Foam {

class mixConstDiffSmoothing: public smoothingModel {

private:

  dictionary propsDict_;

  // @brief 标量场的下限值
  // @note Scalar fields will be bound to this lower value
  scalar lowerLimit_;

  // @brief 标量场的上限值
  // @note Scalar fields will be bound to this upper value
  scalar upperLimit_;

  // @brief 光滑长度尺度
  // @note Length scale over which the exchange fields will be smoothed out
  dimensionedScalar smoothingLength_;

  dimensionedScalar smoothingLengthReferenceField_;

  mutable dimensionedScalar DT_;

  // @brief Flag for debugging output
  bool verbose_;

public:

  // @brief Runtime type information
  TypeName("mixConstDiffSmoothing");

  // @brief Constructors
  mixConstDiffSmoothing(const dictionary& dict, cfdemCloud& sm);

  // @brief Destructor
  ~mixConstDiffSmoothing();

  // Member Functions

  virtual bool doSmoothing() const { return true; }

  void smoothen(volScalarField&) const;

  void smoothen(volVectorField&) const;

  void smoothenReferenceField(volVectorField&) const;

};  // End of class mixConstDiffSmoothing

}  // End of namespace Foam

#endif  // mixConstDiffSmoothing_H