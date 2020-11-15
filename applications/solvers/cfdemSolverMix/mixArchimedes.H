/*---------------------------------------------------------------------------*\
    CFDEMcoupling - Open Source CFD-DEM coupling

    CFDEMcoupling is part of the CFDEMproject
    www.cfdem.com
                                Christoph Goniva, christoph.goniva@cfdem.com
                                Copyright 2009-2012 JKU Linz
                                Copyright 2012-     DCS Computing GmbH, Linz
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

    Archimedes buoyancy force for sphere / cylinder (2d)
    basic "on-off" handling of buoyancy at the interface

Class
    mixArchimedes

SourceFiles
    mixArchimedes.C

\*---------------------------------------------------------------------------*/

#ifndef mixArchimedes_H
#define mixArchimedes_H

#include "forceModel.H"

namespace Foam {

class mixArchimedes: public forceModel {

private:

  dictionary propsDict_;

  bool twoDimensional_;

  word gravityFieldName_;

#ifdef version21
  const uniformDimensionedVectorField& g_;
#elif defined(version16ext) || defined(version15)
  const dimensionedVector& g_;
#endif

public:

  // \brief Runtime type information
  TypeName("mixArchimedes");


  // \brief Constructors
  mixArchimedes(const dictionary& dict, cfdemCloud& sm);

  // \brief Destructor
  ~mixArchimedes();

  // Member Functions
  void calForceKernel(const int& index,
                      const int& cellI,
                      vector& force) const;

  virtual void setForce() const;

  virtual void setMixForce(const std::vector<double>& dimensionRatios) const;

};  // End of class mixArchimedes

}  // End of namespace Foam

#endif  // mixArchimedes_H
