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

    ArchimedesIB buoyancy force for sphere / cylinder (2d)
    basic "on-off" handling of buoyancy at the interface

	Contribution by Alice Hager.

Class
    mixArchimedesIB

SourceFiles
    mixArchimedesIB.C

\*---------------------------------------------------------------------------*/

#ifndef mixArchimedesIB_H
#define mixArchimedesIB_H

#include "forceModel.H"

namespace Foam {

class mixArchimedesIB: public forceModel {

private:

  dictionary propsDict_;

  bool twoDimensional_;

  word voidfractionFieldName_;

  const volScalarField& voidfractions_;

  word volumefractionFieldName_;

  const volScalarField& volumefractions_;

  word gravityFieldName_;

#ifdef version21
  const uniformDimensionedVectorField& g_;
#elif defined(version16ext) || defined(version15)
  const dimensionedVector& g_;
#endif

public:

  //! @brief Runtime type information
  TypeName("mixArchimedesIB");

  //! @brief Constructors
  mixArchimedesIB(const dictionary& dict,
                  cfdemCloud& sm);

  //! @brief Destructor
  ~mixArchimedesIB();

  void calForceKernel(const int& index,
                      vector& force) const;

  virtual void setForce() const;

  virtual void setMixForce(const std::vector<double>& dimensionRatios) const;

};  // End of class mixArchimedesIB

}  // End of namespace Foam

#endif  // mixArchimedesIB_H
