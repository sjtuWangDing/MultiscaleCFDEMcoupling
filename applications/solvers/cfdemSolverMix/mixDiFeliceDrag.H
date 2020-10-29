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

    Di Felice drag law (Zhou et al. (2010), JFM)


Class
    mixDiFeliceDrag

SourceFiles
    mixDiFeliceDrag.C

\*---------------------------------------------------------------------------*/

#ifndef mixDiFeliceDrag_H
#define mixDiFeliceDrag_H

#include "forceModel.H"
#include "interpolationCellPoint.H"

namespace Foam {

class mixDiFeliceDrag: public forceModel {

private:

  dictionary propsDict_;

  //! @brief 速度场名称
  word velFieldName_;

  //! @brief 流体速度场的常引用
  const volVectorField& U_;

  //! @brief 空隙率场的名称
  word voidfractionFieldName_;

  //! @brief 空隙率场的常引用
  const volScalarField& voidfraction_;

  //! @brief 局部平均颗粒速度场名称
  word UsFieldName_;

  //! @brief 局部平均颗粒速度场的常引用
  const volVectorField& UsField_;

public:

  //! @brief Runtime type information
  TypeName("mixDiFeliceDrag");

  //! @brief Constructors
  mixDiFeliceDrag(const dictionary& dict, cfdemCloud& sm);

  //! @brief Destructor
  ~mixDiFeliceDrag();

  // Member Functions
  void calForceKernel(const int& index,
                      const int& cellI,
                      const volScalarField& nufField,
                      const volScalarField& rhoField,
                      vector& drag, vector& Ufluid, scalar& dragCoefficient,
                      scalar& voidfraction, vector& Ur,
                      scalar& Rep, scalar& Cd) const;

  virtual void setForce() const;

  virtual void setMixForce(const std::vector<double>& dimensionRatios) const;

};  // End of class mixDiFeliceDrag

}  // End of namespace Foam

#endif  // mixDiFeliceDrag_H
