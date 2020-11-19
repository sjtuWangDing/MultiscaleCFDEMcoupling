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
  Foam::Particle
\*---------------------------------------------------------------------------*/

#ifndef __PARTICLE_H__
#define __PARTICLE_H__

#include "error.H"
#include "base/tensor/tensor.H"

namespace Foam {

//! \brief 颗粒尺度枚举
enum ScaleType {
  KUnknown = 0, // unknown
  kFine,        // fine
  kMiddle,      // middle
  kCoarse       // coarse
};

/*!
 * \brief base class for particle
 * \tparam PType inheritated class must put their type into this parameter
 */
template<typename PType>
class Particle {
public:
  //! \return  subtype instance of current class
  inline const PType& self() const {
    return *static_cast<const PType*>(this);
  }
  //! \return reference of subtype instance of current class
  inline PType* ptrSelf() {
    return static_cast<PType*>(this);
  }
  //! \brief Constructor
  inline Particle(index_t index, ScaleType scaleType = KUnknown)
      : index_(index), scaleType_(scaleType) {}
  //! \brief Copy constructor
  inline Particle(const Particle<PType>& p)
      : index_(p.index_), scaleType_(p.scaleType_) {}
  //! \brief operator=
  inline Particle& operator= (const Particle<PType>& p) {
    if (*this != p) {
      index_ = p.index_;
      scaleType_ = p.scaleType_;
    }
  }
  /*!
   * \return whether two particle equals
   * \param p the particle to compare against
   * \note 由于同一个颗粒在每一个cpu核心中都有一个，所以这里不使用 this == &p 来比较
   */
  inline bool operator== (const Particle<PType>& p) const {
    return index_ == p.index_;
  }
  inline bool operator!= (const Particle<PType>& p) const {
    return !(*this == p);
  }
protected:
  //! \brief 颗粒索引，颗粒的唯一标识
  index_t index_;
  //! \brief 颗粒尺度类型
  ScaleType scaleType_;
};

class CoarseParticle: public Particle<CoarseParticle> {
public:
  inline CoarseParticle(index_t index)
      : Particle<CoarseParticle>(index, kCoarse),
        radius_(0.0) {}

private:
  double radius_;
  base::CITensor1 cellID_;
  base::CDTensor1 position_;
  base::CDTensor1 voidFraction_;
};

} // namespace Foam

#endif // __PARTICLE_H__
