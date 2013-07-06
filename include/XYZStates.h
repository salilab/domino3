/**
 *  \file IMP/domino3/XYZStates.h
 *  \brief A beyesian infererence-based sampler.
 *
 *  Copyright 2007-2013 IMP Inventors. All rights reserved.
 *
 */

#ifndef IMPDOMINO3_XYZ_STATES_H
#define IMPDOMINO3_XYZ_STATES_H

#include <IMP/domino3/domino3_config.h>
#include <IMP/algebra/Vector3D.h>
#include <IMP/kernel/particle_index.h>
#include <IMP/base/map.h>

IMPDOMINO3_BEGIN_NAMESPACE

/** Store a set of states which explicitly define the XYZ coordinates of
    the particle in question.
*/
class IMPDOMINO3EXPORT XYZStates : public States {
  algebra::Vector3Ds states_;
protected:
 virtual void do_load(unsigned int i, ParticleIndex pi) const IMP_OVERRIDE;
 public:
  XYZStates(Model *m, const algebra::Vector3Ds &states)
      : States(m, "XYZStates %1%"),
        states_(states) {}
  virtual unsigned int get_number() const IMP_OVERRIDE;
  IMP_OBJECT_METHODS(XYZStates);
};

IMPDOMINO3_END_NAMESPACE

#endif /* IMPDOMINO3_XYZ_STATES_H */
