/**
 *  \file domino/DominoSampler.h
 *  \brief A beyesian infererence-based sampler.
 *
 *  Copyright 2007-2013 IMP Inventors. All rights reserved.
 *
 */

#include <IMP/domino3/XYZStates.h>
#include <IMP/core/XYZ.h>

IMPDOMINO3_BEGIN_NAMESPACE

unsigned int XYZStates::get_number() const {
  return states_.size();
}
void XYZStates::do_load(unsigned int i, ParticleIndex pi) const {
  core::XYZ(get_model(), pi).set_coordinates(states_[i]);
}

IMPDOMINO3_END_NAMESPACE
