/**
 *  \file domino/DominoSampler.h
 *  \brief A beyesian infererence-based sampler.
 *
 *  Copyright 2007-2013 IMP Inventors. All rights reserved.
 *
 */

#include <IMP/domino3/IndexStates.h>
#include <RMF/decorators.h>

IMPDOMINO3_BEGIN_NAMESPACE

IndexStates::IndexStates(Model *m, const IMP::kernel::Particles &states)
      : States(m, "IndexStates %1%"),
        states_(states) {
}

unsigned int IndexStates::get_number() const {
  return states_.size();
}

void IndexStates::do_load(unsigned int i, ParticleIndex pi) const {

}


void IndexStates::add_to_rmf(ParticleIndex, RMF::NodeHandle parent) const {

}

void IndexStates::update_rmf(ParticleIndex, RMF::NodeHandle parent,
                           Marginals *m) const {
}

IMPDOMINO3_END_NAMESPACE