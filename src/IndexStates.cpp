/**
 *  \file domino/DominoSampler.h
 *  \brief A beyesian infererence-based sampler.
 *
 *  Copyright 2007-2015 IMP Inventors. All rights reserved.
 *
 */

#include <IMP/domino3/IndexStates.h>
#include <RMF/decorators.h>

IMPDOMINO3_BEGIN_NAMESPACE

IndexStates::IndexStates(Model *m, const std::vector<int> &states)
      : States(m, "IndexStates %1%"),
        states_(states) {
}

unsigned int IndexStates::get_number() const {
  return states_.size();
}

void IndexStates::do_load(unsigned int, ParticleIndex) const {

}


void IndexStates::add_to_rmf(ParticleIndex, RMF::NodeHandle) const {

}

void IndexStates::update_rmf(ParticleIndex, RMF::NodeHandle,
                             Marginals *) const {
}

IMPDOMINO3_END_NAMESPACE
