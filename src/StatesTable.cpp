/**
 *  \file domino/DominoSampler.h
 *  \brief A beyesian infererence-based sampler.
 *
 *  Copyright 2007-2013 IMP Inventors. All rights reserved.
 *
 */

#include <IMP/domino3/StatesTable.h>

IMPDOMINO3_BEGIN_NAMESPACE

kernel::ParticleIndexes StatesTable::get_particle_indexes() const {
  kernel::ParticleIndexes ret;
  ret.reserve(states_.size());
  for (Map::const_iterator it = states_.begin();
       it != states_.end(); ++it) {
    ret.push_back(it->first);
  }
  std::sort(ret.begin(), ret.end());
  return ret;
}

void StatesTable::add(kernel::ParticleIndex pi,
                      States *e,
                      Marginals *m) {
  IMP_USAGE_CHECK(!get_has(pi),
                  "I already know about particle " << pi);
  IMP_USAGE_CHECK(e->get_number() == m->get_number(),
                  "Sizes don't match.");
  IMP_USAGE_CHECK(
                  e->get_number() > 0,
                  "Cannot have 0 states for a particle: \"" << pi << "\"\n");
  states_[pi] = ParticleData(e, m);
}

void StatesTable::set_rmf(RMF::NodeHandle parent) {
  parent_ = parent;
  for (Map::iterator it = states_.begin(); it != states_.end(); ++it) {
    RMF::NodeHandle cur = parent.add_child(m_->get_particle_name(it->first),
                                           RMF::ORGANIZATIONAL);
    it->second.set_node(cur);
    it->second.get_states()->add_to_rmf(it->first, cur);
  }
}

void StatesTable::add_to_frame() {
  for (Map::iterator it = states_.begin(); it != states_.end(); ++it) {
    it->second.get_states()->update_rmf(it->first,
                                        it->second.get_node(),
                                        it->second.get_marginals());
  }
}


IMPDOMINO3_END_NAMESPACE
