/**
 *  \file domino/DominoSampler.h
 *  \brief A beyesian infererence-based sampler.
 *
 *  Copyright 2007-2017 IMP Inventors. All rights reserved.
 *
 */

#include <IMP/domino3/StatesTable.h>

IMPDOMINO3_BEGIN_NAMESPACE

ParticleIndexes StatesTable::get_particle_indexes() const {
  ParticleIndexes ret;
  ret.reserve(states_.size());
  for (Map::const_iterator it = states_.begin();
       it != states_.end(); ++it) {
    ret.push_back(it->first);
  }
  std::sort(ret.begin(), ret.end());
  return ret;
}

void StatesTable::add(ParticleIndex pi,
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
  states_incoming_order_.push_back(pi);
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

void StatesTable::show_marginal(){
      for (unsigned int i = 0; i < states_incoming_order_.size(); ++i) {
          IMP::domino3::Marginals * marg = this->get_marginals(states_incoming_order_[i]);
	  std::cout << states_incoming_order_[i] << ": ";
          for(unsigned int y = 0; y < marg->get_number(); y++){
              std::cout << LogMathFunctions::convert_to_linear(marg->get_current_marginal(y)) << "\t";
          }
          std::cout << std::endl;
      }
}

IMPDOMINO3_END_NAMESPACE
