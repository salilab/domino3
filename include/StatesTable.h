/**
 *  \file IMP/domino3/StatesTable.h
 *  \brief A beyesian infererence-based sampler.
 *
 *  Copyright 2007-2013 IMP Inventors. All rights reserved.
 *
 */

#ifndef IMPDOMINO3_STATES_TABLE_H
#define IMPDOMINO3_STATES_TABLE_H

#include <IMP/domino3/domino3_config.h>
#include "States.h"
#include <IMP/kernel/particle_index.h>
#include <IMP/base/map.h>

IMPDOMINO3_BEGIN_NAMESPACE
/** Store the association between particles and the classes
    which manage their states. I'm not a huge fan of having
    this class, but I haven't thought of a better way to store
    the information that is easily exposed to python
    and gets to all the right places. It is initialized internally
    in the DominoSampler.
 */
class IMPDOMINO3EXPORT StatesTable : public IMP::base::Object {
  typedef IMP::base::map<kernel::ParticleIndex,
  IMP::base::OwnerPointer<States> >
      Map;
  Map states_;

 public:
  StatesTable(kernel::Model *,
                      std::string name = "StatesTable%1%"):
                       Object(name) {}
  // implementation methods use this to get the enumerator
  States *get(kernel::ParticleIndex pi) const {
    IMP_USAGE_CHECK(get_has(pi),
                    "I don't know about particle " << pi);
    return states_.find(pi)->second;
  }
  bool get_has(kernel::ParticleIndex pi) const {
    return states_.find(pi) != states_.end();
  }
  kernel::ParticleIndexes get_particle_indexes() const {
    kernel::ParticleIndexes ret;
    ret.reserve(states_.size());
    for (Map::const_iterator it = states_.begin();
         it != states_.end(); ++it) {
      ret.push_back(it->first);
    }
    std::sort(ret.begin(), ret.end());
    return ret;
  }
  /** One can set the states more than once. If you do that, be
      careful.
  */
  void set_particle_states(kernel::ParticleIndex pi, States *e) {
    IMP_USAGE_CHECK(!get_has(pi),
                    "I already know about particle " << pi);
    IMP_USAGE_CHECK(
        e->get_number() > 0,
        "Cannot have 0 states for a particle: \"" << pi << "\"\n");
    states_[pi] = e;
  }
  IMP_OBJECT_METHODS(StatesTable);
};

IMP_OBJECTS(StatesTable, StatesTables);

IMPDOMINO3_END_NAMESPACE

#endif /* IMPDOMINO3_STATES_TABLE_H */
