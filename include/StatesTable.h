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
  IMP_NAMED_TUPLE_3(ParticleData, ParticleDatas,
                    IMP::base::OwnerPointer<States>, states,
                    IMP::base::OwnerPointer<Marginals>, marginals,
                    RMF::NodeHandle, node,);
  typedef IMP::base::map<kernel::ParticleIndex, ParticleData> Map;
  Map states_;
  base::WeakPointer<kernel::Model> m_;
  RMF::NodeHandle parent_;
 public:
  StatesTable(kernel::Model *m,
              std::string name = "StatesTable%1%"):
    Object(name), m_(m) {}
  // implementation methods use this to get the enumerator
  States *get_states(kernel::ParticleIndex pi) const {
    IMP_USAGE_CHECK(get_has(pi),
                    "I don't know about particle " << pi);
    return states_.find(pi)->second.get_states();
  }
  Marginals *get_marginals(kernel::ParticleIndex pi) const {
    IMP_USAGE_CHECK(get_has(pi),
                    "I don't know about particle " << pi);
    return states_.find(pi)->second.get_marginals();
  }
  bool get_has(kernel::ParticleIndex pi) const {
    return states_.find(pi) != states_.end();
  }
  kernel::ParticleIndexes get_particle_indexes() const;
  void add(kernel::ParticleIndex pi, States *e, Marginals *m);
  void set_rmf(RMF::NodeHandle parent);
  void add_to_frame();
  IMP_OBJECT_METHODS(StatesTable);
};

IMP_OBJECTS(StatesTable, StatesTables);

IMPDOMINO3_END_NAMESPACE

#endif /* IMPDOMINO3_STATES_TABLE_H */
