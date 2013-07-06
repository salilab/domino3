/**
 *  \file IMP/domino3/States.h
 *  \brief A beyesian infererence-based sampler.
 *
 *  Copyright 2007-2013 IMP Inventors. All rights reserved.
 *
 */

#ifndef IMPDOMINO3_STATES_H
#define IMPDOMINO3_STATES_H

#include <IMP/domino3/domino3_config.h>
#include <IMP/kernel/particle_index.h>
#include "Marginals.h"
#include <IMP/base/Pointer.h>

IMPDOMINO3_BEGIN_NAMESPACE
/** Handle the states for a particular particle (or "class" of
    particles. For example a state enumerator class could take
    a bounding box and a number,n, and generate n points in the
    bounding box. Then the get_number function woudld return
    n and update_to_state would modify the particle to have the
    coordiantes for state i.
 */
class IMPDOMINO3EXPORT States : public IMP::base::Object {
  base::WeakPointer<kernel::Model> m_;
protected:
  virtual void do_load(unsigned int, kernel::ParticleIndex pi) const = 0;
 public:
  States(kernel::Model *m, std::string name = "States %1%") :
  Object(name), m_(m) {}
  kernel::Model *get_model() const {return m_;}
  virtual unsigned int get_number() const = 0;
  virtual void load(unsigned int i, kernel::ParticleIndex pi) const {
    IMP_USAGE_CHECK(i < get_number(),
                    "Out of range state requested");
    do_load(i, pi);
  }
  virtual ~States();
};

IMP_OBJECTS(States, StatesList);

IMPDOMINO3_END_NAMESPACE

#endif /* IMPDOMINO3_STATES_H */
