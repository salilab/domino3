/**
 *  \file IMP/domino3/States.h
 *  \brief A beyesian infererence-based sampler.
 *
 *  Copyright 2007-2015 IMP Inventors. All rights reserved.
 *
 */

#ifndef IMPDOMINO3_STATES_H
#define IMPDOMINO3_STATES_H

#include <IMP/domino3/domino3_config.h>
#include <IMP/particle_index.h>
#include <RMF/NodeHandle.h>
#include "Marginals.h"
#include <IMP/base/Pointer.h>
#include <RMF/NodeHandle.h>

IMPDOMINO3_BEGIN_NAMESPACE
/** Handle the states for a particular particle (or "class" of
    particles. For example a state enumerator class could take
    a bounding box and a number,n, and generate n points in the
    bounding box. Then the get_number function woudld return
    n and update_to_state would modify the particle to have the
    coordiantes for state i.
 */
class IMPDOMINO3EXPORT States : public IMP::base::Object {
  base::WeakPointer<Model> m_;
protected:
  virtual void do_load(unsigned int, ParticleIndex pi) const = 0;
 public:
  States(Model *m, std::string name = "States %1%") :
  Object(name), m_(m) {}
  Model *get_model() const {return m_;}
  virtual unsigned int get_number() const = 0;
  virtual void load(unsigned int i, ParticleIndex pi) const {
    IMP_USAGE_CHECK(i < get_number(),
                    "Out of range state requested");
    do_load(i, pi);
  }
  /** Add the marginals to a file for displaying with the passed
      node as a parent.*/
  virtual void add_to_rmf(ParticleIndex pi, RMF::NodeHandle parent) const = 0;
  /** Update the RMF file contents for the current marginals. Assume everything
      under the parent is yours.*/
  virtual void update_rmf(ParticleIndex pi, RMF::NodeHandle parent,
                          Marginals *m) const = 0;
  virtual ~States();
};

IMP_OBJECTS(States, StatesList);

IMPDOMINO3_END_NAMESPACE

#endif /* IMPDOMINO3_STATES_H */
