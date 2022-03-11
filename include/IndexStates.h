/**
 *  \file IMP/domino3/IndexStates.h
 *  \brief A beyesian infererence-based sampler.
 *
 *  Copyright 2007-2017 IMP Inventors. All rights reserved.
 *
 */

#ifndef IMPDOMINO3_INDEX_STATES_H
#define IMPDOMINO3_INDEX_STATES_H

#include <IMP/domino3/domino3_config.h>
#include <IMP/algebra/Vector3D.h>
#include <IMP/particle_index.h>

IMPDOMINO3_BEGIN_NAMESPACE

/** Store a set of states which explicitly define the XYZ coordinates of
    the particle in question.
*/
class IMPDOMINO3EXPORT IndexStates : public States {
 const std::vector<int> states_;
protected:
 virtual void do_load(unsigned int i, ParticleIndex pi) const override;
 public:
 IndexStates(Model *m, const std::vector<int> &states);
  virtual unsigned int get_number() const override;
  virtual void add_to_rmf(ParticleIndex pi, RMF::NodeHandle parent) const override;
  virtual void update_rmf(ParticleIndex pi, RMF::NodeHandle parent,
                          Marginals *m) const override;
  IMP_OBJECT_METHODS(IndexStates);
};

IMPDOMINO3_END_NAMESPACE

#endif /* IMPDOMINO3_INDEX_STATES_H */
