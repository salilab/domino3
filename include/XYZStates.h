/**
 *  \file IMP/domino3/XYZStates.h
 *  \brief A beyesian infererence-based sampler.
 *
 *  Copyright 2007-2017 IMP Inventors. All rights reserved.
 *
 */

#ifndef IMPDOMINO3_XYZ_STATES_H
#define IMPDOMINO3_XYZ_STATES_H

#include <IMP/domino3/domino3_config.h>
#include <IMP/domino3/States.h>
#include <IMP/algebra/Vector3D.h>
#include <IMP/particle_index.h>

IMPDOMINO3_BEGIN_NAMESPACE

/** Store a set of states which explicitly define the XYZ coordinates of
    the particle in question.
*/
class IMPDOMINO3EXPORT XYZStates : public States {
  algebra::Vector3Ds states_;
  double max_radius_;
protected:
 virtual void do_load(unsigned int i, ParticleIndex pi) const override;
 public:
 XYZStates(Model *m, const algebra::Vector3Ds &states);
  virtual unsigned int get_number() const override;
  void set_max_radius(double r) {max_radius_ = r;}
  virtual void add_to_rmf(ParticleIndex pi, RMF::NodeHandle parent) const override;
  virtual void update_rmf(ParticleIndex pi, RMF::NodeHandle parent,
                          Marginals *m) const override;
  IMP_OBJECT_METHODS(XYZStates);
};

IMPDOMINO3_END_NAMESPACE

#endif /* IMPDOMINO3_XYZ_STATES_H */
