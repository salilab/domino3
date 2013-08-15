/**
 *  \file domino/DominoSampler.h
 *  \brief A beyesian infererence-based sampler.
 *
 *  Copyright 2007-2013 IMP Inventors. All rights reserved.
 *
 */

#include <IMP/domino3/XYZStates.h>
#include <IMP/core/XYZ.h>
#include <RMF/decorators.h>

IMPDOMINO3_BEGIN_NAMESPACE

XYZStates::XYZStates(Model *m, const algebra::Vector3Ds &states)
      : States(m, "XYZStates %1%"),
        states_(states), max_radius_(1) {
}

unsigned int XYZStates::get_number() const {
  return states_.size();
}
void XYZStates::do_load(unsigned int i, ParticleIndex pi) const {
  core::XYZ(get_model(), pi).set_coordinates(states_[i]);
}


void XYZStates::add_to_rmf(ParticleIndex, RMF::NodeHandle parent) const {
  RMF::ParticleFactory f(parent.get_file());
  for (unsigned int i = 0; i < states_.size(); ++i) {
    std::ostringstream oss;
    oss << i;
    RMF::NodeHandle n = parent.add_child(oss.str(), RMF::REPRESENTATION);
    RMF::Particle d = f.get(n);
    d.set_mass(1.0);
    d.set_coordinates(RMF::Floats(states_[i].coordinates_begin(),
                                  states_[i].coordinates_end()));
  }
}
void XYZStates::update_rmf(ParticleIndex, RMF::NodeHandle parent,
                           Marginals *m) const {
  RMF::ParticleFactory f(parent.get_file());
  RMF::NodeHandles children = parent.get_children();
  for (unsigned int i = 0; i < states_.size(); ++i) {
    double r = max_radius_ * m->get_current_marginal(i);
    RMF::Particle d = f.get(children[i]);
    d.set_radius(r);
  }
}

IMPDOMINO3_END_NAMESPACE
