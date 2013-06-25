#include <IMP/domino3/Marginals.h>

IMPDOMINO3_BEGIN_NAMESPACE

Marginals::Marginals(Model *m, ParticleIndex pi, unsigned int size):
  base::Object("Marginals"+m->get_particle_name(pi)),
  current_(size, 0.0), next_(size, 0.0), change_(0.0) {
}

double Marginals::get_entropy() const {

}

IMPDOMINO3_END_NAMESPACE
