#include <IMP/domino3/Marginals.h>

IMPDOMINO3_BEGIN_NAMESPACE

Marginals::Marginals(Model *m, ParticleIndex pi, unsigned int size):
  base::Object("Marginals"+m->get_particle_name(pi)),
  current_(new double[size]), next_(new double[size]), change_(0.0) {
  std::fill(current_.get(), current_.get() + size_, 0.0);
  std::fill(next_.get(), next_.get() + size_, 0.0);
}

double Marginals::get_entropy() const {
  double ret = 1.0;
  for (unsigned int i = 0; i < size_; ++i) {
    double cur_prob = current_[i];
    if (cur_prob > 0) {
      double cur = cur_prob * std::log(cur_prob);
      ret -= cur;
    }
  }
  return ret;
}

IMPDOMINO3_END_NAMESPACE
