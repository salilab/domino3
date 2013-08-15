#include <IMP/domino3/Marginals.h>
#include <IMP/domino3/LogMathFunctions.h>
#include <IMP/domino3/LinearMathFunctions.h>

IMPDOMINO3_BEGIN_NAMESPACE

Marginals::Marginals(Model *m, ParticleIndex pi, unsigned int size):
  base::Object("Marginals"+m->get_particle_name(pi)), pi_(pi),
  current_(new double[size]), next_(new double[size]), change_(0.0), size_(size) {
  math = new LogMathFunctions();
  double start_value=1.0/this->get_number();
  start_value = math->convert_to_space(start_value);
  std::fill(current_.get(), current_.get() + size_, start_value);

  math->normalize(current_.get(),size_);
  double init_next = math->convert_to_space(0);

  std::fill(next_.get()   , next_.get()    + size_, init_next);
}

double Marginals::get_entropy() const {
  double ret = 0.0;
  for (unsigned int i = 0; i < this->get_number(); ++i) {
    double cur_prob = current_[i];
    if (cur_prob > 0) {
      double cur = cur_prob * std::log(cur_prob);
      ret += cur;
    }
  }
  return -ret;
}


void Marginals::set_uniform() {
    double start_value=1.0/this->get_number();
    start_value = math->convert_to_space(start_value);
    for (unsigned int i = 0; i < this->get_number(); ++i) {
        this->add_to_next_marginal(i, start_value);
    }
    this->set_current_from_next();
}

void Marginals::set_random() {
    for (unsigned int i = 0; i < this->get_number(); ++i) {
        this->add_to_next_marginal(i, std::rand());
    }
    this->set_current_from_next();
}




IMPDOMINO3_END_NAMESPACE
