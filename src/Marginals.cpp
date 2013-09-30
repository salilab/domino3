#include <IMP/domino3/Marginals.h>
#include <IMP/domino3/LogMathFunctions.h>

IMPDOMINO3_BEGIN_NAMESPACE

Marginals::Marginals(Model *m, ParticleIndex pi, unsigned int size):
  base::Object("Marginals"+m->get_particle_name(pi)), pi_(pi),
  current_(new FP[size]), next_(new FP[size]), change_(0.0), size_(size) {
  FP start_value=1.0/this->get_number();
      start_value = LogMathFunctions::convert_to_space(start_value);
  std::fill(current_.get(), current_.get() + size_, start_value);

  LogMathFunctions::normalize(current_.get(),size_);
      FP init_next = LogMathFunctions::convert_to_space(0);
  this->set_was_used(true);
  std::fill(next_.get()   , next_.get()    + size_, init_next);
}

FP Marginals::get_entropy() const {
  FP ret = 0.0;
  for (unsigned int i = 0; i < this->get_number(); ++i) {
    FP cur_prob = current_[i];
    if (cur_prob > 0) {
        FP cur = cur_prob * LogMathFunctions::convert_to_space(cur_prob);
      ret += cur;
    }
  }
  return -ret;
}

void Marginals::set_init_vector(boost::scoped_array<FP> &array) {
    for (unsigned int i = 0; i < this->get_number(); ++i) {
        this->add_to_next_marginal(i, array[i]);
    }
    this->set_current_from_next();
}


void Marginals::set_uniform() {
    FP start_value=1.0/this->get_number();
    start_value = LogMathFunctions::convert_to_space(start_value);
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
