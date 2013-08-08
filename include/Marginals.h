/** \file domino3/Marginals.h
 *  \brief Store the marginal for a variable.
 */

#ifndef IMPDOMINO3_MARGINAL_H
#define IMPDOMINO3_MARGINAL_H
#include <IMP/domino3/domino3_config.h>
#include <IMP/base/Object.h>
#include <IMP/base/ConstVector.h>
#include <IMP/base/object_macros.h>

IMPDOMINO3_BEGIN_NAMESPACE

class Marginals;

IMP_OBJECTS(Marginals, MarginalsList);

/** Store the marginal for a variable. */
class IMPDOMINO3EXPORT Marginals: public base::Object {
  boost::scoped_array<double> current_, next_;
  unsigned int size_;
  double change_;
  kernel::ParticleIndex pi_;


 public:
  Marginals(kernel::Model *m, kernel::ParticleIndex pi, unsigned int size);

  double get_marginal(unsigned int state) const {
    return current_[state];
  }

  void add_to_marginal(unsigned int state, double value) {
    next_[state] += value;
  }

    
  static void normalize(boost::scoped_array<double> &it,
                                     unsigned int size) {
    double total = std::accumulate(it.get(), it.get() + size, 0.0);
    IMP_USAGE_CHECK(total > .001,
                        "Total is too small to be reliable: " << total);
    for (unsigned int i = 0; i < size; ++i) {
        it[i] /= total;
    }
  }
    
  static void calculate_joint_probability(boost::scoped_array<double> &array,int size, const MarginalsListTemp &others){
    for (unsigned int i = 0; i < others.size(); ++i) {
        for (unsigned int j = 0; j < size; ++j) {
            array[j] *= others[i]->current_[j];
        }
    }
  }
  void merge_probabilities_from_list(const MarginalsListTemp &others) {
        std::copy(current_.get(), current_.get() + size_,next_.get());
        calculate_joint_probability(next_,size_,others);
        set_current_from_next();
  }
  /** Eventually this will be atomic. */
  void set_current_from_next() {
    normalize(next_, size_);
    change_ = 0;
    for (unsigned int i = 0; i < size_; ++i) {
      change_ += std::abs(next_[i] - current_[i]);
    }
    using namespace std;
    swap(current_, next_);
    std::fill(next_.get(), next_.get() + size_, 0.0);
  }


  /** Return a metric on the change (currently L0, could change) */
  double get_change() const {
    return change_;
  }
  
    
  kernel::ParticleIndex get_particle_index() const {
    return pi_;
  }
    
    
  unsigned int get_number() const {
    return size_;
  }

  double get_entropy() const;

  IMP_OBJECT_METHODS(Marginals);
};
 
IMPDOMINO3EXPORT void set_uniform(Marginals *m);
IMPDOMINO3EXPORT void set_random(Marginals *m);
IMPDOMINO3_END_NAMESPACE


#endif // IMPDOMINO3_MARGINAL_H
