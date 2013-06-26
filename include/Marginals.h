/** \file domino3/Marginal.h
 *  \brief Store the marginal for a variable.
 */

#ifndef IMPDOMINO3_MARGINAL_H
#define IMPDOMINO3_MARGINAL_H
#include <IMP/domino3/domino3_config.h>
#include <IMP/base/Object.h>
#include <IMP/base/ConstVector.h>
#include <IMP/base/object_macros.h>

IMPDOMINO3_BEGIN_NAMESPACE

IMP_OBJECTS(Marginals, MarginalsList);

/** Store the marginal for a variable. */
class IMPDOMINO3EXPORT Marginals: public base::Object {
  base::ConstVector<double> current_, next_;
  double change_;

  static void normalize(base::ConstVector<double> &it) {
    double total = std::accumulate(it.begin(), it.end(), 0.0);
    for (unsigned int i = 0; i < it.size(); ++i) {
      it[i] /= total;
    }
  }
 public:
  Marginals(Model *m, ParticleIndex pi, unsigned int size);

  double get_marginal(unsigned int state) const {
    return current_[state];
  }

  double add_to_marginal(unsigned int state, double value) {
    next_[state] += value;
  }

  /** Eventually this will be atomic. */
  void update_current_from_next() {
    normalize(next_);
    change_ = 0;
    for (unsigned int i = 0; i < current_.size(); ++i) {
      change_ += std::abs(next_[i] - current_[i]);
    }
    using namespace std;
    swap(current_, next_);
    std::fill(next_.begin(), next_.end(), 0.0);
  }

  void update_current_from_list(const MarginalsListTemp &others) {
    std::fill(next_.begin(), next_.end(), 0.0);
    for (unsigned int i = 0; i < others.size(); ++i) {
      for (unsigned int j = 0; j < next_.size(); ++j) {
        next_[j] += others[i]->current_[j];
      }
    }
    update_current_from_next();
  }

  /** Return a metric on the change (currently L0, could change) */
  double get_change() const {
    return change_;
  }

  double get_entropy() const;

  IMP_OBJECT_METHODS(Marginals);
};

IMPDOMINO3_END_NAMESPACE


#endif // IMPDOMINO3_MARGINAL_H
