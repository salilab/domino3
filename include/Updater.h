/** \file IMP/domino3/Updater.h
 *  \brief Store the marginal for a variable.
 */

#ifndef IMPDOMINO3_UPDATER_H
#define IMPDOMINO3_UPDATER_H

#include <IMP/domino3/domino3_config.h>
#include <IMP/base/Object.h>
#include "Factor.h"
#include <IMP/base/object_macros.h>
#include <IMP/kernel/particle_index.h>
#include <boost/unordered_set.hpp>
#include "Marginals.h"

IMPDOMINO3_BEGIN_NAMESPACE

/** Updater its marginals based on some criteria. */
class IMPDOMINO3EXPORT Updater: public base::Object {
  Factors factors_;
  typedef std::vector<Factor *> ActiveSet;
  ActiveSet cur_queue_, next_queue_;
  FP change_threshold_;

  void add_factor_to_next_set(Factor *n) {
    next_queue_.push_back(n);
  }
  void swap_active_sets() {
    std::swap(cur_queue_, next_queue_);
  }
  void do_update();
  void fill_queue();
  float change;
  float entropy;

 public:
  Updater(const FactorsTemp &factors,
          std::string name);
 
  void update(unsigned int iterations);

  void set_change_threshold(FP d) {
    change_threshold_ = d;
  }

  IMP_OBJECT_METHODS(Updater);
};

IMP_OBJECTS(Updater, Updaters);

IMPDOMINO3_END_NAMESPACE

#endif /* IMPDOMINO3_UPDATER_H */
