/** \file domino3/Updater.h
 *  \brief Store the marginal for a variable.
 */

#ifndef IMPDOMINO3_UPDATER_H
#define IMPDOMINO3_UPDATER_H
#include <IMP/domino3/domino3_config.h>
#include <IMP/base/Object.h>
#include "Node.h"
#include <IMP/base/object_macros.h>
#include <IMP/base/map.h>
#include <IMP/base/particle_index.h>
#include <IMP/base/set.h>
#include "Marginals.h"

IMPDOMINO3_BEGIN_NAMESPACE

/** Updater its marginals based on some criteria. */
class IMPDOMINO3EXPORT Updater: public base::Object {
  NodeGraph graph_;
  NodeGraphVertexIndex graph_index_;
  Nodes nodes_;
  typedef base::set<unsigned int> ActiveSet;
  ActiveSet cur_queue_, next_queue_;
  double change_threshold_;

  void add_node_to_active_set(unsigned int index) {
    next_queue_.insert(n);
  }
  void swap_active_sets() {
    std::swap(cur_queue_, next_queue_);
    cur_queue_= base::set<unsigned int>();
  }
  void do_update();
  void fill_queue();

 public:
  Updater(const NodeGraph &graph,
          std::string name);

  void update(unsigned int iterations);

  void set_change_threshold(double d) {
    change_threshold_ = d;
  }

  IMP_OBJECT_METHODS(Updater);
};

IMPDOMINO3_END_NAMESPACE

#endif // IMPDOMINO3_UPDATER_H
