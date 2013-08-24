#include <IMP/domino3/Updater.h>

IMPDOMINO3_BEGIN_NAMESPACE

Updater::Updater(const NodesTemp &graph,
                 std::string name):
  base::Object(name), nodes_(graph) {
  for (unsigned int i = 0; i < nodes_.size(); ++i) {
    nodes_[i]->set_index(i);
  }
  set_change_threshold(.01);
  fill_queue();
}

void Updater::do_update() {
  base::set<kernel::ParticleIndex> changed;
  for (ActiveSet::const_iterator it = cur_queue_.begin();
       it != cur_queue_.end(); ++it) {
    IMP_LOG_TERSE("Updating " << (*it)->get_name() << std::endl);
    (*it)->update();
    for (unsigned int i = 0; i < (*it)->get_marginals().size(); ++i) {
      if ((*it)->get_marginals()[i]->get_change() > change_threshold_) {
        changed.insert((*it)->get_particle_indexes()[i]);
      }
    }
  }
  // later search only relevant nodes using the graph
  for (unsigned int i = 0; i< nodes_.size(); ++i) {
    for (unsigned int j = 0; j < nodes_[i]->get_particle_indexes().size(); ++j) {
      if (changed.find(nodes_[i]->get_particle_indexes()[j]) != changed.end()) {
        add_node_to_next_set(nodes_[i]);
        break;
      }
    }
  }
  swap_active_sets();
  next_queue_ = ActiveSet();
}

void Updater::fill_queue() {
  // lazy
  for (unsigned int i = 0; i < nodes_.size(); ++i) {
    add_node_to_next_set(nodes_[i]);
  }
  swap_active_sets();
}

void Updater::update(unsigned int iterations) {
  for (unsigned int i = 0; i < iterations; ++i) {
    do_update();
  }
}

IMPDOMINO3_END_NAMESPACE
