#include <IMP/domino3/Updater.h>

IMPDOMINO3_BEGIN_NAMESPACE

Updater::Updater(const NodeGraph &graph,
                 std::string name):
  base::Object(name), graph_(graph),
  graph_index_(get_vertex_index(graph)) {
  NodeGraphVertexMap vm = boost::get(graph_, boost::vertex_name);
  BOOST_FOREACH(NodeGraphVertex vertex,
                boost::vertices(graph_)) {
    nodes_.push_back(vm[vertex]);
  }
  set_change_threshold(.05);
}

void Updater::do_update() {
  base::set<kernel::ParticleIndex> changed;
  for (ActiveSet::const_iterator it = cur_queue_.begin();
       it != cur_queue_.end(); ++it) {
    unsigned int cur = *it;
    nodes_[cur]->update();
    for (unsigned int i = 0; i < nodes_[cur]->get_marginals().size(); ++i) {
      if (nodes_[cur]->get_marginals()[i]->get_change() > change_threshold_) {
        changed.insert(nodes_[cur]->get_marginals()[i]
                       ->get_particle_indexes()[i]);
      }
    }
  }
  // later search only relevant nodes using the graph
  for (unsigned int i = 0; i< nodes_.size(); ++i) {
    for (unsigned int j = 0; j < nodes_[i]->get_particle_indexes().size(); ++j) {
      if (changed.find(nodes_[i]->get_particle_indexes()[j]) != changed.end()) {
        add_node_to_active_set(i);
        break;
      }
    }
  }
  swap_active_sets();
}

void Updater::fill_queue() {
  // lazy
  for (unsigned int i = 0; i < nodes_.size(); ++i) {
    add_node_to_active_set(i);
  }
  swap_active_sets();
}

void Updater::update(unsigned int iterations) {
  fill_queue();
  for (unsigned int i = 0; i < iterations; ++i) {
    do_update();
  }
}

IMPDOMINO3_END_NAMESPACE
