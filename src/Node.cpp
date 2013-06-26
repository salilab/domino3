#include <IMP/domino3/Node.h>

IMPDOMINO3_BEGIN_NAMESPACE

Node::Node(Model *m,
                 const ParticleIndexes &pis,
                 domino::ParticleStatesTable *pst,
                 std::string name):
  ModelObject(m, name), pis_(pis) {
  std::sort(pis_.begin(), pis_.end());
  mine_.resize(pis_.size());
  inputs_.resize(pis_.size());
  for (unsigned int i = 0; i < pis_.size(); ++i) {
    mine_[i] = new Marginals(m, pis_[i],
                          pst->get_number_of_states(m->get_particle(pis_[i])));
  }
}

void Node::add_input_marginal(ParticleIndex pi, Marginal *m) {
  inputs[pi].push_back(m);
}

void Node::update() {
  for (unsigned int i = 0; i < pis_.size(); ++i) {
    mine_[i]->update_current_from_list(inputs_[i]);
  }

  do_update();

  for (unsigned int i = 0; i < pis_.size(); ++i) {
    mine_[i]->update_current_next();
  }
}

NodeGraph get_node_graph(const NodesTemp &nodes) {
  NodeGraph ret(nodes.size());
  NodeGraphVertexMap vm = boost::get(ret, boost::vertex_name);
  for (unsigned int i = 0; i < nodes.size(); ++i) {
    vm[i] = nodes[i];
    // create self edge too
    for (unsigned int j = 0; j <= i; ++j) {
      kernel::ParticleIndexes intersection;
      std::set_intersection(nodes[i]->get_particle_indexes().begin(),
                            nodes[i]->get_particle_indexes().end(),
                            nodes[j]->get_particle_indexes().begin(),
                            nodes[j]->get_particle_indexes().end(),
                            std::back_inserter(intersection));
      if (!intersection.empty()) {
        boost::add_edge(ret, i, j, intersection);
      }
    }
  }
  return ret;
}

IMPDOMINO3_END_NAMESPACE
