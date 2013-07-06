#include <IMP/domino3/Node.h>

IMPDOMINO3_BEGIN_NAMESPACE

Node::Node(Model *m,
           const ParticleIndexes &pis,
           StatesTable *pst,
           std::string name):
  ModelObject(m, name), pis_(pis) {
  std::sort(pis_.begin(), pis_.end());
  mine_.resize(pis_.size());
  inputs_.resize(pis_.size());
  for (unsigned int i = 0; i < pis_.size(); ++i) {
    mine_[i] = new Marginals(m, pis_[i],
                             pst->get_states(pis_[i])->get_number());
    mine_[i]->update_current_from_list(MarginalsList(1,
                             pst->get_marginals(pis_[i])));
  }
}

void Node::update() {
  for (unsigned int i = 0; i < pis_.size(); ++i) {
    mine_[i]->update_current_from_list(inputs_[i]);
  }

  do_update();

  for (unsigned int i = 0; i < pis_.size(); ++i) {
    mine_[i]->update_current_from_next();
  }
}

void Node::add_neighbor(Node *n) {
  ParticleIndexes pis = n->get_particle_indexes();
  for (unsigned int i = 0; i < pis.size(); ++i) {
    kernel::ParticleIndexes::const_iterator it = std::find(pis_.begin(),
                                                           pis_.end(),
                                                           pis[i]);
    if (it != pis_.end()) {
      unsigned int offset = it - pis_.begin();
      inputs_[offset].push_back(n->get_marginals()[i]);
    }
  }
  neighbors_.push_back(n);
}

NodeGraph get_node_graph(const NodesTemp &nodes) {
  NodeGraph ret(nodes.size());
  NodeGraphVertexName vm = boost::get(boost::vertex_name, ret);
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
        boost::add_edge(i, j, intersection, ret);
      }
    }
  }
  return ret;
}

IMPDOMINO3_END_NAMESPACE
