/** \file domino3/Node.h
 *  \brief Store the marginal for a variable.
 */

#ifndef IMPDOMINO3_NODE_H
#define IMPDOMINO3_NODE_H
#include <IMP/domino3/domino3_config.h>
#include <IMP/kernel/ModelObject.h>
#include <IMP/base/object_macros.h>
#include <IMP/base/graph_macros.h>
#include <IMP/base/map.h>
#include "StatesTable.h"
#include <IMP/kernel/particle_index.h>
#include "Marginals.h"

IMPDOMINO3_BEGIN_NAMESPACE

class Node;
IMP_OBJECTS(Node, Nodes);

/** Node updates its marginals based on some criteria. */
class IMPDOMINO3EXPORT Node: public kernel::ModelObject {
  base::Vector<MarginalsList> inputs_;
  kernel::ParticleIndexes pis_;
  MarginalsList mine_;
  NodesTemp neighbors_;
  unsigned int index_;
 protected:
  /** Node the marginals based on my data. All the marginals have been noded
      from the inputs (via averaging). */
  virtual void do_update() = 0;

  virtual kernel::ModelObjectsTemp do_get_outputs() const IMP_OVERRIDE {
    return kernel::ModelObjectsTemp();
  }

  virtual kernel::ModelObjectsTemp do_get_inputs() const IMP_OVERRIDE {
    return kernel::get_particles(get_model(), get_particle_indexes());
  }

 public:
  Node(kernel::Model *m,
       const kernel::ParticleIndexes &pis,
       StatesTable *pst,
       std::string name);

  /** Node my marginals. */
  void update();

  //! always sorted
  const kernel::ParticleIndexes &get_particle_indexes() const { return pis_; }

  const MarginalsList& get_marginals() const { return mine_; }

  const NodesTemp& get_neighbors() const {return neighbors_;}

  unsigned int get_index() const {return index_;}

  void set_index(unsigned int i) {index_ = i;}

  //! make sure to call it on the neighbor too
  void add_neighbor(Node *n);
};

IMP_GRAPH(NodeGraph, undirected, base::Pointer<Node>,
          kernel::ParticleIndexes,
          out << vertex->get_name() << "\\n"
          << "[" << vertex->get_type_name() << ": "
          << vertex->get_index() << "]");

IMPDOMINO3EXPORT NodeGraph get_node_graph(const NodesTemp &nodes);

IMPDOMINO3_END_NAMESPACE

#endif // IMPDOMINO3_NODE_H
