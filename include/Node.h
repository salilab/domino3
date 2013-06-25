/** \file domino3/Node.h
 *  \brief Store the marginal for a variable.
 */

#ifndef IMPDOMINO3_NODE_H
#define IMPDOMINO3_NODE_H
#include <IMP/domino3/domino3_config.h>
#include <IMP/kernel/ModelObject.h>
#include <IMP/base/object_macros.h>
#include <IMP/base/map.h>
#include <IMP/base/particle_index.h>
#include "Marginals.h"

IMPDOMINO3_BEGIN_NAMESPACE

/** Node its marginals based on some criteria. */
class IMPDOMINO3EXPORT Node: public kernel::ModelObject {
  base::Vector<MarginalsList> inputs_;
  kernel::ParticleIndexes pis_;
  Marginals mine_;

 protected:
  /** Node the marginals based on my data. All the marginals have been noded
      from the inputs (via averaging). */
  virtual void do_update() = 0;

 public:
  Node(Model *m,
       const ParticleIndexes &pis,
       domino::ParticleStatesTable *pst,
       std::string name);

  /** Node my marginals. */
  void update();

  const ParticleIndexes &get_particle_indexes() const { return pis_; }

  const Marginals& get_marginals() const { return mine_; }

  void add_input_marginal(ParticleIndex pi, Marginal *m) {
    inputs[pi].push_back(m);
  }

  IMP_OBJECT_METHODS(Node);
};

IMP_OBJECTS(Node, Nodes);

IMPDOMINO3_END_NAMESPACE

#endif // IMPDOMINO3_NODE_H
