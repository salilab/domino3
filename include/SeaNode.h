#ifndef IMPDOMINO3_SEA_NODE_H
#define IMPDOMINO3_SEA_NODE_H
#include <IMP/domino3/domino3_config.h>
#include <IMP/kernel/ModelObject.h>
#include "Node.h"
#include <IMP/base/object_macros.h>
#include <IMP/base/graph_macros.h>
#include <IMP/base/map.h>
#include <IMP/kernel/particle_index.h>
#include "Marginals.h"

IMPDOMINO3_BEGIN_NAMESPACE

/** A node for a single distance restarint. */
class IMPDOMINO3EXPORT SeaNode: public Node {
  kernel::ParticleIndexPair pis;
  StatesTable *pst;
 protected:
  virtual void do_update() IMP_OVERRIDE;
private:

 public:
  SeaNode(kernel::Model *m,
          const kernel::ParticleIndexPair &pis,
          StatesTable *pst);

  IMP_OBJECT_METHODS(SeaNode);
};

IMPDOMINO3_END_NAMESPACE

#endif // IMPDOMINO3_SEA_NODE_H

