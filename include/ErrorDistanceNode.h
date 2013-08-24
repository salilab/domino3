/** \file domino3/DistanceNode.h
 *  \brief Store the marginal for a variable.
 */

#ifndef IMPDOMINO3_ERROR_DISTANCE_NODE_H
#define IMPDOMINO3_ERROR_DISTANCE_NODE_H
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
class IMPDOMINO3EXPORT ErrorDistanceNode: public Node {
  double distance_, allowed_error_;
  IntPairs allowed_states_;
  kernel::ParticleIndexPair pis;
  StatesTable *pst;
  double distance;
  double allowed_error;
  kernel::Model *m;
 protected:
  virtual void do_update() IMP_OVERRIDE;
private:
  double distance_to_probability(double x);

 public:
  ErrorDistanceNode(kernel::Model *m,
               const kernel::ParticleIndexPair &pis,
               double distance, double allowed_error,
               StatesTable *pst);

  IMP_OBJECT_METHODS(ErrorDistanceNode);
};

IMPDOMINO3_END_NAMESPACE

#endif // IMPDOMINO3_ERROR_DISTANCE_NODE_H
