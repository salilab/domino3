#ifndef IMPDOMINO3_DISTANCE_FACTOR_H
#define IMPDOMINO3_DISTANCE_FACTOR_H
#include <IMP/domino3/domino3_config.h>
#include <IMP/kernel/ModelObject.h>
#include "Factor.h"
#include <IMP/base/object_macros.h>
#include <IMP/base/graph_macros.h>
#include <IMP/base/map.h>
#include <IMP/kernel/particle_index.h>
#include "Marginals.h"

IMPDOMINO3_BEGIN_NAMESPACE

/** A factor for a single distance restarint. */
class IMPDOMINO3EXPORT DistanceFactor: public Factor {
  double distance_, allowed_error_;
  IntPairs allowed_states_;
  kernel::ParticleIndexPair pis;
  StatesTable *pst;
  double distance;
  double allowed_error;
  double * distances;
  kernel::Model *m;
 protected:
  virtual void do_update() IMP_OVERRIDE;
private:
  double distance_to_probability(double x);

 public:
  DistanceFactor(kernel::Model *m,
               const kernel::ParticleIndexPair &pis,
               double distance, double allowed_error,
               StatesTable *pst);

  IMP_OBJECT_METHODS(DistanceFactor);
};

IMPDOMINO3_END_NAMESPACE

#endif // IMPDOMINO3_DISTANCE_FACTOR_H
