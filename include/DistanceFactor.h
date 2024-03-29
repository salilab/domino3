/** \file IMP/domino3/DistanceFactor.h
 */

#ifndef IMPDOMINO3_DISTANCE_FACTOR_H
#define IMPDOMINO3_DISTANCE_FACTOR_H

#include <IMP/domino3/domino3_config.h>
#include <IMP/ModelObject.h>
#include "Factor.h"
#include <IMP/object_macros.h>
#include <IMP/graph_macros.h>
#include <IMP/particle_index.h>
#include "Marginals.h"

IMPDOMINO3_BEGIN_NAMESPACE

/** A factor for a single distance restraint. */
class IMPDOMINO3EXPORT DistanceFactor: public Factor {
  FP distance_, allowed_error_;
  IntPairs allowed_states_;
  ParticleIndexPair pis;
  StatesTable *pst;
  FP distance;
  FP allowed_error;
  FP * distances;
  Model *m;
 protected:
  virtual void do_update() IMP_OVERRIDE;
private:
  FP distance_to_probability(FP x);

 public:
  DistanceFactor(Model *m,
               const ParticleIndexPair &pis,
               FP distance, FP allowed_error,
               StatesTable *pst);

  IMP_OBJECT_METHODS(DistanceFactor);
};

IMPDOMINO3_END_NAMESPACE

#endif /* IMPDOMINO3_DISTANCE_FACTOR_H */
