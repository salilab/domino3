#ifndef IMPDOMINO3_SEA_FACTOR_H
#define IMPDOMINO3_SEA_FACTOR_H
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
class IMPDOMINO3EXPORT ProbabilityFactor: public Factor {
  kernel::ParticleIndexPair pis;
  StatesTable *pst;
  double * probability; 
 protected:
  virtual void do_update() IMP_OVERRIDE;
private:

 public:
  ProbabilityFactor(kernel::Model *m,
          const kernel::ParticleIndexPair &pis,
          StatesTable *pst,double * probability);

  IMP_OBJECT_METHODS(ProbabilityFactor);
};

IMPDOMINO3_END_NAMESPACE

#endif // IMPDOMINO3_SEA_FACTOR_H

