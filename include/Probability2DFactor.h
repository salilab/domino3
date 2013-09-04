#ifndef IMPDOMINO3_PROBABILITY2DFACTOR_FACTOR_H
#define IMPDOMINO3_PROBABILITY2DFACTOR_FACTOR_H
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
class IMPDOMINO3EXPORT Probability2DFactor: public Factor {
  kernel::ParticleIndexPair pis_;
  StatesTable *pst_;
  boost::shared_array<double> probability_;
 protected:
  virtual void do_update() IMP_OVERRIDE;
private:

 public:
  Probability2DFactor(kernel::Model *m,
          const kernel::ParticleIndexPair &pis,
          StatesTable *pst,
          boost::shared_array<double> &probability);

  IMP_OBJECT_METHODS(Probability2DFactor);
};

IMPDOMINO3_END_NAMESPACE

#endif // IMPDOMINO3_SEA_FACTOR_H

