#ifndef IMPDOMINO3_PROBABILITY3D_FACTOR_H
#define IMPDOMINO3_PROBABILITY3D_FACTOR_H
#include <IMP/domino3/domino3_config.h>
#include <IMP/kernel/ModelObject.h>
#include "Factor.h"
#include <IMP/base/object_macros.h>
#include <IMP/base/graph_macros.h>
#include <IMP/base/map.h>
#include <IMP/kernel/particle_index.h>
#include <IMP/domino3/Probability3D.h>

#include "Marginals.h"

IMPDOMINO3_BEGIN_NAMESPACE

/** A factor for a single distance restarint. */
class IMPDOMINO3EXPORT Probability3DFactor: public Factor {
  kernel::ParticleIndexTriplet pis_;
  StatesTable *pst_;
  Probability3D * log_probability;
 protected:
  virtual void do_update() IMP_OVERRIDE;
private:

 public:
  Probability3DFactor(kernel::Model *m,
          const kernel::ParticleIndexTriplet &pis,
          StatesTable *pst,
          Probability3D * log_probability);

  IMP_OBJECT_METHODS(Probability3DFactor);
};

IMPDOMINO3_END_NAMESPACE

#endif // IMPDOMINO3_SEA_FACTOR_H

