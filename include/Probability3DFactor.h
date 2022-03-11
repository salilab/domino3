/** \file IMP/domino3/Probability3DFactor.h
    \brief A factor for a single distance restraint.
 */

#ifndef IMPDOMINO3_PROBABILITY3D_FACTOR_H
#define IMPDOMINO3_PROBABILITY3D_FACTOR_H

#include <IMP/domino3/domino3_config.h>
#include <IMP/ModelObject.h>
#include "Factor.h"
#include <IMP/object_macros.h>
#include <IMP/graph_macros.h>
#include <IMP/particle_index.h>
#include <IMP/domino3/Probability3D.h>

#include "Marginals.h"

IMPDOMINO3_BEGIN_NAMESPACE

/** A factor for a single distance restraint. */
class IMPDOMINO3EXPORT Probability3DFactor: public Factor {
  ParticleIndexTriplet pis_;
  StatesTable *pst_;
  Probability3D * log_probability;
 protected:
  virtual void do_update() override;
private:

 public:
  Probability3DFactor(Model *m,
          const ParticleIndexTriplet &pis,
          StatesTable *pst,
          Probability3D * log_probability);

  IMP_OBJECT_METHODS(Probability3DFactor);
};

IMPDOMINO3_END_NAMESPACE

#endif /* IMPDOMINO3_PROBABILITY3D_FACTOR_H */
