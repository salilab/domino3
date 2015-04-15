/** \file IMP/domino3/Probability2DFactor.h
    \brief A factor for a single distance restraint.
 */

#ifndef IMPDOMINO3_PROBABILITY2DFACTOR_H
#define IMPDOMINO3_PROBABILITY2DFACTOR_H

#include <IMP/domino3/domino3_config.h>
#include <IMP/ModelObject.h>
#include "Factor.h"
#include <IMP/base/object_macros.h>
#include <IMP/base/graph_macros.h>
#include <IMP/particle_index.h>
#include "Marginals.h"

IMPDOMINO3_BEGIN_NAMESPACE

/** A factor for a single distance restraint. */
class IMPDOMINO3EXPORT Probability2DFactor: public Factor {
  ParticleIndexPair pis_;
  StatesTable *pst_;
  boost::shared_array<FP> log_probability_;
 protected:
  virtual void do_update() IMP_OVERRIDE;
private:

 public:
  Probability2DFactor(Model *m,
          const ParticleIndexPair &pis,
          StatesTable *pst,
          boost::shared_array<FP> &log_probability);

  IMP_OBJECT_METHODS(Probability2DFactor);
};

IMPDOMINO3_END_NAMESPACE

#endif /* IMPDOMINO3_PROBABILITY2DFACTOR_H */
