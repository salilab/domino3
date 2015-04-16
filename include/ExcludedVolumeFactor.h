/** \file IMP/domino3/ExcludedVolumeFactor.h
 *  \brief Store the marginal for a variable.
 */

#ifndef IMPDOMINO3_EXCLUDED_VOLUME_FACTOR_H
#define IMPDOMINO3_EXCLUDED_VOLUME_FACTOR_H

#include <IMP/domino3/domino3_config.h>
#include <IMP/ModelObject.h>
#include <IMP/object_macros.h>
#include <IMP/graph_macros.h>
#include <IMP/particle_index.h>
#include "Marginals.h"

IMPDOMINO3_BEGIN_NAMESPACE
/** A factor excluded volume on a pair of factors. They are assumed to share
    the same space of assignments. */
class IMPDOMINO3EXPORT ExcludedVolumeFactor: public Factor {
 protected:
  virtual void do_update() IMP_OVERRIDE;
 public:
  ExcludedVolumeFactor(Model *m,
                     const ParticleIndexPair &pis,
                     StatesTable *pst);

  IMP_OBJECT_METHODS(ExcludedVolumeFactor);
};

IMPDOMINO3_END_NAMESPACE

#endif /* IMPDOMINO3_EXCLUDED_VOLUME_FACTOR_H */
