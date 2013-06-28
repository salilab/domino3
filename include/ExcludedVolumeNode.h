/** \file domino3/ExcludedVolumeNode.h
 *  \brief Store the marginal for a variable.
 */

#ifndef IMPDOMINO3_EXCLUDED_VOLUME_NODE_H
#define IMPDOMINO3_EXCLUDED_VOLUME_NODE_H
#include <IMP/domino3/domino3_config.h>
#include <IMP/kernel/ModelObject.h>
#include <IMP/base/object_macros.h>
#include <IMP/base/graph_macros.h>
#include <IMP/base/map.h>
#include <IMP/kernel/particle_index.h>
#include "Marginals.h"

IMPDOMINO3_BEGIN_NAMESPACE

/** A node excluded volume on a pair of nodes. They are assumed to share
    the same space of assignments. */
class IMPDOMINO3EXPORT ExcludedVolumeNode: public Node {
 protected:
  virtual void do_update() IMP_OVERRIDE;
 public:
  ExcludedVolumeNode(Model *m,
                     const kernel::ParticleIndexPair &pis,
                     domino::ParticleStatesTable *pst);

  IMP_OBJECT_METHODS(ExcludedVolumeNode);
};

IMPDOMINO3_END_NAMESPACE

#endif // IMPDOMINO3_EXCLUDED_VOLUME_NODE_H
