/** \file domino3/Updater.h
 *  \brief Store the marginal for a variable.
 */

#ifndef IMPDOMINO3_UPDATER_H
#define IMPDOMINO3_UPDATER_H
#include <IMP/domino3/domino3_config.h>
#include <IMP/base/Object.h>
#include "Node.h"
#include <IMP/base/object_macros.h>
#include <IMP/base/map.h>
#include <IMP/base/particle_index.h>
#include "Marginals.h"

IMPDOMINO3_BEGIN_NAMESPACE

/** Updater its marginals based on some criteria. */
class IMPDOMINO3EXPORT Updater: public base::Object {
  Nodes nodes_;

 protected:
  /** Updater the marginals based on my data. All the marginals have been updaterd
      from the inputs (via averaging). */
  virtual void do_update() = 0;

 public:
  Updater(const NodesTemp &nodes,
          std::string name);

  void update(unsigned int iterations);

  IMP_OBJECT_METHODS(Updater);
};

IMPDOMINO3_END_NAMESPACE

#endif // IMPDOMINO3_UPDATER_H
