/** \file domino3/Updater.h
 *  \brief Store the marginal for a variable.
 */

#ifndef IMPDOMINO3_UPDATER_H
#define IMPDOMINO3_UPDATER_H
#include <IMP/domino3/domino3_config.h>
#include <IMP/kernel/ModelObject.h>
#include <IMP/base/object_macros.h>
#include <IMP/base/map.h>
#include <IMP/base/particle_index.h>
#include "Marginals.h"

IMPDOMINO3_BEGIN_NAMESPACE

/** Update its marginals based on some criteria. */
class IMPDOMINO3EXPORT Updater: public kernel::ModelObject {
  base::Vector<MarginalsList> inputs_;
  kernel::ParticleIndexes pis_;
  Marginals mine_;
 protected:
  /** Update the marginals based on my data. All the marginals have been updated
      from the inputs (via averaging). */
  virtual void do_update() = 0;
 public:
  Updater(Model *m,
          const ParticleIndexes &pis,
          domino::ParticleStatesTable *pst,
          std::string name);

  /** Update my marginals. */
  void update();

  const ParticleIndexes &get_particle_indexes() const { return pis_; }

  const Marginals& get_marginals() const { return mine_; }

  void add_input_marginal(ParticleIndex pi, Marginal *m) {
    inputs[pi].push_back(m);
  }

  IMP_OBJECT_METHODS(Marginals);
};

IMPDOMINO3_END_NAMESPACE

#endif // IMPDOMINO3_UPDATER_H
