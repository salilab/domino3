/** \file IMP/domino3/Factor.h
 *  \brief Store the marginal for a variable.
 */

#ifndef IMPDOMINO3_FACTOR_H
#define IMPDOMINO3_FACTOR_H

#include <IMP/domino3/domino3_config.h>
#include <IMP/ModelObject.h>
#include <IMP/object_macros.h>
#include <IMP/graph_macros.h>
#include "StatesTable.h"
#include <IMP/particle_index.h>
#include "Marginals.h"

IMPDOMINO3_BEGIN_NAMESPACE

class Factor;
IMP_OBJECTS(Factor, Factors);

struct FactorEdge {
    Factor * from_;
    Factor * to_;
    FactorEdge(Factor *from, Factor *to) : from_(from), to_(to){}
};
typedef std::vector<FactorEdge> FactorEdges;


/** Factor updates its marginals based on some criteria. */
class IMPDOMINO3EXPORT Factor: public ModelObject {
  Vector<MarginalsList> inputs_;
  ParticleIndexes pis_;
  MarginalsList mine_;
  FactorsTemp neighbors_;
  unsigned int index_;
 protected:
  /** Factor the marginals based on my data. All the marginals have been Factord
      from the inputs (via averaging). */
  virtual void do_update() = 0;

  virtual ModelObjectsTemp do_get_outputs() const IMP_OVERRIDE {
    return ModelObjectsTemp();
  }

  virtual ModelObjectsTemp do_get_inputs() const IMP_OVERRIDE {
    return get_particles(get_model(), get_particle_indexes());
  }

 public:
  Factor(Model *m,
       const ParticleIndexes &pis,
       StatesTable *pst,
       std::string name);

  /** Factor my marginals. */
  void update();

  //! always sorted
  const ParticleIndexes &get_particle_indexes() const { return pis_; }

  const MarginalsList& get_marginals() const { return mine_; }

  const FactorsTemp& get_neighbors() const {return neighbors_;}

  unsigned int get_index() const {return index_;}

  void set_index(unsigned int i) {index_ = i;}

  //! make sure to call it on the neighbor too
  void add_neighbor(Factor *n);
  // match the particles and set intput
  void set_matching_inputs(Factor *n);
};

IMPDOMINO3EXPORT void add_neighbors(const FactorsTemp &factors);
IMPDOMINO3EXPORT void add_neighbors_by_factor_edges(const FactorEdges &factor_edges);
IMPDOMINO3EXPORT void print_graph(const FactorsTemp &factors);
IMPDOMINO3EXPORT void update_state_table(const FactorsTemp &factors,const StatesTable *pst);



IMP_GRAPH(FactorGraph, undirected, Pointer<Factor>,
          ParticleIndexes,
          out << vertex->get_name() << "\\n"
          << "[" << vertex->get_type_name() << ": "
          << vertex->get_index() << "]");

IMPDOMINO3_END_NAMESPACE

#endif /* IMPDOMINO3_FACTOR_H */
