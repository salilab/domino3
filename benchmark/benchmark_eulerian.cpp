/**
 * Copyright 2007-2013 IMP Inventors. All rights reserved.
 */

#include <IMP/base/flags.h>
#include <IMP/kernel/Model.h>
#include <IMP/atom/Hierarchy.h>
#include <IMP/atom/hierarchy_tools.h>
#include <IMP/em/DensityMap.h>
#include <IMP/atom/pdb.h>
#include <IMP/base/SetLogState.h>
#include <IMP/core/XYZ.h>
#include <algorithm>

#define MY_LOG_VERBOSE(x) std::cout << x;
#define MY_LOG_TERSE(x) std::cout << x;

namespace {
double cell_size = 1;
typedef IMP::algebra::GridD<3,
                            IMP::algebra::DenseGridStorageD<3, float>,
                            float> Grid;
typedef Grid::Index GridIndex;
typedef IMP::base::Vector<GridIndex> States;

typedef boost::scoped_array<boost::scoped_array<double> >  Probabilities;
typedef IMP::base::Vector<IMP::base::Vector<int> > DistanceIndex;
typedef IMP::base::set<int> DirtyRestraints;

IMP_NAMED_TUPLE_3(Distance, Distances,
                  int, first,
                  int, second,
                  double, distance,);

void fill_one(const Grid &grid,
              const IMP::algebra::BoundingBox3D &bb,
              const IMP::base::map<GridIndex, int> &state_map,
              unsigned int index,
              Probabilities &probabilities) {
  States cur_states(grid.indexes_begin(bb), grid.indexes_end(bb));
  if (cur_states.size() < state_map.size()) {
    MY_LOG_TERSE("Particle " << index << " has initial states " << cur_states
                 << std::endl);
  }
  double cur_prob = 1.0/ cur_states.size();
  for (unsigned int i = 0; i < cur_states.size(); ++i) {
    probabilities[index][state_map.find(cur_states[i])->second]
      = cur_prob;
  }
}

IMP::base::Vector<double> get_entropies(unsigned int num_variables,
                                        unsigned int num_states,
                                        const Probabilities &probabilities) {
  IMP::base::Vector<double> ret(num_variables, 0);
  for (unsigned int i = 0; i < num_variables; ++i) {
    for (unsigned int j = 0; j < num_states; ++j) {
      double cur_prob = probabilities[i][j];
      if (cur_prob > 0) {
        double cur = cur_prob * std::log(cur_prob);
        ret[i] -= cur;
      }
    }
  }
  return ret;
}

/** Be lazy and pin the first 4. I don't think it matters.
*/
void initialize_probabilities(IMP::kernel::Model *m,
                              const IMP::kernel::ParticleIndexes &pis,
                              const Grid &grid,
                              const States &states,
                              Probabilities &probabilities) {
  IMP::base::map<GridIndex, int> state_map;
  for (unsigned int i = 0; i< states.size(); ++i) {
    state_map[states[i]] = i;
  }
  probabilities.reset(new boost::scoped_array<double>[pis.size()]);
  for (unsigned int i = 0; i < pis.size(); ++i) {
    probabilities[i].reset(new double[states.size()]);
    std::fill(probabilities[i].get(), probabilities[i].get() + states.size(),
              0.0);
  }
  for (unsigned int i = 0; i < 4; ++ i) {
    IMP::algebra::Vector3D coords = IMP::core::XYZ(m, pis[i]).get_coordinates();
    fill_one(grid,
             IMP::algebra::BoundingBox3D(coords, coords),
             state_map,
             i, probabilities);
  }
  for (unsigned int i = 4; i < pis.size(); ++i) {
    fill_one(grid,
             grid.get_bounding_box(),
             state_map,
             i, probabilities);

  }
}

void refine_one(const Distance &d, const Grid &grid,
            const States &states, const DistanceIndex &distance_index,
            Probabilities &probabilities,
            DirtyRestraints &dirty_restraints) {
  // cheat on it being symmetric
  IMP::base::Vector<double> accum_first(states.size(), 0);
  IMP::base::Vector<double> accum_second(states.size(), 0);
  for (unsigned int i = 0; i < states.size(); ++i) {
    double pi0 = probabilities[d.get_first()][i];
    if (pi0 < .001) continue;
    for (unsigned int j = 0; j < states.size(); ++j) {
      double pj0 = probabilities[d.get_second()][j];
      if (pj0 < .001) continue;
      double distance = IMP::algebra::get_distance(grid.get_center(states[i]),
                                                   grid.get_center(states[j]));
      if (std::abs(distance - d.get_distance()) < 1.4 * cell_size) {
        accum_first[i] += pi0 * pj0;
        accum_second[j] += pi0 * pj0;
      }
    }
  }
  double total_first = std::accumulate(accum_first.begin(),
                                       accum_first.end(),
                                       0.0);
  double total_second = std::accumulate(accum_second.begin(),
                                        accum_second.end(),
                                        0.0);
  for (unsigned int i = 0; i < accum_first.size(); ++i) {
    probabilities[d.get_first()][i] = accum_first[i]/total_first;
    probabilities[d.get_second()][i] = accum_second[i]/total_second;
  }
  dirty_restraints.insert(distance_index[d.get_first()].begin(),
                          distance_index[d.get_first()].end());
  dirty_restraints.insert(distance_index[d.get_second()].begin(),
                          distance_index[d.get_second()].end());
}

void refine(const Distances &distances, const Grid &grid,
            const States &states, const DistanceIndex &distance_index,
            Probabilities &probabilities,
            DirtyRestraints &dirty_restraints) {
  IMP::base::Vector<int> my_restraints(dirty_restraints.begin(),
                                       dirty_restraints.end());
  dirty_restraints.clear();
  for (unsigned int i = 0; i < my_restraints.size(); ++i) {
    refine_one(distances[my_restraints[i]], grid, states, distance_index,
               probabilities, dirty_restraints);
  }
}

void write_probabilities(unsigned int num_states,
                         const Probabilities &probs, Grid &grid,
                         const States &states, std::string prefix) {
  MY_LOG_TERSE("Writing probabilites to " << prefix << std::endl);
    
  for (unsigned int i = 0; i < num_states; ++i) {
      for (unsigned int j = 0; j < states.size(); ++j) {
          std::cout << probs[i][j] << "\t";
      }
      std::cout << std::endl;
  }
  for (unsigned int i = 0; i < num_states; ++i) {
      
     auto it= std::max_element(probs[i].get(), probs[i].get()+states.size());
      std::cout << *it << std::endl;
      std::cout << std::distance(probs[i].get(), it) << std::endl;
      std::cout << states[std::distance(probs[i].get(), it)] <<std::endl;
  }
    
  for (unsigned int i = 0; i < num_states; ++i) {
    for (unsigned int j = 0; j < states.size(); ++j) {
      grid[states[j]] = probs[i][j];
    }
    std::ostringstream oss;
    oss << prefix << "." << i << ".mrc";
    IMP::base::OwnerPointer<IMP::em::DensityMap> map
    = IMP::em::create_density_map(grid);
    IMP::em::write_map(map, oss.str());
  }
}

boost::int64_t iterations = 100;
IMP::base::AddIntFlag iteration_adder("iterations",
                       "The number of times to iterate through the restraints",
                        &iterations);

IMP::base::AddFloatFlag iteration_cellsize("cellsize",
                       "The cell size of the grid.",
                        &cell_size);


std::string input = IMP::atom::get_example_path("1d3d-protein.pdb");
IMP::base::AddStringFlag input_adder("input",
                       "The input pdb to use",
                        &input);
}


int main(int argc, char **argv) {
  IMP::base::setup_from_argv(argc, argv, "Experiment with loopy domino");
  IMP_NEW(IMP::kernel::Model, m, ());
    IMP::atom::Hierarchy pdb = IMP::atom::read_pdb(input, m, new IMP::atom::CAlphaPDBSelector());
  // add coordinates to residues
  {
    IMP::base::SetLogState sll(IMP::base::SILENT);
    IMP::atom::create_rigid_body(pdb);
  }
    IMP::atom::Atoms atoms
    = IMP::atom::get_by_type(pdb, IMP::atom::ATOM_TYPE);
  // start with all pairs
  Distances distances;
  DistanceIndex distance_index(atoms.size());
  for (unsigned int i = 0; i< atoms.size(); ++i) {
    for (unsigned int j = 0; j < i; ++j) {
      double dist = IMP::core::get_distance(IMP::core::XYZ(atoms[i]),
                                            IMP::core::XYZ(atoms[j]));
      Distance cur(i, j, dist);
      distances.push_back(cur);
      distance_index[i].push_back(distances.size() - 1);
      distance_index[j].push_back(distances.size() - 1);
    }
  }

  IMP::algebra::BoundingBox3D bb = IMP::atom::get_bounding_box(pdb);
  Grid grid(cell_size, bb, 0.0);
  States states;
  std::copy(grid.indexes_begin(grid.get_bounding_box()),
            grid.indexes_end(grid.get_bounding_box()),
            std::back_inserter(states));
  Probabilities probabilities;
  IMP::kernel::ParticleIndexes pis(atoms.size());
  for (unsigned int i = 0; i< atoms.size(); ++i) {
    pis[i] = atoms[i].get_particle_index();
  }
  initialize_probabilities(m, pis, grid, states, probabilities);
  DirtyRestraints dirty_restraints;
  dirty_restraints.insert(distance_index[0].begin(),
                          distance_index[0].end());
  dirty_restraints.insert(distance_index[1].begin(),
                          distance_index[1].end());
  dirty_restraints.insert(distance_index[2].begin(),
                          distance_index[2].end());
  dirty_restraints.insert(distance_index[3].begin(),
                          distance_index[3].end());

    IMP::base::Vector<double> entropies
    = get_entropies(atoms.size(), states.size(), probabilities);
    MY_LOG_VERBOSE("entropies are " << entropies << std::endl);
    MY_LOG_VERBOSE("And the total is " << std::accumulate(entropies.begin(),
                                                          entropies.end(), 0.0)
                   << std::endl);

  for (unsigned int i = 0; i < iterations; ++i) {
    refine(distances, grid, states, distance_index,
           probabilities, dirty_restraints);
    IMP::base::Vector<double> entropies
    = get_entropies(atoms.size(), states.size(), probabilities);
    MY_LOG_VERBOSE("entropies are " << entropies << std::endl);
    MY_LOG_VERBOSE("And the total is " << std::accumulate(entropies.begin(),
                                                          entropies.end(), 0.0)
                   << std::endl);
  }

  write_probabilities(atoms.size(), probabilities, grid,
                      states, "final");
  return 0;
}
