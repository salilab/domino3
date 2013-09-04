/**
 * Copyright 2007-2013 IMP Inventors. All rights reserved.
 */

#include <IMP/base/flags.h>
#include <IMP/kernel/Model.h>
#include <IMP/atom/Hierarchy.h>
#include <IMP/atom/hierarchy_tools.h>
#include <IMP/atom/pdb.h>


#include <IMP/rmf/atom_io.h>
#include <IMP/rmf/restraint_io.h>
#include <IMP/rmf/frames.h>
#include <IMP/rmf/particle_io.h>

#include <IMP/base/SetLogState.h>
#include <IMP/domino3/Factor.h>
#include <IMP/domino3/Probability2DFactor.h>
#include <IMP/domino3/Probability3DFactor.h>
#include <IMP/domino3/ExcludedVolumeFactor.h>


#include <IMP/domino3/Updater.h>
#include <RMF/FileHandle.h>
#include <IMP/domino3/IndexStates.h>
#include <IMP/core/XYZ.h>
#include <algorithm>

namespace {
    static int enzym_id_counter = 0;
    struct Enzyme {
        std::string enzyme_name;
        const unsigned id = enzym_id_counter++;
        Enzyme(std::string name) : enzyme_name(name){}
    };
    
    static int ligand_id_counter = 0;
    struct Ligand {
        std::string ligand_name;
        const unsigned id = ligand_id_counter++;
        Ligand(std::string name) : ligand_name(name){}
    };
    typedef std::pair<Enzyme *,Enzyme *> EnzymeEnzymePair;
    typedef std::pair<Enzyme *,Ligand *> EnzymeLigandPair;
    typedef std::pair<Ligand *,Ligand *> LigandLigandPair;
    typedef std::pair<Enzyme *,LigandLigandPair> EnzymeLigandLigandPair;

    typedef IMP::base::map<EnzymeEnzymePair, double > EnzymeEnzymeScoreLookup;
    typedef IMP::base::map<EnzymeLigandPair, double > EnzymeLigandScoreLookup;
    typedef IMP::base::map<EnzymeLigandLigandPair, double> EnzymeLigandLigandScoreLookup;
    typedef IMP::base::map<std::string, Enzyme *>  StringToEnzyme;
    typedef IMP::base::map<std::string, Ligand *>  StringToLigand;

    
    
    
    std::string input = IMP::atom::get_example_path("1d3d-protein.pdb");
    IMP::base::AddStringFlag asf("input", "Input file name", &input);
    std::string output = "out.rmf";
    IMP::base::AddStringFlag oasf("output", "Output rmf name", &output);
    double dist = std::numeric_limits<double>::infinity();
    IMP::base::AddFloatFlag dasf("distance", "Maximum distance for restraints",
                                 &dist);
    
    boost::int64_t iterations = 1;
    IMP::base::AddIntFlag aif("iterations", "Number of iterations",
                              &iterations);
    
    
    void read_in_sea_file(IMP::kernel::Model *m,
                           std::string path,
                           EnzymeEnzymeScoreLookup &scores,
                           StringToEnzyme &id_to_enzyme ){
        std::ifstream infile(path.c_str());
        std::string line;
        while (std::getline(infile, line))
        {
            std::istringstream iss(line);
            std::string a, b;
            double score;
            iss >> a >> b >> score;
            if(!id_to_enzyme[a])
                id_to_enzyme[a] = new Enzyme(a);
            if(!id_to_enzyme[b])
                id_to_enzyme[b] = new Enzyme(b);
            scores[std::make_pair(id_to_enzyme[a],id_to_enzyme[b])]=score;
            std::cout << a << " " << b << " " << score << std::endl;
        }
    }
    
    
    void read_in_dock_file(IMP::kernel::Model *m,
                           std::string path,
                           EnzymeLigandScoreLookup &scores,
                           StringToEnzyme &id_to_enzyme,
                           StringToLigand &id_to_ligand){
        std::ifstream infile(path.c_str());
        std::string line;
        while (std::getline(infile, line))
        {
            std::istringstream iss(line);
            std::string a, b;
            double score;
            iss >> a >> b >> score;
            if(!id_to_enzyme[a])
                id_to_enzyme[a] = new Enzyme(a);
            if(!id_to_ligand[b])
                id_to_ligand[b] = new Ligand(b);
            scores[std::make_pair(id_to_enzyme[a],id_to_ligand[b])]=score;
            std::cout << a << " " << b << " " << score << std::endl;
        }
    }
    
    
    void read_in_chem_file(IMP::kernel::Model *m,
                           std::string path,
                           EnzymeLigandLigandScoreLookup &scores,
                           StringToEnzyme &id_to_enzyme,
                           StringToLigand &id_to_ligand){
        std::ifstream infile(path.c_str());
        std::string line;
        while (std::getline(infile, line))
        {
            std::istringstream iss(line);
            std::string a, b, c;
            double score;
            iss >> a >> b >> c >> score;
            if(!id_to_enzyme[a])
                id_to_enzyme[a] = new Enzyme(a);
            if(!id_to_ligand[b])
                id_to_ligand[b] = new Ligand(b);
            if(!id_to_ligand[c])
                id_to_ligand[c] = new Ligand(c);
            LigandLigandPair ligand_ligand_pair = std::make_pair(id_to_ligand[b],id_to_ligand[c]);
            scores[std::make_pair(id_to_enzyme[a],ligand_ligand_pair)]=score;
            std::cout << a << " " << b << " " << c << " " <<  score << std::endl;
        }
    }
    
    
    void create_graph(IMP::domino3::Factors &factors,IMP::domino3::StatesTable * st,IMP::kernel::Model *m,
                      boost::shared_array<double> &sea_probability_vector,
                      boost::shared_array<double> &dock_probability_vector,
                      boost::shared_array<double> &chem_probability_vector){
        IMP::domino3::FactorEdgesTemp edges;
        IMP::kernel::ParticleIndexes pis = st->get_particle_indexes();
        IMP::kernel::ParticleIndexPair cur_pair1(pis[0], pis[1]);
        IMP_NEW(IMP::domino3::Probability2DFactor, pf1,(m, cur_pair1, st, sea_probability_vector));
        factors.push_back(pf1);
        pf1->set_was_used(true);
//        
        IMP::kernel::ParticleIndexPair cur_pair2(pis[1], pis[2]);
        IMP_NEW(IMP::domino3::Probability2DFactor, pf2,(m, cur_pair2, st, sea_probability_vector));
        factors.push_back(pf2);
        pf2->set_was_used(true);
//        IMP::domino3::FactorEdge * e1 = new IMP::domino3::FactorEdge(pf1,pf2);
//        IMP::domino3::FactorEdge * e12 = new IMP::domino3::FactorEdge(pf2,pf1);
//        edges.push_back(e1);
//        edges.push_back(e12);
//
        IMP::kernel::ParticleIndexPair cur_pair3(pis[0], pis[3]);
        IMP_NEW(IMP::domino3::Probability2DFactor, pf3,(m, cur_pair3, st, dock_probability_vector));
        factors.push_back(pf3);
        pf3->set_was_used(true);
//
//
//        
        IMP::kernel::ParticleIndexTriplet cur_triplet1(pis[0], pis[3], pis[4]);
        IMP_NEW(IMP::domino3::Probability3DFactor, pf3d1,(m, cur_triplet1, st, chem_probability_vector));
        factors.push_back(pf3d1);
        pf3d1->set_was_used(true);
        
//        IMP::domino3::FactorEdge * e2 = new IMP::domino3::FactorEdge(pf3,pf3d1);
//        IMP::domino3::FactorEdge * e3 = new IMP::domino3::FactorEdge(pf3d1,pf1);
//        IMP::domino3::FactorEdge * e4 = new IMP::domino3::FactorEdge(pf1,pf3d1);
//        edges.push_back(e2); // docking score -> enzyme
//        edges.push_back(e3);
//        edges.push_back(e4);
        
        IMP::kernel::ParticleIndexPair cur_pair5(pis[1],pis[4]);
        IMP_NEW(IMP::domino3::Probability2DFactor, pf5,(m, cur_pair5, st, dock_probability_vector));
        factors.push_back(pf5);
        pf5->set_was_used(true);
//
        IMP::kernel::ParticleIndexTriplet cur_triplet2(pis[1], pis[4], pis[5]);
        IMP_NEW(IMP::domino3::Probability3DFactor, pf3d2,(m, cur_triplet2, st, chem_probability_vector));
        factors.push_back(pf3d2);
        pf3d2->set_was_used(true);
//        IMP::domino3::FactorEdge * e5 = new IMP::domino3::FactorEdge(pf5,pf3d2);
//        IMP::domino3::FactorEdge * e6 = new IMP::domino3::FactorEdge(pf3d2,pf2);
//        IMP::domino3::FactorEdge * e7 = new IMP::domino3::FactorEdge(pf2,pf3d2);
//        IMP::domino3::FactorEdge * e61 = new IMP::domino3::FactorEdge(pf3d2,pf1);
//        IMP::domino3::FactorEdge * e71 = new IMP::domino3::FactorEdge(pf1,pf3d2);
//        edges.push_back(e5);
//        edges.push_back(e6);
//        edges.push_back(e7);
//        edges.push_back(e61);
//        edges.push_back(e71);

        
        IMP::kernel::ParticleIndexPair cur_pair7(pis[2], pis[5]);
        IMP_NEW(IMP::domino3::Probability2DFactor, pf7,(m, cur_pair7, st, dock_probability_vector));
        factors.push_back(pf7);
        pf7->set_was_used(true);
//
//        
        IMP::kernel::ParticleIndexTriplet cur_triplet3(pis[2], pis[5], pis[6]);
        IMP_NEW(IMP::domino3::Probability3DFactor, pf3d3,(m, cur_triplet3, st, chem_probability_vector));
        factors.push_back(pf3d3);
        pf3d3->set_was_used(true);
        for(int x = 0 ; x< 3; x++){
            for(int y = 0; y < x ; y++){
                IMP::kernel::ParticleIndexPair excl_pair(pis[x], pis[y]);
                IMP_NEW(IMP::domino3::ExcludedVolumeFactor, excl_factor,(m, excl_pair, st));
                factors.push_back(excl_factor);
                excl_factor->set_was_used(true);
            }
        }
//        IMP::domino3::FactorEdge * e8 = new IMP::domino3::FactorEdge(pf7,pf2);
//        IMP::domino3::FactorEdge * e9 = new IMP::domino3::FactorEdge(pf3d3,pf2);
//        IMP::domino3::FactorEdge * e10 = new IMP::domino3::FactorEdge(pf2,pf3d3);
//        edges.push_back(e8); // docking score -> enzyme
//        edges.push_back(e9);
//        edges.push_back(e10);
        
        std::cout << "Node size: " << factors.size() << std::endl;
//        IMP::domino3::add_neighbors_by_factor_edges(edges);
        IMP::domino3::add_neighbors(factors);

    }
    
    void normalize(boost::shared_array<double> &probability_vector,
                   int size_x,
                   int size_y,int size_z,
                   bool normalize_square){
        double total = 0;
        for(int x = 0; x < size_x; x++){
            double * it = probability_vector.get()+(x*size_y*size_z);
            if(normalize_square)
                total = std::accumulate(it, it + (size_z*size_y), 0.0);
            std::cout<< "X: " << x << std::endl;

            for(int y = 0; y < size_y; y++){
                std::cout<< "Y: " << y << std::endl;

                double * row = it + (y*size_y);
                if(!normalize_square)
                    total = std::accumulate(row, row + size_z, 0.0);
                for(int z = 0; z < size_z; z++){
                    std::cout << *(row+z) << " ";
                    *(row+z)/=total;
                }
                std::cout << std::endl;
            }
        }
    }
    void run_it(IMP::kernel::Model *m) {
        StringToEnzyme id_to_enzyme;
        EnzymeEnzymeScoreLookup enzyme_enzyme_scores;
        StringToLigand id_to_ligand;
        EnzymeLigandScoreLookup enzyme_ligand_scores;
        EnzymeLigandLigandScoreLookup enzyme_ligand_ligand_scores;

        read_in_sea_file(m,
                          std::string("/Users/aluucard/Documents/workspace/imp/imp/modules/domino3/examples/sea_new.txt"),
                          enzyme_enzyme_scores,
                          id_to_enzyme);
        
        read_in_dock_file(m,
                         std::string("/Users/aluucard/Documents/workspace/imp/imp/modules/domino3/examples/dock.txt"),
                         enzyme_ligand_scores,
                         id_to_enzyme,
                         id_to_ligand);
        read_in_chem_file(m,
                          std::string("/Users/aluucard/Documents/workspace/imp/imp/modules/domino3/examples/chemsim_rxn.txt"),
                          enzyme_ligand_ligand_scores,
                          id_to_enzyme,
                          id_to_ligand);
        std::cout << id_to_enzyme.size() << std::endl;
        std::vector<int> enzyme_ids;
        std::vector<int> ligand_ids;

        
        // read structure to probability_vector and normalize
        boost::shared_array<double> sea_probability_vector(new double[id_to_enzyme.size()*id_to_enzyme.size()]);
        std::fill (sea_probability_vector.get(),sea_probability_vector.get()+(id_to_enzyme.size()*id_to_enzyme.size()),0);
        for (EnzymeEnzymeScoreLookup::iterator i = enzyme_enzyme_scores.begin(); i != enzyme_enzyme_scores.end(); ++i)
        {
            EnzymeEnzymePair enzymePair=i->first;
            double score = i->second;
            Enzyme * e1 = enzymePair.first;
            Enzyme * e2 = enzymePair.second;
            if(e1->id == e2->id){
                sea_probability_vector[e1->id*id_to_enzyme.size()+e2->id] = score/10;
            }else{
                sea_probability_vector[e1->id*id_to_enzyme.size()+e2->id] = score;
            }
        }
        normalize(sea_probability_vector,1,id_to_enzyme.size(),id_to_enzyme.size(),true);

        boost::shared_array<double> dock_probability_vector(new double[id_to_enzyme.size()*id_to_ligand.size()]);
        std::fill (dock_probability_vector.get(),dock_probability_vector.get()+(id_to_enzyme.size()*id_to_ligand.size()),0);
        for (EnzymeLigandScoreLookup::iterator i = enzyme_ligand_scores.begin(); i != enzyme_ligand_scores.end(); ++i)
        {
            EnzymeLigandPair enzymePair=i->first;
            double score = i->second;
            Enzyme * e1 = enzymePair.first;
            Ligand * l1 = enzymePair.second;
            dock_probability_vector[e1->id*id_to_ligand.size()+l1->id] = score;
        }
        normalize(dock_probability_vector,1,id_to_enzyme.size(),id_to_ligand.size(),false);

        const unsigned int chem_size = id_to_enzyme.size()*id_to_ligand.size()*id_to_ligand.size();
        boost::shared_array<double> chem_probability_vector(new double[chem_size]);
        std::fill (chem_probability_vector.get(),chem_probability_vector.get()+chem_size,0);
        for (EnzymeLigandLigandScoreLookup::iterator i = enzyme_ligand_ligand_scores.begin(); i != enzyme_ligand_ligand_scores.end(); ++i)
        {
            EnzymeLigandLigandPair enzymePair=i->first;
            double score = i->second;
            Enzyme * e1 = enzymePair.first;
            LigandLigandPair ligand_pair = enzymePair.second;
            Ligand * l1 = ligand_pair.first;
            Ligand * l2 = ligand_pair.second;
            unsigned int x = e1->id;
            unsigned int y = l1->id;
            unsigned int z = l2->id;

            unsigned int index = x*id_to_ligand.size()*id_to_ligand.size() + y*id_to_ligand.size() + z;
            if(y == z)
                chem_probability_vector[index] = score/100;
            else
                chem_probability_vector[index] = score;
        }
        normalize(chem_probability_vector,id_to_enzyme.size(),id_to_ligand.size(),id_to_ligand.size(),true);


        for (StringToEnzyme::iterator i = id_to_enzyme.begin(); i != id_to_enzyme.end(); ++i){
            Enzyme * p = i->second;
            enzyme_ids.push_back(p->id);
        }
        for (StringToLigand::iterator i = id_to_ligand.begin(); i != id_to_ligand.end(); ++i){
            Ligand * p = i->second;
            ligand_ids.push_back(p->id);
        }
        IMP_NEW(IMP::domino3::StatesTable, st, (m));
        IMP_NEW(IMP::domino3::IndexStates, enzyme_states, (m, enzyme_ids));
        IMP_NEW(IMP::domino3::IndexStates, ligand_states, (m, ligand_ids));

        // three enzymes
        for (int i = 0; i < 3; i++){
            IMP::kernel::Particle * p = new IMP::kernel::Particle(m);
            IMP_NEW(IMP::domino3::Marginals, marginals, (m, p->get_index(), id_to_enzyme.size()));
            marginals->set_uniform();
            st->add(p->get_index(), enzyme_states, marginals);
        }
        // four ligands first knowen
//        IMP::kernel::Particle * p1 = new IMP::kernel::Particle(m);
//        IMP_NEW(IMP::domino3::Marginals, m1, (m, p1->get_index(), id_to_ligand.size()));
//        boost::scoped_array<double> ligand1_probability_vector(new double[id_to_ligand.size()]);
//        std::fill (ligand1_probability_vector.get(),ligand1_probability_vector.get()+id_to_ligand.size(),log(0));
//        ligand1_probability_vector[0]=log(1.0);
//        m1->set_init_vector(ligand1_probability_vector);
//        st->add(p1->get_index(), ligand_states, m1);

        for (int i = 0; i < 4; i++){
            IMP::kernel::Particle * p = new IMP::kernel::Particle(m);
            IMP_NEW(IMP::domino3::Marginals, marginals, (m, p->get_index(), id_to_ligand.size()));
            marginals->set_uniform();
            st->add(p->get_index(), ligand_states, marginals);
        }
        st->print_marginal();
        //      for (StringToEnzyme::iterator i = id_to_enzyme.begin(); i != id_to_enzyme.end(); ++i){
        //          Enzyme * p = i->second;
        //
        //      }
        
        IMP::domino3::Factors factors;
        create_graph(factors,st,m,sea_probability_vector,dock_probability_vector,chem_probability_vector);
        
        //
        IMP_NEW(IMP::domino3::Updater, ud, (factors, "updater"));
        //
        
        //
        for (unsigned int i = 0; i < iterations; ++i) {
            ud->update(100);
        }
        IMP::domino3::update_state_table(factors,st);
        IMP::domino3::print_graph(factors);
        //
        st->print_marginal();
    }
}



int main(int argc, char **argv) {
    IMP::base::setup_from_argv(argc, argv, "Experiment with loopy domino");
    IMP_NEW(IMP::kernel::Model, m, ());
    
    std::vector<IMP::kernel::Particle *> enzyms;
    
    
    run_it(m);
    
    return 0;
}
