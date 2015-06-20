/*
 * debug_header.h
 *
 *  Created on: Jun 1, 2015
 *      Author: zhixuanc
 */
#ifndef DEBUG_HEADER_H
#define DEBUG_HEADER_H

#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <vector>
#include <list>

#ifdef HAVE_MPI_H
#  include <mpi.h>
#endif

#include <cmath>

using namespace std;

#include <hashtab.h>
#include <thashtab.h>
#include <bgmesh.h>
#include <bnd_image.h>
#include <properties.h>
#include <buckhead.h>
#include <hilbert.h>
#include <bucket.h>
#include <particle.h>
#include <multiproc.h>
#include <properties.h>

#include <constant.h>
#include <IndMap.h>
#include <parameters.h>

void
particle_deb (
		//myid
		int
		);
//
///*find particle by the key
// * keyin is input particle key
// * keycheck is the key that is given and all keyin will be compared with keycheck
// * pi is the pointer points to the particle corresponding to keyin
// */
//Particle* find_particle (
//		unsigned*, //keyin
//		unsigned *, //keycheck
//		Particle*  //pi
//		);

/*find particle by the key
 * keyin is input particle key
 * keycheck is the key that is given and all keyin will be compared with keycheck
  */
bool find_particle (
		unsigned*, // keyin
		unsigned* //keycheck
		);

/*find particle by the pos
 *
 */
bool find_particle (
		double*, // kin
		double* //check
		);

/*find particle by the key
 * keyin is input key
 * keycheck is the key that is given and all keyin will be compared with keycheck
  */
bool find_bucket (
		unsigned*, // keyin
		unsigned* //keycheck
		);

//function to find particle with given key, will be useful in debugging.
void check_particle_bykey (
		THashTable *
        );

//function that call output sub_functions
void
write_particles_debug(
		     int , //myid
		     int , //numprocs
             THashTable * , //P_table
             HashTable * ,  //BG_mesh
             vector <BucketHead> &, //partition_table
             TimeProps * , //timeprops
             int,  //format
             char *//prefix
             );

//output certain type of particle
void
write_h5part_bctp(
		int , //myid
		int , //numproc
		THashTable * , //P_table
		TimeProps * ,  //timepros
		int ,//bctp
		char *//prefix
		);

//go through all particle and check their type!
void check_particle_all_type (
		THashTable * //P_table
		);

//check particles in certain region
void check_particle_bypos (
		THashTable * //P_table
		);

#endif /* DEBUG_HEADER_H */
