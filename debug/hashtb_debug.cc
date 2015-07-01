/*
 * hashtb_debug.cc
 *
 *  Created on: Jul 1, 2015
 *      Author: zhixuanc
 */
#include <iostream>
#include <vector>
#include <list>
using namespace std;

#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif

#ifdef MULTI_PROC
#  include <mpi.h>
#  include <multiproc.h>
#endif

#ifdef DEBUG
#  include <debug_header.h>
#endif

#include <hashtab.h>
#include <thashtab.h>
#include <bgmesh.h>
#include <bnd_image.h>
#include <properties.h>
#include <buckhead.h>

#include "sph_header.h"
#include "particler.h"

void hashtb_debug(MatProps* matprops, int myid)
{
	int P_TABLE_SIZE = 10000;
	double mindom[DIMENSION], maxdom[DIMENSION];
	double hvars[6];
	THashTable *P_table;
	char filename[14];
	// Create hash-table for particles
	 // Read Hash table constants
	  // Read Hash-table related constants
	  sprintf(filename, "funky%04d.h5", myid);
	  hid_t fp = GH5_fopen_serial(filename, 'r');
	  GH5_readdata(fp, "/hashtable_constants", hvars);
	  mindom[0] = hvars[0];         // min x
	  maxdom[0] = hvars[1];         // max x
	  mindom[1] = hvars[2];         // min y
	  maxdom[1] = hvars[3];         // max y
	  mindom[2] = hvars[4];         // min z
	  maxdom[2] = hvars[5];         // max z
	  // Create two new Hash-tables
	  for (i = 0; i < DIMENSION; i++)
	  {
	    mindom[i] /= matprops->LENGTH_SCALE;
	    maxdom[i] /= matprops->LENGTH_SCALE;
	  }
	*P_table = new THashTable(P_TABLE_SIZE, 2017, mindom, maxdom);

}
