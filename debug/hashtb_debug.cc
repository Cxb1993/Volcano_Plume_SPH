/*
 * hashtb_debug.cc
 *
 *  Created on: Jul 1, 2015
 *      Author: zhixuanc
 */
#include <iostream>
#include <vector>
#include <list>
#include <cmath>
#include <cassert>
#include <algorithm>

#include <hdf5.h>
#include "hdf5calls.h"

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
#include <bucket.h>
#include <hilbert.h>
#include <particle.h>

#include "sph_header.h"
#include "particler.h"
#include "constant.h"
#include "parameters.h"

/*
 * function for new hashtable debug
 */
void hashtb_debug(MatProps* matprops, int myid, TimeProps* timeprops)
{
    int i, j, k;
    int ii;
	double mindom[DIMENSION], maxdom[DIMENSION];
	double hvars[6];
	char filename[14];

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

	double range_x[2];
	double range_y[2];
	double range_z[2];
	double crd_p[DIMENSION];
	unsigned key[TKEYLENGTH];

	unsigned tempkey[TKEYLENGTH];
	bool erpt;
    unsigned tkeylen = TKEYLENGTH;

    double dist;
    double vel;
    double umax = 2*Vv0_P;
    double rvsq=rv_P*rv_P;
    unsigned add_step = 0;

	THashTable * P_temp = new THashTable(ERUPT_TABLE_SIZE, 2017, mindom, maxdom);; //particle hash table

	//determine the rough range of eruption duck
	range_x[0] = -rv_P;
	range_x[1] = rv_P;
	range_y[0] = -rv_P;
	range_y[1] = rv_P;
    range_z[1] = 0.; // not exact, should use ground height
    range_z[0] = range_z[1]-(matprops->smoothing_length)*1.5*PARTICLE_DENSITY;


	//determine the mass and smooth length of the particle
	double mss, sml;
    double des = rhov_P;
	int np = num_erupt_particles; // total number of erupt ghost particle
	Compute_mass (np, range_x, range_y, range_z, des, &mss, &sml);

	timeprops->mass = mss;
	timeprops->sml = sml;

    double sml2 = 0.5*sml;
    double t_each = sml/Vv0_P;
    double bot = range_z[0];
    timeprops->update_teach(t_each);
    timeprops->update_bot(bot);

	//create new particles and put them into a temporary hash table: coordinate -> key -> new particle

    /*
     * debug constructor and add, bucket_search
     */
    double pi = PI;
    double r = rv_P;
    timeprops->cof = (pi*r*r)/(sml*sml);
    double normc[DIMENSION];
    int num = 0;
//    crd_p[2]=range_z[1]-sml2;
//    while ( crd_p[2] >= range_z[0] )
//    {
//    	crd_p[0]=Pos_v_P[0]+sml2;
//    	while ( crd_p[0] <= range_x[1] )
//    	{
//    		//y+ direction
//    		crd_p[1]=Pos_v_P[1]+sml2;
//    		while ( crd_p[1] <= range_y[1])
//    		{
//    			dist = 0;
//    			for (j=0; j<2; j++)
//    				dist += (crd_p[j]*crd_p[j]);
//    			if (dist <= rvsq) // if particle is in the duct
//    			{
//    			    for (int ii = 0; ii < DIMENSION; ii++)
//    				    normc[ii] = (crd_p[ii] - mindom[ii]) /(maxdom[ii] - mindom[ii]);
//
//    			    THSFC3d (normc, add_step, &tkeylen, key);
//    			    vel = parabolic_vel (rv_P, dist, umax);
////    			    Particle * pnew = new Particle(key, crd_p, mss, sml, myid,
////    			    		     vel, ev0_P,   rhov_P,   pv0_P,   gamma_v_P,
////    			    		    ng0_P,   Cvs_P,   Cvg_P,   Cva_P,   Rg_P,   Ra_P);
//
//    			    int hh = num;
//    			    int * pnew = & hh;
//    			    num++;
//    			    // add to hash-table
//    			    P_temp->add(key, pnew);
//    			}//end of if
//    			crd_p[1] += sml;
//    		}//end of while y+
//    		crd_p[0] += sml;
//    	}//end of while x+
//    	crd_p[2] -= sml;
//    }//end of while z

    int * pnew;
  /*test add*/
    //insert in front
    key = {92851969, 1949954668, 0};
    int hh2 = 102;
    pnew = & hh2;
    P_temp->add(key, pnew);

    key = {92851969, 1949954669, 0};
    int hh3 = 103;
    pnew = & hh3;
    P_temp->add(key, pnew);

//    key = {92851969, 1949954670, 0};
//    int hh4 = 104;
//    pnew = & hh4;
//    P_temp->add(key, pnew);

    key = {92851969, 1949954671, 0};
    int hh5 = 105;
    pnew = & hh5;
    P_temp->add(key, pnew);

    key = {92851969, 1949954672, 0};
    int hh6 = 106;
    pnew = & hh6;
    P_temp->add(key, pnew);

    key = {92851969, 1949954670, 1};
    int hh7 = 107;
    pnew = & hh7;
    P_temp->add(key, pnew);

    key = {92851969, 1949954670, 2};
    int hh8 = 108;
    pnew = & hh8;
    P_temp->add(key, pnew);

    //in total eight

    //The ninth
    key = {92851969, 1949954673, 0};
    int hh9 = 109;
    pnew = & hh9;
    P_temp->add(key, pnew);

    //The tenth
    key = {92851969, 1949954674, 0};
    int hh10 = 110;
    pnew = & hh10;
    P_temp->add(key, pnew);

 /*test destruction*/


 /*test remove*/
    //The tenth
    key = {92851969, 1949954674, 0};
    P_temp->remove(key);

    //The first
    key = {92851969, 1949954668, 0};
    P_temp->remove(key);

    //remove a no-exist one
    key = {92851969, 1949954665, 0};
    P_temp->remove(key);

//    key = {92851969, 1949954669, 0};
//    P_temp->remove(key);
//
//    key = {92851969, 1949954670, 0};
//    P_temp->remove(key);

//    key = {92851969, 1949954671, 0};
//    P_temp->remove(key);

    key = {92851969, 1949954672, 0};
    P_temp->remove(key);

    key = {92851969, 1949954670, 1};
    P_temp->remove(key);

//    key = {92851969, 1949954670, 2};
//    P_temp->remove(key);

 //add and remove and then add
    key = {89851969, 1949954128, 0};
    int mm1 = 201;
    pnew = & mm1;
    P_temp->add(key, pnew);

    key = {89851969, 1949954128, 0};
    P_temp->remove(key);

    key = {89851969, 1949954129, 0};
    int mm2 = 202;
    pnew = & mm2;
    P_temp->add(key, pnew);

    key = {89851969, 1949954130, 0};
    int mm3 = 203;
    pnew = & mm3;
    P_temp->add(key, pnew);


 /*test iterator*/
 THTIterator * itr2 = new THTIterator (P_temp);

 itr2->reset();
 int * Curr_part;
 while ((Curr_part = (int *) itr2->next() ))
 {
	 cout << "The value is: " <<  *Curr_part << endl;
 }


}
