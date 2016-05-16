/*
 * density_massfrac.cc
 *
 *  Created on: Feb 24, 2015
 *      Author: zhixuanc
 */

// Function for density and mass fraction updating

#include <vector>
#include <cassert>
#include <iostream>
using namespace std;

#include <thashtab.h>

#include "particler.h"
#include "constant.h"
#include "sph_header.h"

#ifdef DEBUG
#  include <debug_header.h>
#endif

void
smooth_density(THashTable * P_table)
{
  int i, k;
  double xi[DIMENSION], ds[DIMENSION], s[DIMENSION], hi;
  double wght, density=0.0;
  TKey tmpkey;
  int phs_i;
  double tmprho[PHASE_NUM]={0.0}, phaserho[PHASE_NUM]={0.0}, mssfrc;
  double wnorm[PHASE_NUM];
  // create a Hash Table iterator instance
  THTIterator * itr = new THTIterator (P_table);
  Particle * pi = NULL;

#ifdef DEBUG
   bool check_den = false;
   bool do_search = false;
   bool check_mssfrac = false;
   unsigned keycheck[TKEYLENGTH] =  {69674717, 3383906357, 0};
   unsigned keytemp[TKEYLENGTH] ;
#endif

  // iterate over the table
  while ((pi = (Particle *) itr->next ()))
  {
    if (pi->need_neigh())
    {

#ifdef DEBUG
		if (do_search)
		{
		    for (i = 0; i < TKEYLENGTH; i++)
			    keytemp[i] = pi->getKey ().key[i];

		    if (find_particle (keytemp, keycheck))
			    cout << "The particle found!" << endl;
		}
#endif

      for (i = 0; i < DIMENSION; i++)
        xi[i] = (*(pi->get_coords () + i));

#if DENSITY_UPDATE_SML==0
      hi = pi->get_original_smlen ();
#elif DENSITY_UPDATE_SML==1
      hi = pi->get_smlen ();
#endif

      double supp = 3.0 * hi;
//      TKey ki = pi->getKey ();

      for (i = 0; i < PHASE_NUM; i++)
    	  tmprho[i] = 0.;

      for (i = 0; i < PHASE_NUM; i++)
    	  phaserho[i] = 0.;

      for (i = 0; i < PHASE_NUM; i++)
    	  wnorm [i] = 0.;

      vector <TKey> neighs = pi->get_neighs ();
      vector <TKey> :: iterator p_itr;

      //go through all neighbors
      for (p_itr = neighs.begin (); p_itr != neighs.end (); p_itr++)
      {
          Particle *pj = (Particle *) P_table->lookup (*p_itr);
          assert (pj);

          // guests are included ghosts are not,
          if (pj->contr_dens()) //In my old code, I did not have this so I will get very large density at the very beginning of the eruption.
          //update density for each phase
          //Probably, a switch case is better than for loop in terms of efficiency!
			  for (phs_i = 1; phs_i<= PHASE_NUM; phs_i++)
			  {
				  //for each phase
				  if (pj->which_phase()==phs_i)
				  {

					// the boundary deficiency and interface deficiency is handled by wnorm!
					for (i = 0; i < DIMENSION; i++)
						ds[i] = xi[i] - *(pj->get_coords() + i);

					if (in_support(ds, supp))
					{
						for (k = 0; k < DIMENSION; k++)
						    s[k] = ds[k] / hi;
						wght = weight(s, hi);
						tmprho[phs_i-1] += wght * (pj->get_mass());
#if DENSITY_UPDATE_SPH == 1
						//View multiple phase SPH as one set of discretized point
						wnorm[phs_i-1] += wght * (pj->get_mass()) / pj->get_density();
#elif DENSITY_UPDATE_SPH == 2
						//View multiple phase SPH as two independent sets of discretized point
						double phase_des=*(pj->get_phase_density ()+ phs_i -1);
						if (phase_des>0)
						    wnorm[phs_i-1] += wght * (pj->get_mass())/phase_des;
#endif

					 }
				  }//end of if pj->which_phase()==phs_i

			  }//end of loop for phase
          //end of if particle will contribute to the density.

      }//end of go through all neighbors

/*
 * There are three ways to normalize SPH kernel summation
 * 1)based on the concept that each phase in SPH is essentially independent set of discretized points
 * 2)Another is based on the concept that different types of particles for different phases are essentially the equivalent as discretized points
 * It seems that the former one will smear out density
 * While the second one have some troubles to compute density at the interface.
 */
      for (phs_i = 2; phs_i<= PHASE_NUM; phs_i++) //add wnorm up, so start from second phase.
      {
    	  wnorm[0] += wnorm [phs_i-1];
      }

      density=0.0;
#if DENSITY_UPDATE_SPH == 1
      for (phs_i = 1; phs_i<= PHASE_NUM; phs_i++)
      {
           assert (wnorm[0] > 0);
    	   phaserho[phs_i-1]= tmprho[phs_i-1] / wnorm[0];
    	   density +=phaserho[phs_i-1];
      }
#elif DENSITY_UPDATE_SPH == 2
	 for (phs_i = 1; phs_i<= PHASE_NUM; phs_i++)
	 {
	    if (wnorm[phs_i-1] > 0)
		   phaserho[phs_i-1]= tmprho[phs_i-1] / wnorm[phs_i-1];
	    else
		   phaserho[phs_i-1]=0;

	    density +=phaserho[phs_i-1];
	 }
#endif

#ifdef DEBUG
      if (check_den)
      {
         if (!(density > 0))
        	 cout << "density is invalid" << endl;
      }
#endif

	  assert(density > 0);
      pi->put_new_density(density);

      mssfrc=phaserho[1]/density;

#ifdef DEBUG
      if (check_mssfrac)
      {
         if (!(mssfrc <=1))
        	 cout << "mass fraction larger than 1" << endl;
      }
#endif
      assert(mssfrc <= 1);
      pi->put_new_mass_frac(mssfrc);

//      /*
//       * the bad thing here is: particle do not have information about its bgmesh, the better way is set adapt +=1 only if the bucket where the particle belong to originally does not have involved particles.
//       * Well why not add a new member Key mybucket in particle ---> this is another story, keep track bgmesh for each particle requires additional effort
//       *
//       * Any way, what is the solution? The solution is give up determine the value of adapt(which will determine whether updating of background mesh is necessary!) within density updating function.
//       * update_adpt_domain will be called at more frequent to determine whether updating of background mesh is necessary.
//       */
//
//      if (mssfrc >= MSFRC_THRESH && !pi->is_involved())
//          pi->set_involved_flag(INVOLVED);


    }//end of if --> if the density of that particle need to be updated based summation
  }//end of loop -->go through all particle


  //In the old code, mass fraction is updated while computing, put it separately after density updating done
   THTIterator *it2 = new THTIterator(P_table);
   while ((pi = (Particle *) it2->next()))
     if (pi->need_neigh())
     {
       pi->update_mass_frac();
       mssfrc = pi->get_mass_frac();
       if ( mssfrc>= MSFRC_THRESH && !pi->is_involved())
          pi->set_involved_flag(INVOLVED);

       pi->update_density();
       pi->update_second_var(ng0_P, Cvs_P, Cvg_P, Cva_P, Rg_P, Ra_P, rhoa0_P); //The secondary varible was not updated after updating density and mass fraction.
     }

  // clean up stuff
  delete itr;
  delete it2;

  return;
}

