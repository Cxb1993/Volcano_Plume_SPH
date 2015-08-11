/*******************************************************************
 * Copyright (C) 2003 University at Buffalo
 *
 * This software can be redistributed free of charge.  See COPYING
 * file in the top distribution directory for more details.
 *
 * This software is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *
 * Author: 
 * Description: most of this function is taken from move_data of
 *              TITAN, with functionality for migrating particles
 *              along with elements (called Buckets in SPH code)
 *
 *******************************************************************
 * $Id: move_data.C,v 1.1.1.1 2003/08/13 19:26:11 sorokine Exp $ 
 */

#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif

#include <cassert>
#include <vector>
using namespace std;

#include <constant.h>
#include <hashtab.h>
#include <thashtab.h>
#include <bucket.h>
#include <particle.h>
#include <bnd_image.h>

#include "exvar.h"
#include "multiproc.h"
#include "pack_data.h"
#include "repartition_BSFC.h"

#ifdef DEBUG
#  include <debug_header.h>
#endif

void
move_data (int nump, int myid, int * my_comm, 
           THashTable * P_table, HashTable * BG_mesh)
{
    // if number of procs == 1 , don't waste time here
    if (nump < 2)
        return;

#ifdef DEBUG
   bool do_search = false;
   unsigned keycheck[TKEYLENGTH] = {419430820, 302954358, 0}; //key of its neighbor which is missing
   unsigned keytemp[TKEYLENGTH] ;
#endif


    int i, j, ierr;
    const int NEW = 1;
    const int OLD = -1;

#ifdef DEBUG
    char filename[20];
    sprintf (filename, "move_data%03d.log\0", myid);
    FILE * fp = fopen (filename, "a+");
    for (i = 0; i < nump; i++)
        fprintf (fp, "%8d ", my_comm[i]);
    fprintf (fp, "\n");
#endif

    // send_info array
    int *check_proc = new int[nump];
    int *send_info = new int[2 * nump];
    int *recv_info = new int[2 * nump];

    for (i = 0; i < 2 * nump; i++)
    {
        send_info[i] = 0;
        recv_info[i] = 0;
    }

    MPI_Request * reqinfo = new MPI_Request [2 * nump];
    int tag1 = 5225;
    // post recveives for size info
    for (i = 0; i < nump; i++)
        if ( my_comm[i] )
            ierr = MPI_Irecv ((recv_info + 2*i), 2, MPI_INT, i, tag1,
                              MPI_COMM_WORLD, (reqinfo + i));

    /* count how many buckets we should send and receive from other procs */
    HTIterator *itr = new HTIterator (BG_mesh);
    Bucket *buck;

    while ((buck = (Bucket *) itr->next ()))
        if ((! buck->is_guest ()) &&
            (buck->is_active ()))
        {
            const int *neigh_proc = buck->get_neigh_proc ();
            vector < TKey > plist = buck->get_plist ();
            for (i = 0; i < nump; i++)
                check_proc[i] = 0;

            // find out number of buckets and particles to send-recv
            for (i = 0; i < NEIGH_SIZE; i++)
                if ((neigh_proc[i] > -1) &&    //make sure process id is valid
                    (neigh_proc[i] != myid) && //do not talk to myself
                    (check_proc[neigh_proc[i]] == 0)) //make sure that this process is not checked yet.
                {
                    check_proc[neigh_proc[i]] = 1;  //turn the value corresponding to neigh_proc[i] to 1;
                    send_info[2 * neigh_proc[i]]++; //number of buckets
                    if ( plist.size () > 0 )
                        send_info[2 * neigh_proc[i] + 1] += plist.size (); //number of particles
                }
        }//end of go through all buckets
    send_info[2 * myid] = 0;      // don't need to send info to myself
    send_info[2 * myid + 1] = 0;  // info for # of particles

#ifdef DEBUG
    for (i = 0; i < nump; i++)
        fprintf (fp, "%8d ", send_info[2*i]);
    fprintf (fp, "\n");
    fclose (fp);
#endif
    

    /* send out size information */
    for (i = 0; i < nump; i++)
        if ( my_comm[i] )
            ierr = MPI_Isend ((send_info + 2*i), 2, MPI_INT, i, tag1,
                              MPI_COMM_WORLD, (reqinfo + nump + i));

    /* Mark particles in guest-buckets "old" (if present) */
    itr->reset ();
    while ((buck = (Bucket *) itr->next ()))
        if (buck->is_guest ())
        {
            vector < TKey > plist = buck->get_plist ();
            vector < TKey >::iterator p_itr;
            for (p_itr = plist.begin (); p_itr != plist.end (); p_itr++)
            {
                Particle *p_old = (Particle *) P_table->lookup (*p_itr);
                p_old->put_new_old (OLD);
            }
        }

    // Wait to receive size info from others
    MPI_Status status1;
    for (i = 0; i < nump; i++)
        if ( my_comm[i] )
            MPI_Wait ((reqinfo + i), & status1);

    /* post receives */
    // size of data to be received
    int recv_count[2] = { 0, 0 };
    for (i = 0; i < nump; i++)
    {
        recv_count[0] += recv_info[2 * i];    //number of buckets
        recv_count[1] += recv_info[2 * i + 1];//number of particles
    }

    // allocate space for data to be received
    MPI_Request *recv_req = new MPI_Request[2 * nump];
    BucketPack *recv_buckets = new BucketPack[recv_count[0]];
    ParticlePack *recv_particles = new ParticlePack[recv_count[1]];
    int counter_recv[2] = { 0, 0 };
    int buck_tag = 22674, part_tag = 32491;       // random tags

    for (i = 0; i < nump; i++)
    {
        if (recv_info[2 * i] != 0)
        {
            j = MPI_Irecv ((recv_buckets + counter_recv[0]), recv_info[2 * i],
                   BUCKET_TYPE, i, buck_tag, MPI_COMM_WORLD,
                   (recv_req + 2 * i));
            counter_recv[0] += recv_info[2 * i];
        }
        if (recv_info[2 * i + 1] != 0)
        {
            j = MPI_Irecv ((recv_particles + counter_recv[1]), 
                            recv_info[2 * i + 1], PARTICLE_TYPE, i, 
                            part_tag, MPI_COMM_WORLD, (recv_req + 2 * i + 1));
            counter_recv[1] += recv_info[2 * i + 1];
        }
    }    /* done with receives */

    /* put (GHOST) elements to be moved in the proper arrays */
    // size of data to be sent
    int send_count[2] = { 0, 0 };
    for (i = 0; i < nump; i++)
    {
        send_count[0] += send_info[2 * i];
        send_count[1] += send_info[2 * i + 1];
    }
    // Allocate space for send data
    BucketPack *send_buckets = new BucketPack[send_count[0]];
    ParticlePack *send_particles = new ParticlePack[send_count[1]];
    int *counter_send_proc = new int[2 * nump];

    counter_send_proc[0] = 0;
    counter_send_proc[1] = 0;
    for (i = 1; i < nump; i++)
    {
        counter_send_proc[2 * i] =
        counter_send_proc[2 * (i - 1)] + send_info[2 * (i - 1)];
        counter_send_proc[2 * i + 1] =
            counter_send_proc[2 * (i - 1) + 1] + send_info[2 * (i - 1) + 1];
    }

    /* Pack buckets and particles that are being sent over to other procs */
    itr->reset ();
    while ((buck = (Bucket *) itr->next ()))
        if ((! buck->is_guest ()) &&
            ( buck->is_active ()))
        {
            const int *neigh_proc = buck->get_neigh_proc ();
            vector < TKey > plist = buck->get_plist ();
            vector < TKey >::iterator ip;

            // set check proc to zero
            for (i = 0; i < nump; i++)
                check_proc[i] = 0;

            // pack data to be sent over to neighs
            for (i = 0; i < NEIGH_SIZE; i++)
                if ((neigh_proc[i] > -1) &&
                    (neigh_proc[i] != myid) && 
                    (check_proc[neigh_proc[i]] == 0))//to make sure that data will not be send repeatedly
                {
                    check_proc[neigh_proc[i]] = 1;
                    pack_bucket ((send_buckets + 
                                 counter_send_proc[2 * neigh_proc[i]]),
                                      buck, myid);
                    counter_send_proc[2 * neigh_proc[i]]++;
                    if ( plist.size () > 0 )
                    {
                        for (ip = plist.begin (); ip != plist.end (); ip++)
                        {
                            Particle *psend = (Particle *) P_table->lookup (*ip);
#ifdef DEBUG
                            if ( ! psend )
                            {
                                printf ("attach debugger here\n");
                            }
#endif

#ifdef DEBUG
		                    if (do_search)
		                    {
		                      for (j = 0; j < TKEYLENGTH; j++)
			                      keytemp[j] = psend->getKey().key[j];

		                      if (find_particle (keytemp, keycheck))
			                      cout << "The particle found, will be send out!" << endl;
		                    }
#endif
                            assert (psend);
                            pack_particles (psend,
                                            (send_particles +
                                            counter_send_proc[2 * neigh_proc[i] + 1]));
                            counter_send_proc[2 * neigh_proc[i] + 1]++;
                        }
                    }
                }
        }

    //first nump is for buckets, 2nd nump is for particles
    MPI_Request *send_req = new MPI_Request[2 * nump];
    int counter[2] = { 0, 0 };
    for (i = 0; i < nump; i++)
    {
        if (send_info[2 * i] != 0)
        {
            j = MPI_Isend ((send_buckets + counter[0]), send_info[2 * i], 
                           BUCKET_TYPE, i, buck_tag, MPI_COMM_WORLD, 
                           (send_req + 2 * i));
            counter[0] += send_info[2 * i];
        }
        if (send_info[2 * i + 1] != 0)
        {
            j = MPI_Isend ((send_particles + counter[1]), send_info[2 * i + 1],
                           PARTICLE_TYPE, i, part_tag, MPI_COMM_WORLD,
                           (send_req + 2 * i + 1));
            counter[1] += send_info[2 * i + 1];
        }
    }

    int add_counter = 0, update_counter = 0;
    int count = 0;
    MPI_Status status;
    // unpack particles after receiving them
    for (i = 0; i < nump; i++)
        if (recv_info[2 * i + 1] != 0)
        {
            j = MPI_Wait ((recv_req + 2 * i + 1), &status);
            for (j = 0; j < recv_info[2 * i + 1]; j++)
            {
                Particle *pcurr = ( Particle *) 
                    P_table->lookup ((recv_particles + count)->key);
                // if the particle doesn't exist on this proc, create a new one
                if (!pcurr)
                {
                    Particle *newp = new Particle ();
                    unpack_particle ((recv_particles + count), newp);
                    newp->put_guest_flag (true);
                    newp->put_new_old (NEW);
                    P_table->add (newp->getKey (), newp);
                }
                // if it alread exists, copy the variables
                else
                {
                    unpack_particle ((recv_particles + count), pcurr);
                    pcurr->put_guest_flag (true);
                    pcurr->put_new_old (NEW);
                }
                count++;
            }
        }
    // unpack buckets 
    // All the particle-packets from proc (i) should arrive before we start
    // unpacking buckets-packets from proc (i), otherwise there will an extra 
    // overhead of deleting and re-creating particles that existed already
    count = 0;
    for (i = 0; i < nump; i++)
        if (recv_info[2 * i] != 0)
        {
            ierr = MPI_Wait ((recv_req + 2 * i), &status);
            for (j = 0; j < recv_info[2 * i]; j++)
            {
                Bucket *bcurr =
                    (Bucket *) (BG_mesh->lookup ((recv_buckets + count)->key));
                if (!bcurr)
                {
                    // this bucket doesn't exist on this proc
                    Bucket *new_buck = new Bucket ();
                    unpack_bucket ((recv_buckets + count), new_buck,
                                   (recv_buckets + count)->myprocess);
                    BG_mesh->add (new_buck->getKey (), new_buck);
                    new_buck->put_guest_flag (true);
                    add_counter++;
                }
                else
                {
                    // bucket is present on this proc 
                    // delete old particles 
                    vector < TKey > plist = bcurr->get_plist ();
                    vector < TKey >::iterator ip;
                    for (ip = plist.begin (); ip != plist.end (); ip++)
                    {
                        Particle *p_old = (Particle *) P_table->lookup (*ip);

                        if (p_old->get_new_old () == OLD)
                        {
                            P_table->remove (*ip);
                            delete p_old;
                        }
                    }
                    unpack_bucket ((recv_buckets + count), bcurr,
                                   (recv_buckets + count)->myprocess);
                    update_counter++;
                }
                count++;
            }
        }

    // make sure all the sends are completed
    // first check info sends
    for (i = 0; i < nump; i++)
        if ( my_comm[i] )
            ierr = MPI_Wait ((reqinfo + nump + i), & status1);

    // now check data sends
    for (i = 0; i < 2 * nump; i++)
        if (send_info[i] != 0)
            ierr = MPI_Wait ((send_req + i), &status);

    // clean up
    delete [] reqinfo;
    delete [] check_proc;
    delete [] recv_buckets;
    delete [] recv_particles;
    delete [] recv_info;
    delete [] recv_req;

    delete [] send_req;
    delete [] send_buckets;
    delete [] send_particles;
    delete [] send_info;
    delete [] counter_send_proc;
    delete itr;

    return;
}

/* delete the ghost elements that were put in the element hashtable */
void
delete_guest_buckets (HashTable * BG_mesh, THashTable * P_table)
{
    int delete_counter = 0;
    HTIterator *itr = new HTIterator (BG_mesh);
    Bucket *buck;

    while ((buck = (Bucket *) itr->next ()))
        if (buck->is_guest ())
        {
            vector < TKey > plist = buck->get_plist ();
            vector < TKey >::iterator ip;
            // delete particles in the bucket
            for (ip = plist.begin (); ip != plist.end (); ip++)
            {
                Particle *pdel = (Particle *) P_table->lookup (*ip);
                P_table->remove (*ip);
                delete pdel;
            }

            // now remove the bucket
            BG_mesh->remove (buck->getKey ());
            delete buck;

            delete_counter++;
        }

    delete itr;
    return;
}
