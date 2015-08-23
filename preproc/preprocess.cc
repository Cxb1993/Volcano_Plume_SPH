/*
 * preprocess.cc
 *
 *  Created on: Mar 7, 2015
 *      Author: zhixuanc
 */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <string>
#include <algorithm>
#include <cassert>

#include <sstream>
#include <iomanip>
#include <cctype>

using namespace std;

#include <unistd.h>
#include <hilbert.h>
#include <GisApi.h>

#include "buckstr.h"
#include "preprocess.h"
#include "constant.h"
#include "parameters.h"

//#ifdef DEBUG
bool find_bucket (unsigned* keyin, unsigned* keycheck)
{
	int i;
	for (i=0; i<KEYLENGTH; i++)
		if (keyin[i] != keycheck[i])
			return false;

    return true;
}
//#endif

bool to_bool(std::string str) {
    std::transform(str.begin(), str.end(), str.begin(), ::tolower);
    std::istringstream is(str);
    bool b;
    is >> std::boolalpha >> b;
    return b;
}

void usage()
{
  cerr <<"Real Ground, Usage :"<< "./preprocess <nproc> <smooth-length> <gis-flag> <gis database> <gis location>"
                  << " <gis mapset> <gis map>" << endl;
  exit(0);
}

void usage_flatg()
{
  cerr <<"Ideal Flat Ground, Usage :"<< "./preprocess <nproc> <smooth-length> <gis-flag>" << endl;
  exit(0);
}
void get_plane_coef(double x[], double y[], double z[], double poly[])
{
  double temp1 = z[1] - z[0];
  double temp2 = x[1] - x[0];
  double temp3 = sqrt(pow(temp1,2)+pow(temp2,2));
  poly[0] = -temp1/temp3;
  poly[1] = temp2/temp3;
  poly[2] = (temp1*x[0]-temp2*z[0])/temp3;
  return;
}

bool compare_keys ( const ColumnHead  & buck1, const ColumnHead & buck2)
{

	  if (buck1.key[0] < buck2.key[0])
	    return true;
	  else if (buck1.key[0] > buck2.key[0])
	    return false;
	  else if (buck1.key[1] < buck2.key[1])
	    return true;
	  else
	    return false;

  cerr << "ERROR: Two keys match exactly, expand keylength" << endl;
  cerr << buck1.key[0] <<", " << buck1.key[1] << " == "
          << buck2.key[0] <<", " << buck2.key[1] << endl;
  exit(1);
}

int main(int argc, char *argv[])
{

  int i, j, k, l;
  int np, kk;
  double xcrd[2],ycrd[2], zcrd[2],cntr[DIMENSION],smlen;
  unsigned keylen=KEYLENGTH;
  unsigned key[KEYLENGTH], key2[KEYLENGTH];
  double mindom[DIMENSION], maxdom[DIMENSION], polyconst[DIMENSION+1];
  double mindom_o[DIMENSION], maxdom_o[DIMENSION];
  string gis_db, gis_location, gis_mapname, gis_mapset;
  bool gis_flag = to_bool (argv[3]); //flag that used to indicate whether use real geometry or flat ground : | 1: real ground
 //                                                                                                            | 0: flat ground

#ifdef DEBUG
   bool do_search = false;
   unsigned keycheck[KEYLENGTH] = {92836425, 613566756};
//   unsigned keytemp[KEYLENGTH] ;
#endif

  // hashtable constants
  double htvars[6];

  if (!gis_flag && argc != 4)
        usage_flatg();

  np = atoi(argv[1]);
  smlen = atof(argv[2]);

  double del=PARTICLE_DENSITY*smlen;
  double resolution = del;

  // get domain limits from GIS
	  mindom_o[0]= Lx_P[0];  //original domain
	  mindom_o[1]= Ly_P[0];  //original domain
	  mindom_o[2]= Lz_P[0];  //original domain

	  maxdom_o[0]= Lx_P[1];  //original domain
	  maxdom_o[1]= Ly_P[1];  //original domain
	  maxdom_o[2]= Lz_P[1];  //original domain

  //To impose pressure atmosphere and wall bc, expand the domain
  double extend_cof = EXT_DOM_COF;
  for (i=0; i<DIMENSION; i++)
  {
	  mindom[i] = mindom_o[i] - extend_cof * del;
	  maxdom[i] = maxdom_o[i] + extend_cof * del;
  }

  /*
  if the wall boundary buckets has two ghost buckets
  The lowest bucket will not be able find the bucket which contains
  boundary information
  In addition, the reason to have thicker ghost layer is to avoid
  influence of "near boundary effect" on static pressure ghost,
  It is not necessary for wall ghost to have multiple ghost bucket.
  So, reset value of mindom[2];
  */
  mindom[2] = mindom_o[2] - 1.5*del;

  // max number of buckets along each directions
  int nx = (int) ceil((maxdom[0]-mindom[0])/(del));
  int ny = (int) ceil((maxdom[1]-mindom[1])/(del));
  int nz = (int) ceil((maxdom[2]-mindom[2])/(del));

  maxdom[0] = mindom[0] + nx*del;
  maxdom[1] = mindom[1] + ny*del;
  maxdom[2] = mindom[2] + nz*del;

  // fill hashtable variables
  htvars[0] = mindom[0];
  htvars[1] = maxdom[0];
  htvars[2] = mindom[1];
  htvars[3] = maxdom[1];
  htvars[4] = mindom[2];
  htvars[5] = maxdom[2];

  int Nbucket = nx*ny*nz;
  vector<ColumnHead> partition_table;

  // data-structure for back-ground mesh
  // is a 2-d array of linkded lists
  // using <STL List> for the purpose
  BucketStruct ***bgmesh = new BucketStruct**[nx];
  for (i=0; i<nx; i++)
  {
    bgmesh[i] = new BucketStruct* [ny];

    for (j=0; j<ny; j++)
      bgmesh[i][j] = new BucketStruct [nz];
  }

  double elev[4]; //elev only needed for onground MIXED bucket.
  //In my code, I would like to let all buckets have the information,
  //But only on ground MIXED bucket will need the information.
  for (i=0; i<4; i++)
	  elev[i]=mindom_o[2];

  int type;
  int bk_index[2*DIMENSION];
  double time = 0.;
  for (i = 0; i < nx; i++)
  {
    xcrd[0] = mindom[0] + i*del;
    xcrd[1] = mindom[0] + (i+1)*del;
    for (j=0; j<ny; j++)
    {
      ycrd[0] = mindom[1] + j*del;
      ycrd[1] = mindom[1] + (j+1)*del;

      for (k = 0; k < nz; k++)
      {
        zcrd[0] = mindom[2] + k*del;
        zcrd[1] = mindom[2] + (k+1)*del;

        for (l=0; l<4; l++)
        	bgmesh[i][j][k].elev[l] = elev[l];

        for (l=0; l<2; l++)
        {
          bgmesh[i][j][k].xcoord[l] = xcrd[l];
          bgmesh[i][j][k].ycoord[l] = ycrd[l];
          bgmesh[i][j][k].zcoord[l] = zcrd[l];
        }

        determine_bucket_type (mindom_o, maxdom_o, xcrd, ycrd, zcrd, &type, bk_index);
        bgmesh[i][j][k].buckettype = type;
        for (l=0; l<2*DIMENSION; l++)
        	bgmesh[i][j][k].bucket_index[l]=bk_index[l];

        // generate hash-key for bucket
        double normc[DIMENSION];
        cntr[0] = (xcrd[0]+xcrd[1])*0.5;
        cntr[1] = (ycrd[0]+ycrd[1])*0.5;
        cntr[2] = (zcrd[0]+zcrd[1])*0.5;
        for ( l=0; l<DIMENSION; l++)
            normc[l]=(cntr[l]-mindom[l])/(maxdom[l]-mindom[l]);

        // determine key
        HSFC3d (normc, & keylen, key);

//#ifdef DEBUG
		if (do_search)
		{
		   if (find_bucket (key, keycheck))
			   cout << "The bucket is found!" << endl;
		}
//#endif

        // put key value to BG_mesh structure
        for (l=0; l<KEYLENGTH; l++)
          bgmesh[i][j][k].key[l] = key[l];

        // find key value for partition bucket
        if ( k==0 )
        {
          HSFC2d (normc, & keylen, key2);
          ColumnHead temp_head (i, j, key2);
          partition_table.push_back(temp_head);
        }
      }
    }
  }

  // order buckets according to keys
  sort(partition_table.begin(), partition_table.end());

  int bucks_per_proc = (nx*ny/np);
  int nxny = nx * ny;
  for (i=0; i < nxny; i++ )
  {
    int myid = i/bucks_per_proc;
    if ( myid >= np )  myid = np-1;
    partition_table[i].proc = myid;
    j = partition_table[i].xind;
    k = partition_table[i].yind;
    for (l = 0; l < nz; l++)
      if ( bgmesh[j][k][l].buckettype )
        bgmesh[j][k][l].myproc = myid;
      else
        bgmesh[j][k][l].myproc = -1;
  }

  // determine neighbors
  for (i = 0; i < nx; i++)
    for (j = 0; j < ny; j++)
      for (k = 0; k < nz; k++)
      {
        int ncount=0;
        int npcount=0;
        for (int ii=i-1; ii<i+2; ii++)
          for (int jj=j-1; jj<j+2; jj++)
            for (int kk=k-1; kk<k+2; kk++)
            {
              if ( ! bgmesh[i][j][k].buckettype )//bucket type == 0: means invalid bucket.
                continue;

              if ((ii>-1)&&(ii<nx) &&
                  (jj>-1)&&(jj<ny) &&
                  (kk>-1)&&(kk<nz))
              {
                //The key length of bgmesh should be KEYLENGTH!!!!!
            	for (l=0; l<keylen; l++)
            		bgmesh[i][j][k].neighs[ncount++] = bgmesh[ii][jj][kk].key[l];

                if ((ii==i) && (jj==j) && (kk==k))
                  bgmesh[i][j][k].neigh_proc[npcount++] = -2; // bucket itself can not been its neighbor;
                else
                  bgmesh[i][j][k].neigh_proc[npcount++] = bgmesh[ii][jj][kk].myproc;
              }
              else
              {
                bgmesh[i][j][k].neighs[ncount++]=0;//set key to [0,0]
                bgmesh[i][j][k].neighs[ncount++]=0;
                bgmesh[i][j][k].neigh_proc[npcount++] = -1;
              }
            }
      }

  int kcount;
  vector <ColumnHead> :: iterator c_itr = partition_table.begin ();
  for (int iproc = 0; iproc < np; iproc++)
  {
    vector<BucketStruct> proc_bucks;
    vector <unsigned> partition_keys;
    while ( c_itr->proc == iproc )
    {
      i = c_itr->xind;
      j = c_itr->yind;
      kcount = 0;
      for (k = 0; k < nz; k++)
        if ( bgmesh[i][j][k].buckettype > 0)
        {
          if (bgmesh[i][j][k].myproc == iproc)
          {
            proc_bucks.push_back(bgmesh[i][j][k]);
            kcount++;
          }
          else
          {
            cerr << "Error: proc-ids don't match" << endl;
            exit (1);
          }
        }
        // copy partition table keys
        for (k = 0; k < KEYLENGTH; k++)
          partition_keys.push_back (c_itr->key[k]);

        // search for the first-bucket
        for (l = 0; l < nz; l++)
          if (bgmesh[i][j][l].buckettype == 1)
            break;

        // copy key of the first bucket in the column
        for (k = 0; k < KEYLENGTH; k++)
          partition_keys.push_back (bgmesh[i][j][l].key[k]);

        // advance the iterator
        c_itr++;
        if ( c_itr == partition_table.end () )
          break;
    }
    createfunky (iproc, 6, htvars, proc_bucks, partition_keys);
    proc_bucks.clear();
  }

  // Create write initial data to HDF5 file
  cout << "Total "<< Nbucket <<", "<< bucks_per_proc * nz <<" buckets per proc" << endl;
  cout << "it's ready to run"<<endl;


  // clean up
  for (i=0; i<nx; i++)
  {
    for (j=0; j<ny; j++)
      delete [] bgmesh[i][j];

    delete [] bgmesh[i];
  }
  delete [] bgmesh;
  return 0;
}

void determine_the_key(double norm_coord[], unsigned nkey, unsigned key[],
                             unsigned ma[], unsigned mi[])
{
  // call Hilbert's Space Filling Curve
  HSFC3d (norm_coord, &nkey, key);

    int i;

  if(key[0]>ma[0] || (key[0]==ma[0] && key[1]>ma[1])) {ma[0]=key[0]; ma[1]=key[1];}
  if(key[0]<mi[0] || (key[0]==mi[0] && key[1]<mi[1])) {mi[0]=key[0]; mi[1]=key[1];}

  return;
}

//function that used to determine the value of face by the face's index;
int determine_face_type (double crd, double max, double min)
{
	int flag;

	if (crd<min)
		flag = -1;
	else if (crd > max)
		flag =1;
	else
		flag =0;

	return flag;
}

//function that used to determine the type of bucket
void determine_bucket_type (double *mincrd, double *maxcrd, double *xcrd, double *ycrd, double *zcrd, int* type, int* bt)
{
	int flag[DIMENSION];
	int sum = 0;
	int k;

    bt[0]=determine_face_type(xcrd[0],maxcrd[0],mincrd[0]);
    bt[1]=determine_face_type(xcrd[1],maxcrd[0],mincrd[0]);
    bt[2]=determine_face_type(ycrd[0],maxcrd[1],mincrd[1]);
    bt[3]=determine_face_type(ycrd[1],maxcrd[1],mincrd[1]);
    bt[4]=determine_face_type(zcrd[0],maxcrd[2],mincrd[2]);
    bt[5]=determine_face_type(zcrd[1],maxcrd[2],mincrd[2]);

    for (k=0; k<DIMENSION; k++)
    	flag[k]=abs(bt[2*k]+bt[2*k+1]);

    for (k=0; k<DIMENSION; k++)
    	sum += flag[k];

    if ( (flag[0]==2) || (flag[1]==2) || (flag[2]==2))
    {
    	if ((bt[4] == -1) && (bt[5] == -1)) //only bt[5] == -1 is enough!
    		*type = 1; //underground
    	else
    		*type = 4; //pressure bc
    }

    else if ( (flag[0]==1) || (flag[1]==1) || (flag[2]==1) )
    	*type = 2;     //Mixed
    else if(sum == 0)
    	*type = 3;     //overground
    else
    {
    	*type = 0;     //invalid
    	cout << "Invalid input bucket_index!!!" << endl;
    }

    /*Some buckets which are "underground" but also has ground boundary information are shift to MIXED

        |
        |                                        ORIGINAL DOMAIN
        |
        |                 Left boundary
        |______________|_______!_______
        |Orig:UNDG     |       !       |
        |._._._._._._._|._._._.!_._._._|._._._.  Bottom Boundary
        |Changed to:   |       !       |
        |Mixed___ _____|_______!_______|_______
        |              |       !       |
        |              |       !       |
        |              |       !       |
        |______________|_______!_______|________________________

         The bucket at the left
         corner should still
         be UNDERGROUND!
    */

    if (bt[4]==-1 && bt[5]==0)
    	*type = 2;

    return ;
}
