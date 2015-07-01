/*
 * thashtab.cc
 *
 *  Created on: May 5, 2015
 *      Author: zhixuanc
 */

#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif

#ifdef HAVE_MPI_H
#  include <mpi.h>
#endif

#include <cassert>
#include <climits>
using namespace std;

#include "thashtab.h"

#define THASHTABLE_EXTENDER 1000000
#define MaxBits ( sizeof(unsigned) * CHAR_BIT )
#define IScale  ((unsigned)((MaxBits <= 32) ? ~(0u) : (0xffffffff << (MaxBits - 32))))

THashTable::THashTable (int size, int prime, double minR[], double maxR[])
{
  int i;
  NBUCKETS = size;
  PRIME = prime;
  SIZE02 = NBUCKETS / 10;
  SIZE01 = NBUCKETS - SIZE02;
  umax = (double) IScale;

// allocate table-size, and initialize it
  bucket_vec = new vector<THashEntryPtr> [NBUCKETS];
  for (i = 0; i < NBUCKETS; i++)
	  bucket_vec[i].push_back(NULL);

  for (i = 0; i < DIMENSION; i++)
  {
    minDom[i] = minR[i];
    maxDom[i] = maxR[i];
  }
}

//
////overloading of constructor
//THashTable::THashTable (int size, int prime, double minR[], double maxR[], int maxn)
//{
//  int i;
//  NBUCKETS = size;
//  PRIME = prime;
//  SIZE02 = NBUCKETS / 10;
//  SIZE01 = NBUCKETS - SIZE02;
//  umax = (double) IScale;
//
//// allocate table-size, and initialize it
//  bucket_vec = new vector<THashEntryPtr> [NBUCKETS];
//  for (i = 0; i < NBUCKETS; i++)
//  	 bucket_vec[i].push_back(NULL);
//
//  for (i = 0; i < DIMENSION; i++)
//  {
//    minDom[i] = minR[i];
//    maxDom[i] = maxR[i];
//  }
//}

THashTable::~THashTable ()        //evacuate the table
{
  for (int i = 0; i < NBUCKETS; i++)
  {
    THashEntryPtr p = *(bucket_vec[i].begin());

    while (p)
    {
      THashEntryPtr p_next = p->next;
      delete p;

      p = p_next;
    }
  }
  delete[]bucket_vec;
}

//overloading search bucket -->should be faster!
THashEntryPtr THashTable::searchBucket (int entry, unsigned* keyi)
{
  int i;
  int size = bucket_vec[entry].size();

  if (size == 0)
     return NULL;

  unsigned* keyarr;

/*
 * I think the following is not necessary!
 */

//  unsigned* keymin = (*(bucket_vec[entry].begin()))->key;
//  unsigned* keymax = (*(bucket_vec[entry].end()-1))->key;
//  //if key is smaller than keymin, return null
//  if (*keyi < *keymin)
//	  return NULL;
//  else if (*keyi > *keymin)
//	  break;
//  else if ( (*keyi == *keymin) && (*(keyi+1) < *(keymin+1)) )
//	  return NULL;
//  else if ( (*keyi == *keymin) && (*(keyi+1) > *(keymin+1)) )
//	  break;
//  else if ((*keyi == *keymin) && (*(keyi+1) == *(keymin+1)) && (*(keyi+2) < *(keymin+2)) )
//	  return NULL;
//
//  //if key is greater than keymin, return null
//  if (*keyi > *keymax)
//	  return NULL;
//  else if (*keyi < *keymax)
//	  break;
//  else if ( (*keyi == *keymax) && (*(keyi+1) > *(keymax+1)) )
//	  return NULL;
//  else if ( (*keyi == *keymax) && (*(keyi+1) < *(keymax+1)) )
//	  break;
//  else if ((*keyi == *keymax) && (*(keyi+1) == *(keymax+1)) && (*(keyi+2) > *(keymax+2)) )
//	  return NULL;

  if (size < HASHTABLE_LOOKUP_LINSEARCH)
  {
      for (i = 0; i < size; i++)
      {
    	  keyarr = (*(bucket_vec[entry].begin()+i))->key;
    	  if (*keyi != *keyarr)
    		  break;
    	  else if (*(keyi+1) != *(keyarr+1))
    		  break;
    	  else if (*(keyi+2) != *(keyarr+2))
    		  break;
          else //if ((*keyi == *keyarr) && (*(keyi+1) == *(keyarr+1)) && (*(keyi+2) == *(keyarr+2)))
              return bucket_vec[entry][i];
      }
  }
  else
  {
      int i0, i1, i2;
      i0 = 0;
      i1 = size / 2;
      i2 = size - 1;
      while ((i2 - i0) > HASHTABLE_LOOKUP_LINSEARCH)
      {
    	 keyarr = (*(bucket_vec[entry].begin()+i1))->key;

         if (*keyi > *keyarr)
         {
             i0 = i1 + 1;
             i1 = (i0 + i2) / 2;
         }
         else if (*keyi < *keyarr)
         {
            i2 = i1;
            i1 = (i0 + i2) / 2;
         }
         else if ((*keyi == *keyarr) && (*(keyi+1) > *(keyarr+1)))
         {
             i0 = i1 + 1;
             i1 = (i0 + i2) / 2;
         }
         else if ((*keyi == *keyarr) && (*(keyi+1) < *(keyarr+1)))
         {
            i2 = i1;
            i1 = (i0 + i2) / 2;
         }
         else if ((*keyi == *keyarr) && (*(keyi+1) == *(keyarr+1)) && (*(keyi+2) > *(keyarr+2)))
         {
             i0 = i1 + 1;
             i1 = (i0 + i2) / 2;
         }
         else
         {
            i2 = i1;
            i1 = (i0 + i2) / 2;
         }
      }

      for (i = i0; i <= i2; i++)
      {
    	  keyarr = (*(bucket_vec[entry].begin()+i))->key;
    	  if (*keyi != *keyarr)
    		  break;
    	  else if (*(keyi+1) != *(keyarr+1))
    		  break;
    	  else if (*(keyi+2) != *(keyarr+2))
    		  break;
          else //if ((*keyi == *keyarr) && (*(keyi+1) == *(keyarr+1)) && (*(keyi+2) == *(keyarr+2)))
              return bucket_vec[entry][i];
      }
  }
  return NULL;
}

//overloading search bucket with int i as input.
THashEntryPtr THashTable::searchBucket (int entry, unsigned* keyi, int *index)
{
  int i;
  int size = bucket_vec[entry].size();

  if (size == 0)
     return NULL;

  unsigned* keyarr;

/*
 * I think the following is not necessary!
 */

//  unsigned* keymin = (*(bucket_vec[entry].begin()))->key;
//  unsigned* keymax = (*(bucket_vec[entry].end()-1))->key;
//  //if key is smaller than keymin, return null
//  if (*keyi < *keymin)
//	  return NULL;
//  else if (*keyi > *keymin)
//	  break;
//  else if ( (*keyi == *keymin) && (*(keyi+1) < *(keymin+1)) )
//	  return NULL;
//  else if ( (*keyi == *keymin) && (*(keyi+1) > *(keymin+1)) )
//	  break;
//  else if ((*keyi == *keymin) && (*(keyi+1) == *(keymin+1)) && (*(keyi+2) < *(keymin+2)) )
//	  return NULL;
//
//  //if key is greater than keymin, return null
//  if (*keyi > *keymax)
//	  return NULL;
//  else if (*keyi < *keymax)
//	  break;
//  else if ( (*keyi == *keymax) && (*(keyi+1) > *(keymax+1)) )
//	  return NULL;
//  else if ( (*keyi == *keymax) && (*(keyi+1) < *(keymax+1)) )
//	  break;
//  else if ((*keyi == *keymax) && (*(keyi+1) == *(keymax+1)) && (*(keyi+2) > *(keymax+2)) )
//	  return NULL;


  if (size < HASHTABLE_LOOKUP_LINSEARCH)
  {
      for (i = 0; i < size; i++)
      {
    	  keyarr = (*(bucket_vec[entry].begin()+i))->key;
    	  if (*keyi != *keyarr)
    		  break;
    	  else if (*(keyi+1) != *(keyarr+1))
    		  break;
    	  else if (*(keyi+2) != *(keyarr+2))
    		  break;
          else //if ((*keyi == *keyarr) && (*(keyi+1) == *(keyarr+1)) && (*(keyi+2) == *(keyarr+2)))
          {
              *index = i;
        	  return bucket_vec[entry][i];
          }
      }
  }
  else
  {
      int i0, i1, i2;
      i0 = 0;
      i1 = size / 2;
      i2 = size - 1;
      while ((i2 - i0) > HASHTABLE_LOOKUP_LINSEARCH)
      {
    	 keyarr = (*(bucket_vec[entry].begin()+i1))->key;

         if (*keyi > *keyarr)
         {
             i0 = i1 + 1;
             i1 = (i0 + i2) / 2;
         }
         else if (*keyi < *keyarr)
         {
            i2 = i1;
            i1 = (i0 + i2) / 2;
         }
         else if ((*keyi == *keyarr) && (*(keyi+1) > *(keyarr+1)))
         {
             i0 = i1 + 1;
             i1 = (i0 + i2) / 2;
         }
         else if ((*keyi == *keyarr) && (*(keyi+1) < *(keyarr+1)))
         {
            i2 = i1;
            i1 = (i0 + i2) / 2;
         }
         else if ((*keyi == *keyarr) && (*(keyi+1) == *(keyarr+1)) && (*(keyi+2) > *(keyarr+2)))
         {
             i0 = i1 + 1;
             i1 = (i0 + i2) / 2;
         }
         else
         {
            i2 = i1;
            i1 = (i0 + i2) / 2;
         }
      }

      for (i = i0; i <= i2; i++)
      {
    	  keyarr = (*(bucket_vec[entry].begin()+i))->key;
    	  if (*keyi != *keyarr)
    		  break;
    	  else if (*(keyi+1) != *(keyarr+1))
    		  break;
    	  else if (*(keyi+2) != *(keyarr+2))
    		  break;
          else //if ((*keyi == *keyarr) && (*(keyi+1) == *(keyarr+1)) && (*(keyi+2) == *(keyarr+2)))
          {
              *index = i;
        	  return bucket_vec[entry][i];
          }
      }
  }
  return NULL;
}


/*
 * overload addElement
 * if condition can be optimized! --->overlap judgment is not necessary!
 */
THashEntryPtr THashTable::addElement (int entry, unsigned keyi[])
{
  int i;
  THashEntryPtr p = new THashEntry (keyi);

  int size = bucket_vec[entry].size();

  if (size == 0) //if no element has been added at this row yet
  {
	  bucket_vec[entry].insert(bucket_vec[entry].begin(), p);
	  return p;
  }
  else
  {
	  unsigned* keyarr;
	  if (size < HASHTABLE_LOOKUP_LINSEARCH)
	  {
	      for (i = 0; i < size; i++)
	      {
	    	  keyarr = (*(bucket_vec[entry].begin()+i))->key;
	    	  if (*keyi > *keyarr)
	    		  break;
	    	  else if (*keyi < *keyarr)
	    	  {
	    		  p->next = *(bucket_vec[entry].begin()+i);
	    		  (*(bucket_vec[entry].begin()+i))->pre = p;

	    		  if (i != 0)
	    		  {
		    		  p->pre = *(bucket_vec[entry].begin()+i-1);
		    		  (*(bucket_vec[entry].begin()+i-1))->next = p;
	    		  }
	    		  bucket_vec[entry].insert(bucket_vec[entry].begin()+i,p);
	    		  return p;
	    	  }
	    	  else if ((*keyi == *keyarr) && (*(keyi+1) > *(keyarr+1)))
	    		  break;
	    	  else if ((*keyi == *keyarr) && (*(keyi+1) < *(keyarr+1)))
	    	  {
	    		  p->next = *(bucket_vec[entry].begin()+i);
	    		  (*(bucket_vec[entry].begin()+i))->pre = p;

	    		  if (i != 0)
	    		  {
		    		  p->pre = *(bucket_vec[entry].begin()+i-1);
		    		  (*(bucket_vec[entry].begin()+i-1))->next = p;
	    		  }
	    		  bucket_vec[entry].insert(bucket_vec[entry].begin()+i,p);
	    		  return p;
	    	  }
	    	  else if ((*keyi == *keyarr) && (*(keyi+1) == *(keyarr+1)) && (*(keyi+2) > *(keyarr+2)))
	    		  break;
	          else //if ((*keyi == *keyarr) && (*(keyi+1) == *(keyarr+1)) && (*(keyi+2) == *(keyarr+2)))
	          {
	    		  p->next = *(bucket_vec[entry].begin()+i);
	    		  (*(bucket_vec[entry].begin()+i))->pre = p;

	    		  if (i != 0)
	    		  {
		    		  p->pre = *(bucket_vec[entry].begin()+i-1);
		    		  (*(bucket_vec[entry].begin()+i-1))->next = p;
	    		  }
	    		  bucket_vec[entry].insert(bucket_vec[entry].begin()+i,p);
	    		  return p;
	          }
	      }

	      //if keyi larger than maximum key in the vector, put it at the end of the vector
	      p->pre = *(bucket_vec[entry].end()-1);
	      (*(bucket_vec[entry].end()-1))->next = p;
	      bucket_vec[entry].push_back(p);
		  return p;
	  }
	  else
	  {
	      int i0, i1, i2;
	      i0 = 0;
	      i1 = size / 2;
	      i2 = size - 1;
	      while ((i2 - i0) > HASHTABLE_LOOKUP_LINSEARCH)
	      {
	    	 keyarr = (*(bucket_vec[entry].begin()+i1))->key;

	         if (*keyi > *keyarr)
	         {
	             i0 = i1 + 1;
	             i1 = (i0 + i2) / 2;
	         }
	         else if (*keyi < *keyarr)
	         {
	            i2 = i1;
	            i1 = (i0 + i2) / 2;
	         }
	         else if ((*keyi == *keyarr) && (*(keyi+1) > *(keyarr+1)))
	         {
	             i0 = i1 + 1;
	             i1 = (i0 + i2) / 2;
	         }
	         else if ((*keyi == *keyarr) && (*(keyi+1) < *(keyarr+1)))
	         {
	            i2 = i1;
	            i1 = (i0 + i2) / 2;
	         }
	         else if ((*keyi == *keyarr) && (*(keyi+1) == *(keyarr+1)) && (*(keyi+2) > *(keyarr+2)))
	         {
	             i0 = i1 + 1;
	             i1 = (i0 + i2) / 2;
	         }
	         else
	         {
	            i2 = i1;
	            i1 = (i0 + i2) / 2;
	         }
	      }

	      for (i = i0; i <= i2; i++)
	      {
	    	  keyarr = (*(bucket_vec[entry].begin()+i))->key;
	    	  if (*keyi > *keyarr)
	    		  break;
	    	  else if (*keyi < *keyarr)
	    	  {
	    		  p->next = *(bucket_vec[entry].begin()+i);
	    		  (*(bucket_vec[entry].begin()+i))->pre = p;

	    		  if (i != 0)
	    		  {
		    		  p->pre = *(bucket_vec[entry].begin()+i-1);
		    		  (*(bucket_vec[entry].begin()+i-1))->next = p;
	    		  }
	    		  bucket_vec[entry].insert(bucket_vec[entry].begin()+i,p);
	    		  return p;
	    	  }
	    	  else if ((*keyi == *keyarr) && (*(keyi+1) > *(keyarr+1)))
	    		  break;
	    	  else if ((*keyi == *keyarr) && (*(keyi+1) < *(keyarr+1)))
	    	  {
	    		  p->next = *(bucket_vec[entry].begin()+i);
	    		  (*(bucket_vec[entry].begin()+i))->pre = p;

	    		  if (i != 0)
	    		  {
		    		  p->pre = *(bucket_vec[entry].begin()+i-1);
		    		  (*(bucket_vec[entry].begin()+i-1))->next = p;
	    		  }
	    		  bucket_vec[entry].insert(bucket_vec[entry].begin()+i,p);
	    		  return p;
	    	  }
	    	  else if ((*keyi == *keyarr) && (*(keyi+1) == *(keyarr+1)) && (*(keyi+2) > *(keyarr+2)))
	    		  break;
	          else //if ((*keyi == *keyarr) && (*(keyi+1) == *(keyarr+1)) && (*(keyi+2) == *(keyarr+2)))
	          {
	    		  p->next = *(bucket_vec[entry].begin()+i);
	    		  (*(bucket_vec[entry].begin()+i))->pre = p;

	    		  if (i != 0)
	    		  {
		    		  p->pre = *(bucket_vec[entry].begin()+i-1);
		    		  (*(bucket_vec[entry].begin()+i-1))->next = p;
	    		  }
	    		  bucket_vec[entry].insert(bucket_vec[entry].begin()+i,p);
	    		  return p;
	          }
	      }//end of for

	      //if keyi larger than maximum key in the vector, put it at the end of the vector
	      p->pre = *(bucket_vec[entry].end()-1);
	      (*(bucket_vec[entry].end()-1))->next = p;
	      bucket_vec[entry].push_back(p);
		  return p;
	  }// end else -->when vector size is smaller than HASHTABLE_LOOKUP_LINSEARCH
  }// end else --> when the hashtable row is not occupied!

}

//
//THashEntryPtr THashTable::addElement (int entry, unsigned key[])
//{
//  THashEntryPtr
//    p = new THashEntry (key);
//
//  if ((bucket[entry]))          //this place is already occupied
//  {
//    THashEntryPtr
//      currentPtr = bucket[entry];
//
//    while (currentPtr != 0 && (key[0] > currentPtr->key[0]))
//    {
//      p->pre = currentPtr;
//      currentPtr = currentPtr->next;
//    }
//
//    if (currentPtr != 0 && key[0] == currentPtr->key[0])
//    {
//      while (currentPtr != 0 && (key[1] > currentPtr->key[1]))
//      {
//        p->pre = currentPtr;
//        currentPtr = currentPtr->next;
//      }
//    }
//
//    if (currentPtr != 0 && key[0] == currentPtr->key[0] && key[1] == currentPtr->key[1])
//    {
//      while (currentPtr != 0 && (key[2] > currentPtr->key[2]))
//      {
//        p->pre = currentPtr;
//        currentPtr = currentPtr->next;
//      }
//    }
//
//    if (currentPtr)
//      currentPtr->pre = p;
//    p->next = currentPtr;
//    currentPtr = p->pre;
//    if (currentPtr)
//      currentPtr->next = p;
//    else
//      bucket[entry] = p;//if previous previous does not exist, that means this should be the first element in the hash table
//  }
//
//  //  p->next = *(bucket+entry);        //add the bucket to the head
//  else
//    bucket[entry] = p;          //else eliminate it
//  return p;
//}

void *
THashTable::lookup (unsigned * key)
{
  int entry = hash (key);
  THashEntryPtr p = searchBucket (entry, key);

  if (!p)
    return NULL;                //if not found, return 0
  return p->value;              //if found return a pointer
}

// lookup using key structure
void *
THashTable::lookup (TKey kstr)
{
  unsigned key[TKEYLENGTH];

  kstr.fill_key (key);
  int entry = hash (key);
  THashEntryPtr p = searchBucket (entry, key);

  if (!p)
    return NULL;
  return p->value;
}

void
THashTable::add (unsigned *key, void *value)
{
  int entry = hash (key);
  THashEntryPtr p = searchBucket (entry, key);

  if (p == NULL)    //make sure that the same THashEntryPtr does not exist in the hashtable
  {
    p = addElement (entry, key);
    p->value = value;
  }
  return;
}

void
THashTable::remove (unsigned *key)
{
  int entry = hash (key);
  int i;
  THashEntryPtr p = searchBucket (entry, key, &i);

  if (p == NULL)
    return;
  if ( i == 0)
  {
	  bucket_vec[entry].erase(bucket_vec[entry].begin());
      delete p;
  }
  else
  {
    if (!(p->next))
    {
      delete p;
      bucket_vec[entry].erase(bucket_vec[entry].begin()+i);
    }

    else
    {
      (p->pre)->next = p->next;
      (p->next)->pre = p->pre;
      bucket_vec[entry].erase(bucket_vec[entry].begin()+i);
      delete p;
    }
  }
}

int THashTable::hash (unsigned * key)
{
  int S02 = NBUCKETS/10;
  int S01 = NBUCKETS-S02;
  int i1 = (int) ((double) key[0] / umax * S01);
  int i2 = (int) ((double) key[1] / umax * S02);

  return (i1 + i2);
}

/***************************************************
 *          Hashtable iterator
 ***************************************************/
THTIterator::THTIterator (THashTable * ht)
{
  table = ht;
  current = table->getBucket (0);
  size = table->get_no_of_buckets ();
  index = 0;
}

THashEntryPtr THTIterator::getNextBucket ()
{
  while (table->IsBucketEmpty(++index) && index < size);
  if (index == size)
    return NULL;

  return table->getBucket (index);
}

void *
THTIterator::next ()
{
  void *value;

  if (current)                  // if current pointer has more links
  {
    value = current->value;
    current = current->next;    // move to next link, after getting the value
    return value;
  }
  current = getNextBucket ();   // get next valid Hashtable entry
  if (current == NULL)
    return NULL;
  value = current->value;
  current = current->next;      // move to next link, after getting the value
  return value;
}
