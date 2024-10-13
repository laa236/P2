#include <string.h>

#include "zmorton.hpp"
#include "binhash.hpp"

/*@q
 * ====================================================================
 */

/*@T
 * \subsection{Spatial hashing implementation}
 * 
 * In the current implementation, we assume [[HASH_DIM]] is $2^b$,
 * so that computing a bitwise of an integer with [[HASH_DIM]] extracts
 * the $b$ lowest-order bits.  We could make [[HASH_DIM]] be something
 * other than a power of two, but we would then need to compute an integer
 * modulus or something of that sort.
 * 
 *@c*/

#define HASH_MASK (HASH_DIM-1)

unsigned particle_bucket(particle_t* p, float h, int xd, int yd, int zd)
{
    unsigned ix = p->x[0]/h + xd;
    unsigned iy = p->x[1]/h + yd;
    unsigned iz = p->x[2]/h + zd;
    return zm_encode(ix & HASH_MASK, iy & HASH_MASK, iz & HASH_MASK);
}

unsigned particle_neighborhood(unsigned* buckets, particle_t* p, float h)
{
    /* BEGIN TASK */
    /* END TASK */
    unsigned bin_index = 0;
    for(int i=-1;i<2;++i){
        for(int j=-1;j<2;++j){
            for(int k=-1;k<2;++k){
                buckets[bin_index] = particle_bucket(p,h,i,j,k);
                bin_index++;
            }
        }
    }
    return bin_index;
}

void hash_particles(sim_state_t* s, float h)
{
    /* BEGIN TASK */
    /* END TASK */

    // clear hashmap at every time step
    memset(s->hash, 0, HASH_SIZE * sizeof(particle_t*));
    
    //populating the hashmap
    for (int i = 0; i < s->n; ++i) {
        particle_t *p = s->part + i;
        unsigned b = particle_bucket(p, h, 0, 0, 0);
        p->next = s->hash[b];
        s->hash[b] = p;
    }
}
