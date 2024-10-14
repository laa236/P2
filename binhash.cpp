#include <string.h>

#include "zmorton.hpp"
#include "binhash.hpp"
#include <stdlib.h>

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

// void hash_particles(sim_state_t* s, float h)
// {
//     /* BEGIN TASK */
//     /* END TASK */

//     // clear hashmap at every time step
//     memset(s->hash, 0, HASH_SIZE * sizeof(particle_t*));
    
//     //populating the hashmap
//     for (int i = 0; i < s->n; ++i) {
//         particle_t *p = s->part + i;
//         unsigned b = particle_bucket(p, h, 0, 0, 0);
//         p->next = s->hash[b];
//         s->hash[b] = p;
//     }
// }
void hash_particles(sim_state_t* s, float h) {
    int n = s->n;

    // Clear the global hash map
    memset(s->hash, 0, HASH_SIZE * sizeof(particle_t*));

    // Allocate thread-local hash tables (one per thread)
    #pragma omp parallel
    {
        particle_t** local_hash = (particle_t**) calloc(HASH_SIZE, sizeof(particle_t*));

        #pragma omp for
        for (int i = 0; i < n; ++i) {
            particle_t* p = s->part + i;
            unsigned b = particle_bucket(p, h, 0, 0, 0);

            // Add particle to local hash table
            p->next = local_hash[b];
            local_hash[b] = p;
        }

        // Critical section to merge local hash into global hash table
        #pragma omp critical
        {
            for (unsigned b = 0; b < HASH_SIZE; ++b) {
                if (local_hash[b] != NULL) {
                    particle_t* p = local_hash[b];
                    while (p != NULL) {
                        particle_t* next = p->next;
                        p->next = s->hash[b];
                        s->hash[b] = p;
                        p = next;
                    }
                }
            }
        }

        // Free thread-local hash table
        free(local_hash);
    }
}
