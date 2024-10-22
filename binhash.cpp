#include <string.h>
#include <fstream>
#include <iostream>

#include "zmorton.hpp"
#include "binhash.hpp"
#include "state.hpp"

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

int particle_bucket(particle_t* p, float h, int xd, int yd, int zd)
{
    // Calculate indices
    int ix = static_cast<int>(p->x[0] / h) + xd;
    int iy = static_cast<int>(p->x[1] / h) + yd;
    int iz = static_cast<int>(p->x[2] / h) + zd;

    // Check for negative indices
    if (ix < 0 || iy < 0 || iz < 0) {
        //printf("BAD BUCKET = %f,%f,%f div %f + %d,%d,%d = %d,%d,%d\n",p->x[0],p->x[1],p->x[2],h,xd,yd,zd,ix,iy,iz);
        return -1; // Return an invalid bucket index
    }
    int bound = (1/h)+1;
    // Check for upper boundary conditions
    if (ix >= bound || iy >= bound || iz >= bound) {
        // std::cout << "Upper boundary index detected: ix=" << ix << ", iy=" << iy << ", iz=" << iz << "\n";
        //printf("BAD BUCKET = %f,%f,%f div %f + %d,%d,%d = %d,%d,%d\n",p->x[0],p->x[1],p->x[2],h,xd,yd,zd,ix,iy,iz);
        return -1; // Return an invalid bucket index
    }

    // Mask the indices to ensure they are within bounds
    ix &= HASH_MASK;
    iy &= HASH_MASK;
    iz &= HASH_MASK;

    // Print the masked indices
    // std::cout << "Masked indices: ix=" << ix << ", iy=" << iy << ", iz=" << iz << "\n";

    // Encode the indices
    int bucket = zm_encode(ix, iy, iz);
    // std::cout << "Encoded bucket: " << bucket << "\n";

    return bucket;
}

unsigned particle_neighborhood(int* buckets, particle_t* p, float h)
{
    unsigned bin_index = 0;
    for (int i = -1; i < 2; ++i) {
        for (int j = -1; j < 2; ++j) {
            for (int k = -1; k < 2; ++k) {
                buckets[bin_index] = particle_bucket(p, h, i, j, k);
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

    // Clear hashmap at every time step
    memset(s->hash, 0, HASH_SIZE * sizeof(particle_t*));
    
    // Populating the hashmap
    for (int i = 0; i < s->n; ++i) {
        particle_t *p = s->part + i;
        int b = particle_bucket(p, h, 0, 0, 0);
        if (b == -1 || b >= HASH_SIZE) {
            printf("BAD HASH BUCKET, %f %f %f MATH %d %d %d\n",p->x[0],p->x[1],p->x[2],static_cast<int>(p->x[0] / h),static_cast<int>(p->x[1] / h),static_cast<int>(p->x[2] / h));
            continue; // Skip invalid indices
        }
        p->next = s->hash[b];
        s->hash[b] = p;
    }
}
