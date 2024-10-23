#include <string.h>
#include <fstream>
#include <iostream>

#include "zmorton.hpp"
#include "binhash.hpp"
#include "state.hpp"
#include <cmath>
#include <omp.h>
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

inline int particle_bucket(float fx, float fy, float fz, float h, int xd, int yd, int zd,int bound)
{
    // Fast float-to-int conversion using floorf and casting
    int ix = static_cast<int>(fx) + xd;
    int iy = static_cast<int>(fy) + yd;
    int iz = static_cast<int>(fz) + zd;

    // Branchless check for negative indices and upper boundary
    if ((unsigned int)ix >= (unsigned int)bound ||
        (unsigned int)iy >= (unsigned int)bound ||
        (unsigned int)iz >= (unsigned int)bound) {
        return -1; // Return an invalid bucket index
    }

    // Mask the indices to ensure they are within bounds
    ix &= HASH_MASK;
    iy &= HASH_MASK;
    iz &= HASH_MASK;

    // Encode the indices
    int bucket = zm_encode(ix, iy, iz);

    return bucket;
}

unsigned particle_neighborhood(int* buckets, particle_t* p, float h)
{
    unsigned bin_index = 0;
    // Precompute reciprocal of h
    float inv_h = 1.0f / h;
    // Precompute bound as an integer
    int bound = (int)(inv_h) + 1;
    // Calculate indices using multiplication instead of division
    float fx = p->x[0] * inv_h;
    float fy = p->x[1] * inv_h;
    float fz = p->x[2] * inv_h;

        // Define the 27 neighbor offsets as a compile-time constant
    static const int offsets[27][3] = {
        {-1, -1, -1}, {-1, -1,  0}, {-1, -1,  1},
        {-1,  0, -1}, {-1,  0,  0}, {-1,  0,  1},
        {-1,  1, -1}, {-1,  1,  0}, {-1,  1,  1},
        { 0, -1, -1}, { 0, -1,  0}, { 0, -1,  1},
        { 0,  0, -1}, { 0,  0,  0}, { 0,  0,  1},
        { 0,  1, -1}, { 0,  1,  0}, { 0,  1,  1},
        { 1, -1, -1}, { 1, -1,  0}, { 1, -1,  1},
        { 1,  0, -1}, { 1,  0,  0}, { 1,  0,  1},
        { 1,  1, -1}, { 1,  1,  0}, { 1,  1,  1}
    };
    // Iterate through the offsets and compute bucket indices
    for(int n = 0; n < 27; ++n)
    {
        const int i = offsets[n][0];
        const int j = offsets[n][1];
        const int k = offsets[n][2];
        buckets[bin_index++] = particle_bucket(fx, fy, fz, h, i, j, k, bound);
    }
    return bin_index;
}

void hash_particles(sim_state_t* s, float h)
{
    // Initialize per-bucket locks
    static omp_lock_t locks[HASH_SIZE];
    static bool locks_initialized = false;

    #pragma omp parallel
    {
        #pragma omp single
        {
            if (!locks_initialized) {
                for (int i = 0; i < HASH_SIZE; ++i) {
                    omp_init_lock(&locks[i]);
                }
                locks_initialized = true;
            }
        }
    }

    // Clear hashmap in parallel
    #pragma omp parallel for schedule(static)
    for (int i = 0; i < HASH_SIZE; ++i) {
        s->hash[i] = nullptr;
    }

    // Populating the hashmap in parallel
    #pragma omp parallel for schedule(dynamic, 1000)
    for (int i = 0; i < s->n; ++i) {
        particle_t* p = s->part + i;

        // Precompute reciprocal of h and bound
        float inv_h = 1.0f / h;
        int bound = (int)(inv_h) + 1;

        // Calculate indices
        float fx = p->x[0] * inv_h;
        float fy = p->x[1] * inv_h;
        float fz = p->x[2] * inv_h;
        int b = particle_bucket(fx, fy, fz, h,0,0,0, bound);
        if (b == -1 || b >= HASH_SIZE) {
            printf("BAD HASH BUCKET, %f %f %f MATH %d %d %d\n",p->x[0],p->x[1],p->x[2],static_cast<int>(p->x[0] / h),static_cast<int>(p->x[1] / h),static_cast<int>(p->x[2] / h));
            continue; // Skip invalid indices
        }


        // Lock the bucket, update it, and unlock
        omp_set_lock(&locks[b]);
        p->next = s->hash[b];
        s->hash[b] = p;
        omp_unset_lock(&locks[b]);

    }

    // Locks can be destroyed at the end of the program if necessary
}
// void hash_particles(sim_state_t* s, float h)
// {
//     /* BEGIN TASK */
//     /* END TASK */

//     // Clear hashmap at every time step
//     memset(s->hash, 0, HASH_SIZE * sizeof(particle_t*));
    
//     // Populating the hashmap
//     for (int i = 0; i < s->n; ++i) {
//         particle_t *p = s->part + i;
//         // Precompute reciprocal of h and bound
//         float inv_h = 1.0f / h;
//         int bound = (int)(inv_h) + 1;

//         // Calculate indices
//         float fx = p->x[0] * inv_h;
//         float fy = p->x[1] * inv_h;
//         float fz = p->x[2] * inv_h;
//         int b = particle_bucket(fx, fy, fz, h,0,0,0, bound);
//         if (b == -1 || b >= HASH_SIZE) {
//             printf("BAD HASH BUCKET, %f %f %f MATH %d %d %d\n",p->x[0],p->x[1],p->x[2],static_cast<int>(p->x[0] / h),static_cast<int>(p->x[1] / h),static_cast<int>(p->x[2] / h));
//             continue; // Skip invalid indices
//         }
//         p->next = s->hash[b];
//         s->hash[b] = p;
//     }
// }
