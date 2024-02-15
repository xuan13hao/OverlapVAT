#include <iostream>
#include <stdio.h>
#include "../minimap2.h"
using namespace std;
using std::cout;
using std::endl;
//g++ -g -O2 test.cpp ../minimap2/libminimap2.a -lz -lm -o minimap-lite


int  main(int argc, char *argv[])
{
		mm_idxopt_t iopt;
		mm_mapopt_t mopt;
        int ret;
		mm_set_opt(0, &iopt, &mopt);
		int n_threads = 4;
		mm_verbose = 2;
		mopt.flag |= MM_F_CIGAR; // perform alignment
		mm_idx_reader_t *r = mm_idx_reader_open(argv[1], &iopt, 0);
		mm_idx_t *mi;
		while ((mi = mm_idx_reader_read(r, n_threads)) != 0) { // traverse each part of the index
			mm_mapopt_update(&mopt, mi); // this sets the maximum minimizer occurrence; TODO: set a better default in mm_mapopt_init()!
			ret = 0;
            ret = mm_map_file_frag(mi, 1, (const char**)&argv[2], &mopt, n_threads);
            if (ret < 0) break;
			mm_idx_destroy(mi);
		}
		mm_idx_reader_close(r); // close the index reader
}
