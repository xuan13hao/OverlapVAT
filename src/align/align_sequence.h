

#ifndef ALIGN_SEQUENCE_H_
#define ALIGN_SEQUENCE_H_

#include <vector>
#include "../dp/floating_sw.h"
#include "./pairChimera.h"
#include "./ungappedSeed.h"
#include "findwholeGenoSeeds.h"
#include "splicedunGappedSeed.h"
#include "dp_chain_chimera.h"
#include "chaining.h"
using std::vector;

template<typename _val, typename _locr, typename _locl>
void align_sequence(vector<Segment<_val> > &matches,
		Statistics &stat,
		vector<local_match<_val> > &local,
		unsigned *padding,
		size_t db_letters,
		unsigned dna_len,
		typename Trace_pt_buffer<_locr,_locl>::Vector::iterator &begin,
		typename Trace_pt_buffer<_locr,_locl>::Vector::iterator &end,
		vector<char> &transcript_buf)
{
	std::sort(begin, end, Hits<_locr,_locl>::cmp_normalized_subject);
	const unsigned q_num (begin->query_);
	const sequence<const _val> query (QuerySeqs<_val>::get()[q_num]);
	const unsigned frame = q_num % query_contexts();
	const unsigned query_len = query.length();
	padding[frame] = VATParameters::read_padding<_val>(query_len);
	const SequenceSet<_val> *ref = ReferenceSeqs<_val>::data_;
	vector<DiagonalSeeds<_locr,_locl> > diagonalsegment_;
	for(typename Trace_pt_buffer<_locr,_locl>::Vector::iterator i = begin; i != end; ++i)
	{
		if(i != begin && (i->global_diagonal() - (i-1)->global_diagonal()) <= padding[frame]) 
		{
			stat.inc(Statistics::DUPLICATES);
			continue;
		}
		local_match<_val>  tmp;
		// cout<<i->query_<<"\t"<<i->subject_<<"\t"<<i->seed_offset_<<endl;
		std::pair<size_t, size_t> l = ref->local_position(i->subject_);
		const _val* sbj = ref->data(i->subject_);
		const _val* qry = &query[i->seed_offset_];
		DiagonalSeeds<_locr,_locl> ds = ungappedSeeds<_val, _locr,_locl> (qry, sbj,(int)i->seed_offset_,(int)l.second,*i);
		tmp.query_begin_ = unsigned(ds.i);
		tmp.query_len_ = unsigned(ds.len);
		tmp.subject_begin_ = unsigned(ds.j);
		tmp.subject_len_ = unsigned(ds.len);
		tmp.len_ = unsigned(ds.len);
		tmp.score_ = unsigned(ds.score);
		tmp.gap_openings_ = 0;
		tmp.identities_ = 100;
		tmp.mismatches_ = 0;
		tmp.query_anchor_ = 0;
		const int score = tmp.score_;
		Segment<_val> sw_segment(score, frame, &tmp, l.first);
		matches.push_back(sw_segment);	
	}
}

#endif /* ALIGN_SEQUENCE_H_ */
