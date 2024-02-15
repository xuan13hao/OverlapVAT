


#ifndef ALIGN_READ_H_
#define ALIGN_READ_H_
#include <string>
#include <vector>
#include <assert.h>
#include "../tools/async_buffer.h"
#include "../basic/Hits.h"
#include "../basic/Statistics.h"
#include "../search/XdropUngapped.h"
#include "align_sequence.h"
#include "../tools/text_buffer.h"
#include "../output/output_buffer.h"
#include "link_segments.h"
#include "splicedSeed.h"
using std::vector;
	static bool contain_suffix(const std::string& str, const std::string& suffix) {
		if (str.length() >= suffix.length()) {
			std::string strSuffix = str.substr(str.length() - suffix.length());
			return strSuffix == suffix;
		}
		return false;
	}
	static std::string remove_suffix(const std::string& str, const std::string& suffix) {
		if (str.length() >= suffix.length() && str.substr(str.length() - suffix.length()) == suffix) {
			return str.substr(0, str.length() - suffix.length());
		}
		return str;
	}
template<typename _val, typename _locr, typename _locl>
void align_read(Output_buffer<_val> &buffer,
		Statistics &stat,
		typename Trace_pt_buffer<_locr,_locl>::Vector::iterator &begin,
		typename Trace_pt_buffer<_locr,_locl>::Vector::iterator &end)
{
	static thread_specific_ptr<vector<local_match<_val> > > local_ptr;
	static thread_specific_ptr<vector<Segment<_val> > > matches_ptr;
	static thread_specific_ptr<vector<char> > transcript_ptr;

	Tls<vector<Segment<_val> > > matches (matches_ptr);
	Tls<vector<local_match<_val> > > local (local_ptr);
	Tls<vector<char> > transcript_buf (transcript_ptr);
	local->clear();
	matches->clear();
	transcript_buf->clear();

	assert(end > begin);
	const size_t hit_count = end - begin;
	local->reserve(hit_count);
	const unsigned contexts = query_contexts();

	const unsigned query = begin->query_/contexts;
	const size_t query_len (QuerySeqs<_val>::data_->length(query*contexts));
	const size_t source_query_len = query_len;
	// const size_t source_query_len = query_translated() ? query_seqs<_val>::data_->reverse_translated_len(query*contexts) : query_len;
	const size_t db_letters = ref_header.letters;
	unsigned padding[6];


	typedef Map<typename vector<Hits<_locr,_locl> >::iterator,typename Hits<_locr,_locl>::template Query_id<1> > Map_t;
	Map_t hits (begin, end);

	typename Map_t::Iterator i = hits.begin();
	while(i.valid()) {

		align_sequence<_val,_locr,_locl>(*matches, stat, *local, padding, db_letters, source_query_len, i.begin(), i.end(), *transcript_buf);
		++i;
	}
	if(matches->size() == 0)
		return;
	std::sort(matches->begin(), matches->end());
	unsigned n_hsp = 0, n_target_seq = 0;
	typename vector<Segment<_val> >::iterator it = matches->begin();
	const int min_raw_score = ScoreMatrix::get().rawscore(VATParameters::min_bit_score == 0
			? ScoreMatrix::get().bitscore(VATParameters::max_evalue, ref_header.letters, query_len) : VATParameters::min_bit_score);
			
	const int top_score = matches->operator[](0).score_;
	while(it < matches->end()) 
	{
		// const bool same_subject = it != matches->begin() && (it-1)->subject_id_ == it->subject_id_;
		// if(!same_subject && it->score_ < min_raw_score)
		// 	break;
		// if(same_subject && (it-1)->score_ == it->score_) {
		// 	++it;
		// 	continue;
		// }
		// if(it->subject_id_ == )
		unsigned sid = it->subject_id_;
		const char* q_name = query_ids::get()[query].c_str();
		string qry_name(q_name);
		const char* s_name = ReferenceIds::get()[sid].c_str();
		string sbj_name(s_name);
		string strand = "+";
		if (qry_name != sbj_name)
		{
			if ((contain_suffix(qry_name,"_minus")))
			{
				strand = "-";
				qry_name = remove_suffix(qry_name,"_minus");
			}
						string result = qry_name+"\t" + std::to_string(it->traceback_->query_begin_) + "\t"+std::to_string(it->traceback_->query_begin_+it->traceback_->query_len_)+"\t"+strand+"\t"+sbj_name+"\t"+std::to_string(it->traceback_->subject_begin_)+"\t"+std::to_string(it->traceback_->subject_begin_+it->traceback_->subject_len_)+"\n";
			cout<<result;
			// cout<<qry_name<<"\t"<<it->traceback_->query_begin_<<"\t"<<it->traceback_->query_begin_+it->traceback_->query_len_<<"\t"<<strand<<"\t"<<sbj_name<<"\t"<<it->traceback_->subject_begin_<<"\t"<<it->traceback_->subject_begin_+it->traceback_->subject_len_<<endl;
		}
		++it;
		++n_hsp;


	}
	// if(n_hsp > 0)
	// 	buffer.finish_query_record();
	stat.inc(Statistics::OUT_MATCHES, matches->size());
	if(ref_header.n_blocks == 1) {
		stat.inc(Statistics::MATCHES, n_hsp);
		if(n_hsp > 0)
			stat.inc(Statistics::ALIGNED);
	}
	
}

#endif /* ALIGN_READ_H_ */
// #ifndef ALIGN_READ_H_
// #define ALIGN_READ_H_
// #include <string>
// #include <vector>
// #include <assert.h>
// #include <iostream>
// #include "../tools/async_buffer.h"
// #include "../basic/Hits.h"
// #include "../basic/Statistics.h"
// #include "../search/XdropUngapped.h"
// #include "align_sequence.h"
// #include "../tools/text_buffer.h"
// #include "../output/output_buffer.h"
// #include "link_segments.h"
// #include "splicedSeed.h"
// using std::vector;
// 	static bool contain_suffix(const std::string& str, const std::string& suffix) {
// 		if (str.length() >= suffix.length()) {
// 			std::string strSuffix = str.substr(str.length() - suffix.length());
// 			return strSuffix == suffix;
// 		}
// 		return false;
// 	}
// 	static std::string remove_suffix(const std::string& str, const std::string& suffix) {
// 		if (str.length() >= suffix.length() && str.substr(str.length() - suffix.length()) == suffix) {
// 			return str.substr(0, str.length() - suffix.length());
// 		}
// 		return str;
// 	}
// template<typename _val, typename _locr, typename _locl>
// void align_read(Output_buffer<_val> &buffer,
// 		Statistics &stat,
// 		typename Trace_pt_buffer<_locr,_locl>::Vector::iterator &begin,
// 		typename Trace_pt_buffer<_locr,_locl>::Vector::iterator &end)
// {
// 	static thread_specific_ptr<vector<local_match<_val> > > local_ptr;
// 	static thread_specific_ptr<vector<Segment<_val> > > matches_ptr;
// 	static thread_specific_ptr<vector<char> > transcript_ptr;

// 	Tls<vector<Segment<_val> > > matches (matches_ptr);
// 	Tls<vector<local_match<_val> > > local (local_ptr);
// 	Tls<vector<char> > transcript_buf (transcript_ptr);
// 	local->clear();
// 	matches->clear();
// 	transcript_buf->clear();

// 	assert(end > begin);
// 	const size_t hit_count = end - begin;
// 	local->reserve(hit_count);
// 	const unsigned contexts = query_contexts();

// 	const unsigned query = begin->query_/contexts;
// 	const size_t query_len (QuerySeqs<_val>::data_->length(query*contexts));
// 	const size_t source_query_len = query_len;
// 	// const size_t source_query_len = query_translated() ? query_seqs<_val>::data_->reverse_translated_len(query*contexts) : query_len;
// 	const size_t db_letters = ref_header.letters;
// 	unsigned padding[6];


// 	typedef Map<typename vector<Hits<_locr,_locl> >::iterator,typename Hits<_locr,_locl>::template Query_id<1> > Map_t;
// 	Map_t hits (begin, end);

// 	typename Map_t::Iterator i = hits.begin();
// 	while(i.valid()) {

// 		align_sequence<_val,_locr,_locl>(*matches, stat, *local, padding, db_letters, source_query_len, i.begin(), i.end(), *transcript_buf);
// 		++i;
// 	}
// 	if(matches->size() == 0)
// 		return;
// 	std::sort(matches->begin(), matches->end());
// 	unsigned n_hsp = 0, n_target_seq = 0;
// 	typename vector<Segment<_val> >::iterator it = matches->begin();
// 	const int min_raw_score = ScoreMatrix::get().rawscore(VATParameters::min_bit_score == 0
// 			? ScoreMatrix::get().bitscore(VATParameters::max_evalue, ref_header.letters, query_len) : VATParameters::min_bit_score);
			
// 	const int top_score = matches->operator[](0).score_;
// 	while(it < matches->end()) 
// 	{
// 		const bool same_subject = it != matches->begin() && (it-1)->subject_id_ == it->subject_id_;
// 		if(!same_subject && it->score_ < min_raw_score)
// 			break;
// 		unsigned sid = it->subject_id_;
// 		const char* q_name = query_ids::get()[query].c_str();
// 		string qry_name(q_name);
// 		const char* s_name = ReferenceIds::get()[sid].c_str();
// 		string sbj_name(s_name);
// 		string strand = "+";
// 		// std::lock_guard<std::mutex> lock(fileMutex); 
// 		if (qry_name != sbj_name)
// 		{
// 			if ((contain_suffix(qry_name,"_minus")))
// 			{
// 				strand = "-";
// 				qry_name = remove_suffix(qry_name,"_minus");
// 			}
// 			// printf("%s\t%d\t%d\t%s\t%s\t%d\t%d\n",
// 			// 		qry_name.c_str(),
// 			// 		it->traceback_->query_begin_,
// 			// 		it->traceback_->query_begin_ + it->traceback_->query_len_,
// 			// 		strand.c_str(),
// 			// 		sbj_name.c_str(),
// 			// 		it->traceback_->subject_begin_,
// 			// 		it->traceback_->subject_begin_ + it->traceback_->subject_len_);
// 			// string result = qry_name+"\t" + std::to_string(it->traceback_->query_begin_) + "\t"+std::to_string(it->traceback_->query_begin_+it->traceback_->query_len_)+"\t"+strand+"\t"+sbj_name+"\t"+std::to_string(it->traceback_->subject_begin_)+"\t"+std::to_string(it->traceback_->subject_begin_+it->traceback_->subject_len_)+"\n";
// 			// cout<<result;
// 			cout <<qry_name<<"\t"<<it->traceback_->query_begin_<<"\t"<<it->traceback_->query_begin_+it->traceback_->query_len_<<"\t"<<strand<<"\t"<<sbj_name<<"\t"<<it->traceback_->subject_begin_<<"\t"<<it->traceback_->subject_begin_+it->traceback_->subject_len_<<endl;
// 		}
// 		// if(n_hsp == 0)
// 		// 	buffer.write_query_record(query);

// 		// buffer.print_match(*it, source_query_len, QuerySeqs<_val>::get()[query*contexts + it->frame_], query, *transcript_buf);
// 		// ++it;
// 		// ++n_hsp;


// 	}
// 	// fclose(outputFile);
// 	// if(n_hsp > 0)
// 	// 	buffer.finish_query_record();
// 	// stat.inc(Statistics::OUT_MATCHES, matches->size());
// 	// if(ref_header.n_blocks == 1) {
// 	// 	stat.inc(Statistics::MATCHES, n_hsp);
// 	// 	if(n_hsp > 0)
// 	// 		stat.inc(Statistics::ALIGNED);
// 	// }
	
// }

// #endif /* ALIGN_READ_H_ */
