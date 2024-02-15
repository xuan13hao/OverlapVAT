#ifndef __MMDNA_H__
#define __MMDNA_H__


#include <iostream>
#include <boost/timer/timer.hpp>
#include "../data/Reference.h"
#include "../data/Queries.h"
#include "../basic/Statistics.h"
#include "../basic/ShapeParameter.h"
#include "../output/join_blocks.h"
#include "../align/align_queries.h"
#include "../search/AlignPartition.h"
#include "../basic/ContextSet.h"
#include "SearchContext.h"
#include "minimap2.h"
using std::endl;
using std::cout;
using boost::timer::cpu_timer;
using boost::ptr_vector;
using std::auto_ptr;

template<typename _val, typename _locr>
void MM_DNA( cpu_timer &timer_mapping, cpu_timer &total_timer)
{
	TimerTools timer ("Opening the input file", true);
	timer_mapping.resume();
	timer.go("Opening the output file");
	timer_mapping.stop();
	timer.finish();
	string ref_file = VATParameters::input_subject_file;
	string qry_file  = VATParameters::query_file;

	const char** query_files = stringToConstCharPtrArray(qry_file);
	DNA_MultiThread_MM(ref_file.c_str(),query_files,VATParameters::threads());

	timer.go("Closing the output file");
	timer_mapping.resume();
	timer_mapping.stop();
    free((void*)query_files[0]); // Free the duplicated string
    delete[] query_files; // Free the array itself
	timer.go("Closing the database file");
	timer.finish();
	cout << "Mapping time = " << boost::timer::format(total_timer.elapsed(), 1, "%ws\n");
}

template<typename _val>
void MM_DNA()
{
	cpu_timer timer2, timer_mapping;
	timer_mapping.stop();


	TimerTools timer ("Opening the database", 1);
	timer.finish();
	MM_DNA<_val,uint32_t>( timer_mapping, timer2);
}


#endif // __MMDNA_H__