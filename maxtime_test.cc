///////////////////////////////////////////////////////////////////////////////
// maxtime_test.cc
//
// Unit tests for maxtime.hh
//
///////////////////////////////////////////////////////////////////////////////


#include <cassert>
#include <sstream>


#include "maxtime.hh"
#include "rubrictest.hh"


int main()
{
	Rubric rubric;
	
	RideVector trivial_rides;
	trivial_rides.push_back(std::shared_ptr<RideItem>(new RideItem("test Ferris Wheel", 10.0, 20.0)));
	trivial_rides.push_back(std::shared_ptr<RideItem>(new RideItem("test Speedway", 4.0, 5.0)));
	
	auto all_rides = load_ride_database("ride.csv");
	assert( all_rides );
	
	auto filtered_rides = filter_ride_vector(*all_rides, 1, 2500, all_rides->size());
	
	//
	rubric.criterion(
		"load_ride_database still works", 2,
		[&]()
		{
			TEST_TRUE("non-null", all_rides);
			TEST_EQUAL("size", 8064, all_rides->size());
		}
	);
	
	//
	rubric.criterion(
		"filter_ride_vector", 2,
		[&]()
		{
			auto
				three = filter_ride_vector(*all_rides, 100, 500, 3),
				ten = filter_ride_vector(*all_rides, 100, 500, 10);
			TEST_TRUE("non-null", three);
			TEST_TRUE("non-null", ten);
			TEST_EQUAL("total_size", 3, three->size());
			TEST_EQUAL("total_size", 10, ten->size());
			TEST_EQUAL("contents", "again amazing mystical vertigo", (*ten)[0]->description());
			TEST_EQUAL("contents", "A short enchanted typhoon", (*ten)[9]->description());
			for (int i = 0; i < 3; i++) {
				TEST_EQUAL("contents", (*three)[i]->description(), (*ten)[i]->description());
			}
		}
	);
	
	//
	rubric.criterion(
		"dynamic_max_time trivial cases", 2,
		[&]()
		{
			std::unique_ptr<RideVector> soln;
			
			soln = dynamic_max_time(trivial_rides, 3);
			TEST_TRUE("non-null", soln);
			TEST_TRUE("empty solution", soln->empty());
			
			soln = dynamic_max_time(trivial_rides, 10);
			TEST_TRUE("non-null", soln);
			TEST_EQUAL("Ferris Wheel only", 1, soln->size());
			TEST_EQUAL("Ferris Wheel only", "test Ferris Wheel", (*soln)[0]->description());
			
			soln = dynamic_max_time(trivial_rides, 9);
			TEST_TRUE("non-null", soln);
			TEST_EQUAL("Speedway only", 1, soln->size());
			TEST_EQUAL("Speedway only", "test Speedway", (*soln)[0]->description());
			
			soln = dynamic_max_time(trivial_rides, 14);
			TEST_TRUE("non-null", soln);
			TEST_EQUAL("Ferris Wheel and Speedway", 2, soln->size());
			TEST_EQUAL("Ferris Wheel and Speedway", "test Speedway", (*soln)[0]->description());
			TEST_EQUAL("Ferris Wheel and Speedway", "test Ferris Wheel", (*soln)[1]->description());
		}
	);
	
	//
	rubric.criterion(
		"dynamic_max_time correctness", 4,
		[&]()
		{
			std::unique_ptr<RideVector> soln_small, soln_large;
			
			soln_small = dynamic_max_time(*filtered_rides, 500),
			soln_large = dynamic_max_time(*filtered_rides, 5000);
			
			//print_ride_vector(*soln_small);
			//print_ride_vector(*soln_large);
			
			TEST_TRUE("non-null", soln_small);
			TEST_TRUE("non-null", soln_large);
			
			TEST_FALSE("non-empty", soln_small->empty());
			TEST_FALSE("non-empty", soln_large->empty());
			
			int
				cost_small,
				cost_large
				;
			double
				time_small,
				time_large
				;
			sum_ride_vector(*soln_small, cost_small, time_small);
			sum_ride_vector(*soln_large, cost_large, time_large);
			
			//	Precision
			time_small	= std::round( time_small	* 100 ) / 100;
			time_large	= std::round( time_large	* 100 ) / 100;
			
			TEST_EQUAL("Small solution cost",	500,		cost_small);
			TEST_EQUAL("Small solution time",	9564.92,	time_small);
			TEST_EQUAL("Large solution cost",	5000,		cost_large);
			TEST_EQUAL("Large solution time",	82766.449999999997,	time_large);
		}
	);
	
	//
	rubric.criterion(
		"exhaustive_max_time trivial cases", 2,
		[&]()
		{
			std::unique_ptr<RideVector> soln;
			
			soln = exhaustive_max_time(trivial_rides, 3);
			TEST_TRUE("non-null", soln);
			TEST_TRUE("empty solution", soln->empty());
			
			soln = exhaustive_max_time(trivial_rides, 10);
			TEST_TRUE("non-null", soln);
			TEST_EQUAL("Ferris Wheel only", 1, soln->size());
			TEST_EQUAL("Ferris Wheel only", "test Ferris Wheel", (*soln)[0]->description());
			
			soln = exhaustive_max_time(trivial_rides, 9);
			TEST_TRUE("non-null", soln);
			TEST_EQUAL("Speedway only", 1, soln->size());
			TEST_EQUAL("Speedway only", "test Speedway", (*soln)[0]->description());
			
			soln = exhaustive_max_time(trivial_rides, 14);
			TEST_TRUE("non-null", soln);
			TEST_EQUAL("Ferris Wheel and Speedway", 2, soln->size());
			TEST_EQUAL("Ferris Wheel and Speedway", "test Ferris Wheel", (*soln)[0]->description());
			TEST_EQUAL("Ferris Wheel and Speedway", "test Speedway", (*soln)[1]->description());
		}
	);
	
	//
	rubric.criterion(
		"exhaustive_max_time correctness", 4,
		[&]()
		{
			std::vector<double> optimal_time_totals =
			{
				500,		1033.05,	1500,	2100,	2400,
				2900,		3400,		4200,	4300,	4600,
				5000,		5400,		5800,	6100,	6500,
				7000,		7500,		8100,	8600,	8700
			};
			
			for ( size_t optimal_index = 0; optimal_index < optimal_time_totals.size(); optimal_index++ )
			{
				size_t n = optimal_index + 1;
				double expected_time = optimal_time_totals[optimal_index];
				
				auto small_rides = filter_ride_vector(*filtered_rides, 1, 2000, n);
				TEST_TRUE("non-null", small_rides);
				
				auto solution = exhaustive_max_time(*small_rides, 2000);
				TEST_TRUE("non-null", solution);
				
				int actual_cost;
				double actual_time;
				sum_ride_vector(*solution, actual_cost, actual_time);
				
				// Round
				expected_time	= std::round( expected_time	/ 100.0) * 100;
				actual_time		= std::round( actual_time	/ 100.0) * 100;
				
				std::stringstream ss;
				ss
					<< "exhaustive search n = " << n << " (optimal index = " << optimal_index << ")"
					<< ", expected time = " << expected_time
					<< " but algorithm found = " << actual_time
					;
				TEST_EQUAL(ss.str(), expected_time, actual_time);
				
				auto dynamic_solution = dynamic_max_time(*small_rides, 2000);
				int dynamic_actual_cost;
				double dynamic_actual_time;
				sum_ride_vector(*solution, dynamic_actual_cost, dynamic_actual_time);
				dynamic_actual_time	= std::round( dynamic_actual_time	/ 100.0) * 100;
				TEST_EQUAL("Exhaustive and dynamic get the same answer", actual_time, dynamic_actual_time);
			}
		}
	);
	
	return rubric.run();
}




