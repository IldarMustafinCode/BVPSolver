
#include <gtest/gtest.h>
#include <algorithm> //sort

#include "src/Solver/LinProblem/ConnGraph.h"


using namespace linear_problem;


// Nodes scheme
// |_____|_____|
// 6     7     8
// |_____|_____|
// 3     4     5
// |_____|_____|
// 0     1     2
TEST(ConGraph2D, PointNmbr_3x3)
{

	std::vector<std::vector<int>>
		cntrl{{1, 3}, {0, 2, 4}, {1, 5}, {4, 0, 6}, {3, 5, 1, 7}, {2, 4, 8}, {3, 7}, {6, 8, 4}, {7, 5}};

	auto Nx = cntrl.size();

	ConnGraphCart2D graph2D(Nx, Nx);
	auto graph = graph2D.GetGraph();

	for (int ix = 0; ix < graph.size(); ix++)
	{
		auto cntrlIx = cntrl[ix];
		std::sort(cntrlIx.begin(), cntrlIx.end());
		for (int in = 0; in < cntrlIx.size(); in++)
		{
			ASSERT_EQ(graph[ix][in], cntrlIx[in]);
		}
	}
}

// Nodes scheme

// 2     3  
// |_____|
// 0     1     

TEST(ConGraph2D, PointNmbr_2x2)
{
	std::vector<std::vector<int>>
		cntrl{{1, 2}, {0, 3}, {0, 3}, {2, 1}};

	auto Nx = cntrl.size();

	ConnGraphCart2D graph2D(Nx, Nx);
	auto graph = graph2D.GetGraph();

	for (int ix = 0; ix < graph.size(); ix++)
	{
		auto cntrlIx = cntrl[ix];
		std::sort(cntrlIx.begin(), cntrlIx.end());
		for (int in = 0; in < cntrlIx.size(); in++)
		{
			ASSERT_EQ(graph[ix][in], cntrlIx[in]);
		}
	}
}




int main(int argc, char **argv)
{
	testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}