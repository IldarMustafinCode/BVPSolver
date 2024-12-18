
#include <gtest/gtest.h>
#include <iostream>

#include "src/Solver/LinProblem/ConnGraph.h"


using namespace linear_problem;

TEST(ConGraph1D, PointNmbr_1)
{
	size_t Nx=1;
	ConnGraphCart1D graph1D(Nx);
	auto graph = graph1D.GetGraph();
    
	ASSERT_EQ(graph[0].size(), 0);
}

// Nodes scheme
// |___|
// 0   1

TEST(ConGraph1D, PointNmbr_2)
{

	ConnGraphCart1D graph1D(2);
	const auto& graph = graph1D.GetGraph();

	ASSERT_EQ(graph[0][0], 1);
	ASSERT_EQ(graph[1][0], 0);
}

// Nodes scheme
// |___|___|___|___|___|
// 0   1   2   3   4   5

TEST(ConGraph1D, PointNmbr_6)
{

	std::vector<std::vector<size_t>>
		cntrl{{1}, {0, 2}, {1, 3}, {2, 4}, {3, 5}, {4}};
	auto Nx = cntrl.size();

	ConnGraphCart1D graph1D(Nx);
	const auto& graph = graph1D.GetGraph();

	for (size_t ix = 0; ix < Nx; ++ix)
	{
		for (size_t in = 0; in < cntrl[ix].size(); ++in)
		{
			ASSERT_EQ(graph[ix][in], cntrl[ix][in]);
		}
	}
}

int main(int argc, char **argv)
{
	testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}