

#include <gtest/gtest.h>
#include <iostream>

#include "src/Solver/LinProblem/ConnGraph.h"
#include "src/Solver/LinProblem/SparsityPattern.h"

using namespace linear_problem;

TEST(SparsePattern1D, EqNmbr_1)
{
	unsigned char eqNmbr = 1;
	size_t cellNmbr = 5;

	// Connectivity Graph
	size_t Nx = 1;
	ConnGraphCart1D graph1D(Nx);
	const auto& graph = graph1D.GetGraph();

	// Sparsity Pattern
	SparsityPattern spp(eqNmbr, cellNmbr, graph);





	ASSERT_EQ(1, 1);
}



int main(int argc, char **argv)
{
	testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}