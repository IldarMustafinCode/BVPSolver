#pragma once

#include <fstream>
#include <vector>

namespace linear_problem
{
	class ConnGraph
	{
	protected:
		std::vector<std::vector<int>> connectivityGraph;

	public:
		ConnGraph() = default;
		ConnGraph(ConnGraph &) = default;

		const std::vector<std::vector<int>> &GetGraph() const
		{
			return connectivityGraph;
		}

		// print connectivityGraph to file
		void printConnectivity()
		{
			std::ofstream myfile;
			myfile.open("test_connections.txt");

			for (int l = 0; l < connectivityGraph.size(); l++)
			{
				for (auto neighbourId : connectivityGraph[l])
				{
					myfile << neighbourId << " ";
				}
				myfile << std::endl;
			}
		}
	};

	class ConnGraphCart1D : public ConnGraph
	{
		size_t Nx;

	public:
		ConnGraphCart1D() = default;
		ConnGraphCart1D(size_t Nx_)
			: Nx(Nx_)
		{
			connectivityGraph.reserve(Nx);
			// global cell number
			size_t i_cell = 0;

			// loop through every cell and determine its neighbours
			for (size_t ix = 0; ix < Nx; ++ix)
			{
				// array of neighboring cell numbers
				std::vector<int> neighbours;

				if (ix > 0)
				{
					neighbours.push_back(i_cell - 1);
				}

				if (ix < Nx - 1)
				{ // the cell, following in X direction
					neighbours.push_back(i_cell + 1);
				}

				connectivityGraph.push_back(neighbours);;
				++i_cell;
			}
		}
	};

	class ConnGraphCart2D : public ConnGraph
	{
		size_t Nx, Ny;

	public:
		ConnGraphCart2D() = default;
		ConnGraphCart2D(size_t Nx, size_t Ny) noexcept
			: Nx{Nx}, Ny{Ny}
		{

			connectivityGraph.reserve(Nx * Ny);
			// global cell number
			size_t i_cell = 0;

			// loop through every cell and determine its neighbours

			for (size_t iy = 0; iy < Ny; ++iy)
			{
				for (size_t ix = 0; ix < Nx; ++ix)
				{
					// array of neighboring cell numbers
					std::vector<int> neighbours;

					if (ix > 0)
					{
						neighbours.push_back(i_cell - 1);
					}

					if (ix < Nx - 1)
					{ // the cell, following in X direction
						neighbours.push_back(i_cell + 1);
					}

					if (iy > 0)
					{ // the cell, previous in Y direction
						neighbours.push_back(i_cell - Nx);
					}

					if (iy < Ny - 1)
					{ // the cell, following in Y direction
						neighbours.push_back(i_cell + Nx);
					}

					connectivityGraph.push_back(neighbours);
					++i_cell;
				}
			}
		}
	};

	class ConnGraphCart3D : public ConnGraph
	{
		size_t Nx, Ny, Nz;

	public:
		ConnGraphCart3D() = default;
		ConnGraphCart3D(size_t Nx, size_t Ny, size_t Nz) noexcept
			: Nx{Nx}, Ny{Ny}, Nz{Nz}
		{

			connectivityGraph.reserve(Nx * Ny * Nz);
			// global cell number
			size_t i_cell = 0;

			// loop through every cell and determine its neighbours
			for (size_t iz = 0; iz < Nz; ++iz)
			{
				for (size_t iy = 0; iy < Ny; ++iy)
				{
					for (size_t ix = 0; ix < Nx; ++ix)
					{
						// array of neighboring cell numbers
						std::vector<int> neighbours;

						if (ix > 0)
						{ // the cell, previous in X direction
							neighbours.push_back(i_cell - 1);
						}

						if (ix < Nx - 1)
						{ // the cell, following in X direction
							neighbours.push_back(i_cell + 1);
						}

						if (iy > 0)
						{ // the cell, previous in Y direction
							neighbours.push_back(i_cell - Nx);
						}
						if (iy < Ny - 1)
						{ // the cell, following in Y direction
							neighbours.push_back(i_cell + Nx);
						}

						if (iz > 0)
						{ // the cell, previous in Z direction
							neighbours.push_back(i_cell - Nx * Ny);
						}
						if (iz < Nz - 1)
						{ // the cell, following in Z direction
							neighbours.push_back(i_cell + Nx * Ny);
						}

						connectivityGraph.push_back(neighbours);
						++i_cell;
					}
				}
			}
		}
	};

}