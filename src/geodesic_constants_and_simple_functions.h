// Copyright (C) 2008 Danil Kirsanov, MIT License
#ifndef GEODESIC_CONSTANTS_20071231
#define GEODESIC_CONSTANTS_20071231

// some constants and simple math functions

#include <assert.h>
#include <math.h>
#include <limits>
#include <fstream>

namespace geodesic
{

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

	// double const GEODESIC_INF = std::numeric_limits<double>::max();
	double const GEODESIC_INF = 1e100;

	// in order to avoid numerical problems with "infinitely small" intervals,
	// we drop all the intervals smaller than SMALLEST_INTERVAL_RATIO*edge_length
	double const SMALLEST_INTERVAL_RATIO = 1e-6;
	// double const SMALL_EPSILON = 1e-10;

	inline double cos_from_edges(double const a, // compute the cosine of the angle given the lengths of the edges
								 double const b,
								 double const c)
	{
		assert(a > 1e-50);
		assert(b > 1e-50);
		assert(c > 1e-50);

		double result = (b * b + c * c - a * a) / (2.0 * b * c);
		result = std::max(result, -1.0);
		return std::min(result, 1.0);
	}

	inline double angle_from_edges(double const a, // compute the cosine of the angle given the lengths of the edges
								   double const b,
								   double const c)
	{
		return acos(cos_from_edges(a, b, c));
	}

	template <class Points, class Faces>
	inline bool read_mesh_from_file(std::string filename,
									Points &points,
									Faces &faces)
	{
		std::ifstream file(&filename[0]);
		assert(file.is_open());
		if (!file.is_open())
			return false;

		std::string dum;
		unsigned num_edges;
		file >> dum;
		unsigned num_points;
		file >> num_points;
		assert(num_points >= 3);

		unsigned num_faces;
		file >> num_faces;
		file >> num_edges;

		points.resize(num_points * 3);
		for (typename Points::iterator i = points.begin(); i != points.end(); ++i)
		{
			file >> *i;
		}

		unsigned dummy;
		faces.resize(num_faces * 3);
		for (int i = 0; i < faces.size(); i++)
		{
			if (i % 3 == 0)
			{
				file >> dummy;
				//            std::cout << std::endl << dummy;
			}
			file >> faces[i];
		}
		file.close();
		/*	for(typename Faces::iterator i=faces.begin(); i!=faces.end(); ++i)
			{
				file >> *i;
			}*/

		return true;
	}

	// inline bool read_mesh_from_point_cloud(std::vector<std::string> &terrain,
	// 									   std::vector<double> &points,
	// 									   std::vector<unsigned> &faces)
	// {
	// 	std::string dum;
	// 	unsigned num_points = atoi(terrain[1].c_str());
	// 	unsigned num_faces = atoi(terrain[2].c_str());
	// 	unsigned num_edges = atoi(terrain[3].c_str());
	// 	assert(num_points >= 3);

	// 	// std::cout << num_points << " " << num_faces << " " << num_edges << std::endl;

	// 	points.resize(num_points * 3);
	// 	for (int i = 0; i < points.size(); i++)
	// 	{
	// 		points[i] = atof(terrain[4 + i].c_str());
	// 		// std::cout << points[i] << std::endl;
	// 	}

	// 	faces.resize(num_faces * 3);
	// 	for (int i = 0; i < num_faces; i++)
	// 	{
	// 		faces[3 * i] = atoi(terrain[4 + points.size() + i + 1 + 3 * i].c_str());
	// 		faces[3 * i + 1] = atoi(terrain[4 + points.size() + i + 1 + 3 * i + 1].c_str());
	// 		faces[3 * i + 2] = atoi(terrain[4 + points.size() + i + 1 + 3 * i + 2].c_str());
	// 		// std::cout << faces[3 * i] << " " << faces[3 * i + 1] << " " << faces[3 * i + 2] << std::endl;
	// 	}

	// 	return true;
	// }

} // geodesic

#endif // GEODESIC_CONSTANTS_20071231
