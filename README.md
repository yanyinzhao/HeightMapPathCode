# Efficient Proximity Queries on Simplified Height Maps

## Overview

This project provides the implementation for efficient proximity queries on simplified height maps. We refer the readers to our paper for more details.

We compared 13 algorithms as follows:

- Sur-SimQue-AdpM (baseline simplification algorithm)
- Net-SimQue-AdpM (baseline simplification algorithm)
- Mes-SimQue-AdpM (baseline simplification algorithm)
- Mem-SimQue-LS (variation of our simplification algorithm)
- Mem-SimQue-LST (variation of our simplification algorithm)
- Mem-SimQue-LQT1 (variation of our simplification algorithm)
- Mem-SimQue-LQT2 (variation of our simplification algorithm)
- Mem-SimQue (our simplification algorithm)
- Unf-Que-AdpM (baseline shortest path query algorithm)
- Ste-Que-AdpM (baseline shortest path query algorithm)
- Dij-Que-AdpM (baseline shortest path query algorithm)
- Con-Que-AdpM (baseline shortest path query algorithm)
- Eff-Que (our shortest path query algorithm)

Make sure there is a folder called "input/" and a folder called "output/" under the working directory. They will be used for storing the input/output files.

The source code are stored in "src/" folder.

## Requirment

The project requires the CGAL library. Please download CGAL library at https://www.cgal.org/download.html.

## Dataset

The dataset are stored in "input/" folder.

The datasets are as follows:

- "BH_1014.png" (small version default resolution BH height map dataset with 1014 pixels)
- "BH_10086.png" (small version default resolution BH height map dataset with 10086 pixels)
- "BH_500835.png" (large version default resolution BH height map dataset with 500835 pixels)
- "BH_1000414.png" (large version multiresolution BH height map dataset with 1000414 pixels)
- "BH_1500996.png" (large version multiresolution BH height map dataset with 1500996 pixels)
- "BH_2001610.png" (large version multiresolution BH height map dataset with 2001610 pixels)
- "BH_2502596.png" (large version multiresolution BH height map dataset with 2502596 pixels)
- "EP_10062.png" (small version default resolution EP height map dataset with 10062 pixels)
- "EP_20130.png" (small version multiresolution EP height map dataset with 20130 pixels)
- "EP_30098.png" (small version multiresolution EP height map dataset with 30098 pixels)
- "EP_40076.png" (small version multiresolution EP height map dataset with 40076 pixels)
- "EP_50373.png" (small version multiresolution EP height map dataset with 50373 pixels)
- "EP_500384.png" (large version default resolution EP height map dataset with 500384 pixels)
- "EP_1001040.png" (large version multiresolution EP height map dataset with 1001040 pixels)
- "EP_1501578.png" (large version multiresolution EP height map dataset with 1501578 pixels)
- "EP_2001536.png" (large version multiresolution EP height map dataset with 2001536 pixels)
- "EP_2500560.png" (large version multiresolution EP height map dataset with 2500560 pixels)
- "GF_10092.png" (small version default resolution GF height map dataset with 10092 pixels)
- "GF_500208.png" (large version default resolution GF height map dataset with 500208 pixels)
- "GF_1000518.png" (large version multiresolution GF height map dataset with 1000518 pixels)
- "GF_1501668.png" (large version multiresolution GF height map dataset with 1501668 pixels)
- "GF_2000832.png" (large version multiresolution GF height map dataset with 2000832 pixels)
- "GF_2502075.png" (large version multiresolution GF height map dataset with 2502075 pixels)
- "LM_10092.png" (small version default resolution LM height map dataset with 10092 pixels)
- "LM_500208.png" (large version default resolution LM height map dataset with 500208 pixels)
- "LM_1000518.png" (large version multiresolution LM height map dataset with 1000518 pixels)
- "LM_1501668.png" (large version multiresolution LM height map dataset with 1501668 pixels)
- "LM_2000832.png" (large version multiresolution LM height map dataset with 2000832 pixels)
- "LM_2502075.png" (large version multiresolution LM height map dataset with 2502075 pixels)
- "RM_10092.png" (small version default resolution RM height map dataset with 10092 pixels)
- "RM_500208.png" (large version default resolution RM height map dataset with 500208 pixels)
- "RM_1000518.png" (large version multiresolution RM height map dataset with 1000518 pixels)
- "RM_1501668.png" (large version multiresolution RM height map dataset with 1501668 pixels)
- "RM_2000832.png" (large version multiresolution RM height map dataset with 2000832 pixels)
- "RM_2502075.png" (large version multiresolution RM height map dataset with 2502075 pixels)

## Compile command

```
mkdir build
cd build
cmake ..
make
```

## Run command

```
./main [height_map_data_and_index] [epsilon] [run_knn_and_range_query]
```

The meaning for each parameter is as follows:

- [height_map_data_and_index]: an index for height map data (a integer from 0 to 119)
- [epsilon]: the epsilon value (0 < epsilon <= 1)
- [run_knn_and_range_query]: whether to run knn and range query (0 means not running, 1 means running)

For the [height_map_data_and_index], each index value corresponding to a height map data with different number of pixels, their relationships are as follows:

| Index | Height map data | Pixel number |
| ----------- | ----------- | ----------- |
| 0 | BH | 1014 |
| 1 | BH | 10086 |
| 2 | EP | 10062 |
| 3 | EP | 20130 |
| 4 | EP | 30098 |
| 5 | EP | 40076 |
| 6 | EP | 50373 |
| 7 | GF | 10092 |
| 8 | LM | 10092 |
| 9 | RM | 10092 |
| 10 | BH | 500835 |
| 11 | BH | 1000414 |
| 12 | BH | 1500996 |
| 13 | BH | 2001610 |
| 14 | BH | 2502596 |
| 15 | EP | 500384 |
| 16 | EP | 1001040 |
| 17 | EP | 1501578 |
| 18 | EP | 2001536 |
| 19 | EP | 2500560 |
| 20 | GF | 500208 |
| 21 | GF | 1000518 |
| 22 | GF | 1501668 |
| 23 | GF | 2000832 |
| 24 | GF | 2502075 |
| 25 | LM | 500208 |
| 26 | LM | 1000518 |
| 27 | LM | 1501668 |
| 28 | LM | 2000832 |
| 29 | LM | 2502075 |
| 30 | RM | 500208 |
| 31 | RM | 1000518 |
| 32 | RM | 1501668 |
| 33 | RM | 2000832 |
| 34 | RM | 2502075 |

Since Sur-SimQue-AdpM, Net-SimQue-AdpM, Mem-SimQue-LS, and Mem-SimQue-LST are time consuming or having large memory usage, the project will run all algorithms on small-version dataset ([height_map_data_and_index] <= 9). The project will run all algorithms except these four algorithms on original dataset with default 500 POIs ([height_map_data_and_index] > 9).

In addition, we strongly encourage you to set [run_knn_and_range_query] to 0 if you are not conducting experiments. Otherwise, it will take a very long time to run calculate the knn of all POIs. 

An example:

```
./main 0 0.5 0
```

In this example, [height_map_data_and_index] is 0, [epsilon] is 0.5, [run_knn_and_range_query] is 0. So, it will run BH height map dataset, with pixel number equal to 10086, epsilon is 0.5, it will not run knn and range query. It will run all algorithms.

## Output

The output will be stored in "output/output.txt" file. The format will be as follows:

```
[dataset] [pixel_num] [epsilon] [height_map_to_point_cloud_or_terrain_time (ms)] [height_map_to_point_cloud_or_terrain_memory_usage (MB)] [simplification_time (ms)] [memory_usage (MB)] [output_size (MB)] [query_time (ms)] [distance_error_height_map_or_point_cloud] [distance_error_terrain] [knn_query_time] [knn_error_height_map_or_point_cloud] [knn_error_terrain] [range_query_time] [range_error_height_map_or_point_cloud] [range_error_terrain]
```

These information will also be shown in the terminal. 
