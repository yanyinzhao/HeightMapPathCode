# Efficient Proximity Queries on Simplified Height Maps

## Overview

This project provides the implementation for efficient proximity queries on simplified height maps. We refer the readers to our paper for more details.

We compared 28 algorithms as follows (the algorithms calculate the path passing on a height map as default):

- TIN-SSimplify-Adapt(HM) (baseline simplification algorithm)
- TIN-SSimplify-Adapt(PC) (baseline simplification algorithm that calculates the path passing on a point cloud)
- TIN-SSimplify (baseline simplification algorithm that calculates the path passing on a TIN)
- TIN-NSimplify-Adapt(HM) (baseline simplification algorithm)
- TIN-NSimplify-Adapt(PC) (baseline simplification algorithm that calculates the path passing on a point cloud)
- TIN-NSimplify (baseline simplification algorithm that calculates the path passing on a TIN)
- HM-Simplify-LS (variation of our simplification algorithm)
- HM-Simplify-LST (variation of our simplification algorithm)
- HM-Simplify-LQT1 (variation of our simplification algorithm)
- HM-Simplify-LQT2 (variation of our simplification algorithm)
- HM-Simplify (our simplification algorithm)
- HM-Simplify-Adapt(PC) (our simplification algorithm that calculate the path passing on a point cloud)
- HM-Simplify-Adapt(TIN) (our simplification algorithm that calculate the path passing on a TIN)
- TIN-ESSP-Adapt(HM) (baseline shortest path query algorithm)
- TIN-ESSP-Adapt(PC) (baseline shortest path query algorithm that calculates the path passing on a point cloud)
- TIN-ESSP (baseline shortest path query algorithm that calculates the path passing on a TIN)
- TIN-ASSP-Adapt(HM) (baseline shortest path query algorithm)
- TIN-ASSP-Adapt(PC) (baseline shortest path query algorithm that calculates the path passing on a point cloud)
- TIN-ASSP (baseline shortest path query algorithm that calculates the path passing on a TIN)
- TIN-NSP-Adapt(HM) (baseline shortest path query algorithm)
- TIN-NSP-Adapt(PC) (baseline shortest path query algorithm that calculates the path passing on a point cloud)
- TIN-NSP (baseline shortest path query algorithm that calculates the path passing on a TIN)
- PC-SP-Adapt(HM) (baseline shortest path query algorithm)
- PC-SP (baseline shortest path query algorithm that calculates the path passing on a point cloud)
- PC-SP-Adapt(TIN) (baseline shortest path query algorithm that calculates the path passing on a TIN)
- HM-SP (our shortest path query algorithm)
- HM-SP-Adapt(PC) (our shortest path query algorithm that calculate the path passing on a point cloud)
- HM-SP-Adapt(TIN) (our shortest path query algorithm that calculate the path passing on a TIN)


Make sure there is a folder called "input/" and a folder called "output/" under the working directory. They will be used for storing the input/output files.

The source code are stored in "src/" folder.

## Requirment

The project requires the CGAL library. Please download CGAL library at https://www.cgal.org/download.html.

## Dataset

The dataset are stored in "input/" folder. Before you run, please download some of the datasets from https://drive.google.com/drive/folders/1_Sblnq9XOLFwvPI2lX_l8KJJzVWR828y?usp=drive_link, and save in the "input/" folder, since these datasets are very large, and they exceed the maximum file upload size in Github. 

The datasets are as follows:

- "BH_1014" (small version default resolution BH height map, point cloud or TIN dataset with 1014 pixels, points, or vertices)
- "BH_10086" (small version default resolution BH height map, point cloud or TIN dataset with 10086 pixels, points, or vertices)
- "BH_500835" (large version default resolution BH height map, point cloud or TIN dataset with 500835 pixels, points, or vertices)
- "BH_1000414" (large version multiresolution BH height map, point cloud or TIN dataset with 1000414 pixels, points, or vertices)
- "BH_1500996" (large version multiresolution BH height map, point cloud or TIN dataset with 1500996 pixels, points, or vertices)
- "BH_2001610" (large version multiresolution BH height map, point cloud or TIN dataset with 2001610 pixels, points, or vertices)
- "BH_2502596" (large version multiresolution BH height map, point cloud or TIN dataset with 2502596 pixels, points, or vertices)
- "EP_10062" (small version default resolution EP height map, point cloud or TIN dataset with 10062 pixels, points, or vertices)
- "EP_20130" (small version multiresolution EP height map, point cloud or TIN dataset with 20130 pixels, points, or vertices)
- "EP_30098" (small version multiresolution EP height map, point cloud or TIN dataset with 30098 pixels, points, or vertices)
- "EP_40076" (small version multiresolution EP height map, point cloud or TIN dataset with 40076 pixels, points, or vertices)
- "EP_50373" (small version multiresolution EP height map, point cloud or TIN dataset with 50373 pixels, points, or vertices)
- "EP_500384" (large version default resolution EP height map, point cloud or TIN dataset with 500384 pixels, points, or vertices)
- "EP_1001040" (large version multiresolution EP height map, point cloud or TIN dataset with 1001040 pixels, points, or vertices)
- "EP_1501578" (large version multiresolution EP height map, point cloud or TIN dataset with 1501578 pixels, points, or vertices)
- "EP_2001536" (large version multiresolution EP height map, point cloud or TIN dataset with 2001536 pixels, points, or vertices)
- "EP_2500560" (large version multiresolution EP height map, point cloud or TIN dataset with 2500560 pixels, points, or vertices)
- "GF_10092" (small version default resolution GF height map, point cloud or TIN dataset with 10092 pixels, points, or vertices)
- "GF_500208" (large version default resolution GF height map, point cloud or TIN dataset with 500208 pixels, points, or vertices)
- "GF_1000518" (large version multiresolution GF height map, point cloud or TIN dataset with 1000518 pixels, points, or vertices)
- "GF_1501668" (large version multiresolution GF height map, point cloud or TIN dataset with 1501668 pixels, points, or vertices)
- "GF_2000832" (large version multiresolution GF height map, point cloud or TIN dataset with 2000832 pixels, points, or vertices)
- "GF_2502075" (large version multiresolution GF height map, point cloud or TIN dataset with 2502075 pixels, points, or vertices)
- "LM_10092" (small version default resolution LM height map, point cloud or TIN dataset with 10092 pixels, points, or vertices)
- "LM_500208" (large version default resolution LM height map, point cloud or TIN dataset with 500208 pixels, points, or vertices)
- "LM_1000518" (large version multiresolution LM height map, point cloud or TIN dataset with 1000518 pixels, points, or vertices)
- "LM_1501668" (large version multiresolution LM height map, point cloud or TIN dataset with 1501668 pixels, points, or vertices)
- "LM_2000832" (large version multiresolution LM height map, point cloud or TIN dataset with 2000832 pixels, points, or vertices)
- "LM_2502075" (large version multiresolution LM height map, point cloud or TIN dataset with 2502075 pixels, points, or vertices)
- "RM_10092" (small version default resolution RM height map, point cloud or TIN dataset with 10092 pixels, points, or vertices)
- "RM_500208" (large version default resolution RM height map, point cloud or TIN dataset with 500208 pixels, points, or vertices)
- "RM_1000518" (large version multiresolution RM height map, point cloud or TIN dataset with 1000518 pixels, points, or vertices)
- "RM_1501668" (large version multiresolution RM height map, point cloud or TIN dataset with 1501668 pixels, points, or vertices)
- "RM_2000832" (large version multiresolution RM height map, point cloud or TIN dataset with 2000832 pixels, points, or vertices)
- "RM_2502075" (large version multiresolution RM height map, point cloud or TIN dataset with 2502075 pixels, points, or vertices)

## Compile command

```
mkdir build
cd build
cmake ..
make
```

## Run command

```
./main [data_and_index] [epsilon] [run_knn_and_range_query]
```

The meaning for each parameter is as follows:

- [data_and_index]: an index for height map, point cloud and TIN data (a integer from 0 to 34)
- [epsilon]: the epsilon value (0 < epsilon <= 1)
- [run_knn_and_range_query]: whether to run knn and range query (0 means not running, 1 means running)

For the [data_and_index], each index value corresponding to a height map, point cloud and TIN data with different dataset size, their relationships are as follows:

| Index | Height map, point cloud and TIN data | Dataset size |
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

Data Format:

For the height map dataset, we used the .png format in the experiment.

For the point cloud dataset, we used the .xyz format in the experiment. The content of the .xyz file is as follows:

```
points_num

1st_point_x_coord 1st_point_y_coord 1st_point_z_coord

2nd_point_x_coord 2nd_point_y_coord 2nd_point_z_coord

......

last_point_x_coord last_point_y_coord last_point_z_coord
```

For the TIN dataset, we used the .off format in the experiment. The content of the .off file is as follows:

```
OFF

vertices_num faces_num edges_num

1st_vertex_x_coord 1st_vertex_y_coord 1st_vertex_z_coord

2nd_vertex_x_coord 2nd_vertex_y_coord 2nd_vertex_z_coord

......

last_vertex_x_coord last_vertex_y_coord last_vertex_z_coord

1st_face_1st_vertex_ID 1st_face_2nd_vertex_ID 1st_face_3td_vertex_ID

2nd_face_1st_vertex_ID 2nd_face_2nd_vertex_ID 2nd_face_3td_vertex_ID

......

last_face_1st_vertex_ID last_face_2nd_vertex_ID last_face_3td_vertex_ID
```

For TIN-ESSP-Adapt(HM), TIN-ESSP-Adapt(PC), TIN-ESSP, TIN-SNP-Adapt(HM), TIN-SNP-Adapt(PC), TIN-SNP, HM-SP, HM-SP-Adapt(PC) and HM-SP-Adapt(TIN), we distinguish the case that when epsilon=0 (resp. epsilon>0), we use them on the original (resp. simplified) height maps or TINs.

Since TIN-SSimplify-Adapt(HM), TIN-SSimplify-Adapt(PC), TIN-SSimplify, TIN-NSimplify-Adapt(HM), TIN-NSimplify-Adapt(PC), TIN-NSimplify, HM-Simplify-LS, and HM-Simplify-LST are time consuming or having large memory usage, and TIN-ESSP-Adapt(HM) on the simplified TIN, TIN-ESSP-Adapt(PC) on the simplified TIN, TIN-ESSP on the simplified TIN, TIN-SNP-Adapt(HM) on the simplified TIN, TIN-SNP-Adapt(PC) on the simplified TIN, TIN-SNP on the simplified TIN, HM-SP-LS on the simplified height map, and HM-SP-LST on the simplified height map depends on these algorithms, the project will run all algorithms on small-version dataset ([data_and_index] <= 9). The project will run all algorithms except these mentioned algorithms on original dataset ([data_and_index] > 9).

In addition, we strongly encourage you to set [run_knn_and_range_query] to 0 if you are not conducting experiments. Otherwise, it will take a very long time to calculate them. 

An example:

```
./main 0 0.5 0
```

In this example, [data_and_index] is 0, [epsilon] is 0.5, [run_knn_and_range_query] is 0. So, it will run BH height map, point cloud and TIN dataset, with dataset size equal to 10086, epsilon is 0.5, it will not run knn and range query. It will run all algorithms (for TIN-ESSP-Adapt(HM), TIN-ESSP-Adapt(PC), TIN-ESSP, TIN-SNP-Adapt(HM), TIN-SNP-Adapt(PC), TIN-SNP, HM-SP, HM-SP-Adapt(PC) and HM-SP-Adapt(TIN), they are applied on the simplified height maps or TINs).

Another example:

```
./main 0 0 0
```

In this example, [data_and_index] is 0, [epsilon] is 0, [run_knn_and_range_query] is 0. So, it will run BH height map, point cloud and TIN dataset, with dataset size equal to 10086, epsilon is 0.5, it will not run knn and range query. It will only run TIN-ESSP-Adapt(HM), TIN-ESSP-Adapt(PC), TIN-ESSP, TIN-SNP-Adapt(HM), TIN-SNP-Adapt(PC), TIN-SNP, HM-SP, HM-SP-Adapt(PC) and HM-SP-Adapt(TIN) on the original height maps, point clouds or TINs.

## Output

The output will be stored in "output/output.txt" file. The format will be as follows:

```
[dataset] [dataset_size] [epsilon] [height_map_to_point_cloud_or_terrain_time (ms)] [height_map_to_point_cloud_or_terrain_memory_usage (MB)] [preprocessing_time (ms)] [memory_usage (MB)] [output_size (MB)] [query_time (ms)] [distance_error_height_map_or_point_cloud] [distance_error_terrain] [knn_query_time] [knn_error_height_map_or_point_cloud] [knn_error_terrain] [range_query_time] [range_error_height_map_or_point_cloud] [range_error_terrain]
```

These information will also be shown in the terminal. 
