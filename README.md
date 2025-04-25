# Efficient Proximity Queries on Simplified Height Maps

## Overview

This project provides the implementation for efficient proximity queries on simplified height maps. We refer the readers to our paper for more details.

We compared 37 algorithms as follows (the algorithms calculate the path passing on a height map as default):

- TIN-SSimplify-Adapt(HM) (baseline simplification algorithm)
- TIN-SSimplify-Adapt(PC) (baseline simplification algorithm that calculates the path passing on a point cloud)
- TIN-SSimplify (baseline simplification algorithm that calculates the path passing on a TIN)
- TIN-NSimplify-Adapt(HM) (baseline simplification algorithm)
- TIN-NSimplify-Adapt(PC) (baseline simplification algorithm that calculates the path passing on a point cloud)
- TIN-NSimplify (baseline simplification algorithm that calculates the path passing on a TIN)
- PC-Simplify-Adapt(HM) (baseline simplification algorithm)
- PC-Simplify (baseline simplification algorithm that calculates the path passing on a point cloud)
- PC-Simplify-Adapt(TIN) (baseline simplification algorithm that calculates the path passing on a TIN)
- HM-Simplify-LS (variation of our simplification algorithm)
- HM-Simplify-LST (variation of our simplification algorithm)
- HM-Simplify-LQT1 (variation of our simplification algorithm)
- HM-Simplify-LQT2 (variation of our simplification algorithm)
- HM-Simplify-DS (variation of our simplification algorithm)
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
- HM-SP-LS (variation of our shortest path query algorithm)
- HM-SP-LST (variation of our shortest path query algorithm)
- HM-SP-LQT1 (variation of our shortest path query algorithm)
- HM-SP-LQT2 (variation of our shortest path query algorithm)
- HM-SP-DS (variation of our shortest path query algorithm)

Make sure there is a folder called "input/" and a folder called "output/" under the working directory. They will be used for storing the input/output files.

The source code are stored in "src/" folder.

## Requirment

The project requires the CGAL library. Please download CGAL library at https://www.cgal.org/download.html.

## Dataset

The dataset are stored in "input/" folder. Before you run, please download some of the datasets from https://drive.google.com/file/d/13D3sH0Ryy-GEfxwpua3xA7c8Slb3vZK9/view?usp=sharing (6GB before zip, 34GB after zip), and save in the "input/" folder, since these datasets are very large, and they exceed the maximum file upload size in Github. 

The datasets are as follows. We provide the height map, point cloud and TIN datasets, the number in the file name means the number of cells, points or vertices. 

- "BH_1024" (small version default resolution BH dataset)
- "BH_501264" (large version default resolution BH datasets)
- "BH_5004169" (large version multiresolution BH datasets)
- "BH_10004569" (large version multiresolution BH datasets)
- "BH_15000129" (large version multiresolution BH datasets)
- "BH_20007729" (large version multiresolution BH datasets)
- "BH_25000000" (large version multiresolution BH datasets)
- "EP_1024" (small version default resolution EP dataset)
- "EP_10000" (small version multiresolution resolution EP datasets)
- "EP_20164" (small version multiresolution EP datasets)
- "EP_30276" (small version multiresolution EP datasets)
- "EP_40000" (small version multiresolution EP datasets)
- "EP_50176" (small version multiresolution EP datasets)
- "EP_501264" (large version default resolution EP datasets)
- "EP_5004169" (large version multiresolution EP datasets)
- "EP_10004569" (large version multiresolution EP datasets)
- "EP_15000129" (large version multiresolution EP datasets)
- "EP_20007729" (large version multiresolution EP datasets)
- "EP_25000000" (large version multiresolution EP datasets)
- "GF_1024" (small version default resolution GF dataset)
- "GF_501264" (large version default resolution GF datasets)
- "GF_5004169" (large version multiresolution GF datasets)
- "GF_10004569" (large version multiresolution GF datasets)
- "GF_15000129" (large version multiresolution GF datasets)
- "GF_20007729" (large version multiresolution GF datasets)
- "GF_25000000" (large version multiresolution GF datasets)
- "LM_1024" (small version default resolution LM dataset)
- "LM_501264" (large version default resolution LM datasets)
- "LM_5004169" (large version multiresolution LM datasets)
- "LM_10004569" (large version multiresolution LM datasets)
- "LM_15000129" (large version multiresolution LM datasets)
- "LM_20007729" (large version multiresolution LM datasets)
- "LM_25000000" (large version multiresolution LM datasets)
- "RM_1024" (small version default resolution RM dataset)
- "RM_501264" (large version default resolution RM datasets)
- "RM_5004169" (large version multiresolution RM datasets)
- "RM_10004569" (large version multiresolution RM datasets)
- "RM_15000129" (large version multiresolution RM datasets)
- "RM_20007729" (large version multiresolution RM datasets)
- "RM_25000000" (large version multiresolution RM datasets)

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

- [data_and_index]: an index for height map, point cloud and TIN data (a integer from 0 to 39)
- [epsilon]: the epsilon value (0 < epsilon <= 1)
- [run_knn_and_range_query]: whether to run knn and range query (0 means not running, 1 means running)

For the [data_and_index], each index value corresponding to a height map, point cloud and TIN data with different dataset size, their relationships are as follows:

| Index | Height map, point cloud and TIN data | Dataset size |
| ----------- | ----------- | ----------- |
| 0 | BH | 1024 |
| 1 | EP | 1024 |
| 2 | EP | 10000 |
| 3 | EP | 20164 |
| 4 | EP | 30276 |
| 5 | EP | 40000 |
| 6 | EP | 50176 |
| 7 | GF | 1024 |
| 8 | LM | 1024 |
| 9 | RM | 1024 |
| 10 | BH | 501264 |
| 11 | BH | 5004169 |
| 12 | BH | 10004569 |
| 13 | BH | 15000129 |
| 14 | BH | 20007729 |
| 15 | BH | 25000000 |
| 16 | EP | 501264 |
| 17 | EP | 5004169 |
| 18 | EP | 10004569 |
| 19 | EP | 15000129 |
| 20 | EP | 20007729 |
| 21 | EP | 25000000 |
| 22 | GF | 501264 |
| 23 | GF | 5004169 |
| 24 | GF | 10004569 |
| 25 | GF | 15000129 |
| 26 | GF | 20007729 |
| 27 | GF | 25000000 |
| 28 | LM | 501264 |
| 29 | LM | 5004169 |
| 30 | LM | 10004569 |
| 31 | LM | 15000129 |
| 32 | LM | 20007729 |
| 33 | LM | 25000000 |
| 34 | RM | 501264 |
| 35 | RM | 5004169 |
| 36 | RM | 10004569 |
| 37 | RM | 15000129 |
| 38 | RM | 20007729 |
| 39 | RM | 25000000 |

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

For TIN-ESSP-Adapt(HM), TIN-ESSP-Adapt(PC), TIN-ESSP, TIN-SNP-Adapt(HM), TIN-SNP-Adapt(PC), TIN-SNP, PC-SP-Adapt(HM), PC-SP, PC-SP-Adapt(TIN), HM-SP, HM-SP-Adapt(PC) and HM-SP-Adapt(TIN), we distinguish the case that when epsilon=0 (resp. epsilon>0), we use them on the original (resp. simplified) TINs, point clouds or height maps.

Since TIN-SSimplify-Adapt(HM), TIN-SSimplify-Adapt(PC), TIN-SSimplify, TIN-NSimplify-Adapt(HM), TIN-NSimplify-Adapt(PC), TIN-NSimplify, PC-Simplify-Adapt(HM), PC-Simplify, PC-Simplify-Adapt(TIN), HM-Simplify-LS, and HM-Simplify-LST are time consuming or having large memory usage, and TIN-ESSP-Adapt(HM) on the simplified TIN, TIN-ESSP-Adapt(PC) on the simplified TIN, TIN-ESSP on the simplified TIN, TIN-SNP-Adapt(HM) on the simplified TIN, TIN-SNP-Adapt(PC) on the simplified TIN, TIN-SNP on the simplified TIN, PC-SP-Adapt(HM) on the simplified point cloud, PC-SP on the simplified point cloud, PC-SP-Adapt(TIN) on the simplified point cloud, HM-SP-LS on the simplified height map, and HM-SP-LST on the simplified height map depends on these algorithms, the project will run all algorithms on small-version dataset ([data_and_index] <= 9). The project will run all algorithms except these mentioned algorithms on original dataset ([data_and_index] > 9).

In addition, we strongly encourage you to set [run_knn_and_range_query] to 0 if you are not conducting experiments. Otherwise, it will take a very long time to calculate them. 

An example:

```
./main 0 0.5 0
```

In this example, [data_and_index] is 0, [epsilon] is 0.5, [run_knn_and_range_query] is 0. So, it will run BH height map, point cloud and TIN dataset, with dataset size equal to 10086, epsilon is 0.5, it will not run knn and range query. It will run all algorithms (for TIN-ESSP-Adapt(HM), TIN-ESSP-Adapt(PC), TIN-ESSP, TIN-SNP-Adapt(HM), TIN-SNP-Adapt(PC), TIN-SNP, PC-SP-Adapt(HM), PC-SP, PC-SP-Adapt(TIN), HM-SP, HM-SP-Adapt(PC) and HM-SP-Adapt(TIN), they are applied on the simplified TINs, point clouds or height maps).

Another example:

```
./main 0 0 0
```

In this example, [data_and_index] is 0, [epsilon] is 0, [run_knn_and_range_query] is 0. So, it will run BH height map, point cloud and TIN dataset, with dataset size equal to 10086, epsilon is 0.5, it will not run knn and range query. It will only run TIN-ESSP-Adapt(HM), TIN-ESSP-Adapt(PC), TIN-ESSP, TIN-SNP-Adapt(HM), TIN-SNP-Adapt(PC), TIN-SNP, PC-SP-Adapt(HM), PC-SP, PC-SP-Adapt(TIN), HM-SP, HM-SP-Adapt(PC) and HM-SP-Adapt(TIN) on the original simplified TINs, point clouds or height maps.

## Output

The output will be stored in "output/output.txt" file. The format will be as follows:

```
[dataset] [dataset_size] [epsilon] [height_map_to_point_cloud_or_terrain_time (ms)] [height_map_to_point_cloud_or_terrain_memory_usage (MB)] [preprocessing_time (ms)] [memory_usage (MB)] [output_size (MB)] [num_of_cell_point_vertex] [query_time (ms)] [distance_error_height_map_or_point_cloud] [distance_error_terrain] [knn_query_time] [knn_error_height_map_or_point_cloud] [knn_error_terrain] [range_query_time] [range_error_height_map_or_point_cloud] [range_error_terrain]
```

These information will also be shown in the terminal. 
