#include "helper.h"
#include <chrono>

void height_map_to_terrain_and_initialize_terrain(
    height_map_geodesic::HeightMap *height_map, geodesic::Mesh *mesh,
    std::vector<double> &terrain_points, std::vector<unsigned> &terrain_faces,
    double &height_map_or_point_cloud_to_terrain_time, double &height_map_or_point_cloud_to_terrain_memory_usage)
{
    auto start_height_map_or_point_cloud_to_terrain_time = std::chrono::high_resolution_clock::now();

    std::string terrain_write_path = "temp_terrain.off";
    height_map->height_map_to_terrain(height_map_or_point_cloud_to_terrain_memory_usage);

    geodesic::read_mesh_from_file(&terrain_write_path[0], terrain_points, terrain_faces);
    mesh->initialize_mesh_data(terrain_points, terrain_faces);

    auto stop_height_map_or_point_cloud_to_terrain_time = std::chrono::high_resolution_clock::now();
    auto duration_height_map_or_point_cloud_to_terrain_time = std::chrono::duration_cast<std::chrono::microseconds>(
        stop_height_map_or_point_cloud_to_terrain_time - start_height_map_or_point_cloud_to_terrain_time);
    height_map_or_point_cloud_to_terrain_time = duration_height_map_or_point_cloud_to_terrain_time.count();
    height_map_or_point_cloud_to_terrain_time /= 1000;
}

void calculate_height_map_or_point_cloud_exact_distance(
    height_map_geodesic::HeightMap *org_height_map,
    int source_index, int destination_index, double &height_map_or_point_cloud_exact_distance)
{
    height_map_geodesic::HeightMapGeodesicAlgorithmDijkstra algorithm(org_height_map);
    double const distance_limit = height_map_geodesic::INFIN;
    height_map_geodesic::PathPoint source(&org_height_map->hm_points()[source_index]);
    height_map_geodesic::PathPoint destination(&org_height_map->hm_points()[destination_index]);
    std::vector<height_map_geodesic::PathPoint> one_source_list(1, source);
    std::vector<height_map_geodesic::PathPoint> one_destination_list(1, destination);
    algorithm.propagate(one_source_list, distance_limit, 3);
    algorithm.best_source(destination, height_map_or_point_cloud_exact_distance);
    height_map_or_point_cloud_exact_distance = round(height_map_or_point_cloud_exact_distance * 1000000000.0) / 1000000000.0;
}

void calculate_terrain_exact_distance(
    height_map_geodesic::HeightMap *org_height_map,
    int source_index, int destination_index, double &terrain_exact_distance)
{
    geodesic::Mesh org_mesh;
    std::vector<double> org_terrain_vertex;
    std::vector<unsigned> org_terrain_face;
    org_terrain_vertex.clear();
    org_terrain_face.clear();
    double height_map_or_point_cloud_to_terrain_time = 0;
    double height_map_or_point_cloud_to_terrain_memory_usage = 0;

    height_map_to_terrain_and_initialize_terrain(
        org_height_map, &org_mesh,
        org_terrain_vertex, org_terrain_face,
        height_map_or_point_cloud_to_terrain_time, height_map_or_point_cloud_to_terrain_memory_usage);

    geodesic::GeodesicAlgorithmExact algorithm_face_exact(&org_mesh);
    double const distance_limit = geodesic::GEODESIC_INF;
    geodesic::SurfacePoint source(&org_mesh.vertices()[source_index]);
    geodesic::SurfacePoint destination(&org_mesh.vertices()[destination_index]);
    std::vector<geodesic::SurfacePoint> one_source_list(1, source);
    std::vector<geodesic::SurfacePoint> one_destination_list(1, destination);
    algorithm_face_exact.propagate(one_source_list, distance_limit, 3);
    algorithm_face_exact.best_source(destination, terrain_exact_distance);
    terrain_exact_distance = round(terrain_exact_distance * 1000000000.0) / 1000000000.0;
}

void calculate_height_map_or_point_cloud_exact_knn_and_range_query(
    height_map_geodesic::HeightMap *org_height_map,
    int source_index, int knn_and_range_query_obj_num, int k_value, double range,
    std::vector<int> &knn_list, std::vector<int> &range_list)
{
    std::vector<std::pair<double, int>> obj_to_other_obj_distance_and_index_list;
    obj_to_other_obj_distance_and_index_list.clear();

    height_map_geodesic::HeightMapGeodesicAlgorithmDijkstra algorithm(org_height_map);
    double const distance_limit = height_map_geodesic::INFIN;

    bool contain_source_index = false;
    for (int m = org_height_map->hm_points().size() - 1; m > org_height_map->hm_points().size() - 1 - knn_and_range_query_obj_num; m--)
    {
        if (source_index == m)
        {
            contain_source_index = true;
        }
        int destination_index = m;
        if (contain_source_index)
        {
            destination_index--;
        }
        height_map_geodesic::PathPoint source(&org_height_map->hm_points()[source_index]);
        height_map_geodesic::PathPoint destination(&org_height_map->hm_points()[destination_index]);
        std::vector<height_map_geodesic::PathPoint> one_source_list(1, source);
        std::vector<height_map_geodesic::PathPoint> one_destination_list(1, destination);
        algorithm.propagate(one_source_list, distance_limit, 3);
        double height_map_or_point_cloud_exact_distance;
        algorithm.best_source(destination, height_map_or_point_cloud_exact_distance);
        height_map_or_point_cloud_exact_distance = round(height_map_or_point_cloud_exact_distance * 1000000000.0) / 1000000000.0;
        obj_to_other_obj_distance_and_index_list.push_back(std::make_pair(height_map_or_point_cloud_exact_distance, destination_index));
    }
    knn_or_range_query(1, k_value, range, obj_to_other_obj_distance_and_index_list, knn_list);
    knn_or_range_query(2, k_value, range, obj_to_other_obj_distance_and_index_list, range_list);
}

void calculate_terrain_exact_knn_and_range_query(
    height_map_geodesic::HeightMap *org_height_map,
    int source_index, int knn_and_range_query_obj_num, int k_value, double range,
    std::vector<int> &knn_list, std::vector<int> &range_list)
{
    geodesic::Mesh org_mesh;
    std::vector<double> org_terrain_vertex;
    std::vector<unsigned> org_terrain_face;
    org_terrain_vertex.clear();
    org_terrain_face.clear();
    double height_map_or_point_cloud_to_terrain_time = 0;
    double height_map_or_point_cloud_to_terrain_memory_usage = 0;

    height_map_to_terrain_and_initialize_terrain(
        org_height_map, &org_mesh,
        org_terrain_vertex, org_terrain_face,
        height_map_or_point_cloud_to_terrain_time, height_map_or_point_cloud_to_terrain_memory_usage);

    std::vector<std::pair<double, int>> obj_to_other_obj_distance_and_index_list;
    obj_to_other_obj_distance_and_index_list.clear();

    bool contain_source_index = false;
    for (int m = org_height_map->hm_points().size() - 1; m > org_height_map->hm_points().size() - 1 - knn_and_range_query_obj_num; m--)
    {
        if (source_index == m)
        {
            contain_source_index = true;
        }
        int destination_index = m;
        if (contain_source_index)
        {
            destination_index--;
        }
        geodesic::GeodesicAlgorithmExact algorithm_face_exact(&org_mesh);
        double const distance_limit = geodesic::GEODESIC_INF;
        geodesic::SurfacePoint source(&org_mesh.vertices()[source_index]);
        geodesic::SurfacePoint destination(&org_mesh.vertices()[destination_index]);
        std::vector<geodesic::SurfacePoint> one_source_list(1, source);
        std::vector<geodesic::SurfacePoint> one_destination_list(1, destination);
        algorithm_face_exact.propagate(one_source_list, distance_limit, 3);
        double terrain_exact_distance;
        algorithm_face_exact.best_source(destination, terrain_exact_distance);
        terrain_exact_distance = round(terrain_exact_distance * 1000000000.0) / 1000000000.0;
        obj_to_other_obj_distance_and_index_list.push_back(std::make_pair(terrain_exact_distance, destination_index));
    }
    knn_or_range_query(1, k_value, range, obj_to_other_obj_distance_and_index_list, knn_list);
    knn_or_range_query(2, k_value, range, obj_to_other_obj_distance_and_index_list, range_list);
}

void simplified_height_map_or_point_cloud(height_map_geodesic::HeightMap *org_height_map,
                                          double epsilon,
                                          int source_index, int destination_index,
                                          double &simplification_time, double &simplification_memory_usage, double &output_size,
                                          double &query_time, double &query_memory_usage, double &distance_result,
                                          std::vector<height_map_geodesic::PathPoint> &path_result,
                                          int org_one_lqt1_two_lqt2_three_ls_four_lst_five_point_six,
                                          bool run_knn_and_range_query,
                                          double &knn_query_time, double &range_query_time,
                                          int knn_and_range_query_obj_num,
                                          int k_value, double range,
                                          std::vector<int> &knn_list,
                                          std::vector<int> &range_list,
                                          std::vector<int> &height_map_or_point_cloud_exact_knn_list,
                                          std::vector<int> &terrain_exact_knn_list,
                                          std::vector<int> &height_map_or_point_cloud_exact_range_list,
                                          std::vector<int> &terrain_exact_range_list,
                                          double &height_map_or_point_cloud_knn_query_error,
                                          double &terrain_knn_query_error,
                                          double &height_map_or_point_cloud_range_query_error,
                                          double &terrain_range_query_error)
{
    auto start_simplification_time = std::chrono::high_resolution_clock::now();

    height_map_geodesic::HeightMap new_height_map;
    new_height_map.copy_height_map(org_height_map);

    height_map_geodesic::HeightMapGeodesicAlgorithmDijkstra org_algorithm(org_height_map);
    int org_total_node_size;
    size_t org_node_size;
    org_algorithm.get_memory(org_total_node_size, org_node_size);

    // calculate the distance between each two point on the original height map for naive algorithm
    std::unordered_map<int, double> all_dist_p_to_p_org_height_map;
    all_dist_p_to_p_org_height_map.clear();
    if (org_one_lqt1_two_lqt2_three_ls_four_lst_five_point_six == 5)
    {
        height_map_cal_org_dist_p_to_p(org_height_map, all_dist_p_to_p_org_height_map);
    }

    // use the index information, not the hash information, to store added center and
    // deleted points
    std::unordered_map<int, std::unordered_map<int, int>> dominate_table_map;
    dominate_table_map.clear();
    std::unordered_map<int, std::unordered_map<int, int>> temp_dominate_table_map;
    temp_dominate_table_map.clear();

    // use the index information, not the hash information, to store deleted point as key
    // and deminated by point as value in 1D
    std::unordered_map<int, int> del_p_dom_by_map;
    del_p_dom_by_map.clear();
    std::unordered_map<int, int> temp_del_p_dom_by_map;
    temp_del_p_dom_by_map.clear();

    // store the skipped bottom left point that does not satisfy distance error checking, use index
    std::unordered_map<int, int> skipped_bottom_left_p_index_map;
    skipped_bottom_left_p_index_map.clear();

    // store the information of the height map in previous step
    double added_point_prev_x, added_point_prev_y, added_point_prev_z;
    std::unordered_map<int, bool> prev_deleted;
    std::unordered_map<int, std::unordered_map<int, double>> prev_adjacent_hm_points_and_distance;
    prev_deleted.clear();
    prev_adjacent_hm_points_and_distance.clear();

    // looping for merge four points
    int merge_four_count = 0;
    while (true)
    {
        bool can_merge_four_point = false;
        int merged_bottom_left_index;
        int merged_bottom_left_i;
        int merged_bottom_left_j;
        int added_center_point_index;

        temp_dominate_table_map.clear();
        temp_del_p_dom_by_map.clear();
        temp_dominate_table_map = dominate_table_map;
        temp_del_p_dom_by_map = del_p_dom_by_map;

        height_map_merge_four_point(
            org_height_map, &new_height_map, can_merge_four_point,
            merged_bottom_left_index, merged_bottom_left_i, merged_bottom_left_j,
            added_center_point_index, temp_dominate_table_map,
            temp_del_p_dom_by_map, skipped_bottom_left_p_index_map,
            prev_deleted, prev_adjacent_hm_points_and_distance);

        int merged_top_right_i = merged_bottom_left_i + 1;
        int merged_top_right_j = merged_bottom_left_j + 1;

        if (!can_merge_four_point)
        {
            break;
        }

        merge_four_count++;

        // check whether satisfy distance requirement
        bool merge_four_point_distance_satisfy = true;
        bool merge_four_point_distance_satisfy_naive = true;
        epsilon /= org_one_lqt1_two_lqt2_three_ls_four_lst_five_point_six == 4 ? 2 : 1;
        height_map_cal_simp_dist_and_check(
            org_height_map, &new_height_map, epsilon,
            merge_four_point_distance_satisfy, added_center_point_index,
            temp_dominate_table_map, temp_del_p_dom_by_map);
        if (org_one_lqt1_two_lqt2_three_ls_four_lst_five_point_six == 5)
        {
            height_map_cal_simp_dist_and_check_naive(
                org_height_map, &new_height_map,
                all_dist_p_to_p_org_height_map,
                epsilon, merge_four_point_distance_satisfy_naive,
                temp_dominate_table_map, temp_del_p_dom_by_map);
        }

        if (merge_four_point_distance_satisfy)
        {
            dominate_table_map.clear();
            del_p_dom_by_map.clear();
            dominate_table_map = temp_dominate_table_map;
            del_p_dom_by_map = temp_del_p_dom_by_map;
        }
        else
        {
            skipped_bottom_left_p_index_map[merged_bottom_left_index] = merged_bottom_left_index;
            new_height_map.restore_height_map_merge_four_point(prev_deleted, prev_adjacent_hm_points_and_distance);
            continue;
        }

        if (org_one_lqt1_two_lqt2_three_ls_four_lst_five_point_six == 1 || org_one_lqt1_two_lqt2_three_ls_four_lst_five_point_six == 5)
        {
            // looping for expanding into four directions
            bool keep_exp_four_direct = true;
            int exp_four_direct_count = 0;

            while (keep_exp_four_direct)
            {
                bool can_exp_four_direct = false;

                temp_dominate_table_map.clear();
                temp_del_p_dom_by_map.clear();
                temp_dominate_table_map = dominate_table_map;
                temp_del_p_dom_by_map = del_p_dom_by_map;
                int temp_merged_bottom_left_i = merged_bottom_left_i;
                int temp_merged_bottom_left_j = merged_bottom_left_j;
                int temp_merged_top_right_i = merged_top_right_i;
                int temp_merged_top_right_j = merged_top_right_j;

                std::vector<int> exp_four_direct_v_index_list;
                exp_four_direct_v_index_list.clear();

                height_map_exp_four_direct(
                    org_height_map, &new_height_map, can_exp_four_direct,
                    temp_merged_bottom_left_i, temp_merged_bottom_left_j,
                    temp_merged_top_right_i, temp_merged_top_right_j, added_center_point_index,
                    temp_dominate_table_map, temp_del_p_dom_by_map, exp_four_direct_v_index_list,
                    added_point_prev_x, added_point_prev_y, added_point_prev_z,
                    prev_deleted, prev_adjacent_hm_points_and_distance);

                // check whether satisfy distance requirement,
                // only checking when can expand four direction
                bool exp_four_direct_distance_satisfy = true;
                bool exp_four_direct_distance_satisfy_naive = true;
                if (can_exp_four_direct)
                {
                    height_map_cal_simp_dist_and_check(
                        org_height_map, &new_height_map, epsilon,
                        exp_four_direct_distance_satisfy, added_center_point_index,
                        temp_dominate_table_map, temp_del_p_dom_by_map);
                    if (org_one_lqt1_two_lqt2_three_ls_four_lst_five_point_six == 5)
                    {
                        height_map_cal_simp_dist_and_check_naive(
                            org_height_map, &new_height_map,
                            all_dist_p_to_p_org_height_map,
                            epsilon, exp_four_direct_distance_satisfy_naive,
                            temp_dominate_table_map, temp_del_p_dom_by_map);
                    }
                }

                // if both can expand four direction, and the distance satisfy, then we can
                // further expand four direction
                if (can_exp_four_direct && exp_four_direct_distance_satisfy)
                {
                    exp_four_direct_count++;

                    dominate_table_map.clear();
                    del_p_dom_by_map.clear();
                    dominate_table_map = temp_dominate_table_map;
                    del_p_dom_by_map = temp_del_p_dom_by_map;
                    merged_bottom_left_i = temp_merged_bottom_left_i;
                    merged_bottom_left_j = temp_merged_bottom_left_j;
                    merged_top_right_i = temp_merged_top_right_i;
                    merged_top_right_j = temp_merged_top_right_j;
                }
                else
                {
                    if (can_exp_four_direct && !exp_four_direct_distance_satisfy)
                    {
                        new_height_map.restore_height_map_exp_four_three_two_one_direct(
                            added_center_point_index, added_point_prev_x,
                            added_point_prev_y, added_point_prev_z,
                            prev_deleted, prev_adjacent_hm_points_and_distance);
                    }

                    // looping for 14 cases that expands into three, two, one directions
                    std::vector<move_direction> move_direction_list;
                    move_direction_list.clear();

                    move_direction_list.push_back(move_direction(-1, 1, 1, 0));
                    move_direction_list.push_back(move_direction(-1, 1, 0, -1));
                    move_direction_list.push_back(move_direction(-1, 0, 1, -1));
                    move_direction_list.push_back(move_direction(0, 1, 1, -1));

                    move_direction_list.push_back(move_direction(-1, 1, 0, 0));
                    move_direction_list.push_back(move_direction(-1, 0, 1, 0));
                    move_direction_list.push_back(move_direction(0, 1, 1, 0));
                    move_direction_list.push_back(move_direction(-1, 0, 0, -1));
                    move_direction_list.push_back(move_direction(0, 1, 0, -1));
                    move_direction_list.push_back(move_direction(0, 0, 1, -1));

                    move_direction_list.push_back(move_direction(-1, 0, 0, 0));
                    move_direction_list.push_back(move_direction(0, 1, 0, 0));
                    move_direction_list.push_back(move_direction(0, 0, 1, 0));
                    move_direction_list.push_back(move_direction(0, 0, 0, -1));

                    for (int i = 0; i < move_direction_list.size(); i++)
                    {
                        bool can_exp_three_two_one_direct = false;

                        temp_dominate_table_map.clear();
                        temp_del_p_dom_by_map.clear();
                        temp_dominate_table_map = dominate_table_map;
                        temp_del_p_dom_by_map = del_p_dom_by_map;
                        temp_merged_bottom_left_i = merged_bottom_left_i;
                        temp_merged_bottom_left_j = merged_bottom_left_j;
                        temp_merged_top_right_i = merged_top_right_i;
                        temp_merged_top_right_j = merged_top_right_j;

                        std::vector<int> exp_three_two_one_direct_v_index_list;
                        exp_three_two_one_direct_v_index_list.clear();

                        height_map_exp_three_two_one_direct(
                            org_height_map, &new_height_map, can_exp_three_two_one_direct,
                            temp_merged_bottom_left_i, temp_merged_bottom_left_j,
                            temp_merged_top_right_i, temp_merged_top_right_j, added_center_point_index,
                            temp_dominate_table_map, temp_del_p_dom_by_map,
                            exp_three_two_one_direct_v_index_list, move_direction_list[i],
                            added_point_prev_x, added_point_prev_y, added_point_prev_z,
                            prev_deleted, prev_adjacent_hm_points_and_distance);

                        if (!can_exp_three_two_one_direct)
                        {
                            if (i == move_direction_list.size() - 1)
                            {
                                keep_exp_four_direct = false;
                            }
                            continue;
                        }

                        // check whether satisfy distance requirement
                        bool exp_three_two_one_direct_distance_satisfy = true;
                        bool exp_three_two_one_direct_distance_satisfy_naive = true;
                        height_map_cal_simp_dist_and_check(
                            org_height_map, &new_height_map, epsilon,
                            exp_three_two_one_direct_distance_satisfy, added_center_point_index,
                            temp_dominate_table_map, temp_del_p_dom_by_map);
                        if (org_one_lqt1_two_lqt2_three_ls_four_lst_five_point_six == 5)
                        {
                            height_map_cal_simp_dist_and_check_naive(
                                org_height_map, &new_height_map,
                                all_dist_p_to_p_org_height_map,
                                epsilon, exp_three_two_one_direct_distance_satisfy_naive,
                                temp_dominate_table_map, temp_del_p_dom_by_map);
                        }

                        if (exp_three_two_one_direct_distance_satisfy)
                        {
                            dominate_table_map.clear();
                            del_p_dom_by_map.clear();
                            dominate_table_map = temp_dominate_table_map;
                            del_p_dom_by_map = temp_del_p_dom_by_map;
                            merged_bottom_left_i = temp_merged_bottom_left_i;
                            merged_bottom_left_j = temp_merged_bottom_left_j;
                            merged_top_right_i = temp_merged_top_right_i;
                            merged_top_right_j = temp_merged_top_right_j;
                            break;
                        }
                        else
                        {
                            new_height_map.restore_height_map_exp_four_three_two_one_direct(
                                added_center_point_index, added_point_prev_x,
                                added_point_prev_y, added_point_prev_z,
                                prev_deleted, prev_adjacent_hm_points_and_distance);
                        }
                        keep_exp_four_direct = false;
                    }
                }
            }
        }
    }

    auto stop_simplification_time = std::chrono::high_resolution_clock::now();
    auto duration_simplification_time =
        std::chrono::duration_cast<std::chrono::microseconds>(stop_simplification_time - start_simplification_time);
    simplification_time = duration_simplification_time.count();
    simplification_time /= 1000;

    simp_height_map_query(org_height_map, &new_height_map, dominate_table_map,
                          del_p_dom_by_map, source_index, destination_index,
                          query_time, query_memory_usage, distance_result, path_result,
                          org_one_lqt1_two_lqt2_three_ls_four_lst_five_point_six);

    simplification_memory_usage = org_total_node_size * org_node_size + path_result.size() * sizeof(height_map_geodesic::PathPoint) + sizeof(double) + 4 * del_p_dom_by_map.size() * sizeof(double) + (new_height_map.hm_points().size() - del_p_dom_by_map.size()) * 3 * sizeof(double);
    output_size = (new_height_map.hm_points().size() - del_p_dom_by_map.size()) * 3 * sizeof(double);

    if (run_knn_and_range_query)
    {
        simp_height_map_knn_and_range_query(
            org_height_map, &new_height_map, dominate_table_map,
            del_p_dom_by_map, source_index,
            knn_query_time, range_query_time,
            knn_and_range_query_obj_num, k_value, range, knn_list, range_list,
            org_one_lqt1_two_lqt2_three_ls_four_lst_five_point_six);
        calculate_knn_or_range_query_error(height_map_or_point_cloud_exact_knn_list, knn_list, height_map_or_point_cloud_knn_query_error);
        calculate_knn_or_range_query_error(terrain_exact_knn_list, knn_list, terrain_knn_query_error);
        calculate_knn_or_range_query_error(height_map_or_point_cloud_exact_range_list, range_list, height_map_or_point_cloud_range_query_error);
        calculate_knn_or_range_query_error(terrain_exact_range_list, range_list, terrain_range_query_error);
    }
}

void simplified_terrain_face_exact_and_face_appr_and_vertex(
    height_map_geodesic::HeightMap *org_height_map,
    double epsilon, int source_index, int destination_index,
    int face_exact_one_face_appr_two_vertex_three,
    double &height_map_or_point_cloud_to_terrain_time,
    double &height_map_or_point_cloud_to_terrain_memory_usage,
    double &simplification_time, double &simplification_memory_usage, double &output_size,
    double &query_time, double &query_memory_usage, double &distance_result,
    std::vector<geodesic::SurfacePoint> &path_result,
    bool run_knn_and_range_query,
    double &knn_query_time, double &range_query_time,
    int knn_and_range_query_obj_num,
    int k_value, double range,
    std::vector<int> &knn_list,
    std::vector<int> &range_list,
    std::vector<int> &height_map_or_point_cloud_exact_knn_list,
    std::vector<int> &terrain_exact_knn_list,
    std::vector<int> &height_map_or_point_cloud_exact_range_list,
    std::vector<int> &terrain_exact_range_list,
    double &height_map_or_point_cloud_knn_query_error,
    double &terrain_knn_query_error,
    double &height_map_or_point_cloud_range_query_error,
    double &terrain_range_query_error)
{
    auto start_simplification_time = std::chrono::high_resolution_clock::now();

    geodesic::Mesh org_mesh;
    std::vector<double> org_terrain_vertex;
    std::vector<unsigned> org_terrain_face;
    org_terrain_vertex.clear();
    org_terrain_face.clear();

    std::unordered_map<unsigned, std::unordered_map<unsigned, double>> guest_map;
    std::unordered_map<unsigned, std::unordered_map<unsigned, unsigned>> host_map;
    guest_map.clear();
    host_map.clear();

    height_map_to_terrain_and_initialize_terrain(
        org_height_map, &org_mesh,
        org_terrain_vertex, org_terrain_face,
        height_map_or_point_cloud_to_terrain_time, height_map_or_point_cloud_to_terrain_memory_usage);

    double subdivision_level = 0;
    if (face_exact_one_face_appr_two_vertex_three == 2)
    {
        subdivision_level = epslion_to_subdivision_level(epsilon);
    }

    geodesic::GeodesicAlgorithmExact algorithm_face_exact(&org_mesh);
    geodesic::GeodesicAlgorithmSubdivision algorithm_face_appr(&org_mesh, subdivision_level);
    geodesic::GeodesicAlgorithmSubdivision algorithm_vertex(&org_mesh, 0);

    // simplify the terrain
    std::vector<double> pre_terrain_vertex;
    std::vector<unsigned> pre_terrain_face;
    pre_terrain_vertex.clear();
    pre_terrain_face.clear();
    pre_terrain_vertex = org_terrain_vertex;
    pre_terrain_face = org_terrain_face;

    std::vector<double> post_terrain_vertex;
    std::vector<unsigned> post_terrain_face;

    std::unordered_map<unsigned, std::unordered_map<unsigned, double>> temp_guest_map;
    std::unordered_map<unsigned, std::unordered_map<unsigned, unsigned>> temp_host_map;

    // store the skipped vertex that does not satisfy distance error checking
    std::unordered_map<unsigned, unsigned> skipped_v_hash_map;
    skipped_v_hash_map.clear();

    int count = 0;
    while (true)
    {
        geodesic::Mesh pre_mesh;
        pre_mesh.initialize_mesh_data(pre_terrain_vertex, pre_terrain_face);
        assert(org_mesh.m_x_vertex_num == pre_mesh.m_x_vertex_num &&
               org_mesh.m_y_vertex_num == pre_mesh.m_y_vertex_num &&
               org_mesh.m_xmin == pre_mesh.m_xmin &&
               org_mesh.m_ymin == pre_mesh.m_ymin);

        bool can_simplify;
        unsigned del_vertex_x_y_hash;
        std::vector<std::pair<unsigned, double>> del_v_neig_hash_and_dist_list;
        del_v_neig_hash_and_dist_list.clear();

        post_terrain_vertex.clear();
        post_terrain_face.clear();

        auto start_time2 = std::chrono::high_resolution_clock::now();

        terrain_simplify_one_iteration(
            &org_mesh, &pre_mesh, epsilon, post_terrain_vertex, post_terrain_face, can_simplify,
            del_vertex_x_y_hash, del_v_neig_hash_and_dist_list, skipped_v_hash_map,
            subdivision_level, face_exact_one_face_appr_two_vertex_three);

        if (!can_simplify)
        {
            break;
        }

        temp_guest_map.clear();
        temp_host_map.clear();
        temp_guest_map = guest_map;
        temp_host_map = host_map;

        update_guest_host_map(temp_guest_map, temp_host_map, del_vertex_x_y_hash, del_v_neig_hash_and_dist_list);
        for (auto ite : temp_guest_map)
        {
            unsigned a = ite.first;
            unsigned x, y;
            hash_function_one_key_to_two_keys(pre_mesh.m_hash_max, a, x, y);
            for (auto ite2 : temp_guest_map[a])
            {
                unsigned b = ite2.first;
                double c = ite2.second;
                unsigned x1, y1;
                hash_function_one_key_to_two_keys(pre_mesh.m_hash_max, b, x1, y1);
            }
        }
        for (auto ite : temp_host_map)
        {
            unsigned a = ite.first;
            unsigned x, y;
            hash_function_one_key_to_two_keys(pre_mesh.m_hash_max, a, x, y);
            for (auto ite2 : temp_host_map[a])
            {
                unsigned b = ite2.first;
                unsigned c = ite2.second;
                unsigned x1, y1;
                hash_function_one_key_to_two_keys(pre_mesh.m_hash_max, b, x1, y1);
            }
        }

        geodesic::Mesh post_mesh;
        post_mesh.initialize_mesh_data(post_terrain_vertex, post_terrain_face);
        assert(org_mesh.m_x_vertex_num == post_mesh.m_x_vertex_num &&
               org_mesh.m_y_vertex_num == post_mesh.m_y_vertex_num &&
               org_mesh.m_xmin == post_mesh.m_xmin &&
               org_mesh.m_ymin == post_mesh.m_ymin);

        bool satisfy;
        terrain_face_exact_and_face_appr_and_vertex_cal_simp_dist_and_check(
            &org_mesh, &post_mesh, epsilon, subdivision_level,
            face_exact_one_face_appr_two_vertex_three, satisfy,
            temp_guest_map, temp_host_map,
            del_v_neig_hash_and_dist_list);

        if (satisfy)
        {
            pre_terrain_vertex.clear();
            pre_terrain_face.clear();
            pre_terrain_vertex = post_terrain_vertex;
            pre_terrain_face = post_terrain_face;

            guest_map = temp_guest_map;
            host_map = temp_host_map;
        }
        else
        {
            skipped_v_hash_map[del_vertex_x_y_hash] = del_vertex_x_y_hash;
        }
        count++;
    }

    geodesic::Mesh post_mesh;
    post_mesh.initialize_mesh_data(pre_terrain_vertex, pre_terrain_face);
    assert(org_mesh.m_x_vertex_num == post_mesh.m_x_vertex_num &&
           org_mesh.m_y_vertex_num == post_mesh.m_y_vertex_num &&
           org_mesh.m_xmin == post_mesh.m_xmin &&
           org_mesh.m_ymin == post_mesh.m_ymin);

    auto stop_simplification_time = std::chrono::high_resolution_clock::now();
    auto duration_simplification_time =
        std::chrono::duration_cast<std::chrono::microseconds>(stop_simplification_time - start_simplification_time);
    simplification_time = duration_simplification_time.count();
    simplification_time /= 1000;

    write_terrain_off_file(pre_terrain_vertex, pre_terrain_face);

    simp_terrain_face_exact_face_appr_vertex_query(
        &org_mesh, &post_mesh, guest_map, host_map, subdivision_level,
        face_exact_one_face_appr_two_vertex_three,
        source_index, destination_index, query_time,
        query_memory_usage, distance_result, path_result);

    double map_size = 0;
    for (auto ite : guest_map)
    {
        map_size += sizeof(unsigned);
        for (auto ite2 : ite.second)
        {
            map_size += sizeof(unsigned) + sizeof(double);
        }
    }
    for (auto ite : host_map)
    {
        map_size += sizeof(unsigned);
        for (auto ite2 : ite.second)
        {
            map_size += sizeof(unsigned) + sizeof(unsigned);
        }
    }
    if (face_exact_one_face_appr_two_vertex_three == 1)
    {
        simplification_memory_usage = algorithm_face_exact.get_memory();
    }
    else if (face_exact_one_face_appr_two_vertex_three == 2)
    {
        simplification_memory_usage = algorithm_face_appr.get_memory();
    }
    else if (face_exact_one_face_appr_two_vertex_three == 3)
    {
        simplification_memory_usage = algorithm_vertex.get_memory();
    }

    simplification_memory_usage += path_result.size() * sizeof(geodesic::SurfacePoint) + sizeof(double) + map_size + post_mesh.vertices().size() * 3 * sizeof(double) + post_mesh.faces().size() * 3 * sizeof(int);
    output_size = post_mesh.vertices().size() * 3 * sizeof(double) + post_mesh.faces().size() * 3 * sizeof(int);

    if (run_knn_and_range_query)
    {
        simp_terrain_face_exact_face_appr_vertex_knn_and_range_query(
            &org_mesh, &post_mesh, guest_map, host_map, subdivision_level,
            face_exact_one_face_appr_two_vertex_three,
            source_index,
            knn_query_time, range_query_time,
            knn_and_range_query_obj_num, k_value, range, knn_list, range_list);
        calculate_knn_or_range_query_error(height_map_or_point_cloud_exact_knn_list, knn_list, height_map_or_point_cloud_knn_query_error);
        calculate_knn_or_range_query_error(terrain_exact_knn_list, knn_list, terrain_knn_query_error);
        calculate_knn_or_range_query_error(height_map_or_point_cloud_exact_range_list, range_list, height_map_or_point_cloud_range_query_error);
        calculate_knn_or_range_query_error(terrain_exact_range_list, range_list, terrain_range_query_error);
    }
}

void height_map_or_point_cloud(height_map_geodesic::HeightMap *org_height_map,
                               int source_index, int destination_index,
                               double &query_time, double &query_memory_usage, double &distance_result,
                               std::vector<height_map_geodesic::PathPoint> &path_result)
{
    auto start_query_time = std::chrono::high_resolution_clock::now();

    height_map_geodesic::HeightMapGeodesicAlgorithmDijkstra algorithm(org_height_map);
    double const distance_limit = height_map_geodesic::INFIN;
    height_map_geodesic::PathPoint source(&org_height_map->hm_points()[source_index]);
    height_map_geodesic::PathPoint destination(&org_height_map->hm_points()[destination_index]);
    std::vector<height_map_geodesic::PathPoint> one_source_list(1, source);
    std::vector<height_map_geodesic::PathPoint> one_destination_list(1, destination);
    algorithm.propagate(one_source_list, distance_limit, 3);
    algorithm.trace_back(destination, path_result);
    distance_result = length(path_result);
    distance_result = round(distance_result * 1000000000.0) / 1000000000.0;
    int total_node_size;
    size_t node_size;
    algorithm.get_memory(total_node_size, node_size);
    query_memory_usage += total_node_size * node_size + path_result.size() * sizeof(height_map_geodesic::PathPoint) + sizeof(double);

    auto stop_query_time = std::chrono::high_resolution_clock::now();
    auto duration_query_time = std::chrono::duration_cast<std::chrono::microseconds>(stop_query_time - start_query_time);
    query_time = duration_query_time.count();
    query_time /= 1000;
}

void height_map_or_point_cloud_knn_and_range_query(
    height_map_geodesic::HeightMap *org_height_map,
    int source_index,
    double &knn_query_time, double &range_query_time,
    int knn_and_range_query_obj_num,
    int k_value, double range,
    std::vector<int> &knn_list,
    std::vector<int> &range_list,
    std::vector<int> &height_map_or_point_cloud_exact_knn_list,
    std::vector<int> &terrain_exact_knn_list,
    std::vector<int> &height_map_or_point_cloud_exact_range_list,
    std::vector<int> &terrain_exact_range_list,
    double &height_map_or_point_cloud_knn_query_error,
    double &terrain_knn_query_error,
    double &height_map_or_point_cloud_range_query_error,
    double &terrain_range_query_error)
{
    auto start_knn_or_range_query_time = std::chrono::high_resolution_clock::now();

    std::vector<std::pair<double, int>> obj_to_other_obj_distance_and_index_list;
    obj_to_other_obj_distance_and_index_list.clear();

    bool contain_source_index = false;
    for (int m = org_height_map->hm_points().size() - 1; m > org_height_map->hm_points().size() - 1 - knn_and_range_query_obj_num; m--)
    {
        if (source_index == m)
        {
            contain_source_index = true;
        }
        int destination_index = m;
        if (contain_source_index)
        {
            destination_index--;
        }

        height_map_geodesic::HeightMapGeodesicAlgorithmDijkstra algorithm(org_height_map);
        double const distance_limit = height_map_geodesic::INFIN;
        double distance_result;
        std::vector<height_map_geodesic::PathPoint> path_result;
        path_result.clear();
        height_map_geodesic::PathPoint source(&org_height_map->hm_points()[source_index]);
        height_map_geodesic::PathPoint destination(&org_height_map->hm_points()[destination_index]);
        std::vector<height_map_geodesic::PathPoint> one_source_list(1, source);
        std::vector<height_map_geodesic::PathPoint> one_destination_list(1, destination);
        algorithm.propagate(one_source_list, distance_limit, 3);
        algorithm.trace_back(destination, path_result);
        distance_result = length(path_result);
        distance_result = round(distance_result * 1000000000.0) / 1000000000.0;
        obj_to_other_obj_distance_and_index_list.push_back(std::make_pair(distance_result, destination_index));
    }
    std::sort(obj_to_other_obj_distance_and_index_list.begin(), obj_to_other_obj_distance_and_index_list.end());

    auto stop_knn_or_range_query_time = std::chrono::high_resolution_clock::now();
    auto duration_knn_or_range_query_time = std::chrono::duration_cast<std::chrono::microseconds>(stop_knn_or_range_query_time - start_knn_or_range_query_time);

    auto start_knn_query_time = std::chrono::high_resolution_clock::now();
    knn_or_range_query(1, k_value, range, obj_to_other_obj_distance_and_index_list, knn_list);
    auto stop_knn_query_time = std::chrono::high_resolution_clock::now();
    auto duration_knn_query_time = std::chrono::duration_cast<std::chrono::microseconds>(stop_knn_query_time - start_knn_query_time);
    knn_query_time = duration_knn_or_range_query_time.count() + duration_knn_query_time.count();
    knn_query_time /= 1000;

    auto start_range_query_time = std::chrono::high_resolution_clock::now();
    knn_or_range_query(2, k_value, range, obj_to_other_obj_distance_and_index_list, range_list);
    auto stop_range_query_time = std::chrono::high_resolution_clock::now();
    auto duration_range_query_time = std::chrono::duration_cast<std::chrono::microseconds>(stop_range_query_time - start_range_query_time);
    range_query_time = duration_knn_or_range_query_time.count() + duration_range_query_time.count();
    range_query_time /= 1000;

    calculate_knn_or_range_query_error(height_map_or_point_cloud_exact_knn_list, knn_list, height_map_or_point_cloud_knn_query_error);
    calculate_knn_or_range_query_error(terrain_exact_knn_list, knn_list, terrain_knn_query_error);
    calculate_knn_or_range_query_error(height_map_or_point_cloud_exact_range_list, range_list, height_map_or_point_cloud_range_query_error);
    calculate_knn_or_range_query_error(terrain_exact_range_list, range_list, terrain_range_query_error);
}

void terrain_face_exact_and_face_appr_and_vertex(
    height_map_geodesic::HeightMap *org_height_map,
    double epsilon, int source_index, int destination_index,
    int face_exact_one_face_appr_two_vertex_three,
    double &height_map_or_point_cloud_to_terrain_time,
    double &height_map_or_point_cloud_to_terrain_memory_usage,
    double &query_time, double &query_memory_usage, double &distance_result,
    std::vector<geodesic::SurfacePoint> &path_result)
{
    geodesic::Mesh org_mesh;
    std::vector<double> org_terrain_vertex;
    std::vector<unsigned> org_terrain_face;
    org_terrain_vertex.clear();
    org_terrain_face.clear();

    height_map_to_terrain_and_initialize_terrain(
        org_height_map, &org_mesh,
        org_terrain_vertex, org_terrain_face,
        height_map_or_point_cloud_to_terrain_time, height_map_or_point_cloud_to_terrain_memory_usage);

    auto start_query_time = std::chrono::high_resolution_clock::now();

    double subdivision_level = 0;
    if (face_exact_one_face_appr_two_vertex_three == 2)
    {
        subdivision_level = epslion_to_subdivision_level(epsilon);
    }

    geodesic::GeodesicAlgorithmExact algorithm_face_exact(&org_mesh);
    geodesic::GeodesicAlgorithmSubdivision algorithm_face_appr(&org_mesh, subdivision_level);
    geodesic::GeodesicAlgorithmSubdivision algorithm_vertex(&org_mesh, 0);

    double const distance_limit = geodesic::GEODESIC_INF;
    geodesic::SurfacePoint source(&org_mesh.vertices()[source_index]);
    geodesic::SurfacePoint destination(&org_mesh.vertices()[destination_index]);
    std::vector<geodesic::SurfacePoint> one_source_list(1, source);
    std::vector<geodesic::SurfacePoint> one_destination_list(1, destination);
    if (face_exact_one_face_appr_two_vertex_three == 1)
    {
        algorithm_face_exact.propagate(one_source_list, distance_limit, 3);
        algorithm_face_exact.trace_back(destination, path_result);
    }
    else if (face_exact_one_face_appr_two_vertex_three == 2)
    {
        algorithm_face_appr.propagate(one_source_list, distance_limit, 3);
        algorithm_face_appr.trace_back(destination, path_result);
    }
    else if (face_exact_one_face_appr_two_vertex_three == 3)
    {
        algorithm_vertex.propagate(one_source_list, distance_limit, 3);
        algorithm_vertex.trace_back(destination, path_result);
        modify_path(path_result);
    }

    distance_result = length(path_result);
    distance_result = round(distance_result * 1000000000.0) / 1000000000.0;
    if (face_exact_one_face_appr_two_vertex_three == 1)
    {
        query_memory_usage += algorithm_face_exact.get_memory() + path_result.size() * sizeof(geodesic::SurfacePoint) + sizeof(double);
    }
    else if (face_exact_one_face_appr_two_vertex_three == 2)
    {
        query_memory_usage += algorithm_face_appr.get_memory() + path_result.size() * sizeof(geodesic::SurfacePoint) + sizeof(double);
    }
    else if (face_exact_one_face_appr_two_vertex_three == 3)
    {
        query_memory_usage += algorithm_vertex.get_memory() + path_result.size() * sizeof(geodesic::SurfacePoint) + sizeof(double);
    }

    auto stop_query_time = std::chrono::high_resolution_clock::now();
    auto duration_query_time = std::chrono::duration_cast<std::chrono::microseconds>(stop_query_time - start_query_time);
    query_time = duration_query_time.count();
    query_time /= 1000;
}

void terrain_face_exact_and_face_appr_and_vertex_knn_and_range_query(
    height_map_geodesic::HeightMap *org_height_map,
    double epsilon, int source_index,
    int face_exact_one_face_appr_two_vertex_three,
    bool run_knn_and_range_query,
    double &knn_query_time, double &range_query_time,
    int knn_and_range_query_obj_num,
    int k_value, double range,
    std::vector<int> &knn_list,
    std::vector<int> &range_list,
    std::vector<int> &height_map_or_point_cloud_exact_knn_list,
    std::vector<int> &terrain_exact_knn_list,
    std::vector<int> &height_map_or_point_cloud_exact_range_list,
    std::vector<int> &terrain_exact_range_list,
    double &height_map_or_point_cloud_knn_query_error,
    double &terrain_knn_query_error,
    double &height_map_or_point_cloud_range_query_error,
    double &terrain_range_query_error)
{
    geodesic::Mesh org_mesh;
    std::vector<double> org_terrain_vertex;
    std::vector<unsigned> org_terrain_face;
    org_terrain_vertex.clear();
    org_terrain_face.clear();

    double height_map_or_point_cloud_to_terrain_time;
    double height_map_or_point_cloud_to_terrain_memory_usage;

    height_map_to_terrain_and_initialize_terrain(
        org_height_map, &org_mesh,
        org_terrain_vertex, org_terrain_face,
        height_map_or_point_cloud_to_terrain_time, height_map_or_point_cloud_to_terrain_memory_usage);

    auto start_knn_or_range_query_time = std::chrono::high_resolution_clock::now();

    double subdivision_level = 0;
    if (face_exact_one_face_appr_two_vertex_three == 2)
    {
        subdivision_level = epslion_to_subdivision_level(epsilon);
    }

    std::vector<std::pair<double, int>> obj_to_other_obj_distance_and_index_list;
    obj_to_other_obj_distance_and_index_list.clear();

    bool contain_source_index = false;
    for (int m = org_height_map->hm_points().size() - 1; m > org_height_map->hm_points().size() - 1 - knn_and_range_query_obj_num; m--)
    {
        if (source_index == m)
        {
            contain_source_index = true;
        }
        int destination_index = m;
        if (contain_source_index)
        {
            destination_index--;
        }
        geodesic::GeodesicAlgorithmExact algorithm_face_exact(&org_mesh);
        geodesic::GeodesicAlgorithmSubdivision algorithm_face_appr(&org_mesh, subdivision_level);
        geodesic::GeodesicAlgorithmSubdivision algorithm_vertex(&org_mesh, 0);

        double const distance_limit = geodesic::GEODESIC_INF;
        double distance_result;
        std::vector<geodesic::SurfacePoint> path_result;
        path_result.clear();
        geodesic::SurfacePoint source(&org_mesh.vertices()[source_index]);
        geodesic::SurfacePoint destination(&org_mesh.vertices()[destination_index]);
        std::vector<geodesic::SurfacePoint> one_source_list(1, source);
        std::vector<geodesic::SurfacePoint> one_destination_list(1, destination);
        if (face_exact_one_face_appr_two_vertex_three == 1)
        {
            algorithm_face_exact.propagate(one_source_list, distance_limit, 3);
            algorithm_face_exact.trace_back(destination, path_result);
        }
        else if (face_exact_one_face_appr_two_vertex_three == 2)
        {
            algorithm_face_appr.propagate(one_source_list, distance_limit, 3);
            algorithm_face_appr.trace_back(destination, path_result);
        }
        else if (face_exact_one_face_appr_two_vertex_three == 3)
        {
            algorithm_vertex.propagate(one_source_list, distance_limit, 3);
            algorithm_vertex.trace_back(destination, path_result);
            modify_path(path_result);
        }

        distance_result = length(path_result);
        distance_result = round(distance_result * 1000000000.0) / 1000000000.0;
        obj_to_other_obj_distance_and_index_list.push_back(std::make_pair(distance_result, destination_index));
    }

    std::sort(obj_to_other_obj_distance_and_index_list.begin(), obj_to_other_obj_distance_and_index_list.end());

    auto stop_knn_or_range_query_time = std::chrono::high_resolution_clock::now();
    auto duration_knn_or_range_query_time = std::chrono::duration_cast<std::chrono::microseconds>(stop_knn_or_range_query_time - start_knn_or_range_query_time);

    auto start_knn_query_time = std::chrono::high_resolution_clock::now();
    knn_or_range_query(1, k_value, range, obj_to_other_obj_distance_and_index_list, knn_list);
    auto stop_knn_query_time = std::chrono::high_resolution_clock::now();
    auto duration_knn_query_time = std::chrono::duration_cast<std::chrono::microseconds>(stop_knn_query_time - start_knn_query_time);
    knn_query_time = duration_knn_or_range_query_time.count() + duration_knn_query_time.count();
    knn_query_time /= 1000;

    auto start_range_query_time = std::chrono::high_resolution_clock::now();
    knn_or_range_query(2, k_value, range, obj_to_other_obj_distance_and_index_list, range_list);
    auto stop_range_query_time = std::chrono::high_resolution_clock::now();
    auto duration_range_query_time = std::chrono::duration_cast<std::chrono::microseconds>(stop_range_query_time - start_range_query_time);
    range_query_time = duration_knn_or_range_query_time.count() + duration_range_query_time.count();
    range_query_time /= 1000;

    calculate_knn_or_range_query_error(height_map_or_point_cloud_exact_knn_list, knn_list, height_map_or_point_cloud_knn_query_error);
    calculate_knn_or_range_query_error(terrain_exact_knn_list, knn_list, terrain_knn_query_error);
    calculate_knn_or_range_query_error(height_map_or_point_cloud_exact_range_list, range_list, height_map_or_point_cloud_range_query_error);
    calculate_knn_or_range_query_error(terrain_exact_range_list, range_list, terrain_range_query_error);
}

void simplified_height_map_or_point_cloud_with_output(
    std::string output_file,
    height_map_geodesic::HeightMap *org_height_map,
    double epsilon, int source_index, int destination_index,
    int org_one_lqt1_two_lqt2_three_ls_four_lst_five_point_six,
    double height_map_or_point_cloud_exact_distance,
    double terrain_exact_distance, bool run_knn_and_range_query,
    int knn_and_range_query_obj_num, int k_value, double range,
    std::vector<int> &height_map_or_point_cloud_exact_knn_list,
    std::vector<int> &terrain_exact_knn_list,
    std::vector<int> &height_map_or_point_cloud_exact_range_list,
    std::vector<int> &terrain_exact_range_list,
    std::string write_file_header,
    int height_map_one_point_cloud_two_terrain_three)
{
    double simplification_time = 0;
    double simplification_memory_usage = 0;
    double output_size = 0;
    double query_time = 0;
    double query_memory_usage = 0;
    double distance_result = 0;
    double knn_query_time = 0;
    double range_query_time = 0;
    double height_map_or_point_cloud_knn_query_error = 0;
    double terrain_knn_query_error = 0;
    double height_map_or_point_cloud_range_query_error = 0;
    double terrain_range_query_error = 0;
    std::vector<height_map_geodesic::PathPoint> path_result;
    std::vector<int> knn_list;
    std::vector<int> range_list;
    path_result.clear();
    knn_list.clear();
    range_list.clear();
    simplified_height_map_or_point_cloud(
        org_height_map, epsilon, source_index, destination_index,
        simplification_time, simplification_memory_usage, output_size,
        query_time, query_memory_usage, distance_result, path_result,
        org_one_lqt1_two_lqt2_three_ls_four_lst_five_point_six, run_knn_and_range_query,
        knn_query_time, range_query_time, knn_and_range_query_obj_num,
        k_value, range, knn_list, range_list,
        height_map_or_point_cloud_exact_knn_list, terrain_exact_knn_list,
        height_map_or_point_cloud_exact_range_list, terrain_exact_range_list,
        height_map_or_point_cloud_knn_query_error, terrain_knn_query_error,
        height_map_or_point_cloud_range_query_error, terrain_range_query_error);

    std::cout << "Simplification time: " << simplification_time << " ms" << std::endl;
    std::cout << "Simplification memory usage: " << simplification_memory_usage / 1e6 << " MB" << std::endl;
    std::cout << "Output size: " << output_size / 1e6 << " MB" << std::endl;
    std::cout << "Query time: " << query_time << " ms" << std::endl;
    std::cout << "Calculated distance: " << distance_result
              << ", height map or point cloud exact distance: " << height_map_or_point_cloud_exact_distance
              << ", height map or point cloud distance error: " << std::abs(distance_result / height_map_or_point_cloud_exact_distance - 1)
              << ", terrain exact distance: " << terrain_exact_distance
              << ", terrain distance error: " << std::abs(distance_result / terrain_exact_distance - 1) << std::endl;
    if (run_knn_and_range_query)
    {
        std::cout << "Knn query time: " << knn_query_time << " ms" << std::endl;
        std::cout << "Height map or point cloud knn error: " << height_map_or_point_cloud_knn_query_error << ", terrain knn error: " << terrain_knn_query_error << std::endl;
        std::cout << "Range query time: " << range_query_time << " ms" << std::endl;
        std::cout << "Height map or point cloud range error: " << height_map_or_point_cloud_range_query_error << ", terrain range error: " << terrain_range_query_error << std::endl;
    }

    std::ofstream ofs(output_file, std::ios_base::app);
    if (output_file == "../output/output.txt")
    {
        if (org_one_lqt1_two_lqt2_three_ls_four_lst_five_point_six == 1)
        {
            if (height_map_one_point_cloud_two_terrain_three == 1)
            {
                ofs << "\n== Mem_SimQue ==\n";
            }
            else if (height_map_one_point_cloud_two_terrain_three == 2)
            {
                ofs << "\n== Mem_SimQue_AdpC ==\n";
            }
            else if (height_map_one_point_cloud_two_terrain_three == 3)
            {
                ofs << "\n== Mem_SimQue_AdpT ==\n";
            }
        }
        else if (org_one_lqt1_two_lqt2_three_ls_four_lst_five_point_six == 2)
        {
            ofs << "\n== Mem_SimQue_LQT1 ==\n";
        }
        else if (org_one_lqt1_two_lqt2_three_ls_four_lst_five_point_six == 3)
        {
            ofs << "\n== Mem_SimQue_LQT2 ==\n";
        }
        else if (org_one_lqt1_two_lqt2_three_ls_four_lst_five_point_six == 4)
        {
            ofs << "\n== Mem_SimQue_LS ==\n";
        }
        else if (org_one_lqt1_two_lqt2_three_ls_four_lst_five_point_six == 5)
        {
            ofs << "\n== Mem_SimQue_LST ==\n";
        }
        else if (org_one_lqt1_two_lqt2_three_ls_four_lst_five_point_six == 6)
        {
            if (height_map_one_point_cloud_two_terrain_three == 1)
            {
                ofs << "\n== Mes_SimQue_AdpM ==\n";
            }
            else if (height_map_one_point_cloud_two_terrain_three == 2)
            {
                ofs << "\n== Mes_SimQue ==\n";
            }
            else if (height_map_one_point_cloud_two_terrain_three == 3)
            {
                ofs << "\n== Mes_SimQue_AdpT ==\n";
            }
        }
    }
    ofs << write_file_header << "\t"
        << 0 << "\t"
        << 0 / 1e6 << "\t"
        << simplification_time << "\t"
        << simplification_memory_usage / 1e6 << "\t"
        << output_size / 1e6 << "\t"
        << query_time << "\t"
        << std::abs(distance_result / height_map_or_point_cloud_exact_distance - 1) << "\t"
        << std::abs(distance_result / terrain_exact_distance - 1) << "\t"
        << knn_query_time << "\t"
        << height_map_or_point_cloud_knn_query_error << "\t"
        << terrain_knn_query_error << "\t"
        << range_query_time << "\t"
        << height_map_or_point_cloud_range_query_error << "\t"
        << terrain_range_query_error << "\n";
    ofs.close();
}

void simplified_terrain_face_exact_and_face_appr_and_vertex_with_output(
    std::string output_file,
    height_map_geodesic::HeightMap *org_height_map,
    double epsilon, int source_index, int destination_index,
    int face_exact_one_face_appr_two_vertex_three,
    double height_map_or_point_cloud_exact_distance, double terrain_exact_distance,
    bool run_knn_and_range_query, int knn_and_range_query_obj_num,
    int k_value, double range,
    std::vector<int> &height_map_or_point_cloud_exact_knn_list,
    std::vector<int> &terrain_exact_knn_list,
    std::vector<int> &height_map_or_point_cloud_exact_range_list,
    std::vector<int> &terrain_exact_range_list,
    std::string write_file_header,
    int height_map_one_point_cloud_two_terrain_three)
{
    double height_map_or_point_cloud_to_terrain_time = 0;
    double height_map_or_point_cloud_to_terrain_memory_usage = 0;
    double simplification_time = 0;
    double simplification_memory_usage = 0;
    double output_size = 0;
    double query_time = 0;
    double query_memory_usage = 0;
    double distance_result = 0;
    double knn_query_time = 0;
    double range_query_time = 0;
    double height_map_or_point_cloud_knn_query_error = 0;
    double terrain_knn_query_error = 0;
    double height_map_or_point_cloud_range_query_error = 0;
    double terrain_range_query_error = 0;
    std::vector<geodesic::SurfacePoint> path_result;
    std::vector<int> knn_list;
    std::vector<int> range_list;
    path_result.clear();
    knn_list.clear();
    range_list.clear();
    simplified_terrain_face_exact_and_face_appr_and_vertex(
        org_height_map, epsilon, source_index, destination_index,
        face_exact_one_face_appr_two_vertex_three,
        height_map_or_point_cloud_to_terrain_time, height_map_or_point_cloud_to_terrain_memory_usage,
        simplification_time, simplification_memory_usage, output_size,
        query_time, query_memory_usage, distance_result, path_result,
        run_knn_and_range_query, knn_query_time, range_query_time,
        knn_and_range_query_obj_num, k_value, range,
        knn_list, range_list,
        height_map_or_point_cloud_exact_knn_list, terrain_exact_knn_list,
        height_map_or_point_cloud_exact_range_list, terrain_exact_range_list,
        height_map_or_point_cloud_knn_query_error, terrain_knn_query_error,
        height_map_or_point_cloud_range_query_error, terrain_range_query_error);

    if (height_map_one_point_cloud_two_terrain_three == 3)
    {
        height_map_or_point_cloud_to_terrain_time = 0;
        height_map_or_point_cloud_to_terrain_memory_usage = 0;
        std::cout << "Height map or point cloud to terrain time: " << height_map_or_point_cloud_to_terrain_time << " ms" << std::endl;
        std::cout << "Height map or point cloud to terrain memory usage: " << height_map_or_point_cloud_to_terrain_memory_usage / 1e6 << " MB" << std::endl;
    }
    std::cout << "Simplification time: " << simplification_time << " ms" << std::endl;
    std::cout << "Simplification memory usage: " << simplification_memory_usage / 1e6 << " MB" << std::endl;
    std::cout << "Output size: " << output_size / 1e6 << " MB" << std::endl;
    std::cout << "Query time: " << query_time << " ms" << std::endl;
    std::cout << "Calculated distance: " << distance_result
              << ", height map or point cloud exact distance: " << height_map_or_point_cloud_exact_distance
              << ", height map or point cloud distance error: " << std::abs(distance_result / height_map_or_point_cloud_exact_distance - 1)
              << ", terrain exact distance: " << terrain_exact_distance
              << ", terrain distance error: " << std::abs(distance_result / terrain_exact_distance - 1) << std::endl;
    if (run_knn_and_range_query)
    {
        std::cout << "Knn query time: " << knn_query_time << " ms" << std::endl;
        std::cout << "Height map or point cloud knn error: " << height_map_or_point_cloud_knn_query_error << ", terrain knn error: " << terrain_knn_query_error << std::endl;
        std::cout << "Range query time: " << range_query_time << " ms" << std::endl;
        std::cout << "Height map or point cloud range error: " << height_map_or_point_cloud_range_query_error << ", terrain range error: " << terrain_range_query_error << std::endl;
    }

    std::ofstream ofs(output_file, std::ios_base::app);
    if (output_file == "../output/output.txt")
    {
        if (face_exact_one_face_appr_two_vertex_three == 1)
        {
            if (height_map_one_point_cloud_two_terrain_three == 1)
            {
                ofs << "\n== Sur_SimQue_AdpM ==\n";
            }
            else if (height_map_one_point_cloud_two_terrain_three == 2)
            {
                ofs << "\n== Sur_SimQue_AdpC ==\n";
            }
            else if (height_map_one_point_cloud_two_terrain_three == 3)
            {
                ofs << "\n== Sur_SimQue ==\n";
            }
        }
        else if (face_exact_one_face_appr_two_vertex_three == 2)
        {
            ofs << "\n== Simplified_Face_Approximate ==\n";
        }
        else if (face_exact_one_face_appr_two_vertex_three == 3)
        {
            if (height_map_one_point_cloud_two_terrain_three == 1)
            {
                ofs << "\n== Net_SimQue_AdpM ==\n";
            }
            else if (height_map_one_point_cloud_two_terrain_three == 2)
            {
                ofs << "\n== Net_SimQue_AdpC ==\n";
            }
            else if (height_map_one_point_cloud_two_terrain_three == 3)
            {
                ofs << "\n== Net_SimQue ==\n";
            }
        }
    }
    ofs << write_file_header << "\t"
        << height_map_or_point_cloud_to_terrain_time << "\t"
        << height_map_or_point_cloud_to_terrain_memory_usage / 1e6 << "\t"
        << simplification_time << "\t"
        << simplification_memory_usage / 1e6 << "\t"
        << output_size / 1e6 << "\t"
        << query_time << "\t"
        << std::abs(distance_result / height_map_or_point_cloud_exact_distance - 1) << "\t"
        << std::abs(distance_result / terrain_exact_distance - 1) << "\t"
        << knn_query_time << "\t"
        << height_map_or_point_cloud_knn_query_error << "\t"
        << terrain_knn_query_error << "\t"
        << range_query_time << "\t"
        << height_map_or_point_cloud_range_query_error << "\t"
        << terrain_range_query_error << "\n";
    ofs.close();
}

void height_map_or_point_cloud_with_output(
    std::string output_file,
    height_map_geodesic::HeightMap *org_height_map,
    int source_index, int destination_index,
    double height_map_or_point_cloud_exact_distance,
    double terrain_exact_distance, bool run_knn_and_range_query,
    int knn_and_range_query_obj_num, int k_value, double range,
    std::vector<int> &height_map_or_point_cloud_exact_knn_list,
    std::vector<int> &terrain_exact_knn_list,
    std::vector<int> &height_map_or_point_cloud_exact_range_list,
    std::vector<int> &terrain_exact_range_list,
    int height_map_one_point_cloud_two,
    std::string write_file_header,
    int height_map_one_point_cloud_two_terrain_three)
{
    double query_time = 0;
    double query_memory_usage = 0;
    double distance_result = 0;
    double knn_query_time = 0;
    double range_query_time = 0;
    double height_map_or_point_cloud_knn_query_error = 0;
    double terrain_knn_query_error = 0;
    double height_map_or_point_cloud_range_query_error = 0;
    double terrain_range_query_error = 0;
    std::vector<height_map_geodesic::PathPoint> path_result;
    std::vector<int> knn_list;
    std::vector<int> range_list;
    path_result.clear();
    knn_list.clear();
    range_list.clear();
    height_map_or_point_cloud(
        org_height_map, source_index, destination_index,
        query_time, query_memory_usage, distance_result, path_result);

    if (run_knn_and_range_query)
    {
        height_map_or_point_cloud_knn_and_range_query(
            org_height_map, source_index,
            knn_query_time, range_query_time,
            knn_and_range_query_obj_num, k_value, range,
            knn_list, range_list,
            height_map_or_point_cloud_exact_knn_list, terrain_exact_knn_list,
            height_map_or_point_cloud_exact_range_list, terrain_exact_range_list,
            height_map_or_point_cloud_knn_query_error, terrain_knn_query_error,
            height_map_or_point_cloud_range_query_error, terrain_range_query_error);
    }

    std::cout << "Query time: " << query_time << " ms" << std::endl;
    std::cout << "Query memory usage: " << query_memory_usage / 1e6 << " MB" << std::endl;
    std::cout << "Calculated distance: " << distance_result
              << ", height map exact distance: " << height_map_or_point_cloud_exact_distance
              << ", height map distance error: " << distance_result / height_map_or_point_cloud_exact_distance - 1
              << ", terrain exact distance: " << terrain_exact_distance
              << ", terrain distance error: " << distance_result / terrain_exact_distance - 1 << std::endl;
    if (run_knn_and_range_query)
    {
        std::cout << "Knn query time: " << knn_query_time << " ms" << std::endl;
        std::cout << "Height map or point cloud knn error: " << height_map_or_point_cloud_knn_query_error << ", terrain knn error: " << terrain_knn_query_error << std::endl;
        std::cout << "Range query time: " << range_query_time << " ms" << std::endl;
        std::cout << "Height map or point cloud range error: " << height_map_or_point_cloud_range_query_error << ", terrain range error: " << terrain_range_query_error << std::endl;
    }

    std::ofstream ofs(output_file, std::ios_base::app);
    if (output_file == "../output/output.txt")
    {
        if (height_map_one_point_cloud_two == 1)
        {
            if (height_map_one_point_cloud_two_terrain_three == 1)
            {
                ofs << "\n== Eff_Que ==\n";
            }
            else if (height_map_one_point_cloud_two_terrain_three == 2)
            {
                ofs << "\n== Eff_Que_AdpC ==\n";
            }
            else if (height_map_one_point_cloud_two_terrain_three == 3)
            {
                ofs << "\n== Eff_Que_AdpT ==\n";
            }
        }
        else if (height_map_one_point_cloud_two == 2)
        {
            if (height_map_one_point_cloud_two_terrain_three == 1)
            {
                ofs << "\n== Con_Que_AdpM ==\n";
            }
            else if (height_map_one_point_cloud_two_terrain_three == 2)
            {
                ofs << "\n== Con_Que ==\n";
            }
            else if (height_map_one_point_cloud_two_terrain_three == 3)
            {
                ofs << "\n== Con_Que_AdpT ==\n";
            }
        }
    }
    ofs << write_file_header << "\t"
        << 0 << "\t"
        << 0 / 1e6 << "\t"
        << 0 << "\t"
        << query_memory_usage / 1e6 << "\t"
        << 0 / 1e6 << "\t"
        << query_time << "\t"
        << std::abs(distance_result / height_map_or_point_cloud_exact_distance - 1) << "\t"
        << std::abs(distance_result / terrain_exact_distance - 1) << "\t"
        << knn_query_time << "\t"
        << height_map_or_point_cloud_knn_query_error << "\t"
        << terrain_knn_query_error << "\t"
        << range_query_time << "\t"
        << height_map_or_point_cloud_range_query_error << "\t"
        << terrain_range_query_error << "\n";
    ofs.close();
}

void terrain_face_exact_and_face_appr_and_vertex_with_output(
    std::string output_file,
    height_map_geodesic::HeightMap *org_height_map,
    double epsilon, int source_index, int destination_index,
    int face_exact_one_face_appr_two_vertex_three,
    double height_map_or_point_cloud_exact_distance, double terrain_exact_distance,
    bool run_knn_and_range_query, int knn_and_range_query_obj_num,
    int k_value, double range,
    std::vector<int> &height_map_or_point_cloud_exact_knn_list,
    std::vector<int> &terrain_exact_knn_list,
    std::vector<int> &height_map_or_point_cloud_exact_range_list,
    std::vector<int> &terrain_exact_range_list,
    std::string write_file_header,
    int height_map_one_point_cloud_two_terrain_three)
{
    double height_map_or_point_cloud_to_terrain_time = 0;
    double height_map_or_point_cloud_to_terrain_memory_usage = 0;
    double query_time = 0;
    double query_memory_usage = 0;
    double distance_result = 0;
    double knn_query_time = 0;
    double range_query_time = 0;
    double height_map_or_point_cloud_knn_query_error = 0;
    double terrain_knn_query_error = 0;
    double height_map_or_point_cloud_range_query_error = 0;
    double terrain_range_query_error = 0;
    std::vector<geodesic::SurfacePoint> path_result;
    std::vector<int> knn_list;
    std::vector<int> range_list;
    path_result.clear();
    knn_list.clear();
    range_list.clear();
    terrain_face_exact_and_face_appr_and_vertex(
        org_height_map, epsilon,
        source_index, destination_index,
        face_exact_one_face_appr_two_vertex_three,
        height_map_or_point_cloud_to_terrain_time, height_map_or_point_cloud_to_terrain_memory_usage,
        query_time, query_memory_usage, distance_result, path_result);

    if (run_knn_and_range_query)
    {
        terrain_face_exact_and_face_appr_and_vertex_knn_and_range_query(
            org_height_map, epsilon,
            source_index,
            face_exact_one_face_appr_two_vertex_three,
            run_knn_and_range_query, knn_query_time, range_query_time,
            knn_and_range_query_obj_num, k_value, range,
            knn_list, range_list,
            height_map_or_point_cloud_exact_knn_list, terrain_exact_knn_list,
            height_map_or_point_cloud_exact_range_list, terrain_exact_range_list,
            height_map_or_point_cloud_knn_query_error, terrain_knn_query_error,
            height_map_or_point_cloud_range_query_error, terrain_range_query_error);
    }

    if (height_map_one_point_cloud_two_terrain_three == 3)
    {
        height_map_or_point_cloud_to_terrain_time = 0;
        height_map_or_point_cloud_to_terrain_memory_usage = 0;
        std::cout << "Height map or point cloud to terrain time: " << height_map_or_point_cloud_to_terrain_time << " ms" << std::endl;
        std::cout << "Height map or point cloud to terrain memory usage: " << height_map_or_point_cloud_to_terrain_memory_usage / 1e6 << " MB" << std::endl;
    }
    std::cout << "Query time: " << query_time << " ms" << std::endl;
    std::cout << "Query memory usage: " << query_memory_usage / 1e6 << " MB" << std::endl;
    std::cout << "Calculated distance: " << distance_result
              << ", height map exact distance: " << height_map_or_point_cloud_exact_distance
              << ", height map distance error: " << distance_result / height_map_or_point_cloud_exact_distance - 1
              << ", terrain exact distance: " << terrain_exact_distance
              << ", terrain distance error: " << distance_result / terrain_exact_distance - 1 << std::endl;
    if (run_knn_and_range_query)
    {
        std::cout << "Knn query time: " << knn_query_time << " ms" << std::endl;
        std::cout << "Height map or point cloud knn error: " << height_map_or_point_cloud_knn_query_error << ", terrain knn error: " << terrain_knn_query_error << std::endl;
        std::cout << "Range query time: " << range_query_time << " ms" << std::endl;
        std::cout << "Height map or point cloud range error: " << height_map_or_point_cloud_range_query_error << ", terrain range error: " << terrain_range_query_error << std::endl;
    }

    std::ofstream ofs(output_file, std::ios_base::app);
    if (output_file == "../output/output.txt")
    {
        if (face_exact_one_face_appr_two_vertex_three == 1)
        {
            if (height_map_one_point_cloud_two_terrain_three == 1)
            {
                ofs << "\n== Unf_Que_AdpM ==\n";
            }
            else if (height_map_one_point_cloud_two_terrain_three == 2)
            {
                ofs << "\n== Unf_Que_AdpC ==\n";
            }
            else if (height_map_one_point_cloud_two_terrain_three == 3)
            {
                ofs << "\n== Unf_Que ==\n";
            }
        }
        else if (face_exact_one_face_appr_two_vertex_three == 2)
        {
            if (height_map_one_point_cloud_two_terrain_three == 1)
            {
                ofs << "\n== Ste_Que_AdpM ==\n";
            }
            else if (height_map_one_point_cloud_two_terrain_three == 2)
            {
                ofs << "\n== Ste_Que_AdpC ==\n";
            }
            else if (height_map_one_point_cloud_two_terrain_three == 3)
            {
                ofs << "\n== Ste_Que ==\n";
            }
        }
        else if (face_exact_one_face_appr_two_vertex_three == 3)
        {
            if (height_map_one_point_cloud_two_terrain_three == 1)
            {
                ofs << "\n== Dij_Que_AdpM ==\n";
            }
            else if (height_map_one_point_cloud_two_terrain_three == 2)
            {
                ofs << "\n== Dij_Que_AdpC ==\n";
            }
            else if (height_map_one_point_cloud_two_terrain_three == 3)
            {
                ofs << "\n== Dij_Que ==\n";
            }
        }
    }
    ofs << write_file_header << "\t"
        << height_map_or_point_cloud_to_terrain_time << "\t"
        << height_map_or_point_cloud_to_terrain_memory_usage / 1e6 << "\t"
        << 0 << "\t"
        << query_memory_usage / 1e6 << "\t"
        << 0 / 1e6 << "\t"
        << query_time << "\t"
        << std::abs(distance_result / height_map_or_point_cloud_exact_distance - 1) << "\t"
        << std::abs(distance_result / terrain_exact_distance - 1) << "\t"
        << knn_query_time << "\t"
        << height_map_or_point_cloud_knn_query_error << "\t"
        << terrain_knn_query_error << "\t"
        << range_query_time << "\t"
        << height_map_or_point_cloud_range_query_error << "\t"
        << terrain_range_query_error << "\n";
    ofs.close();
}