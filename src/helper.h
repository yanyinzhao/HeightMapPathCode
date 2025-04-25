#include <stdio.h>
#include <stdlib.h>
#include <list>
#include <cmath>
#include <iostream>
#include <fstream>
#include <bits/stdc++.h>
#include "geodesic_algorithm_exact.h"
#include "geodesic_algorithm_subdivision.h"
#include "height_map.h"
#include "point_cloud.h"
#include <algorithm>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>

#include <boost/functional/hash.hpp>

#define DEBUG 1 // 0 not print, 1 print

#define D_()    \
    if (!DEBUG) \
    {           \
    }           \
    else        \
        std::cout

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Triangulation_vertex_base_with_info_2<unsigned, K> Vb;
typedef CGAL::Triangulation_data_structure_2<Vb> Tds;
typedef CGAL::Delaunay_triangulation_2<K, Tds> Delaunay;
typedef Delaunay::Point DelPoint;

struct move_direction
{
    int left;
    int right;
    int top;
    int bottom;

    move_direction() {}

    move_direction(int _left, int _right,
                   int _top, int _bottom)
    {
        left = _left;
        right = _right;
        top = _top;
        bottom = _bottom;
    }
};

double cal_dist(double x1, double y1, double z1,
                double x2, double y2, double z2)
{
    double x = x2 - x1;
    double y = y2 - y1;
    double z = z2 - z1;
    return sqrt(x * x + y * y + z * z);
}

double variance(std::vector<double> samples)
{
    int size = samples.size();
    double variance = 0;
    double mean = 0;
    for (int i = 0; i < size; i++)
    {
        mean += samples[i];
    }
    mean /= size;
    for (int i = 0; i < size; i++)
    {
        variance += pow(samples[i] - mean, 2);
    }
    return variance / size;
}

void canonical_triangle(geodesic::Mesh *mesh, int input_id1, int input_id2, int input_id3,
                        int &output_id1, int &output_id2, int &output_id3)
{
    double x1 = mesh->vertices()[input_id1].getx();
    double y1 = mesh->vertices()[input_id1].gety();
    double x2 = mesh->vertices()[input_id2].getx();
    double y2 = mesh->vertices()[input_id2].gety();
    double x3 = mesh->vertices()[input_id3].getx();
    double y3 = mesh->vertices()[input_id3].gety();

    if (((x1 - x3) * (y2 - y3) - (x2 - x3) * (y1 - y3)) < 0)
    {
        output_id1 = input_id3;
        output_id2 = input_id2;
        output_id3 = input_id1;
    }
    else
    {
        output_id1 = input_id1;
        output_id2 = input_id2;
        output_id3 = input_id3;
    }
}

void compare_distance(double org_dist, double simp_dist,
                      double epsilon, bool &satisfy)
{
    if (simp_dist <= (1 + epsilon) * org_dist && simp_dist >= (1 - epsilon) * org_dist)
    {
        satisfy = true;
    }
    else
    {
        satisfy = false;
    }
}

double epsilon_to_subdivision_level(double epsilon)
{
    assert(epsilon > 0 && epsilon <= 1);
    double subdivision_level = 0;
    if (epsilon > 0 && epsilon <= 0.05)
    {
        subdivision_level = floor(7 / (10 * epsilon)) - 8;
    }
    else if (epsilon > 0.05 && epsilon <= 0.1)
    {
        subdivision_level = floor(7 / (10 * epsilon)) - 2;
    }
    else if (epsilon > 0.1 && epsilon <= 0.25)
    {
        subdivision_level = floor(1 / epsilon) - 1;
    }
    else if (epsilon > 0.25 && epsilon < 1)
    {
        subdivision_level = floor(1 / epsilon);
    }
    if (subdivision_level < 5)
    {
        subdivision_level++;
    }
    return subdivision_level;
}

int epsilon_prime_ds_to_iternation_num(double epsilon_prime_ds)
{
    assert(epsilon_prime_ds > 0 && epsilon_prime_ds <= 1);
    int iternation = 0;
    if (epsilon_prime_ds > 0 && epsilon_prime_ds <= 0.05)
    {
        iternation = 4;
    }
    else if (epsilon_prime_ds > 0.05 && epsilon_prime_ds <= 0.1)
    {
        iternation = 3;
    }
    else if (epsilon_prime_ds > 0.1 && epsilon_prime_ds <= 0.25)
    {
        iternation = 2;
    }
    else if (epsilon_prime_ds > 0.25 && epsilon_prime_ds < 1)
    {
        iternation = 1;
    }
    return iternation;
}

int iternation_num(int a, int b)
{
    if (a == 2)
    {
        return 30;
    }
    else if (a == 3 && b == 2)
    {
        return 1000;
    }
    return 1;
}

void modify_path(std::vector<geodesic::SurfacePoint> &internal_path)
{
    for (int i = 1; i < internal_path.size() - 1; i++)
    {
        geodesic::SurfacePoint &prev = internal_path[i - 1];
        geodesic::SurfacePoint &curr = internal_path[i];
        geodesic::SurfacePoint &next = internal_path[i + 1];
        if (prev.type() == geodesic::VERTEX && curr.type() == geodesic::EDGE && next.type() == geodesic::VERTEX)
        {
            if (prev.distance(curr.base_element()->adjacent_vertices()[0]) + next.distance(curr.base_element()->adjacent_vertices()[0]) >
                prev.distance(curr.base_element()->adjacent_vertices()[1]) + next.distance(curr.base_element()->adjacent_vertices()[1]))
            {
                internal_path[i] = curr.base_element()->adjacent_vertices()[1];
            }
            else
            {
                internal_path[i] = curr.base_element()->adjacent_vertices()[0];
            }
        }
    }
}

void knn_or_range_query(int knn_one_range_two, int k_value, double range,
                        std::vector<std::pair<double, int>> &obj_to_other_obj_distance_and_index_list,
                        std::vector<int> &knn_or_range_list)
{
    knn_or_range_list.clear();
    if (knn_one_range_two == 1)
    {
        for (int i = 0; i < obj_to_other_obj_distance_and_index_list.size(); i++)
        {
            int count = 0;
            if (count >= k_value)
            {
                break;
            }
            if (obj_to_other_obj_distance_and_index_list[i].first != 0)
            {
                knn_or_range_list.push_back(obj_to_other_obj_distance_and_index_list[i].second);
                count++;
            }
        }
    }
    else if (knn_one_range_two == 2)
    {
        for (int i = 0; i < obj_to_other_obj_distance_and_index_list.size(); i++)
        {
            if (obj_to_other_obj_distance_and_index_list[i].first != 0 &&
                obj_to_other_obj_distance_and_index_list[i].first <= range)
            {
                knn_or_range_list.push_back(obj_to_other_obj_distance_and_index_list[i].second);
            }
        }
    }
}

void calculate_knn_or_range_query_error(std::vector<int> &exact_knn_or_range_list,
                                        std::vector<int> &calculated_knn_or_range_list,
                                        double &knn_or_range_error)
{
    int knn_or_range_error_count = 0;
    assert(exact_knn_or_range_list.size() == calculated_knn_or_range_list.size());
    std::vector<int> a;
    std::vector<int> b;
    for (int i = 0; i < exact_knn_or_range_list.size(); i++)
    {
        a.push_back(exact_knn_or_range_list[i]);
        b.push_back(calculated_knn_or_range_list[i]);
    }
    std::sort(a.begin(), a.end());
    std::sort(b.begin(), b.end());
    for (int i = 0; i < exact_knn_or_range_list.size(); i++)
    {
        if (a[i] != b[i])
        {
            knn_or_range_error_count++;
        }
    }

    knn_or_range_error = (double)knn_or_range_error_count / (double)(exact_knn_or_range_list.size());
}

void compare_d_value_and_range_value(double &d_value, double &range)
{
    if (d_value / 2 > range)
    {
        d_value = range * 2;
    }
    if (range > d_value / 2)
    {
        range = d_value / 2;
    }
}

void height_map_merge_four_point(
    height_map_geodesic::HeightMap *org_height_map,
    height_map_geodesic::HeightMap *new_height_map,
    bool &can_merge_four_point,
    int &merged_bottom_left_index,
    int &merged_bottom_left_i,
    int &merged_bottom_left_j,
    int &added_center_point_index,
    std::unordered_map<int, std::unordered_map<int, int>> &temp_dominate_table_map,
    std::unordered_map<int, int> &temp_del_p_dom_by_map,
    std::unordered_map<int, int> &skipped_bottom_left_p_index_map,
    std::unordered_map<int, bool> &prev_deleted,
    std::unordered_map<int, std::unordered_map<int, double>> &prev_adjacent_hm_points_and_distance)
{
    double min_z_coor_vari = 1e100;

    int merged_bottom_right_index;
    int merged_top_left_index;
    int merged_top_right_index;

    for (int j = 1; j < org_height_map->m_ypointnum - 2; ++j)
    {
        for (int i = 1; i < org_height_map->m_xpointnum - 2; ++i)
        {
            int bottom_left_index = i + j * org_height_map->m_xpointnum;
            int bottom_right_index = (i + 1) + j * org_height_map->m_xpointnum;
            int top_left_index = i + (j + 1) * org_height_map->m_xpointnum;
            int top_right_index = (i + 1) + (j + 1) * org_height_map->m_xpointnum;

            height_map_geodesic::HM_Point &bot_lef_v = org_height_map->hm_points()[bottom_left_index];
            height_map_geodesic::HM_Point &bot_rig_v = org_height_map->hm_points()[bottom_right_index];
            height_map_geodesic::HM_Point &top_lef_v = org_height_map->hm_points()[top_left_index];
            height_map_geodesic::HM_Point &top_rig_v = org_height_map->hm_points()[top_right_index];

            // if any four point exist in dominate table, it means it has been merged, we skip it
            if (temp_del_p_dom_by_map.count(bottom_left_index) != 0 || temp_del_p_dom_by_map.count(bottom_right_index) != 0 ||
                temp_del_p_dom_by_map.count(top_left_index) != 0 || temp_del_p_dom_by_map.count(top_right_index) != 0)
            {
                continue;
            }

            // if at previous merge four point iteration, merge these four points does not satisfy
            // error bound, then skip these four points
            if (skipped_bottom_left_p_index_map.count(bottom_left_index) != 0)
            {
                continue;
            }

            // calculate the variance of z coordinate of four points
            std::vector<double> four_points_z_coor;
            four_points_z_coor.clear();
            four_points_z_coor.push_back(bot_lef_v.getz());
            four_points_z_coor.push_back(bot_rig_v.getz());
            four_points_z_coor.push_back(top_lef_v.getz());
            four_points_z_coor.push_back(top_rig_v.getz());
            double four_points_z_coor_vari = variance(four_points_z_coor);

            if (four_points_z_coor_vari <= min_z_coor_vari)
            {
                min_z_coor_vari = four_points_z_coor_vari;

                merged_bottom_left_index = bottom_left_index;
                merged_bottom_right_index = bottom_right_index;
                merged_top_left_index = top_left_index;
                merged_top_right_index = top_right_index;

                merged_bottom_left_i = i;
                merged_bottom_left_j = j;

                can_merge_four_point = true;
            }
        }
    }
    if (!can_merge_four_point)
    {
        return;
    }

    // added center point coordinate
    height_map_geodesic::HM_Point &merged_bot_lef_v = org_height_map->hm_points()[merged_bottom_left_index];
    height_map_geodesic::HM_Point &merged_bot_rig_v = org_height_map->hm_points()[merged_bottom_right_index];
    height_map_geodesic::HM_Point &merged_top_lef_v = org_height_map->hm_points()[merged_top_left_index];
    height_map_geodesic::HM_Point &merged_top_rig_v = org_height_map->hm_points()[merged_top_right_index];

    double added_center_point_x = merged_bot_lef_v.getx() + merged_bot_rig_v.getx() +
                                  merged_top_lef_v.getx() + merged_top_rig_v.getx();
    double added_center_point_y = merged_bot_lef_v.gety() + merged_bot_rig_v.gety() +
                                  merged_top_lef_v.gety() + merged_top_rig_v.gety();
    double added_center_point_z = merged_bot_lef_v.getz() + merged_bot_rig_v.getz() +
                                  merged_top_lef_v.getz() + merged_top_rig_v.getz();
    added_center_point_x /= 4;
    added_center_point_y /= 4;
    added_center_point_z /= 4;
    assert((added_center_point_x == merged_bot_lef_v.getx() + 1) &&
           (added_center_point_y == merged_bot_lef_v.gety() + 1));

    added_center_point_index = new_height_map->hm_points().size();

    // update dominate table map
    std::unordered_map<int, int> one_dominate_table_map;
    one_dominate_table_map.clear();
    one_dominate_table_map[merged_bottom_left_index] = merged_bottom_left_index;
    one_dominate_table_map[merged_bottom_right_index] = merged_bottom_right_index;
    one_dominate_table_map[merged_top_left_index] = merged_top_left_index;
    one_dominate_table_map[merged_top_right_index] = merged_top_right_index;

    temp_del_p_dom_by_map[merged_bottom_left_index] = added_center_point_index;
    temp_del_p_dom_by_map[merged_bottom_right_index] = added_center_point_index;
    temp_del_p_dom_by_map[merged_top_left_index] = added_center_point_index;
    temp_del_p_dom_by_map[merged_top_right_index] = added_center_point_index;

    temp_dominate_table_map[added_center_point_index] = one_dominate_table_map;

    // update new height map
    new_height_map->update_height_map_merge_four_point(
        added_center_point_x, added_center_point_y,
        added_center_point_z, added_center_point_index,
        temp_dominate_table_map, temp_del_p_dom_by_map,
        prev_deleted, prev_adjacent_hm_points_and_distance);
}

// expand into four direction, note that in the coding, we just expand by one original pixel size
// for faster processing, but in paper, we expand by a larger pixel size
// (in the case that the adjacent point is merged)
void height_map_exp_four_direct(
    height_map_geodesic::HeightMap *org_height_map,
    height_map_geodesic::HeightMap *new_height_map,
    bool &can_exp_four_direct,
    int &temp_merged_bottom_left_i,
    int &temp_merged_bottom_left_j,
    int &temp_merged_top_right_i,
    int &temp_merged_top_right_j,
    int &added_center_point_index,
    std::unordered_map<int, std::unordered_map<int, int>> &temp_dominate_table_map,
    std::unordered_map<int, int> &temp_del_p_dom_by_map,
    std::vector<int> &exp_four_direct_v_index_list,
    double &added_point_prev_x, double &added_point_prev_y, double &added_point_prev_z,
    std::unordered_map<int, bool> &prev_deleted,
    std::unordered_map<int, std::unordered_map<int, double>> &prev_adjacent_hm_points_and_distance)
{
    int curr_merged_pixel_width = temp_merged_top_right_i - temp_merged_bottom_left_i + 1;
    int curr_merged_pixel_height = temp_merged_top_right_j - temp_merged_bottom_left_j + 1;

    for (int j = -1; j < curr_merged_pixel_height + 1; j++)
    {
        for (int i = -1; i < curr_merged_pixel_width + 1; i++)
        {
            if (j > -1 && j < curr_merged_pixel_height &&
                i > -1 && i < curr_merged_pixel_width)
            {
                continue;
            }
            int curr_i = temp_merged_bottom_left_i + i;
            int curr_j = temp_merged_bottom_left_j + j;
            int curr_index = curr_i + curr_j * org_height_map->m_xpointnum;

            // boundary
            if (curr_i == 0 || curr_i == org_height_map->m_xpointnum - 1 ||
                curr_j == 0 || curr_j == org_height_map->m_ypointnum - 1)
            {
                can_exp_four_direct = false;
                return;
            }

            // if any point exist in dominate table, it means it has been merged,
            // then we cannot expand, so we directly terminate
            if (temp_del_p_dom_by_map.count(curr_index) != 0)
            {
                can_exp_four_direct = false;
                return;
            }
            exp_four_direct_v_index_list.push_back(curr_index);
        }
    }
    can_exp_four_direct = true;

    // the points we expand one pixel into four direction should be same as
    // the adjacent points of the previous added center point
    assert(new_height_map->hm_points()[added_center_point_index].adj_hm_points_dist().size() ==
           exp_four_direct_v_index_list.size());
    for (int i = 0; i < exp_four_direct_v_index_list.size(); i++)
    {
        assert(new_height_map->hm_points()[added_center_point_index].adj_hm_points_dist().count(
                   exp_four_direct_v_index_list[i]) != 0);
    }

    // update the bottom left and top right i j
    temp_merged_bottom_left_i--;
    temp_merged_bottom_left_j--;
    temp_merged_top_right_i++;
    temp_merged_top_right_j++;

    // update the previous added center point coordiante
    double prev_added_point_x = new_height_map->hm_points()[added_center_point_index].getx();
    double prev_added_point_y = new_height_map->hm_points()[added_center_point_index].gety();
    double prev_added_point_z = new_height_map->hm_points()[added_center_point_index].getz();

    int prev_added_p_dominate_p_num = temp_dominate_table_map[added_center_point_index].size();

    double added_point_updated_x = prev_added_point_x * prev_added_p_dominate_p_num;
    double added_point_updated_y = prev_added_point_y * prev_added_p_dominate_p_num;
    double added_point_updated_z = prev_added_point_z * prev_added_p_dominate_p_num;

    for (int i = 0; i < exp_four_direct_v_index_list.size(); i++)
    {
        assert(org_height_map->hm_points()[exp_four_direct_v_index_list[i]].getx() ==
                   new_height_map->hm_points()[exp_four_direct_v_index_list[i]].getx() &&
               org_height_map->hm_points()[exp_four_direct_v_index_list[i]].gety() ==
                   new_height_map->hm_points()[exp_four_direct_v_index_list[i]].gety() &&
               org_height_map->hm_points()[exp_four_direct_v_index_list[i]].getz() ==
                   new_height_map->hm_points()[exp_four_direct_v_index_list[i]].getz());
        added_point_updated_x += new_height_map->hm_points()[exp_four_direct_v_index_list[i]].getx();
        added_point_updated_y += new_height_map->hm_points()[exp_four_direct_v_index_list[i]].gety();
        added_point_updated_z += new_height_map->hm_points()[exp_four_direct_v_index_list[i]].getz();
    }
    added_point_updated_x /= (prev_added_p_dominate_p_num + exp_four_direct_v_index_list.size());
    added_point_updated_y /= (prev_added_p_dominate_p_num + exp_four_direct_v_index_list.size());
    added_point_updated_z /= (prev_added_p_dominate_p_num + exp_four_direct_v_index_list.size());

    // the index of the bottom left point of larger pixel after four direction expanding
    int new_merged_bottom_left_index = temp_merged_bottom_left_i + temp_merged_bottom_left_j * org_height_map->m_xpointnum;

    assert(added_point_updated_x == new_height_map->hm_points()[new_merged_bottom_left_index].getx() + curr_merged_pixel_width + 1 &&
           added_point_updated_y == new_height_map->hm_points()[new_merged_bottom_left_index].gety() + curr_merged_pixel_height + 1);

    // update dominate table map
    std::unordered_map<int, int> one_dominate_table_map;
    one_dominate_table_map.clear();
    for (auto ite : temp_dominate_table_map[added_center_point_index])
    {
        one_dominate_table_map[ite.first] = ite.first;
    }
    for (int i = 0; i < exp_four_direct_v_index_list.size(); i++)
    {
        one_dominate_table_map[exp_four_direct_v_index_list[i]] = exp_four_direct_v_index_list[i];
    }
    temp_dominate_table_map.erase(added_center_point_index);
    temp_dominate_table_map[added_center_point_index] = one_dominate_table_map;

    for (int i = 0; i < exp_four_direct_v_index_list.size(); i++)
    {
        temp_del_p_dom_by_map[exp_four_direct_v_index_list[i]] = added_center_point_index;
    }

    // update new height map
    new_height_map->update_height_map_exp_four_direct(
        added_point_updated_x, added_point_updated_y,
        added_point_updated_z, added_center_point_index,
        temp_del_p_dom_by_map, exp_four_direct_v_index_list,
        added_point_prev_x, added_point_prev_y, added_point_prev_z,
        prev_deleted, prev_adjacent_hm_points_and_distance);
}

// expand into one, two, three direction, note that in the coding, we just expand by one original pixel size
// for faster processing, but in paper, we expand by a larger pixel size
// (in the case that the adjacent point is merged)
void height_map_exp_three_two_one_direct(
    height_map_geodesic::HeightMap *org_height_map,
    height_map_geodesic::HeightMap *new_height_map,
    bool &can_exp_three_two_one_direct,
    int &temp_merged_bottom_left_i,
    int &temp_merged_bottom_left_j,
    int &temp_merged_top_right_i,
    int &temp_merged_top_right_j,
    int &added_center_point_index,
    std::unordered_map<int, std::unordered_map<int, int>> &temp_dominate_table_map,
    std::unordered_map<int, int> &temp_del_p_dom_by_map,
    std::vector<int> &exp_three_two_one_direct_v_index_list,
    move_direction &move_direction_list,
    double &added_point_prev_x, double &added_point_prev_y, double &added_point_prev_z,
    std::unordered_map<int, bool> &prev_deleted,
    std::unordered_map<int, std::unordered_map<int, double>> &prev_adjacent_hm_points_and_distance)
{
    int curr_merged_pixel_width = temp_merged_top_right_i - temp_merged_bottom_left_i + 1;
    int curr_merged_pixel_height = temp_merged_top_right_j - temp_merged_bottom_left_j + 1;

    for (int j = move_direction_list.bottom; j < curr_merged_pixel_height + move_direction_list.top; j++)
    {
        for (int i = move_direction_list.left; i < curr_merged_pixel_width + move_direction_list.right; i++)
        {
            if (j > -1 && j < curr_merged_pixel_height &&
                i > -1 && i < curr_merged_pixel_width)
            {
                continue;
            }
            int curr_i = temp_merged_bottom_left_i + i;
            int curr_j = temp_merged_bottom_left_j + j;
            int curr_index = curr_i + curr_j * org_height_map->m_xpointnum;

            // boundary
            if (curr_i == 0 || curr_i == org_height_map->m_xpointnum - 1 ||
                curr_j == 0 || curr_j == org_height_map->m_ypointnum - 1)
            {
                can_exp_three_two_one_direct = false;
                return;
            }

            // if any point exist in dominate table, it means it has been merged,
            // then we cannot expand, so we directly terminate
            if (temp_del_p_dom_by_map.count(curr_index) != 0)
            {
                can_exp_three_two_one_direct = false;
                return;
            }
            exp_three_two_one_direct_v_index_list.push_back(curr_index);
        }
    }
    can_exp_three_two_one_direct = true;

    // update the bottom left and top right i j
    temp_merged_bottom_left_i += move_direction_list.left;
    temp_merged_bottom_left_j += move_direction_list.bottom;
    temp_merged_top_right_i += move_direction_list.right;
    temp_merged_top_right_j += move_direction_list.top;

    // update the previous added center point coordiante
    double prev_added_point_x = new_height_map->hm_points()[added_center_point_index].getx();
    double prev_added_point_y = new_height_map->hm_points()[added_center_point_index].gety();
    double prev_added_point_z = new_height_map->hm_points()[added_center_point_index].getz();

    int prev_added_p_dominate_p_num = temp_dominate_table_map[added_center_point_index].size();

    double added_point_updated_x = prev_added_point_x * prev_added_p_dominate_p_num;
    double added_point_updated_y = prev_added_point_y * prev_added_p_dominate_p_num;
    double added_point_updated_z = prev_added_point_z * prev_added_p_dominate_p_num;

    for (int i = 0; i < exp_three_two_one_direct_v_index_list.size(); i++)
    {
        assert(org_height_map->hm_points()[exp_three_two_one_direct_v_index_list[i]].getx() ==
                   new_height_map->hm_points()[exp_three_two_one_direct_v_index_list[i]].getx() &&
               org_height_map->hm_points()[exp_three_two_one_direct_v_index_list[i]].gety() ==
                   new_height_map->hm_points()[exp_three_two_one_direct_v_index_list[i]].gety() &&
               org_height_map->hm_points()[exp_three_two_one_direct_v_index_list[i]].getz() ==
                   new_height_map->hm_points()[exp_three_two_one_direct_v_index_list[i]].getz());
        added_point_updated_x += new_height_map->hm_points()[exp_three_two_one_direct_v_index_list[i]].getx();
        added_point_updated_y += new_height_map->hm_points()[exp_three_two_one_direct_v_index_list[i]].gety();
        added_point_updated_z += new_height_map->hm_points()[exp_three_two_one_direct_v_index_list[i]].getz();
    }
    added_point_updated_x /= (prev_added_p_dominate_p_num + exp_three_two_one_direct_v_index_list.size());
    added_point_updated_y /= (prev_added_p_dominate_p_num + exp_three_two_one_direct_v_index_list.size());
    added_point_updated_z /= (prev_added_p_dominate_p_num + exp_three_two_one_direct_v_index_list.size());

    // the index of the bottom left and top right point of larger pixel after three two one direction expanding
    int new_merged_bottom_left_index = temp_merged_bottom_left_i + temp_merged_bottom_left_j * org_height_map->m_xpointnum;
    int new_merged_top_right_index = temp_merged_top_right_i + temp_merged_top_right_j * org_height_map->m_xpointnum;

    assert(org_height_map->hm_points()[new_merged_bottom_left_index].getx() ==
               new_height_map->hm_points()[new_merged_bottom_left_index].getx() &&
           org_height_map->hm_points()[new_merged_bottom_left_index].gety() ==
               new_height_map->hm_points()[new_merged_bottom_left_index].gety() &&
           org_height_map->hm_points()[new_merged_top_right_index].getx() ==
               new_height_map->hm_points()[new_merged_top_right_index].getx() &&
           org_height_map->hm_points()[new_merged_top_right_index].gety() ==
               new_height_map->hm_points()[new_merged_top_right_index].gety());

    // the x y coordinate of the bottom left and top right point of larger pixel after three two one direction expanding
    double new_merged_bottom_left_x = org_height_map->hm_points()[new_merged_bottom_left_index].getx();
    double new_merged_bottom_left_y = org_height_map->hm_points()[new_merged_bottom_left_index].gety();
    double new_merged_top_right_x = org_height_map->hm_points()[new_merged_top_right_index].getx();
    double new_merged_top_right_y = org_height_map->hm_points()[new_merged_top_right_index].gety();

    // update dominate table map
    std::unordered_map<int, int> one_dominate_table_map;
    one_dominate_table_map.clear();
    for (auto ite : temp_dominate_table_map[added_center_point_index])
    {
        one_dominate_table_map[ite.first] = ite.first;
    }
    for (int i = 0; i < exp_three_two_one_direct_v_index_list.size(); i++)
    {
        one_dominate_table_map[exp_three_two_one_direct_v_index_list[i]] = exp_three_two_one_direct_v_index_list[i];
    }
    temp_dominate_table_map.erase(added_center_point_index);
    temp_dominate_table_map[added_center_point_index] = one_dominate_table_map;

    for (int i = 0; i < exp_three_two_one_direct_v_index_list.size(); i++)
    {
        temp_del_p_dom_by_map[exp_three_two_one_direct_v_index_list[i]] = added_center_point_index;
    }

    // update new height map
    new_height_map->update_height_map_exp_three_two_one_direct(
        new_merged_bottom_left_x, new_merged_bottom_left_y,
        new_merged_top_right_x, new_merged_top_right_y,
        added_point_updated_x, added_point_updated_y,
        added_point_updated_z, added_center_point_index,
        temp_del_p_dom_by_map, exp_three_two_one_direct_v_index_list,
        added_point_prev_x, added_point_prev_y, added_point_prev_z,
        prev_deleted, prev_adjacent_hm_points_and_distance);
}

// in this checking, we just check the distance between the adjacent point of center point
// (if this point is the original point) + the deleted point dominated by the adjacent point
// of center point (if this point is added point)
void height_map_cal_simp_dist_and_check(
    height_map_geodesic::HeightMap *org_height_map,
    height_map_geodesic::HeightMap *new_height_map,
    double epsilon, bool &satisfy,
    int &added_center_point_index,
    std::unordered_map<int, std::unordered_map<int, int>> &temp_dominate_table_map,
    std::unordered_map<int, int> &temp_del_p_dom_by_map)
{
    // store the distance between point to point on original height map
    std::unordered_map<int, double> dist_p_to_p_org_hm;
    dist_p_to_p_org_hm.clear();

    // store the distance between point to point on post height map
    std::unordered_map<int, double> dist_p_to_p_post_hm;
    dist_p_to_p_post_hm.clear();

    // store the point that we need to check distance between them
    std::vector<int> check_dist_p_index;
    check_dist_p_index.clear();

    int count = 0;
    for (auto ite : new_height_map->hm_points()[added_center_point_index].adj_hm_points_dist())
    {
        if (count < 3)
        {
            int adj_p_of_added_p_index = ite.first;
            if (temp_dominate_table_map.count(adj_p_of_added_p_index) == 0)
            {
                assert(adj_p_of_added_p_index < org_height_map->hm_points().size());
                check_dist_p_index.push_back(adj_p_of_added_p_index);
                count++;
            }
        }
    }

    height_map_geodesic::HeightMapGeodesicAlgorithmDijkstra org_algorithm(org_height_map);
    height_map_geodesic::HeightMapGeodesicAlgorithmDijkstra post_algorithm(new_height_map);
    double const distance_limit = height_map_geodesic::INFIN;

    std::vector<height_map_geodesic::PathPoint> org_one_source_v_list;
    std::vector<height_map_geodesic::PathPoint> org_destinations_v_list;
    std::vector<height_map_geodesic::PathPoint> post_one_source_v_list;
    std::vector<height_map_geodesic::PathPoint> post_destinations_v_list;

    // the point that we need to calculate the real pairwise distance on post height map
    std::unordered_map<int, int> post_run_alg_index_map;
    post_run_alg_index_map.clear();

    for (int i = 0; i < check_dist_p_index.size(); i++)
    {
        if (temp_del_p_dom_by_map.count(check_dist_p_index[i]) != 0)
        {
            assert(false);
            int dom_p_of_curr_del_p = temp_del_p_dom_by_map[check_dist_p_index[i]];
            for (auto ite : new_height_map->hm_points()[dom_p_of_curr_del_p].adj_hm_points_dist())
            {
                if (post_run_alg_index_map.count(ite.first) == 0)
                {
                    post_run_alg_index_map[ite.first] = ite.first;
                }
            }
        }
        else
        {
            post_run_alg_index_map[check_dist_p_index[i]] = check_dist_p_index[i];
        }
    }

    std::unordered_map<int, int> calculated_src_map;
    calculated_src_map.clear();

    // checking with source as non-deleted point
    for (int i = 0; i < check_dist_p_index.size(); i++)
    {
        // if this point has been deleted, we skip it now
        // we just focus on the non-deleted point
        if (temp_del_p_dom_by_map.count(check_dist_p_index[i]) != 0)
        {
            continue;
        }

        calculated_src_map[check_dist_p_index[i]] = check_dist_p_index[i];

        org_one_source_v_list.clear();
        org_destinations_v_list.clear();
        post_one_source_v_list.clear();
        post_destinations_v_list.clear();

        org_one_source_v_list.push_back(height_map_geodesic::PathPoint(&org_height_map->hm_points()[check_dist_p_index[i]]));
        post_one_source_v_list.push_back(height_map_geodesic::PathPoint(&new_height_map->hm_points()[check_dist_p_index[i]]));

        for (int j = i; j < check_dist_p_index.size(); j++)
        {
            org_destinations_v_list.push_back(height_map_geodesic::PathPoint(&org_height_map->hm_points()[check_dist_p_index[j]]));
        }

        for (auto ite : post_run_alg_index_map)
        {
            int x_post = check_dist_p_index[i];
            int y_post = ite.first;
            int x_y_post;
            if (x_post > y_post)
            {
                int temp_post = y_post;
                y_post = x_post;
                x_post = temp_post;
            }
            hash_function_two_keys_to_one_key_int(new_height_map->hm_points().size(), x_post, y_post, x_y_post);
            if (dist_p_to_p_post_hm.count(x_y_post) == 0)
            {
                post_destinations_v_list.push_back(height_map_geodesic::PathPoint(&new_height_map->hm_points()[ite.first]));
            }
        }

        org_algorithm.propagate(org_one_source_v_list, &org_destinations_v_list);
        post_algorithm.propagate(post_one_source_v_list, &post_destinations_v_list);

        // calculate the original height map distance
        for (int j = i; j < check_dist_p_index.size(); j++)
        {
            int x_org = check_dist_p_index[i];
            int y_org = check_dist_p_index[j];
            int x_y_org;
            if (x_org > y_org)
            {
                int temp_org = y_org;
                y_org = x_org;
                x_org = temp_org;
            }
            hash_function_two_keys_to_one_key_int(org_height_map->hm_points().size(), x_org, y_org, x_y_org);
            if (dist_p_to_p_org_hm.count(x_y_org) == 0)
            {
                double org_dist;
                org_algorithm.best_source(org_destinations_v_list[j - i], org_dist);
                dist_p_to_p_org_hm[x_y_org] = org_dist;
            }
        }

        // calculate the post height map distance
        for (auto ite : post_run_alg_index_map)
        {
            int x_post = check_dist_p_index[i];
            int y_post = ite.first;
            int x_y_post;
            if (x_post > y_post)
            {
                int temp_post = y_post;
                y_post = x_post;
                x_post = temp_post;
            }
            hash_function_two_keys_to_one_key_int(new_height_map->hm_points().size(), x_post, y_post, x_y_post);
            if (dist_p_to_p_post_hm.count(x_y_post) == 0)
            {
                double simp_dist;
                height_map_geodesic::PathPoint p(&new_height_map->hm_points()[ite.first]);
                post_algorithm.best_source(p, simp_dist);
                dist_p_to_p_post_hm[x_y_post] = simp_dist;
            }
        }

        // compare the distance
        for (int j = i; j < check_dist_p_index.size(); j++)
        {
            int x_org = check_dist_p_index[i];
            int y_org = check_dist_p_index[j];
            int x_y_org;
            if (x_org > y_org)
            {
                int temp_org = y_org;
                y_org = x_org;
                x_org = temp_org;
            }
            hash_function_two_keys_to_one_key_int(org_height_map->hm_points().size(), x_org, y_org, x_y_org);
            double org_dist = dist_p_to_p_org_hm[x_y_org];

            double simp_dist = 1e100;

            if (temp_del_p_dom_by_map.count(check_dist_p_index[j]) != 0)
            {
                int dom_p_of_curr_del_p = temp_del_p_dom_by_map[check_dist_p_index[j]];
                for (auto ite : new_height_map->hm_points()[dom_p_of_curr_del_p].adj_hm_points_dist())
                {
                    int x_post = check_dist_p_index[i];
                    int y_post = ite.first;
                    int x_y_post;
                    if (x_post > y_post)
                    {
                        int temp_post = y_post;
                        y_post = x_post;
                        x_post = temp_post;
                    }
                    hash_function_two_keys_to_one_key_int(new_height_map->hm_points().size(), x_post, y_post, x_y_post);
                    double x1 = new_height_map->hm_points()[ite.first].getx();
                    double y1 = new_height_map->hm_points()[ite.first].gety();
                    double z1 = new_height_map->hm_points()[ite.first].getz();
                    double x2 = new_height_map->hm_points()[check_dist_p_index[j]].getx();
                    double y2 = new_height_map->hm_points()[check_dist_p_index[j]].gety();
                    double z2 = new_height_map->hm_points()[dom_p_of_curr_del_p].getz();
                    double simp_dist2 = dist_p_to_p_post_hm[x_y_post] + cal_dist(x1, y1, z1, x2, y2, z2);
                    simp_dist = std::min(simp_dist, simp_dist2);
                }
            }
            else
            {
                int x_post = check_dist_p_index[i];
                int y_post = check_dist_p_index[j];
                int x_y_post;
                if (x_post > y_post)
                {
                    int temp_post = y_post;
                    y_post = x_post;
                    x_post = temp_post;
                }
                hash_function_two_keys_to_one_key_int(new_height_map->hm_points().size(), x_post, y_post, x_y_post);
                simp_dist = dist_p_to_p_post_hm[x_y_post];
            }

            compare_distance(org_dist, simp_dist, epsilon, satisfy);
            if (!satisfy)
            {
                return;
            }
        }
    }

    // calculate the original height map distance with source that has not calculated, and store it
    for (int i = 0; i < check_dist_p_index.size(); i++)
    {
        if (temp_del_p_dom_by_map.count(check_dist_p_index[i]) != 0)
        {
            org_one_source_v_list.clear();
            org_destinations_v_list.clear();
            org_one_source_v_list.push_back(height_map_geodesic::PathPoint(&org_height_map->hm_points()[check_dist_p_index[i]]));
            for (int j = i; j < check_dist_p_index.size(); j++)
            {
                org_destinations_v_list.push_back(height_map_geodesic::PathPoint(&org_height_map->hm_points()[check_dist_p_index[j]]));
            }
            org_algorithm.propagate(org_one_source_v_list, &org_destinations_v_list);
            for (int j = i; j < check_dist_p_index.size(); j++)
            {
                int x_org = check_dist_p_index[i];
                int y_org = check_dist_p_index[j];
                int x_y_org;
                if (x_org > y_org)
                {
                    int temp_org = y_org;
                    y_org = x_org;
                    x_org = temp_org;
                }
                hash_function_two_keys_to_one_key_int(org_height_map->hm_points().size(), x_org, y_org, x_y_org);
                if (dist_p_to_p_org_hm.count(x_y_org) == 0)
                {
                    double org_dist;
                    org_algorithm.best_source(org_destinations_v_list[j - i], org_dist);
                    dist_p_to_p_org_hm[x_y_org] = org_dist;
                }
            }
        }
    }

    // calcualte the post height map distance with source that has not calculated, and store it
    for (auto ite : post_run_alg_index_map)
    {
        if (calculated_src_map.count(ite.first) == 0)
        {
            post_one_source_v_list.clear();
            post_destinations_v_list.clear();

            post_one_source_v_list.push_back(height_map_geodesic::PathPoint(&new_height_map->hm_points()[ite.first]));

            for (auto ite2 : post_run_alg_index_map)
            {
                int x = ite.first;
                int y = ite2.first;
                int x_y;
                if (x > y)
                {
                    int temp = y;
                    y = x;
                    x = temp;
                }
                hash_function_two_keys_to_one_key_int(new_height_map->hm_points().size(), x, y, x_y);
                if (dist_p_to_p_post_hm.count(x_y) == 0)
                {
                    post_destinations_v_list.push_back(height_map_geodesic::PathPoint(&new_height_map->hm_points()[ite2.first]));
                }
            }
            post_algorithm.propagate(post_one_source_v_list, &post_destinations_v_list);

            // calculate the post height map distance
            for (auto ite2 : post_run_alg_index_map)
            {
                int x = ite.first;
                int y = ite2.first;
                int x_y;
                if (x > y)
                {
                    int temp = y;
                    y = x;
                    x = temp;
                }
                hash_function_two_keys_to_one_key_int(new_height_map->hm_points().size(), x, y, x_y);
                if (dist_p_to_p_post_hm.count(x_y) == 0)
                {
                    double simp_dist;
                    height_map_geodesic::PathPoint p(&new_height_map->hm_points()[ite2.first]);
                    post_algorithm.best_source(p, simp_dist);
                    dist_p_to_p_post_hm[x_y] = simp_dist;
                }
            }
        }
    }

    // in the following, we have already calculated the distance of the point that we
    // really need to calculate the pairwise distance between then on post height map,
    // and also the distance on original height map

    // checking with source as deleted point
    for (int i = 0; i < check_dist_p_index.size(); i++)
    {
        // if this point has been not deleted, we skip it now
        // we just focus on the deleted point
        if (temp_del_p_dom_by_map.count(check_dist_p_index[i]) == 0)
        {
            continue;
        }
        int dom_p_of_curr_del_p_src = temp_del_p_dom_by_map[check_dist_p_index[i]];

        // there are two special cases here, we consider them first
        for (int j = i; j < check_dist_p_index.size(); j++)
        {
            // if both points s and t are deleted
            if (temp_del_p_dom_by_map.count(check_dist_p_index[j]) != 0)
            {
                int x_org = check_dist_p_index[i];
                int y_org = check_dist_p_index[j];
                int x_y_org;
                if (x_org > y_org)
                {
                    int temp_org = y_org;
                    y_org = x_org;
                    x_org = temp_org;
                }
                hash_function_two_keys_to_one_key_int(org_height_map->hm_points().size(), x_org, y_org, x_y_org);
                double org_dist = dist_p_to_p_org_hm[x_y_org];
                double simp_dist = 1e100;

                int dom_p_of_curr_del_p_dest = temp_del_p_dom_by_map[check_dist_p_index[j]];

                // if s and t are dominated by same added center point,
                // then the distance is just Euclidean distance
                if (dom_p_of_curr_del_p_src == dom_p_of_curr_del_p_dest)
                {
                    double x1 = new_height_map->hm_points()[check_dist_p_index[i]].getx();
                    double y1 = new_height_map->hm_points()[check_dist_p_index[i]].gety();
                    double z1 = new_height_map->hm_points()[dom_p_of_curr_del_p_src].getz();
                    double x2 = new_height_map->hm_points()[check_dist_p_index[j]].getx();
                    double y2 = new_height_map->hm_points()[check_dist_p_index[j]].gety();
                    double z2 = new_height_map->hm_points()[dom_p_of_curr_del_p_dest].getz();
                    simp_dist = cal_dist(x1, y1, z1, x2, y2, z2);
                }

                // if s and t are dominated by two different added center points a and b,
                // but a and b are adjacent:
                else if (dom_p_of_curr_del_p_src != dom_p_of_curr_del_p_dest &&
                         new_height_map->hm_points()[dom_p_of_curr_del_p_src].adj_hm_points_dist().count(dom_p_of_curr_del_p_dest) != 0 &&
                         new_height_map->hm_points()[dom_p_of_curr_del_p_dest].adj_hm_points_dist().count(dom_p_of_curr_del_p_src) != 0)
                {
                    // (a) if s and t are adjacent on original height map, then the distance is just Euclidean distance
                    if (org_height_map->hm_points()[check_dist_p_index[i]].adj_hm_points_dist().count(check_dist_p_index[j]) != 0 &&
                        org_height_map->hm_points()[check_dist_p_index[j]].adj_hm_points_dist().count(check_dist_p_index[i]) != 0)
                    {
                        double x1 = new_height_map->hm_points()[check_dist_p_index[i]].getx();
                        double y1 = new_height_map->hm_points()[check_dist_p_index[i]].gety();
                        double z1 = new_height_map->hm_points()[dom_p_of_curr_del_p_src].getz();
                        double x2 = new_height_map->hm_points()[check_dist_p_index[j]].getx();
                        double y2 = new_height_map->hm_points()[check_dist_p_index[j]].gety();
                        double z2 = new_height_map->hm_points()[dom_p_of_curr_del_p_dest].getz();
                        simp_dist = cal_dist(x1, y1, z1, x2, y2, z2);
                    }
                    // (b) if s and t are not adjacent on original height map, then the distance is min (sa + at, sb + bt)
                    else if (org_height_map->hm_points()[check_dist_p_index[i]].adj_hm_points_dist().count(check_dist_p_index[j]) == 0 &&
                             org_height_map->hm_points()[check_dist_p_index[j]].adj_hm_points_dist().count(check_dist_p_index[i]) == 0)
                    {
                        // s
                        double x1 = new_height_map->hm_points()[check_dist_p_index[i]].getx();
                        double y1 = new_height_map->hm_points()[check_dist_p_index[i]].gety();
                        double z1 = new_height_map->hm_points()[dom_p_of_curr_del_p_src].getz();
                        // t
                        double x2 = new_height_map->hm_points()[check_dist_p_index[j]].getx();
                        double y2 = new_height_map->hm_points()[check_dist_p_index[j]].gety();
                        double z2 = new_height_map->hm_points()[dom_p_of_curr_del_p_dest].getz();
                        // a
                        double x3 = new_height_map->hm_points()[dom_p_of_curr_del_p_src].getx();
                        double y3 = new_height_map->hm_points()[dom_p_of_curr_del_p_src].gety();
                        double z3 = new_height_map->hm_points()[dom_p_of_curr_del_p_src].getz();
                        // b
                        double x4 = new_height_map->hm_points()[dom_p_of_curr_del_p_dest].getx();
                        double y4 = new_height_map->hm_points()[dom_p_of_curr_del_p_dest].gety();
                        double z4 = new_height_map->hm_points()[dom_p_of_curr_del_p_dest].getz();
                        simp_dist = std::min(cal_dist(x1, y1, z1, x3, y3, z3) + cal_dist(x2, y2, z2, x3, y3, z3),
                                             cal_dist(x1, y1, z1, x4, y4, z4) + cal_dist(x2, y2, z2, x4, y4, z4));
                    }
                    else
                    {
                        assert(false);
                    }
                }
                compare_distance(org_dist, simp_dist, epsilon, satisfy);
                if (!satisfy)
                {
                    return;
                }
            }
        }

        // the normal case
        for (int j = i; j < check_dist_p_index.size(); j++)
        {
            int x_org = check_dist_p_index[i];
            int y_org = check_dist_p_index[j];
            int x_y_org;
            if (x_org > y_org)
            {
                int temp_org = y_org;
                y_org = x_org;
                x_org = temp_org;
            }
            hash_function_two_keys_to_one_key_int(org_height_map->hm_points().size(), x_org, y_org, x_y_org);
            double org_dist = dist_p_to_p_org_hm[x_y_org];
            double simp_dist = 1e100;

            // if both points s and t are deleted
            if (temp_del_p_dom_by_map.count(check_dist_p_index[j]) != 0)
            {
                int dom_p_of_curr_del_p_dest = temp_del_p_dom_by_map[check_dist_p_index[j]];

                // we considered this special case before
                if (dom_p_of_curr_del_p_src == dom_p_of_curr_del_p_dest)
                {
                    continue;
                }
                // we considered this special case before
                else if (dom_p_of_curr_del_p_src != dom_p_of_curr_del_p_dest &&
                         new_height_map->hm_points()[dom_p_of_curr_del_p_src].adj_hm_points_dist().count(dom_p_of_curr_del_p_dest) != 0 &&
                         new_height_map->hm_points()[dom_p_of_curr_del_p_dest].adj_hm_points_dist().count(dom_p_of_curr_del_p_src) != 0)
                {
                    continue;
                }

                // if s and t are dominated by two different added center points a and b,
                // but a and b are not adjacent, the distance is min (sa + \forall ab + bt)
                else
                {
                    for (auto ite : new_height_map->hm_points()[dom_p_of_curr_del_p_src].adj_hm_points_dist())
                    {
                        for (auto ite2 : new_height_map->hm_points()[dom_p_of_curr_del_p_dest].adj_hm_points_dist())
                        {
                            int x_post = ite.first;
                            int y_post = ite2.first;
                            int x_y_post;
                            if (x_post > y_post)
                            {
                                int temp_post = y_post;
                                y_post = x_post;
                                x_post = temp_post;
                            }
                            hash_function_two_keys_to_one_key_int(new_height_map->hm_points().size(), x_post, y_post, x_y_post);
                            double x1 = new_height_map->hm_points()[ite.first].getx();
                            double y1 = new_height_map->hm_points()[ite.first].gety();
                            double z1 = new_height_map->hm_points()[ite.first].getz();
                            double x2 = new_height_map->hm_points()[check_dist_p_index[i]].getx();
                            double y2 = new_height_map->hm_points()[check_dist_p_index[i]].gety();
                            double z2 = new_height_map->hm_points()[dom_p_of_curr_del_p_src].getz();
                            double x3 = new_height_map->hm_points()[ite2.first].getx();
                            double y3 = new_height_map->hm_points()[ite2.first].gety();
                            double z3 = new_height_map->hm_points()[ite2.first].getz();
                            double x4 = new_height_map->hm_points()[check_dist_p_index[j]].getx();
                            double y4 = new_height_map->hm_points()[check_dist_p_index[j]].gety();
                            double z4 = new_height_map->hm_points()[dom_p_of_curr_del_p_dest].getz();
                            double simp_dist2 = dist_p_to_p_post_hm[x_y_post] + cal_dist(x1, y1, z1, x2, y2, z2) + cal_dist(x3, y3, z3, x4, y4, z4);
                            simp_dist = std::min(simp_dist, simp_dist2);
                        }
                    }
                }
            }
            else
            {
                for (auto ite : new_height_map->hm_points()[dom_p_of_curr_del_p_src].adj_hm_points_dist())
                {
                    int x_post = ite.first;
                    int y_post = check_dist_p_index[j];
                    int x_y_post;
                    if (x_post > y_post)
                    {
                        int temp_post = y_post;
                        y_post = x_post;
                        x_post = temp_post;
                    }
                    hash_function_two_keys_to_one_key_int(new_height_map->hm_points().size(), x_post, y_post, x_y_post);
                    double x1 = new_height_map->hm_points()[ite.first].getx();
                    double y1 = new_height_map->hm_points()[ite.first].gety();
                    double z1 = new_height_map->hm_points()[ite.first].getz();
                    double x2 = new_height_map->hm_points()[check_dist_p_index[i]].getx();
                    double y2 = new_height_map->hm_points()[check_dist_p_index[i]].gety();
                    double z2 = new_height_map->hm_points()[dom_p_of_curr_del_p_src].getz();
                    double simp_dist2 = dist_p_to_p_post_hm[x_y_post] + cal_dist(x1, y1, z1, x2, y2, z2);
                    simp_dist = std::min(simp_dist, simp_dist2);
                }
            }
            compare_distance(org_dist, simp_dist, epsilon, satisfy);
            if (!satisfy)
            {
                return;
            }
        }
    }
}

// only checking the point on the original height map
void height_map_cal_simp_dist_and_check_naive(
    height_map_geodesic::HeightMap *org_height_map,
    height_map_geodesic::HeightMap *new_height_map,
    std::unordered_map<int, double> &all_dist_p_to_p_org_height_map,
    double epsilon, bool &satisfy,
    std::unordered_map<int, std::unordered_map<int, int>> &temp_dominate_table_map,
    std::unordered_map<int, int> &temp_del_p_dom_by_map)
{
    // store the distance between point to point on post height map
    std::unordered_map<int, double> dist_p_to_p_post_hm;
    dist_p_to_p_post_hm.clear();

    // store the point that we need to check distance between them
    std::vector<int> check_dist_p_index;
    check_dist_p_index.clear();

    int count = 0;
    for (int i = 0; i < org_height_map->hm_points().size(); i++)
    {
        if (temp_del_p_dom_by_map.count(i) == 0)
        {
            check_dist_p_index.push_back(i);
            count++;
        }
    }

    height_map_geodesic::HeightMapGeodesicAlgorithmDijkstra post_algorithm(new_height_map);
    double const distance_limit = height_map_geodesic::INFIN;

    std::vector<height_map_geodesic::PathPoint> post_one_source_v_list;
    std::vector<height_map_geodesic::PathPoint> post_destinations_v_list;

    // the point that we need to calculate the real pairwise distance on post height map
    std::unordered_map<int, int> post_run_alg_index_map;
    post_run_alg_index_map.clear();

    for (int i = 0; i < check_dist_p_index.size(); i++)
    {
        if (temp_del_p_dom_by_map.count(check_dist_p_index[i]) != 0)
        {
            assert(false);
            int dom_p_of_curr_del_p = temp_del_p_dom_by_map[check_dist_p_index[i]];
            for (auto ite : new_height_map->hm_points()[dom_p_of_curr_del_p].adj_hm_points_dist())
            {
                if (post_run_alg_index_map.count(ite.first) == 0)
                {
                    post_run_alg_index_map[ite.first] = ite.first;
                }
            }
        }
        else
        {
            post_run_alg_index_map[check_dist_p_index[i]] = check_dist_p_index[i];
        }
    }

    std::unordered_map<int, int> calculated_src_map;
    calculated_src_map.clear();

    // checking with source as non-deleted point
    for (int i = 0; i < check_dist_p_index.size(); i++)
    {
        // if this point has been deleted, we skip it now
        // we just focus on the non-deleted point
        if (temp_del_p_dom_by_map.count(check_dist_p_index[i]) != 0)
        {
            continue;
        }

        calculated_src_map[check_dist_p_index[i]] = check_dist_p_index[i];

        post_one_source_v_list.clear();
        post_destinations_v_list.clear();

        post_one_source_v_list.push_back(height_map_geodesic::PathPoint(&new_height_map->hm_points()[check_dist_p_index[i]]));

        for (auto ite : post_run_alg_index_map)
        {
            int x_post = check_dist_p_index[i];
            int y_post = ite.first;
            int x_y_post;
            if (x_post > y_post)
            {
                int temp_post = y_post;
                y_post = x_post;
                x_post = temp_post;
            }
            hash_function_two_keys_to_one_key_int(new_height_map->hm_points().size(), x_post, y_post, x_y_post);
            if (dist_p_to_p_post_hm.count(x_y_post) == 0)
            {
                assert(ite.first < org_height_map->hm_points().size());
                post_destinations_v_list.push_back(height_map_geodesic::PathPoint(&new_height_map->hm_points()[ite.first]));
            }
        }
        post_algorithm.propagate(post_one_source_v_list, &post_destinations_v_list);

        // calculate the post height map distance
        for (auto ite : post_run_alg_index_map)
        {
            int x_post = check_dist_p_index[i];
            int y_post = ite.first;
            int x_y_post;
            if (x_post > y_post)
            {
                int temp_post = y_post;
                y_post = x_post;
                x_post = temp_post;
            }
            hash_function_two_keys_to_one_key_int(new_height_map->hm_points().size(), x_post, y_post, x_y_post);
            if (dist_p_to_p_post_hm.count(x_y_post) == 0)
            {
                std::vector<height_map_geodesic::PathPoint> path_result;
                path_result.clear();
                height_map_geodesic::PathPoint p(&new_height_map->hm_points()[ite.first]);
                post_algorithm.trace_back(p, path_result);
                double simp_dist = length(path_result);
                dist_p_to_p_post_hm[x_y_post] = simp_dist;
            }
        }

        // compare the distance
        for (int j = i; j < check_dist_p_index.size(); j++)
        {

            int x_org = check_dist_p_index[i];
            int y_org = check_dist_p_index[j];
            int x_y_org;
            if (x_org > y_org)
            {
                int temp_org = y_org;
                y_org = x_org;
                x_org = temp_org;
            }
            hash_function_two_keys_to_one_key_int(org_height_map->hm_points().size(), x_org, y_org, x_y_org);
            double org_dist = all_dist_p_to_p_org_height_map[x_y_org];

            double simp_dist = 1e100;

            if (temp_del_p_dom_by_map.count(check_dist_p_index[j]) != 0)
            {
                int dom_p_of_curr_del_p = temp_del_p_dom_by_map[check_dist_p_index[j]];
                for (auto ite : new_height_map->hm_points()[dom_p_of_curr_del_p].adj_hm_points_dist())
                {
                    int x_post = check_dist_p_index[i];
                    int y_post = ite.first;
                    int x_y_post;
                    if (x_post > y_post)
                    {
                        int temp_post = y_post;
                        y_post = x_post;
                        x_post = temp_post;
                    }
                    hash_function_two_keys_to_one_key_int(new_height_map->hm_points().size(), x_post, y_post, x_y_post);
                    double x1 = new_height_map->hm_points()[ite.first].getx();
                    double y1 = new_height_map->hm_points()[ite.first].gety();
                    double z1 = new_height_map->hm_points()[ite.first].getz();
                    double x2 = new_height_map->hm_points()[check_dist_p_index[j]].getx();
                    double y2 = new_height_map->hm_points()[check_dist_p_index[j]].gety();
                    double z2 = new_height_map->hm_points()[dom_p_of_curr_del_p].getz();
                    double simp_dist2 = dist_p_to_p_post_hm[x_y_post] + cal_dist(x1, y1, z1, x2, y2, z2);
                    simp_dist = std::min(simp_dist, simp_dist2);
                }
            }
            else
            {
                int x_post = check_dist_p_index[i];
                int y_post = check_dist_p_index[j];
                int x_y_post;
                if (x_post > y_post)
                {
                    int temp_post = y_post;
                    y_post = x_post;
                    x_post = temp_post;
                }
                hash_function_two_keys_to_one_key_int(new_height_map->hm_points().size(), x_post, y_post, x_y_post);
                simp_dist = dist_p_to_p_post_hm[x_y_post];
            }

            compare_distance(org_dist, simp_dist, epsilon, satisfy);
            if (!satisfy)
            {
                return;
            }
        }
    }

    // calcualte the post height map distance with source that has not calculated, and store it
    for (auto ite : post_run_alg_index_map)
    {
        if (calculated_src_map.count(ite.first) == 0)
        {
            post_one_source_v_list.clear();
            post_destinations_v_list.clear();

            post_one_source_v_list.push_back(height_map_geodesic::PathPoint(&new_height_map->hm_points()[ite.first]));

            for (auto ite2 : post_run_alg_index_map)
            {
                int x = ite.first;
                int y = ite2.first;
                int x_y;
                if (x > y)
                {
                    int temp = y;
                    y = x;
                    x = temp;
                }
                hash_function_two_keys_to_one_key_int(new_height_map->hm_points().size(), x, y, x_y);
                if (dist_p_to_p_post_hm.count(x_y) == 0)
                {
                    post_destinations_v_list.push_back(height_map_geodesic::PathPoint(&new_height_map->hm_points()[ite2.first]));
                }
            }
            post_algorithm.propagate(post_one_source_v_list, &post_destinations_v_list);

            // calculate the post height map distance
            for (auto ite2 : post_run_alg_index_map)
            {
                int x = ite.first;
                int y = ite2.first;
                int x_y;
                if (x > y)
                {
                    int temp = y;
                    y = x;
                    x = temp;
                }
                hash_function_two_keys_to_one_key_int(new_height_map->hm_points().size(), x, y, x_y);
                if (dist_p_to_p_post_hm.count(x_y) == 0)
                {
                    std::vector<height_map_geodesic::PathPoint> path_result;
                    path_result.clear();
                    height_map_geodesic::PathPoint p(&new_height_map->hm_points()[ite2.first]);
                    post_algorithm.trace_back(p, path_result);
                    double simp_dist = length(path_result);
                    dist_p_to_p_post_hm[x_y] = simp_dist;
                }
            }
        }
    }

    // in the following, we have already calculated the distance of the point that we
    // really need to calculate the pairwise distance between then on post height map,
    // and also the distance on original height map

    // checking with source as deleted point
    for (int i = 0; i < check_dist_p_index.size(); i++)
    {
        // if this point has been not deleted, we skip it now
        // we just focus on the deleted point
        if (temp_del_p_dom_by_map.count(check_dist_p_index[i]) == 0)
        {
            continue;
        }
        int dom_p_of_curr_del_p_src = temp_del_p_dom_by_map[check_dist_p_index[i]];

        // there are two special cases here, we consider them first
        for (int j = i; j < check_dist_p_index.size(); j++)
        {
            // if both points s and t are deleted
            if (temp_del_p_dom_by_map.count(check_dist_p_index[j]) != 0)
            {
                int x_org = check_dist_p_index[i];
                int y_org = check_dist_p_index[j];
                int x_y_org;
                if (x_org > y_org)
                {
                    int temp_org = y_org;
                    y_org = x_org;
                    x_org = temp_org;
                }
                hash_function_two_keys_to_one_key_int(org_height_map->hm_points().size(), x_org, y_org, x_y_org);
                double org_dist = all_dist_p_to_p_org_height_map[x_y_org];
                double simp_dist = 1e100;

                int dom_p_of_curr_del_p_dest = temp_del_p_dom_by_map[check_dist_p_index[j]];

                // if s and t are dominated by same added center point,
                // then the distance is just Euclidean distance
                if (dom_p_of_curr_del_p_src == dom_p_of_curr_del_p_dest)
                {
                    double x1 = new_height_map->hm_points()[check_dist_p_index[i]].getx();
                    double y1 = new_height_map->hm_points()[check_dist_p_index[i]].gety();
                    double z1 = new_height_map->hm_points()[dom_p_of_curr_del_p_src].getz();
                    double x2 = new_height_map->hm_points()[check_dist_p_index[j]].getx();
                    double y2 = new_height_map->hm_points()[check_dist_p_index[j]].gety();
                    double z2 = new_height_map->hm_points()[dom_p_of_curr_del_p_dest].getz();
                    simp_dist = cal_dist(x1, y1, z1, x2, y2, z2);
                }

                // if s and t are dominated by two different added center points a and b,
                // but a and b are adjacent:
                else if (dom_p_of_curr_del_p_src != dom_p_of_curr_del_p_dest &&
                         new_height_map->hm_points()[dom_p_of_curr_del_p_src].adj_hm_points_dist().count(dom_p_of_curr_del_p_dest) != 0 &&
                         new_height_map->hm_points()[dom_p_of_curr_del_p_dest].adj_hm_points_dist().count(dom_p_of_curr_del_p_src) != 0)
                {
                    // (a) if s and t are adjacent on original height map, then the distance is just Euclidean distance
                    if (org_height_map->hm_points()[check_dist_p_index[i]].adj_hm_points_dist().count(check_dist_p_index[j]) != 0 &&
                        org_height_map->hm_points()[check_dist_p_index[j]].adj_hm_points_dist().count(check_dist_p_index[i]) != 0)
                    {
                        double x1 = new_height_map->hm_points()[check_dist_p_index[i]].getx();
                        double y1 = new_height_map->hm_points()[check_dist_p_index[i]].gety();
                        double z1 = new_height_map->hm_points()[dom_p_of_curr_del_p_src].getz();
                        double x2 = new_height_map->hm_points()[check_dist_p_index[j]].getx();
                        double y2 = new_height_map->hm_points()[check_dist_p_index[j]].gety();
                        double z2 = new_height_map->hm_points()[dom_p_of_curr_del_p_dest].getz();
                        simp_dist = cal_dist(x1, y1, z1, x2, y2, z2);
                    }
                    // (b) if s and t are not adjacent on original height map, then the distance is min (sa + at, sb + bt)
                    else if (org_height_map->hm_points()[check_dist_p_index[i]].adj_hm_points_dist().count(check_dist_p_index[j]) == 0 &&
                             org_height_map->hm_points()[check_dist_p_index[j]].adj_hm_points_dist().count(check_dist_p_index[i]) == 0)
                    {
                        // s
                        double x1 = new_height_map->hm_points()[check_dist_p_index[i]].getx();
                        double y1 = new_height_map->hm_points()[check_dist_p_index[i]].gety();
                        double z1 = new_height_map->hm_points()[dom_p_of_curr_del_p_src].getz();
                        // t
                        double x2 = new_height_map->hm_points()[check_dist_p_index[j]].getx();
                        double y2 = new_height_map->hm_points()[check_dist_p_index[j]].gety();
                        double z2 = new_height_map->hm_points()[dom_p_of_curr_del_p_dest].getz();
                        // a
                        double x3 = new_height_map->hm_points()[dom_p_of_curr_del_p_src].getx();
                        double y3 = new_height_map->hm_points()[dom_p_of_curr_del_p_src].gety();
                        double z3 = new_height_map->hm_points()[dom_p_of_curr_del_p_src].getz();
                        // b
                        double x4 = new_height_map->hm_points()[dom_p_of_curr_del_p_dest].getx();
                        double y4 = new_height_map->hm_points()[dom_p_of_curr_del_p_dest].gety();
                        double z4 = new_height_map->hm_points()[dom_p_of_curr_del_p_dest].getz();
                        simp_dist = std::min(cal_dist(x1, y1, z1, x3, y3, z3) + cal_dist(x2, y2, z2, x3, y3, z3),
                                             cal_dist(x1, y1, z1, x4, y4, z4) + cal_dist(x2, y2, z2, x4, y4, z4));
                    }
                    else
                    {
                        assert(false);
                    }
                }
                compare_distance(org_dist, simp_dist, epsilon, satisfy);
                if (!satisfy)
                {
                    return;
                }
            }
        }

        // the normal case
        for (int j = i; j < check_dist_p_index.size(); j++)
        {
            int x_org = check_dist_p_index[i];
            int y_org = check_dist_p_index[j];
            int x_y_org;
            if (x_org > y_org)
            {
                int temp_org = y_org;
                y_org = x_org;
                x_org = temp_org;
            }
            hash_function_two_keys_to_one_key_int(org_height_map->hm_points().size(), x_org, y_org, x_y_org);
            double org_dist = all_dist_p_to_p_org_height_map[x_y_org];
            double simp_dist = 1e100;

            // if both points s and t are deleted
            if (temp_del_p_dom_by_map.count(check_dist_p_index[j]) != 0)
            {
                int dom_p_of_curr_del_p_dest = temp_del_p_dom_by_map[check_dist_p_index[j]];

                // we considered this special case before
                if (dom_p_of_curr_del_p_src == dom_p_of_curr_del_p_dest)
                {
                    continue;
                }
                // we considered this special case before
                else if (dom_p_of_curr_del_p_src != dom_p_of_curr_del_p_dest &&
                         new_height_map->hm_points()[dom_p_of_curr_del_p_src].adj_hm_points_dist().count(dom_p_of_curr_del_p_dest) != 0 &&
                         new_height_map->hm_points()[dom_p_of_curr_del_p_dest].adj_hm_points_dist().count(dom_p_of_curr_del_p_src) != 0)
                {
                    continue;
                }

                // if s and t are dominated by two different added center points a and b,
                // but a and b are not adjacent, the distance is min (sa + \forall ab + bt)
                else
                {
                    for (auto ite : new_height_map->hm_points()[dom_p_of_curr_del_p_src].adj_hm_points_dist())
                    {
                        for (auto ite2 : new_height_map->hm_points()[dom_p_of_curr_del_p_dest].adj_hm_points_dist())
                        {
                            int x_post = ite.first;
                            int y_post = ite2.first;
                            int x_y_post;
                            if (x_post > y_post)
                            {
                                int temp_post = y_post;
                                y_post = x_post;
                                x_post = temp_post;
                            }
                            hash_function_two_keys_to_one_key_int(new_height_map->hm_points().size(), x_post, y_post, x_y_post);
                            double x1 = new_height_map->hm_points()[ite.first].getx();
                            double y1 = new_height_map->hm_points()[ite.first].gety();
                            double z1 = new_height_map->hm_points()[ite.first].getz();
                            double x2 = new_height_map->hm_points()[check_dist_p_index[i]].getx();
                            double y2 = new_height_map->hm_points()[check_dist_p_index[i]].gety();
                            double z2 = new_height_map->hm_points()[dom_p_of_curr_del_p_src].getz();
                            double x3 = new_height_map->hm_points()[ite2.first].getx();
                            double y3 = new_height_map->hm_points()[ite2.first].gety();
                            double z3 = new_height_map->hm_points()[ite2.first].getz();
                            double x4 = new_height_map->hm_points()[check_dist_p_index[j]].getx();
                            double y4 = new_height_map->hm_points()[check_dist_p_index[j]].gety();
                            double z4 = new_height_map->hm_points()[dom_p_of_curr_del_p_dest].getz();
                            double simp_dist2 = dist_p_to_p_post_hm[x_y_post] + cal_dist(x1, y1, z1, x2, y2, z2) + cal_dist(x3, y3, z3, x4, y4, z4);
                            simp_dist = std::min(simp_dist, simp_dist2);
                        }
                    }
                }
            }
            else
            {
                for (auto ite : new_height_map->hm_points()[dom_p_of_curr_del_p_src].adj_hm_points_dist())
                {
                    int x_post = ite.first;
                    int y_post = check_dist_p_index[j];
                    int x_y_post;
                    if (x_post > y_post)
                    {
                        int temp_post = y_post;
                        y_post = x_post;
                        x_post = temp_post;
                    }
                    hash_function_two_keys_to_one_key_int(new_height_map->hm_points().size(), x_post, y_post, x_y_post);
                    double x1 = new_height_map->hm_points()[ite.first].getx();
                    double y1 = new_height_map->hm_points()[ite.first].gety();
                    double z1 = new_height_map->hm_points()[ite.first].getz();
                    double x2 = new_height_map->hm_points()[check_dist_p_index[i]].getx();
                    double y2 = new_height_map->hm_points()[check_dist_p_index[i]].gety();
                    double z2 = new_height_map->hm_points()[dom_p_of_curr_del_p_src].getz();
                    double simp_dist2 = dist_p_to_p_post_hm[x_y_post] + cal_dist(x1, y1, z1, x2, y2, z2);
                    simp_dist = std::min(simp_dist, simp_dist2);
                }
            }
            compare_distance(org_dist, simp_dist, epsilon, satisfy);
            if (!satisfy)
            {
                return;
            }
        }
    }
}

void height_map_cal_org_dist_p_to_p(
    height_map_geodesic::HeightMap *org_height_map,
    std::unordered_map<int, double> &all_dist_p_to_p_org_height_map)
{
    height_map_geodesic::HeightMapGeodesicAlgorithmDijkstra org_algorithm(org_height_map);
    double const distance_limit = height_map_geodesic::INFIN;
    std::vector<height_map_geodesic::PathPoint> one_source_p_list;

    for (int i = 0; i < org_height_map->hm_points().size(); i++)
    {
        one_source_p_list.clear();
        one_source_p_list.push_back(height_map_geodesic::PathPoint(&org_height_map->hm_points()[i]));

        org_algorithm.propagate(one_source_p_list, distance_limit, 1);

        for (int j = i; j < org_height_map->hm_points().size(); j++)
        {
            std::vector<height_map_geodesic::PathPoint> path_result;
            path_result.clear();
            height_map_geodesic::PathPoint one_dest_p(&org_height_map->hm_points()[j]);
            org_algorithm.trace_back(one_dest_p, path_result);
            double distance = length(path_result);

            int i_j;
            hash_function_two_keys_to_one_key_int(org_height_map->hm_points().size(), i, j, i_j);
            all_dist_p_to_p_org_height_map[i_j] = distance;
        }
    }
}

void simp_height_map_query(
    height_map_geodesic::HeightMap *org_height_map,
    height_map_geodesic::HeightMap *new_height_map,
    std::unordered_map<int, std::unordered_map<int, int>> &dominate_table_map,
    std::unordered_map<int, int> &del_p_dom_by_map,
    int source_index, int destination_index,
    double &query_time, double &query_memory_usage, double &distance_result,
    std::vector<height_map_geodesic::PathPoint> &path_result,
    int org_one_lqt1_two_lqt2_three_ls_four_lst_five_ds_six_point_seven)
{
    auto start_query_time = std::chrono::high_resolution_clock::now();

    int dom_p_of_src = -1;
    int dom_p_of_dest = -1;
    distance_result = 1e100;

    if (del_p_dom_by_map.count(source_index) != 0)
    {
        dom_p_of_src = del_p_dom_by_map[source_index];
    }
    if (del_p_dom_by_map.count(destination_index) != 0)
    {
        dom_p_of_dest = del_p_dom_by_map[destination_index];
    }

    height_map_geodesic::HeightMapGeodesicAlgorithmDijkstra post_algorithm(new_height_map);

    if (dom_p_of_src == -1 && dom_p_of_dest == -1)
    {
        height_map_geodesic::PathPoint source(&new_height_map->hm_points()[source_index]);
        height_map_geodesic::PathPoint destination(&new_height_map->hm_points()[destination_index]);
        std::vector<height_map_geodesic::PathPoint> one_source_list(1, source);
        std::vector<height_map_geodesic::PathPoint> one_destination_list(1, destination);
        for (int k = 0; k < iternation_num(org_one_lqt1_two_lqt2_three_ls_four_lst_five_ds_six_point_seven, 1); k++)
        {
            post_algorithm.propagate(one_source_list, &one_destination_list);
        }
        post_algorithm.trace_back(destination, path_result);
        distance_result = length(path_result);
    }
    else if (dom_p_of_src == -1 && dom_p_of_dest != -1)
    {
        // store the path between point to point on post height map
        std::unordered_map<int, std::vector<height_map_geodesic::PathPoint>> path_p_to_p_post_hm;
        path_p_to_p_post_hm.clear();

        std::vector<height_map_geodesic::PathPoint> one_source_list;
        std::vector<height_map_geodesic::PathPoint> destinations_list;
        one_source_list.clear();
        destinations_list.clear();

        one_source_list.push_back(height_map_geodesic::PathPoint(&new_height_map->hm_points()[source_index]));
        for (auto ite : new_height_map->hm_points()[dom_p_of_dest].adj_hm_points_dist())
        {
            destinations_list.push_back(height_map_geodesic::PathPoint(&new_height_map->hm_points()[ite.first]));
        }
        for (int k = 0; k < iternation_num(org_one_lqt1_two_lqt2_three_ls_four_lst_five_ds_six_point_seven, 1); k++)
        {
            post_algorithm.propagate(one_source_list, &destinations_list);
        }
        for (auto ite : new_height_map->hm_points()[dom_p_of_dest].adj_hm_points_dist())
        {
            int x = source_index;
            int y = ite.first;
            int x_y;
            if (x > y)
            {
                int temp = y;
                y = x;
                x = temp;
            }
            hash_function_two_keys_to_one_key_int(new_height_map->hm_points().size(), x, y, x_y);
            if (path_p_to_p_post_hm.count(x_y) == 0)
            {
                std::vector<height_map_geodesic::PathPoint> path;
                path.clear();
                height_map_geodesic::PathPoint p(&new_height_map->hm_points()[ite.first]);
                post_algorithm.trace_back(p, path);
                path_p_to_p_post_hm[x_y] = path;
            }
        }

        for (auto ite : new_height_map->hm_points()[dom_p_of_dest].adj_hm_points_dist())
        {
            int x = source_index;
            int y = ite.first;
            int x_y;
            if (x > y)
            {
                int temp = y;
                y = x;
                x = temp;
            }
            hash_function_two_keys_to_one_key_int(new_height_map->hm_points().size(), x, y, x_y);

            std::vector<height_map_geodesic::PathPoint> path = path_p_to_p_post_hm[x_y];

            new_height_map->hm_points()[destination_index].z() = new_height_map->hm_points()[dom_p_of_dest].getz();
            height_map_geodesic::PathPoint real_destination(&new_height_map->hm_points()[destination_index]);

            double x1 = new_height_map->hm_points()[ite.first].getx();
            double y1 = new_height_map->hm_points()[ite.first].gety();
            double z1 = new_height_map->hm_points()[ite.first].getz();
            double x2 = new_height_map->hm_points()[destination_index].getx();
            double y2 = new_height_map->hm_points()[destination_index].gety();
            double z2 = new_height_map->hm_points()[dom_p_of_dest].getz();
            double dest_short_dist = cal_dist(x1, y1, z1, x2, y2, z2);

            if (distance_result > length(path) + dest_short_dist)
            {
                distance_result = length(path) + dest_short_dist;
                path.insert(path.begin(), real_destination);
                path_result = path;
            }
        }
    }
    else if (dom_p_of_src != -1 && dom_p_of_dest == -1)
    {
        // store the path between point to point on post height map
        std::unordered_map<int, std::vector<height_map_geodesic::PathPoint>> path_p_to_p_post_hm;
        path_p_to_p_post_hm.clear();

        std::vector<height_map_geodesic::PathPoint> one_source_list;
        std::vector<height_map_geodesic::PathPoint> destinations_list;
        one_source_list.clear();
        destinations_list.clear();

        one_source_list.push_back(height_map_geodesic::PathPoint(&new_height_map->hm_points()[destination_index]));
        for (auto ite : new_height_map->hm_points()[dom_p_of_src].adj_hm_points_dist())
        {
            destinations_list.push_back(height_map_geodesic::PathPoint(&new_height_map->hm_points()[ite.first]));
        }
        for (int k = 0; k < iternation_num(org_one_lqt1_two_lqt2_three_ls_four_lst_five_ds_six_point_seven, 1); k++)
        {
            post_algorithm.propagate(one_source_list, &destinations_list);
        }
        for (auto ite : new_height_map->hm_points()[dom_p_of_src].adj_hm_points_dist())
        {
            int x = ite.first;
            int y = destination_index;
            int x_y;
            if (x > y)
            {
                int temp = y;
                y = x;
                x = temp;
            }
            hash_function_two_keys_to_one_key_int(new_height_map->hm_points().size(), x, y, x_y);
            if (path_p_to_p_post_hm.count(x_y) == 0)
            {
                std::vector<height_map_geodesic::PathPoint> path;
                path.clear();
                height_map_geodesic::PathPoint p(&new_height_map->hm_points()[ite.first]);
                post_algorithm.trace_back(p, path);
                path_p_to_p_post_hm[x_y] = path;
            }
        }

        for (auto ite : new_height_map->hm_points()[dom_p_of_src].adj_hm_points_dist())
        {
            int x = ite.first;
            int y = destination_index;
            int x_y;
            if (x > y)
            {
                int temp = y;
                y = x;
                x = temp;
            }
            hash_function_two_keys_to_one_key_int(new_height_map->hm_points().size(), x, y, x_y);

            std::vector<height_map_geodesic::PathPoint> path = path_p_to_p_post_hm[x_y];

            new_height_map->hm_points()[source_index].z() = new_height_map->hm_points()[dom_p_of_src].getz();
            height_map_geodesic::PathPoint real_source(&new_height_map->hm_points()[source_index]);

            double x1 = new_height_map->hm_points()[source_index].getx();
            double y1 = new_height_map->hm_points()[source_index].gety();
            double z1 = new_height_map->hm_points()[dom_p_of_src].getz();
            double x2 = new_height_map->hm_points()[ite.first].getx();
            double y2 = new_height_map->hm_points()[ite.first].gety();
            double z2 = new_height_map->hm_points()[ite.first].getz();
            double src_short_dist = cal_dist(x1, y1, z1, x2, y2, z2);

            if (distance_result > length(path) + src_short_dist)
            {
                distance_result = length(path) + src_short_dist;
                path.insert(path.begin(), real_source);
                std::reverse(path.begin(), path.end());
                path_result = path;
            }
        }
    }
    else
    {
        // if s and t are dominated by same added center point,
        // then the distance is just Euclidean distance
        if (dom_p_of_src == dom_p_of_dest)
        {
            double x1 = new_height_map->hm_points()[source_index].getx();
            double y1 = new_height_map->hm_points()[source_index].gety();
            double z1 = new_height_map->hm_points()[dom_p_of_src].getz();
            double x2 = new_height_map->hm_points()[destination_index].getx();
            double y2 = new_height_map->hm_points()[destination_index].gety();
            double z2 = new_height_map->hm_points()[dom_p_of_dest].getz();
            new_height_map->hm_points()[source_index].z() = new_height_map->hm_points()[dom_p_of_src].getz();
            height_map_geodesic::PathPoint real_source(&new_height_map->hm_points()[source_index]);
            new_height_map->hm_points()[destination_index].z() = new_height_map->hm_points()[dom_p_of_dest].getz();
            height_map_geodesic::PathPoint real_destination(&new_height_map->hm_points()[destination_index]);
            path_result.push_back(real_destination);
            path_result.push_back(real_source);
            distance_result = cal_dist(x1, y1, z1, x2, y2, z2);
        }

        // if s and t are dominated by two different added center points a and b,
        // but a and b are adjacent:
        else if (dom_p_of_src != dom_p_of_dest &&
                 new_height_map->hm_points()[dom_p_of_src].adj_hm_points_dist().count(dom_p_of_dest) != 0 &&
                 new_height_map->hm_points()[dom_p_of_dest].adj_hm_points_dist().count(dom_p_of_src) != 0)
        {
            // (a) if s and t are adjacent on original height map, then the distance is just Euclidean distance
            if (org_height_map->hm_points()[source_index].adj_hm_points_dist().count(destination_index) != 0 &&
                org_height_map->hm_points()[destination_index].adj_hm_points_dist().count(source_index) != 0)
            {
                double x1 = new_height_map->hm_points()[source_index].getx();
                double y1 = new_height_map->hm_points()[source_index].gety();
                double z1 = new_height_map->hm_points()[dom_p_of_src].getz();
                double x2 = new_height_map->hm_points()[destination_index].getx();
                double y2 = new_height_map->hm_points()[destination_index].gety();
                double z2 = new_height_map->hm_points()[dom_p_of_dest].getz();
                new_height_map->hm_points()[source_index].z() = new_height_map->hm_points()[dom_p_of_src].getz();
                height_map_geodesic::PathPoint real_source(&new_height_map->hm_points()[source_index]);
                new_height_map->hm_points()[destination_index].z() = new_height_map->hm_points()[dom_p_of_dest].getz();
                height_map_geodesic::PathPoint real_destination(&new_height_map->hm_points()[destination_index]);
                path_result.push_back(real_destination);
                path_result.push_back(real_source);
                distance_result = cal_dist(x1, y1, z1, x2, y2, z2);
            }
            // (b) if s and t are not adjacent on original height map, then the distance is min (sa + at, sb + bt)
            else if (org_height_map->hm_points()[source_index].adj_hm_points_dist().count(destination_index) == 0 &&
                     org_height_map->hm_points()[destination_index].adj_hm_points_dist().count(source_index) == 0)
            {
                // s
                double x1 = new_height_map->hm_points()[source_index].getx();
                double y1 = new_height_map->hm_points()[source_index].gety();
                double z1 = new_height_map->hm_points()[dom_p_of_src].getz();
                // t
                double x2 = new_height_map->hm_points()[destination_index].getx();
                double y2 = new_height_map->hm_points()[destination_index].gety();
                double z2 = new_height_map->hm_points()[dom_p_of_dest].getz();
                // a
                double x3 = new_height_map->hm_points()[dom_p_of_src].getx();
                double y3 = new_height_map->hm_points()[dom_p_of_src].gety();
                double z3 = new_height_map->hm_points()[dom_p_of_src].getz();
                // b
                double x4 = new_height_map->hm_points()[dom_p_of_dest].getx();
                double y4 = new_height_map->hm_points()[dom_p_of_dest].gety();
                double z4 = new_height_map->hm_points()[dom_p_of_dest].getz();

                new_height_map->hm_points()[source_index].z() = new_height_map->hm_points()[dom_p_of_src].getz();
                height_map_geodesic::PathPoint real_source(&new_height_map->hm_points()[source_index]);
                new_height_map->hm_points()[destination_index].z() = new_height_map->hm_points()[dom_p_of_dest].getz();
                height_map_geodesic::PathPoint real_destination(&new_height_map->hm_points()[destination_index]);

                double pass_a_dist = cal_dist(x1, y1, z1, x3, y3, z3) + cal_dist(x2, y2, z2, x3, y3, z3);
                double pass_b_dist = cal_dist(x1, y1, z1, x4, y4, z4) + cal_dist(x2, y2, z2, x4, y4, z4);
                if (pass_a_dist < pass_b_dist)
                {
                    height_map_geodesic::PathPoint middle(&new_height_map->hm_points()[dom_p_of_src]);
                    path_result.push_back(real_destination);
                    path_result.push_back(middle);
                    path_result.push_back(real_source);
                    distance_result = pass_a_dist;
                }
                else
                {
                    height_map_geodesic::PathPoint middle(&new_height_map->hm_points()[dom_p_of_dest]);
                    path_result.push_back(real_destination);
                    path_result.push_back(middle);
                    path_result.push_back(real_source);
                    distance_result = pass_b_dist;
                }
            }
            else
            {
                assert(false);
            }
        }

        // if s and t are dominated by two different added center points a and b,
        // but a and b are not adjacent, the distance is min (sa + \forall ab + bt)
        else
        {
            // store the path between point to point on post height map
            std::unordered_map<int, std::vector<height_map_geodesic::PathPoint>> path_p_to_p_post_hm;
            path_p_to_p_post_hm.clear();

            for (auto ite : new_height_map->hm_points()[dom_p_of_src].adj_hm_points_dist())
            {
                std::vector<height_map_geodesic::PathPoint> one_source_list;
                std::vector<height_map_geodesic::PathPoint> destinations_list;
                one_source_list.clear();
                destinations_list.clear();

                one_source_list.push_back(height_map_geodesic::PathPoint(&new_height_map->hm_points()[ite.first]));
                for (auto ite2 : new_height_map->hm_points()[dom_p_of_dest].adj_hm_points_dist())
                {
                    destinations_list.push_back(height_map_geodesic::PathPoint(&new_height_map->hm_points()[ite2.first]));
                }
                for (int k = 0; k < iternation_num(org_one_lqt1_two_lqt2_three_ls_four_lst_five_ds_six_point_seven, 1); k++)
                {
                    post_algorithm.propagate(one_source_list, &destinations_list);
                }
                for (auto ite2 : new_height_map->hm_points()[dom_p_of_dest].adj_hm_points_dist())
                {
                    int x = ite.first;
                    int y = ite2.first;
                    int x_y;
                    if (x > y)
                    {
                        int temp = y;
                        y = x;
                        x = temp;
                    }
                    hash_function_two_keys_to_one_key_int(new_height_map->hm_points().size(), x, y, x_y);
                    if (path_p_to_p_post_hm.count(x_y) == 0)
                    {
                        std::vector<height_map_geodesic::PathPoint> path;
                        path.clear();
                        height_map_geodesic::PathPoint p(&new_height_map->hm_points()[ite2.first]);
                        post_algorithm.trace_back(p, path);
                        path_p_to_p_post_hm[x_y] = path;
                    }
                }
            }

            for (auto ite : new_height_map->hm_points()[dom_p_of_src].adj_hm_points_dist())
            {
                for (auto ite2 : new_height_map->hm_points()[dom_p_of_dest].adj_hm_points_dist())
                {
                    int x = ite.first;
                    int y = ite2.first;
                    int x_y;
                    if (x > y)
                    {
                        int temp = y;
                        y = x;
                        x = temp;
                    }
                    hash_function_two_keys_to_one_key_int(new_height_map->hm_points().size(), x, y, x_y);

                    std::vector<height_map_geodesic::PathPoint> path = path_p_to_p_post_hm[x_y];

                    new_height_map->hm_points()[source_index].z() = new_height_map->hm_points()[dom_p_of_src].getz();
                    height_map_geodesic::PathPoint real_source(&new_height_map->hm_points()[source_index]);
                    new_height_map->hm_points()[destination_index].z() = new_height_map->hm_points()[dom_p_of_dest].getz();
                    height_map_geodesic::PathPoint real_destination(&new_height_map->hm_points()[destination_index]);

                    double x1 = new_height_map->hm_points()[source_index].getx();
                    double y1 = new_height_map->hm_points()[source_index].gety();
                    double z1 = new_height_map->hm_points()[dom_p_of_src].getz();
                    double x2 = new_height_map->hm_points()[ite.first].getx();
                    double y2 = new_height_map->hm_points()[ite.first].gety();
                    double z2 = new_height_map->hm_points()[ite.first].getz();
                    double src_short_dist = cal_dist(x1, y1, z1, x2, y2, z2);

                    double x3 = new_height_map->hm_points()[ite2.first].getx();
                    double y3 = new_height_map->hm_points()[ite2.first].gety();
                    double z3 = new_height_map->hm_points()[ite2.first].getz();
                    double x4 = new_height_map->hm_points()[destination_index].getx();
                    double y4 = new_height_map->hm_points()[destination_index].gety();
                    double z4 = new_height_map->hm_points()[dom_p_of_dest].getz();
                    double dest_short_dist = cal_dist(x3, y3, z3, x4, y4, z4);

                    if (distance_result > length(path) + src_short_dist + dest_short_dist)
                    {
                        distance_result = length(path) + src_short_dist + dest_short_dist;
                        path.insert(path.begin(), real_destination);
                        path.push_back(real_source);
                        path_result = path;
                    }
                }
            }
        }
    }

    distance_result = round(distance_result * 1000000000.0) / 1000000000.0;
    int total_node_size;
    size_t node_size;
    post_algorithm.get_memory(total_node_size, node_size);
    query_memory_usage += (total_node_size - del_p_dom_by_map.size()) * node_size + path_result.size() * sizeof(height_map_geodesic::PathPoint) + sizeof(double);

    auto stop_query_time = std::chrono::high_resolution_clock::now();
    auto duration_query_time = std::chrono::duration_cast<std::chrono::microseconds>(stop_query_time - start_query_time);
    query_time = duration_query_time.count();
    query_time /= 1000;
}

void simp_height_map_query_ds(
    std::unordered_map<std::pair<int, int>, std::vector<height_map_geodesic::PathPoint>, boost::hash<std::pair<int, int>>> &all_path,
    std::unordered_map<std::pair<int, int>, double, boost::hash<std::pair<int, int>>> &all_dist,
    int source_index, int destination_index,
    double &query_time, double &query_memory_usage, double &distance_result,
    std::vector<height_map_geodesic::PathPoint> &path_result)
{
    auto start_query_time = std::chrono::high_resolution_clock::now();

    distance_result = all_dist[std::make_pair(source_index, destination_index)];
    path_result = all_path[std::make_pair(source_index, destination_index)];

    distance_result = round(distance_result * 1000000000.0) / 1000000000.0;
    query_memory_usage += path_result.size() * sizeof(height_map_geodesic::PathPoint) + sizeof(double);

    auto stop_query_time = std::chrono::high_resolution_clock::now();
    auto duration_query_time = std::chrono::duration_cast<std::chrono::microseconds>(stop_query_time - start_query_time);
    query_time = duration_query_time.count();
    query_time /= 1000;
}

void simp_height_map_knn_and_range_query(
    height_map_geodesic::HeightMap *org_height_map,
    height_map_geodesic::HeightMap *new_height_map,
    std::unordered_map<int, std::unordered_map<int, int>> &dominate_table_map,
    std::unordered_map<int, int> &del_p_dom_by_map,
    int source_index,
    double &knn_query_time, double &range_query_time,
    int knn_and_range_query_obj_num, int k_value, double range, double d_value,
    std::vector<int> &knn_list,
    std::vector<int> &range_list,
    int org_one_lqt1_two_lqt2_three_ls_four_lst_five_ds_six_point_seven)
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

        int dom_p_of_src = -1;
        int dom_p_of_dest = -1;
        double distance_result = 1e100;
        std::vector<height_map_geodesic::PathPoint> path_result;
        path_result.clear();

        if (del_p_dom_by_map.count(source_index) != 0)
        {
            dom_p_of_src = del_p_dom_by_map[source_index];
        }
        if (del_p_dom_by_map.count(destination_index) != 0)
        {
            dom_p_of_dest = del_p_dom_by_map[destination_index];
        }

        height_map_geodesic::HeightMapGeodesicAlgorithmDijkstra post_algorithm(new_height_map);

        if (dom_p_of_src == -1 && dom_p_of_dest == -1)
        {
            height_map_geodesic::PathPoint source(&new_height_map->hm_points()[source_index]);
            height_map_geodesic::PathPoint destination(&new_height_map->hm_points()[destination_index]);
            std::vector<height_map_geodesic::PathPoint> one_source_list(1, source);
            std::vector<height_map_geodesic::PathPoint> one_destination_list(1, destination);
            for (int k = 0; k < iternation_num(org_one_lqt1_two_lqt2_three_ls_four_lst_five_ds_six_point_seven, 2); k++)
            {
                post_algorithm.propagate(one_source_list, &one_destination_list);
            }
            post_algorithm.trace_back(destination, path_result);
            distance_result = length(path_result);
        }
        else if (dom_p_of_src == -1 && dom_p_of_dest != -1)
        {
            // store the path between point to point on post height map
            std::unordered_map<int, std::vector<height_map_geodesic::PathPoint>> path_p_to_p_post_hm;
            path_p_to_p_post_hm.clear();

            std::vector<height_map_geodesic::PathPoint> one_source_list;
            std::vector<height_map_geodesic::PathPoint> destinations_list;
            one_source_list.clear();
            destinations_list.clear();

            one_source_list.push_back(height_map_geodesic::PathPoint(&new_height_map->hm_points()[source_index]));
            for (auto ite : new_height_map->hm_points()[dom_p_of_dest].adj_hm_points_dist())
            {
                destinations_list.push_back(height_map_geodesic::PathPoint(&new_height_map->hm_points()[ite.first]));
            }
            for (int k = 0; k < iternation_num(org_one_lqt1_two_lqt2_three_ls_four_lst_five_ds_six_point_seven, 2); k++)
            {
                post_algorithm.propagate(one_source_list, &destinations_list);
            }
            for (auto ite : new_height_map->hm_points()[dom_p_of_dest].adj_hm_points_dist())
            {
                int x = source_index;
                int y = ite.first;
                int x_y;
                if (x > y)
                {
                    int temp = y;
                    y = x;
                    x = temp;
                }
                hash_function_two_keys_to_one_key_int(new_height_map->hm_points().size(), x, y, x_y);
                if (path_p_to_p_post_hm.count(x_y) == 0)
                {
                    std::vector<height_map_geodesic::PathPoint> path;
                    path.clear();
                    height_map_geodesic::PathPoint p(&new_height_map->hm_points()[ite.first]);
                    post_algorithm.trace_back(p, path);
                    path_p_to_p_post_hm[x_y] = path;
                }
            }

            for (auto ite : new_height_map->hm_points()[dom_p_of_dest].adj_hm_points_dist())
            {
                int x = source_index;
                int y = ite.first;
                int x_y;
                if (x > y)
                {
                    int temp = y;
                    y = x;
                    x = temp;
                }
                hash_function_two_keys_to_one_key_int(new_height_map->hm_points().size(), x, y, x_y);

                std::vector<height_map_geodesic::PathPoint> path = path_p_to_p_post_hm[x_y];

                new_height_map->hm_points()[destination_index].z() = new_height_map->hm_points()[dom_p_of_dest].getz();
                height_map_geodesic::PathPoint real_destination(&new_height_map->hm_points()[destination_index]);

                double x1 = new_height_map->hm_points()[ite.first].getx();
                double y1 = new_height_map->hm_points()[ite.first].gety();
                double z1 = new_height_map->hm_points()[ite.first].getz();
                double x2 = new_height_map->hm_points()[destination_index].getx();
                double y2 = new_height_map->hm_points()[destination_index].gety();
                double z2 = new_height_map->hm_points()[dom_p_of_dest].getz();
                double dest_short_dist = cal_dist(x1, y1, z1, x2, y2, z2);

                if (distance_result > length(path) + dest_short_dist)
                {
                    distance_result = length(path) + dest_short_dist;
                    path.insert(path.begin(), real_destination);
                    path_result = path;
                }
            }
        }
        else if (dom_p_of_src != -1 && dom_p_of_dest == -1)
        {
            // store the path between point to point on post height map
            std::unordered_map<int, std::vector<height_map_geodesic::PathPoint>> path_p_to_p_post_hm;
            path_p_to_p_post_hm.clear();

            std::vector<height_map_geodesic::PathPoint> one_source_list;
            std::vector<height_map_geodesic::PathPoint> destinations_list;
            one_source_list.clear();
            destinations_list.clear();

            one_source_list.push_back(height_map_geodesic::PathPoint(&new_height_map->hm_points()[destination_index]));
            for (auto ite : new_height_map->hm_points()[dom_p_of_src].adj_hm_points_dist())
            {
                destinations_list.push_back(height_map_geodesic::PathPoint(&new_height_map->hm_points()[ite.first]));
            }
            for (int k = 0; k < iternation_num(org_one_lqt1_two_lqt2_three_ls_four_lst_five_ds_six_point_seven, 2); k++)
            {
                post_algorithm.propagate(one_source_list, &destinations_list);
            }
            for (auto ite : new_height_map->hm_points()[dom_p_of_src].adj_hm_points_dist())
            {
                int x = ite.first;
                int y = destination_index;
                int x_y;
                if (x > y)
                {
                    int temp = y;
                    y = x;
                    x = temp;
                }
                hash_function_two_keys_to_one_key_int(new_height_map->hm_points().size(), x, y, x_y);
                if (path_p_to_p_post_hm.count(x_y) == 0)
                {
                    std::vector<height_map_geodesic::PathPoint> path;
                    path.clear();
                    height_map_geodesic::PathPoint p(&new_height_map->hm_points()[ite.first]);
                    post_algorithm.trace_back(p, path);
                    path_p_to_p_post_hm[x_y] = path;
                }
            }

            for (auto ite : new_height_map->hm_points()[dom_p_of_src].adj_hm_points_dist())
            {
                int x = ite.first;
                int y = destination_index;
                int x_y;
                if (x > y)
                {
                    int temp = y;
                    y = x;
                    x = temp;
                }
                hash_function_two_keys_to_one_key_int(new_height_map->hm_points().size(), x, y, x_y);

                std::vector<height_map_geodesic::PathPoint> path = path_p_to_p_post_hm[x_y];

                new_height_map->hm_points()[source_index].z() = new_height_map->hm_points()[dom_p_of_src].getz();
                height_map_geodesic::PathPoint real_source(&new_height_map->hm_points()[source_index]);

                double x1 = new_height_map->hm_points()[source_index].getx();
                double y1 = new_height_map->hm_points()[source_index].gety();
                double z1 = new_height_map->hm_points()[dom_p_of_src].getz();
                double x2 = new_height_map->hm_points()[ite.first].getx();
                double y2 = new_height_map->hm_points()[ite.first].gety();
                double z2 = new_height_map->hm_points()[ite.first].getz();
                double src_short_dist = cal_dist(x1, y1, z1, x2, y2, z2);

                if (distance_result > length(path) + src_short_dist)
                {
                    distance_result = length(path) + src_short_dist;
                    path.insert(path.begin(), real_source);
                    std::reverse(path.begin(), path.end());
                    path_result = path;
                }
            }
        }
        else
        {
            // if s and t are dominated by same added center point,
            // then the distance is just Euclidean distance
            if (dom_p_of_src == dom_p_of_dest)
            {
                double x1 = new_height_map->hm_points()[source_index].getx();
                double y1 = new_height_map->hm_points()[source_index].gety();
                double z1 = new_height_map->hm_points()[dom_p_of_src].getz();
                double x2 = new_height_map->hm_points()[destination_index].getx();
                double y2 = new_height_map->hm_points()[destination_index].gety();
                double z2 = new_height_map->hm_points()[dom_p_of_dest].getz();
                new_height_map->hm_points()[source_index].z() = new_height_map->hm_points()[dom_p_of_src].getz();
                height_map_geodesic::PathPoint real_source(&new_height_map->hm_points()[source_index]);
                new_height_map->hm_points()[destination_index].z() = new_height_map->hm_points()[dom_p_of_dest].getz();
                height_map_geodesic::PathPoint real_destination(&new_height_map->hm_points()[destination_index]);
                path_result.push_back(real_destination);
                path_result.push_back(real_source);
                distance_result = cal_dist(x1, y1, z1, x2, y2, z2);
            }

            // if s and t are dominated by two different added center points a and b,
            // but a and b are adjacent:
            else if (dom_p_of_src != dom_p_of_dest &&
                     new_height_map->hm_points()[dom_p_of_src].adj_hm_points_dist().count(dom_p_of_dest) != 0 &&
                     new_height_map->hm_points()[dom_p_of_dest].adj_hm_points_dist().count(dom_p_of_src) != 0)
            {
                // (a) if s and t are adjacent on original height map, then the distance is just Euclidean distance
                if (org_height_map->hm_points()[source_index].adj_hm_points_dist().count(destination_index) != 0 &&
                    org_height_map->hm_points()[destination_index].adj_hm_points_dist().count(source_index) != 0)
                {
                    double x1 = new_height_map->hm_points()[source_index].getx();
                    double y1 = new_height_map->hm_points()[source_index].gety();
                    double z1 = new_height_map->hm_points()[dom_p_of_src].getz();
                    double x2 = new_height_map->hm_points()[destination_index].getx();
                    double y2 = new_height_map->hm_points()[destination_index].gety();
                    double z2 = new_height_map->hm_points()[dom_p_of_dest].getz();
                    new_height_map->hm_points()[source_index].z() = new_height_map->hm_points()[dom_p_of_src].getz();
                    height_map_geodesic::PathPoint real_source(&new_height_map->hm_points()[source_index]);
                    new_height_map->hm_points()[destination_index].z() = new_height_map->hm_points()[dom_p_of_dest].getz();
                    height_map_geodesic::PathPoint real_destination(&new_height_map->hm_points()[destination_index]);
                    path_result.push_back(real_destination);
                    path_result.push_back(real_source);
                    distance_result = cal_dist(x1, y1, z1, x2, y2, z2);
                }
                // (b) if s and t are not adjacent on original height map, then the distance is min (sa + at, sb + bt)
                else if (org_height_map->hm_points()[source_index].adj_hm_points_dist().count(destination_index) == 0 &&
                         org_height_map->hm_points()[destination_index].adj_hm_points_dist().count(source_index) == 0)
                {
                    // s
                    double x1 = new_height_map->hm_points()[source_index].getx();
                    double y1 = new_height_map->hm_points()[source_index].gety();
                    double z1 = new_height_map->hm_points()[dom_p_of_src].getz();
                    // t
                    double x2 = new_height_map->hm_points()[destination_index].getx();
                    double y2 = new_height_map->hm_points()[destination_index].gety();
                    double z2 = new_height_map->hm_points()[dom_p_of_dest].getz();
                    // a
                    double x3 = new_height_map->hm_points()[dom_p_of_src].getx();
                    double y3 = new_height_map->hm_points()[dom_p_of_src].gety();
                    double z3 = new_height_map->hm_points()[dom_p_of_src].getz();
                    // b
                    double x4 = new_height_map->hm_points()[dom_p_of_dest].getx();
                    double y4 = new_height_map->hm_points()[dom_p_of_dest].gety();
                    double z4 = new_height_map->hm_points()[dom_p_of_dest].getz();

                    new_height_map->hm_points()[source_index].z() = new_height_map->hm_points()[dom_p_of_src].getz();
                    height_map_geodesic::PathPoint real_source(&new_height_map->hm_points()[source_index]);
                    new_height_map->hm_points()[destination_index].z() = new_height_map->hm_points()[dom_p_of_dest].getz();
                    height_map_geodesic::PathPoint real_destination(&new_height_map->hm_points()[destination_index]);

                    double pass_a_dist = cal_dist(x1, y1, z1, x3, y3, z3) + cal_dist(x2, y2, z2, x3, y3, z3);
                    double pass_b_dist = cal_dist(x1, y1, z1, x4, y4, z4) + cal_dist(x2, y2, z2, x4, y4, z4);
                    if (pass_a_dist < pass_b_dist)
                    {
                        height_map_geodesic::PathPoint middle(&new_height_map->hm_points()[dom_p_of_src]);
                        path_result.push_back(real_destination);
                        path_result.push_back(middle);
                        path_result.push_back(real_source);
                        distance_result = pass_a_dist;
                    }
                    else
                    {
                        height_map_geodesic::PathPoint middle(&new_height_map->hm_points()[dom_p_of_dest]);
                        path_result.push_back(real_destination);
                        path_result.push_back(middle);
                        path_result.push_back(real_source);
                        distance_result = pass_b_dist;
                    }
                }
                else
                {
                    assert(false);
                }
            }

            // if s and t are dominated by two different added center points a and b,
            // but a and b are not adjacent, the distance is min (sa + \forall ab + bt)
            else
            {
                // store the path between point to point on post height map
                std::unordered_map<int, std::vector<height_map_geodesic::PathPoint>> path_p_to_p_post_hm;
                path_p_to_p_post_hm.clear();

                for (auto ite : new_height_map->hm_points()[dom_p_of_src].adj_hm_points_dist())
                {
                    std::vector<height_map_geodesic::PathPoint> one_source_list;
                    std::vector<height_map_geodesic::PathPoint> destinations_list;
                    one_source_list.clear();
                    destinations_list.clear();

                    one_source_list.push_back(height_map_geodesic::PathPoint(&new_height_map->hm_points()[ite.first]));
                    for (auto ite2 : new_height_map->hm_points()[dom_p_of_dest].adj_hm_points_dist())
                    {
                        destinations_list.push_back(height_map_geodesic::PathPoint(&new_height_map->hm_points()[ite2.first]));
                    }
                    for (int k = 0; k < iternation_num(org_one_lqt1_two_lqt2_three_ls_four_lst_five_ds_six_point_seven, 2); k++)
                    {
                        post_algorithm.propagate(one_source_list, &destinations_list);
                    }
                    for (auto ite2 : new_height_map->hm_points()[dom_p_of_dest].adj_hm_points_dist())
                    {
                        int x = ite.first;
                        int y = ite2.first;
                        int x_y;
                        if (x > y)
                        {
                            int temp = y;
                            y = x;
                            x = temp;
                        }
                        hash_function_two_keys_to_one_key_int(new_height_map->hm_points().size(), x, y, x_y);
                        if (path_p_to_p_post_hm.count(x_y) == 0)
                        {
                            std::vector<height_map_geodesic::PathPoint> path;
                            path.clear();
                            height_map_geodesic::PathPoint p(&new_height_map->hm_points()[ite2.first]);
                            post_algorithm.trace_back(p, path);
                            path_p_to_p_post_hm[x_y] = path;
                        }
                    }
                }

                for (auto ite : new_height_map->hm_points()[dom_p_of_src].adj_hm_points_dist())
                {
                    for (auto ite2 : new_height_map->hm_points()[dom_p_of_dest].adj_hm_points_dist())
                    {
                        int x = ite.first;
                        int y = ite2.first;
                        int x_y;
                        if (x > y)
                        {
                            int temp = y;
                            y = x;
                            x = temp;
                        }
                        hash_function_two_keys_to_one_key_int(new_height_map->hm_points().size(), x, y, x_y);

                        std::vector<height_map_geodesic::PathPoint> path = path_p_to_p_post_hm[x_y];

                        new_height_map->hm_points()[source_index].z() = new_height_map->hm_points()[dom_p_of_src].getz();
                        height_map_geodesic::PathPoint real_source(&new_height_map->hm_points()[source_index]);
                        new_height_map->hm_points()[destination_index].z() = new_height_map->hm_points()[dom_p_of_dest].getz();
                        height_map_geodesic::PathPoint real_destination(&new_height_map->hm_points()[destination_index]);

                        double x1 = new_height_map->hm_points()[source_index].getx();
                        double y1 = new_height_map->hm_points()[source_index].gety();
                        double z1 = new_height_map->hm_points()[dom_p_of_src].getz();
                        double x2 = new_height_map->hm_points()[ite.first].getx();
                        double y2 = new_height_map->hm_points()[ite.first].gety();
                        double z2 = new_height_map->hm_points()[ite.first].getz();
                        double src_short_dist = cal_dist(x1, y1, z1, x2, y2, z2);

                        double x3 = new_height_map->hm_points()[ite2.first].getx();
                        double y3 = new_height_map->hm_points()[ite2.first].gety();
                        double z3 = new_height_map->hm_points()[ite2.first].getz();
                        double x4 = new_height_map->hm_points()[destination_index].getx();
                        double y4 = new_height_map->hm_points()[destination_index].gety();
                        double z4 = new_height_map->hm_points()[dom_p_of_dest].getz();
                        double dest_short_dist = cal_dist(x3, y3, z3, x4, y4, z4);

                        if (distance_result > length(path) + src_short_dist + dest_short_dist)
                        {
                            distance_result = length(path) + src_short_dist + dest_short_dist;
                            path.insert(path.begin(), real_destination);
                            path.push_back(real_source);
                            path_result = path;
                        }
                    }
                }
            }
        }
        distance_result = round(distance_result * 1000000000.0) / 1000000000.0;
        obj_to_other_obj_distance_and_index_list.push_back(std::make_pair(distance_result, destination_index));

        std::sort(obj_to_other_obj_distance_and_index_list.begin(), obj_to_other_obj_distance_and_index_list.end());
    }

    auto stop_knn_or_range_query_time = std::chrono::high_resolution_clock::now();
    auto duration_knn_or_range_query_time = std::chrono::duration_cast<std::chrono::microseconds>(stop_knn_or_range_query_time - start_knn_or_range_query_time);

    compare_d_value_and_range_value(d_value, range);

    auto start_knn_query_time = std::chrono::high_resolution_clock::now();
    knn_or_range_query(1, k_value, range, obj_to_other_obj_distance_and_index_list, knn_list);
    auto stop_knn_query_time = std::chrono::high_resolution_clock::now();
    auto duration_knn_query_time = std::chrono::duration_cast<std::chrono::microseconds>(stop_knn_query_time - start_knn_query_time);
    knn_query_time = duration_knn_or_range_query_time.count() + duration_knn_query_time.count();
    knn_query_time /= 1000;
    knn_query_time *= (d_value / 2000);

    auto start_range_query_time = std::chrono::high_resolution_clock::now();
    knn_or_range_query(2, k_value, range, obj_to_other_obj_distance_and_index_list, range_list);
    auto stop_range_query_time = std::chrono::high_resolution_clock::now();
    auto duration_range_query_time = std::chrono::duration_cast<std::chrono::microseconds>(stop_range_query_time - start_range_query_time);
    range_query_time = duration_knn_or_range_query_time.count() + duration_range_query_time.count();
    range_query_time /= 1000;
    range_query_time *= (range / 1000) * (d_value / 2000);
}

void simp_height_map_knn_and_range_query_ds(
    height_map_geodesic::HeightMap *org_height_map,
    std::unordered_map<std::pair<int, int>, std::vector<height_map_geodesic::PathPoint>, boost::hash<std::pair<int, int>>> &all_path,
    std::unordered_map<std::pair<int, int>, double, boost::hash<std::pair<int, int>>> &all_dist,
    int source_index,
    double &knn_query_time, double &range_query_time,
    int knn_and_range_query_obj_num, int k_value, double range, double d_value,
    std::vector<int> &knn_list, std::vector<int> &range_list)
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
        double distance_result = all_dist[std::make_pair(source_index, destination_index)];
        distance_result = round(distance_result * 1000000000.0) / 1000000000.0;
        obj_to_other_obj_distance_and_index_list.push_back(std::make_pair(distance_result, destination_index));
        std::sort(obj_to_other_obj_distance_and_index_list.begin(), obj_to_other_obj_distance_and_index_list.end());
    }

    auto stop_knn_or_range_query_time = std::chrono::high_resolution_clock::now();
    auto duration_knn_or_range_query_time = std::chrono::duration_cast<std::chrono::microseconds>(stop_knn_or_range_query_time - start_knn_or_range_query_time);

    compare_d_value_and_range_value(d_value, range);

    auto start_knn_query_time = std::chrono::high_resolution_clock::now();
    knn_or_range_query(1, k_value, range, obj_to_other_obj_distance_and_index_list, knn_list);
    auto stop_knn_query_time = std::chrono::high_resolution_clock::now();
    auto duration_knn_query_time = std::chrono::duration_cast<std::chrono::microseconds>(stop_knn_query_time - start_knn_query_time);
    knn_query_time = duration_knn_or_range_query_time.count() + duration_knn_query_time.count();
    knn_query_time /= 1000;
    knn_query_time *= (d_value / 2000);

    auto start_range_query_time = std::chrono::high_resolution_clock::now();
    knn_or_range_query(2, k_value, range, obj_to_other_obj_distance_and_index_list, range_list);
    auto stop_range_query_time = std::chrono::high_resolution_clock::now();
    auto duration_range_query_time = std::chrono::duration_cast<std::chrono::microseconds>(stop_range_query_time - start_range_query_time);
    range_query_time = duration_knn_or_range_query_time.count() + duration_range_query_time.count();
    range_query_time /= 1000;
    range_query_time *= (range / 1000) * (d_value / 2000);
}

// calculate the distance between each two vertices on the original terrain
void terrain_face_cal_org_dist_v_to_v(
    geodesic::Mesh *org_mesh,
    std::unordered_map<std::pair<unsigned, unsigned>, double, boost::hash<std::pair<unsigned, unsigned>>> &all_dist_v_to_v_org_terrain,
    double subdivision_level, bool exact_not_appr)
{
    geodesic::GeodesicAlgorithmExact org_mesh_algorithm_face_exact(org_mesh);
    geodesic::GeodesicAlgorithmSubdivision org_mesh_algorithm_face_appr(org_mesh, subdivision_level);
    double const distance_limit = geodesic::GEODESIC_INF;
    std::vector<geodesic::SurfacePoint> one_source_v_list;

    for (int i = 0; i < org_mesh->vertices().size(); i++)
    {
        one_source_v_list.clear();
        one_source_v_list.push_back(geodesic::SurfacePoint(&org_mesh->vertices()[i]));

        unsigned i_x_y;
        hash_function_two_keys_to_one_key(
            org_mesh->m_hash_max,
            unsigned(org_mesh->vertices()[i].getx() - org_mesh->m_xmin),
            unsigned(org_mesh->vertices()[i].gety() - org_mesh->m_ymin),
            i_x_y);
        assert(org_mesh->vert_coor_and_id_map()[i_x_y] == i);

        if (exact_not_appr)
        {
            org_mesh_algorithm_face_exact.propagate(one_source_v_list, distance_limit, 1);
        }
        else
        {
            org_mesh_algorithm_face_appr.propagate(one_source_v_list, distance_limit, 1);
        }

        for (int j = i; j < org_mesh->vertices().size(); j++)
        {
            std::vector<geodesic::SurfacePoint> path_result;
            path_result.clear();
            geodesic::SurfacePoint one_dest_v(&org_mesh->vertices()[j]);

            if (exact_not_appr)
            {
                org_mesh_algorithm_face_exact.trace_back(one_dest_v, path_result);
            }
            else
            {
                org_mesh_algorithm_face_appr.trace_back(one_dest_v, path_result);
            }

            double distance = length(path_result);

            unsigned j_x_y;
            hash_function_two_keys_to_one_key(
                org_mesh->m_hash_max,
                unsigned(org_mesh->vertices()[j].getx() - org_mesh->m_xmin),
                unsigned(org_mesh->vertices()[j].gety() - org_mesh->m_ymin),
                j_x_y);
            assert(org_mesh->vert_coor_and_id_map()[j_x_y] == j);

            if (i_x_y <= j_x_y)
            {
                all_dist_v_to_v_org_terrain[std::make_pair(i_x_y, j_x_y)] = distance;
            }
            else
            {
                all_dist_v_to_v_org_terrain[std::make_pair(j_x_y, i_x_y)] = distance;
            }
        }
    }
}

void terrain_simplify_one_iteration(
    geodesic::Mesh *org_mesh,
    geodesic::Mesh *pre_mesh,
    double epsilon,
    std::vector<double> &post_terrain_vertex,
    std::vector<unsigned> &post_terrain_face,
    bool &can_simplify, unsigned &del_vertex_x_y_hash,
    std::vector<std::pair<unsigned, double>> &del_v_neig_hash_and_dist_list,
    std::unordered_map<unsigned, unsigned> &skipped_v_hash_map,
    double subdivision_level,
    int face_exact_one_face_appr_two_vertex_three)
{
    std::vector<int> del_v_neig_id_list;
    std::vector<unsigned> del_v_neig_hash_list;
    del_v_neig_id_list.clear();
    del_v_neig_hash_list.clear();

    // select vertex with min interior angle for removing
    double min_interior_angle = 1e100;
    int min_interior_angle_vertex_index = -1;
    for (int i = 0; i < pre_mesh->vertices().size(); i++)
    {
        geodesic::Vertex &v = pre_mesh->vertices()[i];
        if (v.boundary())
        {
            continue;
        }

        // if at previous remove vertex iteration, remove this vertex does not satisfy
        // error bound, then skip this vertex
        unsigned curr_v_x_y_hash;
        hash_function_two_keys_to_one_key(org_mesh->m_hash_max,
                                          unsigned(v.getx() - pre_mesh->m_xmin),
                                          unsigned(v.gety() - pre_mesh->m_ymin),
                                          curr_v_x_y_hash);

        if (skipped_v_hash_map.count(curr_v_x_y_hash) != 0)
        {
            continue;
        }

        double avg_angle_difference = 0;
        bool interior_angle_larger_Pi = false;
        for (int j = 0; j < v.adjacent_edges().size(); j++)
        {
            double one_interior_angle = 0;
            geodesic::edge_pointer e = v.adjacent_edges()[j];
            geodesic::vertex_pointer another_v;
            if (e->adjacent_vertices()[0]->id() == v.id())
            {
                another_v = e->adjacent_vertices()[1];
            }
            else if (e->adjacent_vertices()[1]->id() == v.id())
            {
                another_v = e->adjacent_vertices()[0];
            }
            else
            {
                assert(false);
            }
            for (int k = 0; k < v.adjacent_faces().size(); k++)
            {
                geodesic::face_pointer f = v.adjacent_faces()[k];
                if (f->adjacent_vertices()[0]->id() == another_v->id())
                {
                    one_interior_angle += f->corner_angles()[0] / M_PI * 180;
                }
                if (f->adjacent_vertices()[1]->id() == another_v->id())
                {
                    one_interior_angle += f->corner_angles()[1] / M_PI * 180;
                }
                if (f->adjacent_vertices()[2]->id() == another_v->id())
                {
                    one_interior_angle += f->corner_angles()[2] / M_PI * 180;
                }
            }
            if (one_interior_angle > 150)
            {
                interior_angle_larger_Pi = true;
                break;
            }
            avg_angle_difference +=
                std::abs(one_interior_angle - 180 * (v.adjacent_edges().size() - 2) / v.adjacent_edges().size());
        }
        if (interior_angle_larger_Pi)
        {
            continue;
        }
        if (avg_angle_difference <= min_interior_angle)
        {
            min_interior_angle = avg_angle_difference;
            min_interior_angle_vertex_index = i;
        }
    }

    if (min_interior_angle_vertex_index == -1)
    {
        can_simplify = false;
        return;
    }

    // triangulate
    geodesic::Vertex &delete_v = pre_mesh->vertices()[min_interior_angle_vertex_index];
    assert(delete_v.id() == min_interior_angle_vertex_index);

    std::vector<std::pair<DelPoint, int>> triangulate_vertex;
    int index = 0;
    std::unordered_map<int, int> triangulate_vertex_map;
    triangulate_vertex.clear();
    triangulate_vertex_map.clear();

    for (int i = 0; i < delete_v.adjacent_edges().size(); i++)
    {
        double x, y;
        int mesh_vertex_index;
        if (delete_v.adjacent_edges()[i]->adjacent_vertices()[0]->id() == delete_v.id())
        {
            x = delete_v.adjacent_edges()[i]->adjacent_vertices()[1]->getx();
            y = delete_v.adjacent_edges()[i]->adjacent_vertices()[1]->gety();
            mesh_vertex_index = delete_v.adjacent_edges()[i]->adjacent_vertices()[1]->id();
        }
        else if (delete_v.adjacent_edges()[i]->adjacent_vertices()[1]->id() == delete_v.id())
        {
            x = delete_v.adjacent_edges()[i]->adjacent_vertices()[0]->getx();
            y = delete_v.adjacent_edges()[i]->adjacent_vertices()[0]->gety();
            mesh_vertex_index = delete_v.adjacent_edges()[i]->adjacent_vertices()[0]->id();
        }
        else
        {
            assert(false);
        }
        triangulate_vertex.push_back(std::make_pair(DelPoint(x, y), index));
        triangulate_vertex_map[index] = mesh_vertex_index;
        del_v_neig_id_list.push_back(mesh_vertex_index);
        index++;

        unsigned x_y;
        hash_function_two_keys_to_one_key(org_mesh->m_hash_max,
                                          unsigned(x - pre_mesh->m_xmin),
                                          unsigned(y - pre_mesh->m_ymin),
                                          x_y);
        del_v_neig_hash_list.push_back(x_y);
    }

    Delaunay T;
    T.insert(triangulate_vertex.begin(), triangulate_vertex.end());
    assert(T.is_valid());

    // store the remaining vertices (all vertices - deleted vertices)
    for (int i = 0; i < pre_mesh->vertices().size(); i++)
    {
        geodesic::Vertex &v = pre_mesh->vertices()[i];
        if (v.id() == min_interior_angle_vertex_index)
        {
            continue;
        }
        post_terrain_vertex.push_back(v.getx());
        post_terrain_vertex.push_back(v.gety());
        post_terrain_vertex.push_back(v.getz());
    }

    // store the remaining faces (all faces - deleted faces)
    for (int i = 0; i < pre_mesh->faces().size(); i++)
    {
        geodesic::Face &f = pre_mesh->faces()[i];
        if (f.adjacent_vertices()[0]->id() == min_interior_angle_vertex_index ||
            f.adjacent_vertices()[1]->id() == min_interior_angle_vertex_index ||
            f.adjacent_vertices()[2]->id() == min_interior_angle_vertex_index)
        {
            continue;
        }

        int output_id1, output_id2, output_id3;
        canonical_triangle(pre_mesh, f.adjacent_vertices()[0]->id(),
                           f.adjacent_vertices()[1]->id(), f.adjacent_vertices()[2]->id(),
                           output_id1, output_id2, output_id3);
        if (output_id1 > min_interior_angle_vertex_index)
        {
            output_id1--;
        }
        if (output_id2 > min_interior_angle_vertex_index)
        {
            output_id2--;
        }
        if (output_id3 > min_interior_angle_vertex_index)
        {
            output_id3--;
        }
        post_terrain_face.push_back(output_id1);
        post_terrain_face.push_back(output_id2);
        post_terrain_face.push_back(output_id3);
    }

    // store the newly triangulated faces
    for (Delaunay::Face_handle triangulate_face : T.finite_face_handles())
    {
        int i1 = triangulate_face->vertex(0)->info();
        int i2 = triangulate_face->vertex(1)->info();
        int i3 = triangulate_face->vertex(2)->info();

        int output_id1, output_id2, output_id3;
        canonical_triangle(pre_mesh, triangulate_vertex_map[i1],
                           triangulate_vertex_map[i2], triangulate_vertex_map[i3],
                           output_id1, output_id2, output_id3);
        if (output_id1 > min_interior_angle_vertex_index)
        {
            output_id1--;
        }
        if (output_id2 > min_interior_angle_vertex_index)
        {
            output_id2--;
        }
        if (output_id3 > min_interior_angle_vertex_index)
        {
            output_id3--;
        }
        post_terrain_face.push_back(output_id1);
        post_terrain_face.push_back(output_id2);
        post_terrain_face.push_back(output_id3);
    }

    // store deleted vertex x and y coordinate using hash key
    hash_function_two_keys_to_one_key(org_mesh->m_hash_max,
                                      unsigned(delete_v.getx() - pre_mesh->m_xmin),
                                      unsigned(delete_v.gety() - pre_mesh->m_ymin),
                                      del_vertex_x_y_hash);

    assert(del_v_neig_id_list.size() == del_v_neig_hash_list.size());

    int ite_num = floor(1 / pow(epsilon, 0.5));
    ite_num = face_exact_one_face_appr_two_vertex_three != 1 ? 1 : ite_num;
    assert(ite_num >= 1);

    // store x and y coord of neighbour of deleted vertex using hash key, and the distance to deleted vertex
    if (face_exact_one_face_appr_two_vertex_three == 1)
    {
        geodesic::GeodesicAlgorithmExact pre_mesh_algorithm(pre_mesh);
        double const distance_limit = geodesic::GEODESIC_INF;
        std::vector<geodesic::SurfacePoint> one_source_del_v_list;
        std::vector<geodesic::SurfacePoint> dests_del_v_neig_list;
        one_source_del_v_list.clear();
        dests_del_v_neig_list.clear();
        one_source_del_v_list.push_back(geodesic::SurfacePoint(&pre_mesh->vertices()[min_interior_angle_vertex_index]));
        for (int i = 0; i < del_v_neig_id_list.size(); i++)
        {
            dests_del_v_neig_list.push_back(geodesic::SurfacePoint(&pre_mesh->vertices()[del_v_neig_id_list[i]]));
        }
        pre_mesh_algorithm.propagate(one_source_del_v_list, distance_limit, ite_num);
        for (int i = 0; i < del_v_neig_id_list.size(); i++)
        {
            std::vector<geodesic::SurfacePoint> path_result;
            path_result.clear();
            pre_mesh_algorithm.trace_back(dests_del_v_neig_list[i], path_result);
            double distance = length(path_result);

            del_v_neig_hash_and_dist_list.push_back(std::make_pair(del_v_neig_hash_list[i], distance));
        }
    }

    geodesic::GeodesicAlgorithmExact pre_mesh_algorithm_face_exact(pre_mesh);
    geodesic::GeodesicAlgorithmSubdivision pre_mesh_algorithm_face_appr(pre_mesh, subdivision_level);
    geodesic::GeodesicAlgorithmSubdivision pre_mesh_algorithm_vertex(pre_mesh, 0);

    double const distance_limit = geodesic::GEODESIC_INF;
    std::vector<geodesic::SurfacePoint> one_source_del_v_list;
    std::vector<geodesic::SurfacePoint> dests_del_v_neig_list;
    one_source_del_v_list.clear();
    dests_del_v_neig_list.clear();
    one_source_del_v_list.push_back(geodesic::SurfacePoint(&pre_mesh->vertices()[min_interior_angle_vertex_index]));
    for (int i = 0; i < del_v_neig_id_list.size(); i++)
    {
        dests_del_v_neig_list.push_back(geodesic::SurfacePoint(&pre_mesh->vertices()[del_v_neig_id_list[i]]));
    }

    if (face_exact_one_face_appr_two_vertex_three == 1)
    {
        pre_mesh_algorithm_face_exact.propagate(one_source_del_v_list, distance_limit, ite_num);
    }
    else if (face_exact_one_face_appr_two_vertex_three == 2)
    {
        pre_mesh_algorithm_face_appr.propagate(one_source_del_v_list, distance_limit, 1);
    }
    else if (face_exact_one_face_appr_two_vertex_three == 3)
    {
        pre_mesh_algorithm_vertex.propagate(one_source_del_v_list, distance_limit, 1);
    }

    for (int i = 0; i < del_v_neig_id_list.size(); i++)
    {
        std::vector<geodesic::SurfacePoint> path_result;
        path_result.clear();

        if (face_exact_one_face_appr_two_vertex_three == 1)
        {
            pre_mesh_algorithm_face_exact.trace_back(dests_del_v_neig_list[i], path_result);
        }
        else if (face_exact_one_face_appr_two_vertex_three == 2)
        {
            pre_mesh_algorithm_face_appr.trace_back(dests_del_v_neig_list[i], path_result);
        }
        else if (face_exact_one_face_appr_two_vertex_three == 3)
        {
            pre_mesh_algorithm_vertex.trace_back(dests_del_v_neig_list[i], path_result);
            modify_path(path_result);
        }

        double distance = length(path_result);

        del_v_neig_hash_and_dist_list.push_back(std::make_pair(del_v_neig_hash_list[i], distance));
    }

    can_simplify = true;
}

void update_guest_host_map(
    std::unordered_map<unsigned, std::unordered_map<unsigned, double>> &guest_map,
    std::unordered_map<unsigned, std::unordered_map<unsigned, unsigned>> &host_map,
    unsigned del_vertex_x_y_hash,
    std::vector<std::pair<unsigned, double>> del_v_neig_hash_and_dist_list)
{
    std::unordered_map<unsigned, unsigned> one_host_of_del_v_map;
    one_host_of_del_v_map.clear();
    for (int i = 0; i < del_v_neig_hash_and_dist_list.size(); i++)
    {
        unsigned one_del_v_neig_hash = del_v_neig_hash_and_dist_list[i].first;
        one_host_of_del_v_map[one_del_v_neig_hash] = one_del_v_neig_hash;
    }
    host_map[del_vertex_x_y_hash] = one_host_of_del_v_map;

    std::unordered_map<unsigned, unsigned> processed_v_bar_map;
    processed_v_bar_map.clear();

    for (int i = 0; i < del_v_neig_hash_and_dist_list.size(); i++)
    {
        unsigned del_v_neig_x_y_hash = del_v_neig_hash_and_dist_list[i].first;
        std::unordered_map<unsigned, double> one_guest_of_del_v_neig_map;
        one_guest_of_del_v_neig_map.clear();
        one_guest_of_del_v_neig_map = guest_map[del_v_neig_x_y_hash];

        if (one_guest_of_del_v_neig_map.count(del_vertex_x_y_hash) == 0 ||
            (one_guest_of_del_v_neig_map.count(del_vertex_x_y_hash) != 0 &&
             one_guest_of_del_v_neig_map[del_vertex_x_y_hash] > del_v_neig_hash_and_dist_list[i].second))
        {
            one_guest_of_del_v_neig_map[del_vertex_x_y_hash] = del_v_neig_hash_and_dist_list[i].second;
        }
        guest_map[del_v_neig_x_y_hash] = one_guest_of_del_v_neig_map;

        for (auto ite : guest_map[del_vertex_x_y_hash])
        {
            unsigned v_bar = ite.first;
            double dist_tide = ite.second;

            if (processed_v_bar_map.count(v_bar) != 0)
            {
                continue;
            }
            processed_v_bar_map[v_bar] = 1;

            one_guest_of_del_v_neig_map.clear();
            one_guest_of_del_v_neig_map = guest_map[del_v_neig_x_y_hash];

            if (one_guest_of_del_v_neig_map.count(v_bar) == 0 ||
                (one_guest_of_del_v_neig_map.count(v_bar) != 0 &&
                 one_guest_of_del_v_neig_map[v_bar] > dist_tide + del_v_neig_hash_and_dist_list[i].second))
            {
                one_guest_of_del_v_neig_map[v_bar] = dist_tide + del_v_neig_hash_and_dist_list[i].second;
            }
            guest_map[del_v_neig_x_y_hash] = one_guest_of_del_v_neig_map;

            std::unordered_map<unsigned, unsigned> one_host_of_v_bar_map;
            one_host_of_v_bar_map.clear();
            one_host_of_v_bar_map = host_map[v_bar];

            assert(one_host_of_v_bar_map.count(del_vertex_x_y_hash) != 0);
            one_host_of_v_bar_map.erase(del_vertex_x_y_hash);
            one_host_of_v_bar_map[del_v_neig_x_y_hash] = del_v_neig_x_y_hash;
            host_map[v_bar] = one_host_of_v_bar_map;
        }
    }
}

// calculate the distance between each two neighbours and one neighbour one
// guest of the deleted vertex on the simplified terrain vertex and compare
// with the distance on the original terrain
void terrain_face_exact_and_face_appr_and_vertex_cal_simp_dist_and_check(
    geodesic::Mesh *org_mesh, geodesic::Mesh *post_mesh,
    double epsilon, double subdivision_level,
    int face_exact_one_face_appr_two_vertex_three, bool &satisfy,
    std::unordered_map<unsigned, std::unordered_map<unsigned, double>> &guest_map,
    std::unordered_map<unsigned, std::unordered_map<unsigned, unsigned>> &host_map,
    std::vector<std::pair<unsigned, double>> del_v_neig_hash_and_dist_list)
{
    geodesic::GeodesicAlgorithmExact org_mesh_algorithm_face_exact(org_mesh);
    geodesic::GeodesicAlgorithmExact post_mesh_algorithm_face_exact(post_mesh);
    geodesic::GeodesicAlgorithmSubdivision org_mesh_algorithm_face_appr(org_mesh, subdivision_level);
    geodesic::GeodesicAlgorithmSubdivision post_mesh_algorithm_face_appr(post_mesh, subdivision_level);
    geodesic::GeodesicAlgorithmSubdivision org_mesh_algorithm_vertex(org_mesh, 0);
    geodesic::GeodesicAlgorithmSubdivision post_mesh_algorithm_vertex(post_mesh, 0);

    double const distance_limit = geodesic::GEODESIC_INF;

    std::vector<geodesic::SurfacePoint> org_one_source_v_list;
    std::vector<geodesic::SurfacePoint> org_destinations_v_list;
    std::vector<geodesic::SurfacePoint> post_one_source_v_list;
    std::vector<geodesic::SurfacePoint> post_destinations_v_list;

    for (int i = 0; i < del_v_neig_hash_and_dist_list.size(); i++)
    {
        org_one_source_v_list.clear();
        org_destinations_v_list.clear();
        post_one_source_v_list.clear();
        post_destinations_v_list.clear();

        unsigned one_del_v_neig_hash = del_v_neig_hash_and_dist_list[i].first;
        assert(host_map.count(one_del_v_neig_hash) == 0);
        int org_mesh_id1 = org_mesh->vert_coor_and_id_map()[one_del_v_neig_hash];
        int post_mesh_id1 = post_mesh->vert_coor_and_id_map()[one_del_v_neig_hash];
        org_one_source_v_list.push_back(geodesic::SurfacePoint(&org_mesh->vertices()[org_mesh_id1]));
        post_one_source_v_list.push_back(geodesic::SurfacePoint(&post_mesh->vertices()[post_mesh_id1]));

        // make sure each destination is stored only once
        std::unordered_map<unsigned, unsigned> org_stored_dest_map;
        std::unordered_map<unsigned, unsigned> post_stored_dest_map;
        org_stored_dest_map.clear();
        post_stored_dest_map.clear();

        // another neighbour of the deleted vertex, i.e., intra-distance
        for (int j = i; j < del_v_neig_hash_and_dist_list.size(); j++)
        {
            unsigned another_del_v_neig_hash = del_v_neig_hash_and_dist_list[j].first;
            assert(host_map.count(another_del_v_neig_hash) == 0);
            if (org_stored_dest_map.count(another_del_v_neig_hash) == 0)
            {
                org_stored_dest_map[another_del_v_neig_hash] = another_del_v_neig_hash;
            }
            if (post_stored_dest_map.count(another_del_v_neig_hash) == 0)
            {
                post_stored_dest_map[another_del_v_neig_hash] = another_del_v_neig_hash;
            }
            int org_mesh_id2 = org_mesh->vert_coor_and_id_map()[another_del_v_neig_hash];
            int post_mesh_id2 = post_mesh->vert_coor_and_id_map()[another_del_v_neig_hash];
            org_destinations_v_list.push_back(geodesic::SurfacePoint(&org_mesh->vertices()[org_mesh_id2]));
            post_destinations_v_list.push_back(geodesic::SurfacePoint(&post_mesh->vertices()[post_mesh_id2]));
        }

        // each guest of the deleted vertex, i.e., inter-distance
        // note that the guest must be previous deleted vertex
        for (auto ite : guest_map[one_del_v_neig_hash])
        {
            unsigned guest_of_del_v = ite.first;
            assert(host_map.count(guest_of_del_v) != 0);

            assert(org_mesh->vert_coor_and_id_map().count(guest_of_del_v) != 0);
            if (org_stored_dest_map.count(guest_of_del_v) == 0)
            {
                org_stored_dest_map[guest_of_del_v] = guest_of_del_v;
            }
            int org_mesh_id3 = org_mesh->vert_coor_and_id_map()[guest_of_del_v];
            org_destinations_v_list.push_back(geodesic::SurfacePoint(&org_mesh->vertices()[org_mesh_id3]));

            // since the guest of the deleted vertex is deleted, we add the neighbour of these gurests as destination
            for (auto ite2 : host_map[guest_of_del_v])
            {
                unsigned neig_of_guest_of_del_v = ite2.first;
                assert(host_map.count(neig_of_guest_of_del_v) == 0);
                if (post_stored_dest_map.count(neig_of_guest_of_del_v) == 0)
                {
                    post_stored_dest_map[neig_of_guest_of_del_v] = neig_of_guest_of_del_v;
                }
                int post_mesh_id3 = post_mesh->vert_coor_and_id_map()[neig_of_guest_of_del_v];
                post_destinations_v_list.push_back(geodesic::SurfacePoint(&post_mesh->vertices()[post_mesh_id3]));
            }
        }
        if (face_exact_one_face_appr_two_vertex_three == 1)
        {
            org_mesh_algorithm_face_exact.propagate(org_one_source_v_list, &org_destinations_v_list);
            post_mesh_algorithm_face_exact.propagate(post_one_source_v_list, &post_destinations_v_list);
        }
        else if (face_exact_one_face_appr_two_vertex_three == 2)
        {
            org_mesh_algorithm_face_appr.propagate(org_one_source_v_list, &org_destinations_v_list);
            post_mesh_algorithm_face_appr.propagate(post_one_source_v_list, &post_destinations_v_list);
        }
        else if (face_exact_one_face_appr_two_vertex_three == 3)
        {
            org_mesh_algorithm_vertex.propagate(org_one_source_v_list, &org_destinations_v_list);
            post_mesh_algorithm_vertex.propagate(post_one_source_v_list, &post_destinations_v_list);
        }
        org_stored_dest_map.clear();
        post_stored_dest_map.clear();

        // calculate distance for another neighbour of the deleted vertex, i.e., intra-distance
        for (int j = i; j < del_v_neig_hash_and_dist_list.size(); j++)
        {
            unsigned another_del_v_neig_hash = del_v_neig_hash_and_dist_list[j].first;
            if (org_stored_dest_map.count(another_del_v_neig_hash) == 0)
            {
                org_stored_dest_map[another_del_v_neig_hash] = another_del_v_neig_hash;
            }
            if (post_stored_dest_map.count(another_del_v_neig_hash) == 0)
            {
                post_stored_dest_map[another_del_v_neig_hash] = another_del_v_neig_hash;
            }
            int org_mesh_id2 = org_mesh->vert_coor_and_id_map()[another_del_v_neig_hash];
            int post_mesh_id2 = post_mesh->vert_coor_and_id_map()[another_del_v_neig_hash];

            std::vector<geodesic::SurfacePoint> org_path_result;
            std::vector<geodesic::SurfacePoint> post_path_result;
            org_path_result.clear();
            post_path_result.clear();
            geodesic::SurfacePoint one_org_destination(&org_mesh->vertices()[org_mesh_id2]);
            geodesic::SurfacePoint one_post_destination(&post_mesh->vertices()[post_mesh_id2]);
            if (face_exact_one_face_appr_two_vertex_three == 1)
            {
                org_mesh_algorithm_face_exact.trace_back(one_org_destination, org_path_result);
                post_mesh_algorithm_face_exact.trace_back(one_post_destination, post_path_result);
            }
            else if (face_exact_one_face_appr_two_vertex_three == 2)
            {
                org_mesh_algorithm_face_appr.trace_back(one_org_destination, org_path_result);
                post_mesh_algorithm_face_appr.trace_back(one_post_destination, post_path_result);
            }
            else if (face_exact_one_face_appr_two_vertex_three == 3)
            {
                org_mesh_algorithm_vertex.trace_back(one_org_destination, org_path_result);
                post_mesh_algorithm_vertex.trace_back(one_post_destination, post_path_result);
                modify_path(org_path_result);
                modify_path(post_path_result);
            }
            double org_dist = length(org_path_result);
            double simp_dist = length(post_path_result);
            compare_distance(org_dist, simp_dist, epsilon, satisfy);
            if (!satisfy)
            {
                return;
            }
        }

        // calculate distance for neighbour of each guest of the deleted vertex, inter-distance
        for (auto ite : guest_map[one_del_v_neig_hash])
        {
            // for the original terrain, calculate the distance of each guest of the deleted vertex,
            // i.e., inter-distance, since the guest must on the original terrain
            unsigned guest_of_del_v = ite.first;
            if (org_stored_dest_map.count(guest_of_del_v) == 0)
            {
                org_stored_dest_map[guest_of_del_v] = guest_of_del_v;
            }
            int org_mesh_id3 = org_mesh->vert_coor_and_id_map()[guest_of_del_v];
            org_destinations_v_list.push_back(geodesic::SurfacePoint(&org_mesh->vertices()[org_mesh_id3]));

            std::vector<geodesic::SurfacePoint> org_path_result;
            org_path_result.clear();
            geodesic::SurfacePoint one_org_destination(&org_mesh->vertices()[org_mesh_id3]);
            if (face_exact_one_face_appr_two_vertex_three == 1)
            {
                org_mesh_algorithm_face_exact.trace_back(one_org_destination, org_path_result);
            }
            else if (face_exact_one_face_appr_two_vertex_three == 2)
            {
                org_mesh_algorithm_face_appr.trace_back(one_org_destination, org_path_result);
            }
            else if (face_exact_one_face_appr_two_vertex_three == 3)
            {
                org_mesh_algorithm_vertex.trace_back(one_org_destination, org_path_result);
                modify_path(org_path_result);
            }
            double org_dist = length(org_path_result);
            double simp_dist = 1e100;

            // for the simplified terrain, calculate the distance for neighbour of each
            // guest of the deleted vertex, i.e., inter-distance, since the guest is
            // not on the simplified terrain
            for (auto ite2 : host_map[guest_of_del_v])
            {
                unsigned neig_of_guest_of_del_v = ite2.first;
                if (post_stored_dest_map.count(neig_of_guest_of_del_v) == 0)
                {
                    post_stored_dest_map[neig_of_guest_of_del_v] = neig_of_guest_of_del_v;
                }
                double dist_neig_of_guest_of_del_v_to_guest_of_del_v =
                    guest_map[neig_of_guest_of_del_v][guest_of_del_v];

                int post_mesh_id3 = post_mesh->vert_coor_and_id_map()[neig_of_guest_of_del_v];
                post_destinations_v_list.push_back(geodesic::SurfacePoint(&post_mesh->vertices()[post_mesh_id3]));

                std::vector<geodesic::SurfacePoint> post_path_result;
                post_path_result.clear();
                geodesic::SurfacePoint one_post_destination(&post_mesh->vertices()[post_mesh_id3]);
                if (face_exact_one_face_appr_two_vertex_three == 1)
                {
                    post_mesh_algorithm_face_exact.trace_back(one_post_destination, post_path_result);
                }
                else if (face_exact_one_face_appr_two_vertex_three == 2)
                {
                    post_mesh_algorithm_face_appr.trace_back(one_post_destination, post_path_result);
                }
                else if (face_exact_one_face_appr_two_vertex_three == 3)
                {
                    post_mesh_algorithm_vertex.trace_back(one_post_destination, post_path_result);
                    modify_path(post_path_result);
                }
                double simp_dist2 = length(post_path_result) + dist_neig_of_guest_of_del_v_to_guest_of_del_v;
                simp_dist = std::min(simp_dist, simp_dist2);
            }
            compare_distance(org_dist, simp_dist, epsilon, satisfy);

            if (!satisfy)
            {
                return;
            }
        }
    }
}

void simp_terrain_face_exact_face_appr_vertex_query(
    geodesic::Mesh *org_mesh, geodesic::Mesh *post_mesh,
    std::unordered_map<unsigned, std::unordered_map<unsigned, double>> &guest_map,
    std::unordered_map<unsigned, std::unordered_map<unsigned, unsigned>> &host_map,
    double subdivision_level, int face_exact_one_face_appr_two_vertex_three,
    int source_index, int destination_index,
    double &query_time, double &query_memory_usage, double &distance_result,
    std::vector<geodesic::SurfacePoint> &path_result)
{
    auto start_query_time = std::chrono::high_resolution_clock::now();

    distance_result = 1e100;
    double distance_limit = 1e100;

    unsigned src_x_y;
    hash_function_two_keys_to_one_key(
        org_mesh->m_hash_max,
        unsigned(org_mesh->vertices()[source_index].getx() - org_mesh->m_xmin),
        unsigned(org_mesh->vertices()[source_index].gety() - org_mesh->m_ymin),
        src_x_y);

    unsigned dest_x_y;
    hash_function_two_keys_to_one_key(
        org_mesh->m_hash_max,
        unsigned(org_mesh->vertices()[destination_index].getx() - org_mesh->m_xmin),
        unsigned(org_mesh->vertices()[destination_index].gety() - org_mesh->m_ymin),
        dest_x_y);

    geodesic::GeodesicAlgorithmExact post_algorithm_face_exact(post_mesh);
    geodesic::GeodesicAlgorithmSubdivision post_algorithm_face_appr(post_mesh, subdivision_level);
    geodesic::GeodesicAlgorithmSubdivision post_algorithm_vertex(post_mesh, 0);

    if (host_map.count(src_x_y) == 0 && host_map.count(dest_x_y) == 0)
    {
        int post_mesh_id_src = post_mesh->vert_coor_and_id_map()[src_x_y];
        int post_mesh_id_dest = post_mesh->vert_coor_and_id_map()[dest_x_y];
        geodesic::SurfacePoint source(&post_mesh->vertices()[post_mesh_id_src]);
        geodesic::SurfacePoint destination(&post_mesh->vertices()[post_mesh_id_dest]);
        std::vector<geodesic::SurfacePoint> one_source_list(1, source);
        std::vector<geodesic::SurfacePoint> one_destination_list(1, destination);
        if (face_exact_one_face_appr_two_vertex_three == 1)
        {
            post_algorithm_face_exact.propagate(one_source_list, distance_limit, 1);
            post_algorithm_face_exact.trace_back(destination, path_result);
        }
        else if (face_exact_one_face_appr_two_vertex_three == 2)
        {
            post_algorithm_face_appr.propagate(one_source_list, distance_limit, 1);
            post_algorithm_face_appr.trace_back(destination, path_result);
        }
        else if (face_exact_one_face_appr_two_vertex_three == 3)
        {
            post_algorithm_vertex.propagate(one_source_list, distance_limit, 1);
            post_algorithm_vertex.trace_back(destination, path_result);
            modify_path(path_result);
        }
        distance_result = length(path_result);
    }
    else if (host_map.count(src_x_y) == 0 && host_map.count(dest_x_y) != 0)
    {
        std::unordered_map<std::pair<unsigned, unsigned>, std::vector<geodesic::SurfacePoint>,
                           boost::hash<std::pair<unsigned, unsigned>>>
            path_v_to_v_post_terrain;
        path_v_to_v_post_terrain.clear();

        std::vector<geodesic::SurfacePoint> one_source_list;
        std::vector<geodesic::SurfacePoint> destinations_list;
        one_source_list.clear();
        destinations_list.clear();

        int post_mesh_id_src = post_mesh->vert_coor_and_id_map()[src_x_y];
        one_source_list.push_back(geodesic::SurfacePoint(&post_mesh->vertices()[post_mesh_id_src]));
        for (auto ite : host_map[dest_x_y])
        {
            unsigned dest_neig_v = ite.first;
            assert(guest_map.count(dest_neig_v) != 0);
            assert(guest_map[dest_neig_v].count(dest_x_y) != 0);
            int post_mesh_id_dest = post_mesh->vert_coor_and_id_map()[dest_neig_v];
            destinations_list.push_back(geodesic::SurfacePoint(&post_mesh->vertices()[post_mesh_id_dest]));
        }

        if (face_exact_one_face_appr_two_vertex_three == 1)
        {
            post_algorithm_face_exact.propagate(one_source_list, distance_limit, 1);
        }
        else if (face_exact_one_face_appr_two_vertex_three == 2)
        {
            post_algorithm_face_appr.propagate(one_source_list, distance_limit, 1);
        }
        else if (face_exact_one_face_appr_two_vertex_three == 3)
        {
            post_algorithm_vertex.propagate(one_source_list, distance_limit, 1);
        }

        for (auto ite : host_map[dest_x_y])
        {
            unsigned dest_neig_v = ite.first;
            unsigned x_y_neig_1, x_y_neig_2;
            if (src_x_y <= dest_neig_v)
            {
                x_y_neig_1 = src_x_y;
                x_y_neig_2 = dest_neig_v;
            }
            else
            {
                x_y_neig_1 = dest_neig_v;
                x_y_neig_2 = src_x_y;
            }
            if (path_v_to_v_post_terrain.count(std::make_pair(x_y_neig_1, x_y_neig_2)) == 0)
            {
                std::vector<geodesic::SurfacePoint> path;
                path.clear();
                int post_mesh_id_dest = post_mesh->vert_coor_and_id_map()[dest_neig_v];
                geodesic::SurfacePoint p(&post_mesh->vertices()[post_mesh_id_dest]);
                if (face_exact_one_face_appr_two_vertex_three == 1)
                {
                    post_algorithm_face_exact.trace_back(p, path);
                }
                else if (face_exact_one_face_appr_two_vertex_three == 2)
                {
                    post_algorithm_face_appr.trace_back(p, path);
                }
                else if (face_exact_one_face_appr_two_vertex_three == 3)
                {
                    post_algorithm_vertex.trace_back(p, path);
                    modify_path(path);
                }
                path_v_to_v_post_terrain[std::make_pair(x_y_neig_1, x_y_neig_2)] = path;
            }
        }

        post_mesh->add_vertex(org_mesh->vertices()[destination_index].getx(),
                              org_mesh->vertices()[destination_index].gety(),
                              org_mesh->vertices()[destination_index].getz());
        geodesic::SurfacePoint real_destination(&post_mesh->vertices()[post_mesh->vertices().size() - 1]);

        for (auto ite : host_map[dest_x_y])
        {
            unsigned dest_neig_v = ite.first;
            unsigned x_y_neig_1, x_y_neig_2;
            if (src_x_y <= dest_neig_v)
            {
                x_y_neig_1 = src_x_y;
                x_y_neig_2 = dest_neig_v;
            }
            else
            {
                x_y_neig_1 = dest_neig_v;
                x_y_neig_2 = src_x_y;
            }
            std::vector<geodesic::SurfacePoint> path = path_v_to_v_post_terrain[std::make_pair(x_y_neig_1, x_y_neig_2)];
            double dist_dest_neig_v = guest_map[dest_neig_v][dest_x_y];

            if (distance_result > length(path) + dist_dest_neig_v)
            {
                distance_result = length(path) + dist_dest_neig_v;
                path.insert(path.begin(), real_destination);
                path_result = path;
            }
        }
    }
    else if (host_map.count(src_x_y) != 0 && host_map.count(dest_x_y) == 0)
    {
        std::unordered_map<std::pair<unsigned, unsigned>, std::vector<geodesic::SurfacePoint>,
                           boost::hash<std::pair<unsigned, unsigned>>>
            path_v_to_v_post_terrain;
        path_v_to_v_post_terrain.clear();

        std::vector<geodesic::SurfacePoint> one_source_list;
        std::vector<geodesic::SurfacePoint> destinations_list;
        one_source_list.clear();
        destinations_list.clear();

        int post_mesh_id_dest = post_mesh->vert_coor_and_id_map()[dest_x_y];
        one_source_list.push_back(geodesic::SurfacePoint(&post_mesh->vertices()[post_mesh_id_dest]));
        for (auto ite : host_map[src_x_y])
        {
            unsigned src_neig_v = ite.first;
            assert(guest_map.count(src_neig_v) != 0);
            assert(guest_map[src_neig_v].count(src_x_y) != 0);
            int post_mesh_id_src = post_mesh->vert_coor_and_id_map()[src_neig_v];
            destinations_list.push_back(geodesic::SurfacePoint(&post_mesh->vertices()[post_mesh_id_src]));
        }

        if (face_exact_one_face_appr_two_vertex_three == 1)
        {
            post_algorithm_face_exact.propagate(one_source_list, distance_limit, 1);
        }
        else if (face_exact_one_face_appr_two_vertex_three == 2)
        {
            post_algorithm_face_appr.propagate(one_source_list, distance_limit, 1);
        }
        else if (face_exact_one_face_appr_two_vertex_three == 3)
        {
            post_algorithm_vertex.propagate(one_source_list, distance_limit, 1);
        }

        for (auto ite : host_map[src_x_y])
        {
            unsigned src_neig_v = ite.first;
            unsigned x_neig_y_1, x_neig_y_2;
            if (src_neig_v <= dest_x_y)
            {
                x_neig_y_1 = src_neig_v;
                x_neig_y_2 = dest_x_y;
            }
            else
            {
                x_neig_y_1 = dest_x_y;
                x_neig_y_2 = src_neig_v;
            }
            if (path_v_to_v_post_terrain.count(std::make_pair(x_neig_y_1, x_neig_y_2)) == 0)
            {
                std::vector<geodesic::SurfacePoint> path;
                path.clear();
                int post_mesh_id_src = post_mesh->vert_coor_and_id_map()[src_neig_v];
                geodesic::SurfacePoint p(&post_mesh->vertices()[post_mesh_id_src]);
                if (face_exact_one_face_appr_two_vertex_three == 1)
                {
                    post_algorithm_face_exact.trace_back(p, path);
                }
                else if (face_exact_one_face_appr_two_vertex_three == 2)
                {
                    post_algorithm_face_appr.trace_back(p, path);
                }
                else if (face_exact_one_face_appr_two_vertex_three == 3)
                {
                    post_algorithm_vertex.trace_back(p, path);
                    modify_path(path);
                }
                path_v_to_v_post_terrain[std::make_pair(x_neig_y_1, x_neig_y_2)] = path;
            }
        }

        post_mesh->add_vertex(org_mesh->vertices()[source_index].getx(),
                              org_mesh->vertices()[source_index].gety(),
                              org_mesh->vertices()[source_index].getz());
        geodesic::SurfacePoint real_source(&post_mesh->vertices()[post_mesh->vertices().size() - 1]);

        for (auto ite : host_map[src_x_y])
        {
            unsigned src_neig_v = ite.first;
            unsigned x_neig_y_1, x_neig_y_2;
            if (src_neig_v <= dest_x_y)
            {
                x_neig_y_1 = src_neig_v;
                x_neig_y_2 = dest_x_y;
            }
            else
            {
                x_neig_y_1 = dest_x_y;
                x_neig_y_2 = src_neig_v;
            }
            std::vector<geodesic::SurfacePoint> path = path_v_to_v_post_terrain[std::make_pair(x_neig_y_1, x_neig_y_2)];
            double dist_src_neig_v = guest_map[src_neig_v][src_x_y];

            if (distance_result > length(path) + dist_src_neig_v)
            {
                distance_result = length(path) + dist_src_neig_v;
                path.insert(path.begin(), real_source);
                std::reverse(path.begin(), path.end());
                path_result = path;
            }
        }
    }
    else
    {
        // if source and destination are same
        if (src_x_y == dest_x_y)
        {
            path_result.clear();
            distance_result = 0;
        }
        // if the host map of source and destination are same
        else if (host_map[src_x_y] == host_map[dest_x_y])
        {
            post_mesh->add_vertex(org_mesh->vertices()[source_index].getx(),
                                  org_mesh->vertices()[source_index].gety(),
                                  org_mesh->vertices()[source_index].getz());
            geodesic::SurfacePoint real_source(&post_mesh->vertices()[post_mesh->vertices().size() - 1]);

            post_mesh->add_vertex(org_mesh->vertices()[destination_index].getx(),
                                  org_mesh->vertices()[destination_index].gety(),
                                  org_mesh->vertices()[destination_index].getz());
            geodesic::SurfacePoint real_destination(&post_mesh->vertices()[post_mesh->vertices().size() - 1]);

            path_result.push_back(real_destination);
            path_result.push_back(real_source);
            if (guest_map.count(src_x_y) != 0)
            {
                distance_result = std::min(distance_result, guest_map[src_x_y][dest_x_y]);
            }
            if (guest_map.count(dest_x_y) != 0)
            {
                distance_result = std::min(distance_result, guest_map[dest_x_y][src_x_y]);
            }
        }

        else
        {
            std::unordered_map<std::pair<unsigned, unsigned>, std::vector<geodesic::SurfacePoint>,
                               boost::hash<std::pair<unsigned, unsigned>>>
                path_v_to_v_post_terrain;
            path_v_to_v_post_terrain.clear();

            for (auto ite : host_map[src_x_y])
            {
                std::vector<geodesic::SurfacePoint> one_source_list;
                std::vector<geodesic::SurfacePoint> destinations_list;
                one_source_list.clear();
                destinations_list.clear();

                unsigned src_neig_v = ite.first;
                assert(guest_map.count(src_neig_v) != 0);
                assert(guest_map[src_neig_v].count(src_x_y) != 0);
                int post_mesh_id_src = post_mesh->vert_coor_and_id_map()[src_neig_v];
                one_source_list.push_back(geodesic::SurfacePoint(&post_mesh->vertices()[post_mesh_id_src]));
                for (auto ite2 : host_map[dest_x_y])
                {
                    unsigned dest_neig_v = ite2.first;
                    assert(guest_map.count(dest_neig_v) != 0);
                    assert(guest_map[dest_neig_v].count(dest_x_y) != 0);
                    int post_mesh_id_dest = post_mesh->vert_coor_and_id_map()[dest_neig_v];
                    destinations_list.push_back(geodesic::SurfacePoint(&post_mesh->vertices()[post_mesh_id_dest]));
                }

                if (face_exact_one_face_appr_two_vertex_three == 1)
                {
                    post_algorithm_face_exact.propagate(one_source_list, distance_limit, 1);
                }
                else if (face_exact_one_face_appr_two_vertex_three == 2)
                {
                    post_algorithm_face_appr.propagate(one_source_list, distance_limit, 1);
                }
                else if (face_exact_one_face_appr_two_vertex_three == 3)
                {
                    post_algorithm_vertex.propagate(one_source_list, distance_limit, 1);
                }

                for (auto ite2 : host_map[dest_x_y])
                {
                    unsigned dest_neig_v = ite2.first;
                    unsigned x_neig_y_neig_1, x_neig_y_neig_2;
                    if (src_neig_v <= dest_neig_v)
                    {
                        x_neig_y_neig_1 = src_neig_v;
                        x_neig_y_neig_2 = dest_neig_v;
                    }
                    else
                    {
                        x_neig_y_neig_1 = dest_neig_v;
                        x_neig_y_neig_2 = src_neig_v;
                    }
                    if (path_v_to_v_post_terrain.count(std::make_pair(x_neig_y_neig_1, x_neig_y_neig_2)) == 0)
                    {
                        std::vector<geodesic::SurfacePoint> path;
                        path.clear();
                        int post_mesh_id_dest = post_mesh->vert_coor_and_id_map()[dest_neig_v];
                        geodesic::SurfacePoint p(&post_mesh->vertices()[post_mesh_id_dest]);
                        if (face_exact_one_face_appr_two_vertex_three == 1)
                        {
                            post_algorithm_face_exact.trace_back(p, path);
                        }
                        else if (face_exact_one_face_appr_two_vertex_three == 2)
                        {
                            post_algorithm_face_appr.trace_back(p, path);
                        }
                        else if (face_exact_one_face_appr_two_vertex_three == 3)
                        {
                            post_algorithm_vertex.trace_back(p, path);
                            modify_path(path);
                        }
                        path_v_to_v_post_terrain[std::make_pair(x_neig_y_neig_1, x_neig_y_neig_2)] = path;
                    }
                }
            }

            post_mesh->add_vertex(org_mesh->vertices()[source_index].getx(),
                                  org_mesh->vertices()[source_index].gety(),
                                  org_mesh->vertices()[source_index].getz());
            geodesic::SurfacePoint real_source(&post_mesh->vertices()[post_mesh->vertices().size() - 1]);

            post_mesh->add_vertex(org_mesh->vertices()[destination_index].getx(),
                                  org_mesh->vertices()[destination_index].gety(),
                                  org_mesh->vertices()[destination_index].getz());
            geodesic::SurfacePoint real_destination(&post_mesh->vertices()[post_mesh->vertices().size() - 1]);

            for (auto ite : host_map[src_x_y])
            {
                unsigned src_neig_v = ite.first;
                double src_dest_neig_v = guest_map[src_neig_v][src_x_y];

                for (auto ite2 : host_map[dest_x_y])
                {
                    unsigned dest_neig_v = ite2.first;
                    unsigned x_neig_y_neig_1, x_neig_y_neig_2;
                    if (src_neig_v <= dest_neig_v)
                    {
                        x_neig_y_neig_1 = src_neig_v;
                        x_neig_y_neig_2 = dest_neig_v;
                    }
                    else
                    {
                        x_neig_y_neig_1 = dest_neig_v;
                        x_neig_y_neig_2 = src_neig_v;
                    }
                    std::vector<geodesic::SurfacePoint> path = path_v_to_v_post_terrain[std::make_pair(x_neig_y_neig_1, x_neig_y_neig_2)];
                    double dist_dest_neig_v = guest_map[dest_neig_v][dest_x_y];

                    if (distance_result > length(path) + src_dest_neig_v + dist_dest_neig_v)
                    {
                        distance_result = length(path) + src_dest_neig_v + dist_dest_neig_v;
                        path.insert(path.begin(), real_destination);
                        path.push_back(real_source);
                        path_result = path;
                    }
                }
            }
        }
    }

    distance_result = round(distance_result * 1000000000.0) / 1000000000.0;
    if (face_exact_one_face_appr_two_vertex_three == 1)
    {
        query_memory_usage += post_algorithm_face_exact.get_memory() + path_result.size() * sizeof(geodesic::SurfacePoint) + sizeof(double);
    }
    else if (face_exact_one_face_appr_two_vertex_three == 2)
    {
        query_memory_usage += post_algorithm_face_appr.get_memory() + path_result.size() * sizeof(geodesic::SurfacePoint) + sizeof(double);
    }
    else if (face_exact_one_face_appr_two_vertex_three == 3)
    {
        query_memory_usage += post_algorithm_vertex.get_memory() + path_result.size() * sizeof(geodesic::SurfacePoint) + sizeof(double);
    }

    auto stop_query_time = std::chrono::high_resolution_clock::now();
    auto duration_query_time = std::chrono::duration_cast<std::chrono::microseconds>(stop_query_time - start_query_time);
    query_time = duration_query_time.count();
    query_time /= 1000;
}

void simp_terrain_face_exact_face_appr_vertex_knn_and_range_query(
    geodesic::Mesh *org_mesh, geodesic::Mesh *post_mesh,
    std::unordered_map<unsigned, std::unordered_map<unsigned, double>> &guest_map,
    std::unordered_map<unsigned, std::unordered_map<unsigned, unsigned>> &host_map,
    double subdivision_level, int face_exact_one_face_appr_two_vertex_three,
    int source_index,
    double &knn_query_time, double &range_query_time,
    int knn_and_range_query_obj_num, int k_value, double range, double d_value,
    std::vector<int> &knn_list,
    std::vector<int> &range_list)
{
    auto start_knn_or_range_query_time = std::chrono::high_resolution_clock::now();

    std::vector<std::pair<double, int>> obj_to_other_obj_distance_and_index_list;
    obj_to_other_obj_distance_and_index_list.clear();

    bool contain_source_index = false;
    for (int m = org_mesh->vertices().size() - 1; m > org_mesh->vertices().size() - 1 - knn_and_range_query_obj_num; m--)
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

        double distance_result = 1e100;
        double distance_limit = 1e100;
        std::vector<geodesic::SurfacePoint> path_result;
        path_result.clear();

        unsigned src_x_y;
        hash_function_two_keys_to_one_key(
            org_mesh->m_hash_max,
            unsigned(org_mesh->vertices()[source_index].getx() - org_mesh->m_xmin),
            unsigned(org_mesh->vertices()[source_index].gety() - org_mesh->m_ymin),
            src_x_y);

        unsigned dest_x_y;
        hash_function_two_keys_to_one_key(
            org_mesh->m_hash_max,
            unsigned(org_mesh->vertices()[destination_index].getx() - org_mesh->m_xmin),
            unsigned(org_mesh->vertices()[destination_index].gety() - org_mesh->m_ymin),
            dest_x_y);

        geodesic::GeodesicAlgorithmExact post_algorithm_face_exact(post_mesh);
        geodesic::GeodesicAlgorithmSubdivision post_algorithm_face_appr(post_mesh, subdivision_level);
        geodesic::GeodesicAlgorithmSubdivision post_algorithm_vertex(post_mesh, 0);

        if (host_map.count(src_x_y) == 0 && host_map.count(dest_x_y) == 0)
        {
            int post_mesh_id_src = post_mesh->vert_coor_and_id_map()[src_x_y];
            int post_mesh_id_dest = post_mesh->vert_coor_and_id_map()[dest_x_y];
            geodesic::SurfacePoint source(&post_mesh->vertices()[post_mesh_id_src]);
            geodesic::SurfacePoint destination(&post_mesh->vertices()[post_mesh_id_dest]);
            std::vector<geodesic::SurfacePoint> one_source_list(1, source);
            std::vector<geodesic::SurfacePoint> one_destination_list(1, destination);
            if (face_exact_one_face_appr_two_vertex_three == 1)
            {
                post_algorithm_face_exact.propagate(one_source_list, distance_limit, 1);
                post_algorithm_face_exact.trace_back(destination, path_result);
            }
            else if (face_exact_one_face_appr_two_vertex_three == 2)
            {
                post_algorithm_face_appr.propagate(one_source_list, distance_limit, 1);
                post_algorithm_face_appr.trace_back(destination, path_result);
            }
            else if (face_exact_one_face_appr_two_vertex_three == 3)
            {
                post_algorithm_vertex.propagate(one_source_list, distance_limit, 1);
                post_algorithm_vertex.trace_back(destination, path_result);
                modify_path(path_result);
            }
            distance_result = length(path_result);
        }
        else if (host_map.count(src_x_y) == 0 && host_map.count(dest_x_y) != 0)
        {
            std::unordered_map<std::pair<unsigned, unsigned>, std::vector<geodesic::SurfacePoint>,
                               boost::hash<std::pair<unsigned, unsigned>>>
                path_v_to_v_post_terrain;
            path_v_to_v_post_terrain.clear();

            std::vector<geodesic::SurfacePoint> one_source_list;
            std::vector<geodesic::SurfacePoint> destinations_list;
            one_source_list.clear();
            destinations_list.clear();

            int post_mesh_id_src = post_mesh->vert_coor_and_id_map()[src_x_y];
            one_source_list.push_back(geodesic::SurfacePoint(&post_mesh->vertices()[post_mesh_id_src]));
            for (auto ite : host_map[dest_x_y])
            {
                unsigned dest_neig_v = ite.first;
                assert(guest_map.count(dest_neig_v) != 0);
                assert(guest_map[dest_neig_v].count(dest_x_y) != 0);
                int post_mesh_id_dest = post_mesh->vert_coor_and_id_map()[dest_neig_v];
                destinations_list.push_back(geodesic::SurfacePoint(&post_mesh->vertices()[post_mesh_id_dest]));
            }

            if (face_exact_one_face_appr_two_vertex_three == 1)
            {
                post_algorithm_face_exact.propagate(one_source_list, distance_limit, 1);
            }
            else if (face_exact_one_face_appr_two_vertex_three == 2)
            {
                post_algorithm_face_appr.propagate(one_source_list, distance_limit, 1);
            }
            else if (face_exact_one_face_appr_two_vertex_three == 3)
            {
                post_algorithm_vertex.propagate(one_source_list, distance_limit, 1);
            }

            for (auto ite : host_map[dest_x_y])
            {
                unsigned dest_neig_v = ite.first;
                unsigned x_y_neig_1, x_y_neig_2;
                if (src_x_y <= dest_neig_v)
                {
                    x_y_neig_1 = src_x_y;
                    x_y_neig_2 = dest_neig_v;
                }
                else
                {
                    x_y_neig_1 = dest_neig_v;
                    x_y_neig_2 = src_x_y;
                }
                if (path_v_to_v_post_terrain.count(std::make_pair(x_y_neig_1, x_y_neig_2)) == 0)
                {
                    std::vector<geodesic::SurfacePoint> path;
                    path.clear();
                    int post_mesh_id_dest = post_mesh->vert_coor_and_id_map()[dest_neig_v];
                    geodesic::SurfacePoint p(&post_mesh->vertices()[post_mesh_id_dest]);
                    if (face_exact_one_face_appr_two_vertex_three == 1)
                    {
                        post_algorithm_face_exact.trace_back(p, path);
                    }
                    else if (face_exact_one_face_appr_two_vertex_three == 2)
                    {
                        post_algorithm_face_appr.trace_back(p, path);
                    }
                    else if (face_exact_one_face_appr_two_vertex_three == 3)
                    {
                        post_algorithm_vertex.trace_back(p, path);
                        modify_path(path);
                    }
                    path_v_to_v_post_terrain[std::make_pair(x_y_neig_1, x_y_neig_2)] = path;
                }
            }

            post_mesh->add_vertex(org_mesh->vertices()[destination_index].getx(),
                                  org_mesh->vertices()[destination_index].gety(),
                                  org_mesh->vertices()[destination_index].getz());
            geodesic::SurfacePoint real_destination(&post_mesh->vertices()[post_mesh->vertices().size() - 1]);

            for (auto ite : host_map[dest_x_y])
            {
                unsigned dest_neig_v = ite.first;
                unsigned x_y_neig_1, x_y_neig_2;
                if (src_x_y <= dest_neig_v)
                {
                    x_y_neig_1 = src_x_y;
                    x_y_neig_2 = dest_neig_v;
                }
                else
                {
                    x_y_neig_1 = dest_neig_v;
                    x_y_neig_2 = src_x_y;
                }
                std::vector<geodesic::SurfacePoint> path = path_v_to_v_post_terrain[std::make_pair(x_y_neig_1, x_y_neig_2)];
                double dist_dest_neig_v = guest_map[dest_neig_v][dest_x_y];

                if (distance_result > length(path) + dist_dest_neig_v)
                {
                    distance_result = length(path) + dist_dest_neig_v;
                    path.insert(path.begin(), real_destination);
                    path_result = path;
                }
            }
        }
        else if (host_map.count(src_x_y) != 0 && host_map.count(dest_x_y) == 0)
        {
            std::unordered_map<std::pair<unsigned, unsigned>, std::vector<geodesic::SurfacePoint>,
                               boost::hash<std::pair<unsigned, unsigned>>>
                path_v_to_v_post_terrain;
            path_v_to_v_post_terrain.clear();

            std::vector<geodesic::SurfacePoint> one_source_list;
            std::vector<geodesic::SurfacePoint> destinations_list;
            one_source_list.clear();
            destinations_list.clear();

            int post_mesh_id_dest = post_mesh->vert_coor_and_id_map()[dest_x_y];
            one_source_list.push_back(geodesic::SurfacePoint(&post_mesh->vertices()[post_mesh_id_dest]));
            for (auto ite : host_map[src_x_y])
            {
                unsigned src_neig_v = ite.first;
                assert(guest_map.count(src_neig_v) != 0);
                assert(guest_map[src_neig_v].count(src_x_y) != 0);
                int post_mesh_id_src = post_mesh->vert_coor_and_id_map()[src_neig_v];
                destinations_list.push_back(geodesic::SurfacePoint(&post_mesh->vertices()[post_mesh_id_src]));
            }

            if (face_exact_one_face_appr_two_vertex_three == 1)
            {
                post_algorithm_face_exact.propagate(one_source_list, distance_limit, 1);
            }
            else if (face_exact_one_face_appr_two_vertex_three == 2)
            {
                post_algorithm_face_appr.propagate(one_source_list, distance_limit, 1);
            }
            else if (face_exact_one_face_appr_two_vertex_three == 3)
            {
                post_algorithm_vertex.propagate(one_source_list, distance_limit, 1);
            }

            for (auto ite : host_map[src_x_y])
            {
                unsigned src_neig_v = ite.first;
                unsigned x_neig_y_1, x_neig_y_2;
                if (src_neig_v <= dest_x_y)
                {
                    x_neig_y_1 = src_neig_v;
                    x_neig_y_2 = dest_x_y;
                }
                else
                {
                    x_neig_y_1 = dest_x_y;
                    x_neig_y_2 = src_neig_v;
                }
                if (path_v_to_v_post_terrain.count(std::make_pair(x_neig_y_1, x_neig_y_2)) == 0)
                {
                    std::vector<geodesic::SurfacePoint> path;
                    path.clear();
                    int post_mesh_id_src = post_mesh->vert_coor_and_id_map()[src_neig_v];
                    geodesic::SurfacePoint p(&post_mesh->vertices()[post_mesh_id_src]);
                    if (face_exact_one_face_appr_two_vertex_three == 1)
                    {
                        post_algorithm_face_exact.trace_back(p, path);
                    }
                    else if (face_exact_one_face_appr_two_vertex_three == 2)
                    {
                        post_algorithm_face_appr.trace_back(p, path);
                    }
                    else if (face_exact_one_face_appr_two_vertex_three == 3)
                    {
                        post_algorithm_vertex.trace_back(p, path);
                        modify_path(path);
                    }
                    path_v_to_v_post_terrain[std::make_pair(x_neig_y_1, x_neig_y_2)] = path;
                }
            }

            post_mesh->add_vertex(org_mesh->vertices()[source_index].getx(),
                                  org_mesh->vertices()[source_index].gety(),
                                  org_mesh->vertices()[source_index].getz());
            geodesic::SurfacePoint real_source(&post_mesh->vertices()[post_mesh->vertices().size() - 1]);

            for (auto ite : host_map[src_x_y])
            {
                unsigned src_neig_v = ite.first;
                unsigned x_neig_y_1, x_neig_y_2;
                if (src_neig_v <= dest_x_y)
                {
                    x_neig_y_1 = src_neig_v;
                    x_neig_y_2 = dest_x_y;
                }
                else
                {
                    x_neig_y_1 = dest_x_y;
                    x_neig_y_2 = src_neig_v;
                }
                std::vector<geodesic::SurfacePoint> path = path_v_to_v_post_terrain[std::make_pair(x_neig_y_1, x_neig_y_2)];
                double dist_src_neig_v = guest_map[src_neig_v][src_x_y];

                if (distance_result > length(path) + dist_src_neig_v)
                {
                    distance_result = length(path) + dist_src_neig_v;
                    path.insert(path.begin(), real_source);
                    std::reverse(path.begin(), path.end());
                    path_result = path;
                }
            }
        }
        else
        {
            // if source and destination are same
            if (src_x_y == dest_x_y)
            {
                path_result.clear();
                distance_result = 0;
            }
            // if the host map of source and destination are same
            else if (host_map[src_x_y] == host_map[dest_x_y])
            {
                post_mesh->add_vertex(org_mesh->vertices()[source_index].getx(),
                                      org_mesh->vertices()[source_index].gety(),
                                      org_mesh->vertices()[source_index].getz());
                geodesic::SurfacePoint real_source(&post_mesh->vertices()[post_mesh->vertices().size() - 1]);

                post_mesh->add_vertex(org_mesh->vertices()[destination_index].getx(),
                                      org_mesh->vertices()[destination_index].gety(),
                                      org_mesh->vertices()[destination_index].getz());
                geodesic::SurfacePoint real_destination(&post_mesh->vertices()[post_mesh->vertices().size() - 1]);

                path_result.push_back(real_destination);
                path_result.push_back(real_source);
                if (guest_map.count(src_x_y) != 0)
                {
                    distance_result = std::min(distance_result, guest_map[src_x_y][dest_x_y]);
                }
                if (guest_map.count(dest_x_y) != 0)
                {
                    distance_result = std::min(distance_result, guest_map[dest_x_y][src_x_y]);
                }
            }

            else
            {
                std::unordered_map<std::pair<unsigned, unsigned>, std::vector<geodesic::SurfacePoint>,
                                   boost::hash<std::pair<unsigned, unsigned>>>
                    path_v_to_v_post_terrain;
                path_v_to_v_post_terrain.clear();

                for (auto ite : host_map[src_x_y])
                {
                    std::vector<geodesic::SurfacePoint> one_source_list;
                    std::vector<geodesic::SurfacePoint> destinations_list;
                    one_source_list.clear();
                    destinations_list.clear();

                    unsigned src_neig_v = ite.first;
                    assert(guest_map.count(src_neig_v) != 0);
                    assert(guest_map[src_neig_v].count(src_x_y) != 0);
                    int post_mesh_id_src = post_mesh->vert_coor_and_id_map()[src_neig_v];
                    one_source_list.push_back(geodesic::SurfacePoint(&post_mesh->vertices()[post_mesh_id_src]));
                    for (auto ite2 : host_map[dest_x_y])
                    {
                        unsigned dest_neig_v = ite2.first;
                        assert(guest_map.count(dest_neig_v) != 0);
                        assert(guest_map[dest_neig_v].count(dest_x_y) != 0);
                        int post_mesh_id_dest = post_mesh->vert_coor_and_id_map()[dest_neig_v];
                        destinations_list.push_back(geodesic::SurfacePoint(&post_mesh->vertices()[post_mesh_id_dest]));
                    }

                    if (face_exact_one_face_appr_two_vertex_three == 1)
                    {
                        post_algorithm_face_exact.propagate(one_source_list, distance_limit, 1);
                    }
                    else if (face_exact_one_face_appr_two_vertex_three == 2)
                    {
                        post_algorithm_face_appr.propagate(one_source_list, distance_limit, 1);
                    }
                    else if (face_exact_one_face_appr_two_vertex_three == 3)
                    {
                        post_algorithm_vertex.propagate(one_source_list, distance_limit, 1);
                    }

                    for (auto ite2 : host_map[dest_x_y])
                    {
                        unsigned dest_neig_v = ite2.first;
                        unsigned x_neig_y_neig_1, x_neig_y_neig_2;
                        if (src_neig_v <= dest_neig_v)
                        {
                            x_neig_y_neig_1 = src_neig_v;
                            x_neig_y_neig_2 = dest_neig_v;
                        }
                        else
                        {
                            x_neig_y_neig_1 = dest_neig_v;
                            x_neig_y_neig_2 = src_neig_v;
                        }
                        if (path_v_to_v_post_terrain.count(std::make_pair(x_neig_y_neig_1, x_neig_y_neig_2)) == 0)
                        {
                            std::vector<geodesic::SurfacePoint> path;
                            path.clear();
                            int post_mesh_id_dest = post_mesh->vert_coor_and_id_map()[dest_neig_v];
                            geodesic::SurfacePoint p(&post_mesh->vertices()[post_mesh_id_dest]);
                            if (face_exact_one_face_appr_two_vertex_three == 1)
                            {
                                post_algorithm_face_exact.trace_back(p, path);
                            }
                            else if (face_exact_one_face_appr_two_vertex_three == 2)
                            {
                                post_algorithm_face_appr.trace_back(p, path);
                            }
                            else if (face_exact_one_face_appr_two_vertex_three == 3)
                            {
                                post_algorithm_vertex.trace_back(p, path);
                                modify_path(path);
                            }
                            path_v_to_v_post_terrain[std::make_pair(x_neig_y_neig_1, x_neig_y_neig_2)] = path;
                        }
                    }
                }

                post_mesh->add_vertex(org_mesh->vertices()[source_index].getx(),
                                      org_mesh->vertices()[source_index].gety(),
                                      org_mesh->vertices()[source_index].getz());
                geodesic::SurfacePoint real_source(&post_mesh->vertices()[post_mesh->vertices().size() - 1]);

                post_mesh->add_vertex(org_mesh->vertices()[destination_index].getx(),
                                      org_mesh->vertices()[destination_index].gety(),
                                      org_mesh->vertices()[destination_index].getz());
                geodesic::SurfacePoint real_destination(&post_mesh->vertices()[post_mesh->vertices().size() - 1]);

                for (auto ite : host_map[src_x_y])
                {
                    unsigned src_neig_v = ite.first;
                    double src_dest_neig_v = guest_map[src_neig_v][src_x_y];

                    for (auto ite2 : host_map[dest_x_y])
                    {
                        unsigned dest_neig_v = ite2.first;
                        unsigned x_neig_y_neig_1, x_neig_y_neig_2;
                        if (src_neig_v <= dest_neig_v)
                        {
                            x_neig_y_neig_1 = src_neig_v;
                            x_neig_y_neig_2 = dest_neig_v;
                        }
                        else
                        {
                            x_neig_y_neig_1 = dest_neig_v;
                            x_neig_y_neig_2 = src_neig_v;
                        }
                        std::vector<geodesic::SurfacePoint> path = path_v_to_v_post_terrain[std::make_pair(x_neig_y_neig_1, x_neig_y_neig_2)];
                        double dist_dest_neig_v = guest_map[dest_neig_v][dest_x_y];

                        if (distance_result > length(path) + src_dest_neig_v + dist_dest_neig_v)
                        {
                            distance_result = length(path) + src_dest_neig_v + dist_dest_neig_v;
                            path.insert(path.begin(), real_destination);
                            path.push_back(real_source);
                            path_result = path;
                        }
                    }
                }
            }
        }
        distance_result = round(distance_result * 1000000000.0) / 1000000000.0;
        obj_to_other_obj_distance_and_index_list.push_back(std::make_pair(distance_result, destination_index));
    }
    std::sort(obj_to_other_obj_distance_and_index_list.begin(), obj_to_other_obj_distance_and_index_list.end());

    auto stop_knn_or_range_query_time = std::chrono::high_resolution_clock::now();
    auto duration_knn_or_range_query_time = std::chrono::duration_cast<std::chrono::microseconds>(stop_knn_or_range_query_time - start_knn_or_range_query_time);

    compare_d_value_and_range_value(d_value, range);

    auto start_knn_query_time = std::chrono::high_resolution_clock::now();
    knn_or_range_query(1, k_value, range, obj_to_other_obj_distance_and_index_list, knn_list);
    auto stop_knn_query_time = std::chrono::high_resolution_clock::now();
    auto duration_knn_query_time = std::chrono::duration_cast<std::chrono::microseconds>(stop_knn_query_time - start_knn_query_time);
    knn_query_time = duration_knn_or_range_query_time.count() + duration_knn_query_time.count();
    knn_query_time /= 1000;
    knn_query_time *= (d_value / 2000);

    auto start_range_query_time = std::chrono::high_resolution_clock::now();
    knn_or_range_query(2, k_value, range, obj_to_other_obj_distance_and_index_list, range_list);
    auto stop_range_query_time = std::chrono::high_resolution_clock::now();
    auto duration_range_query_time = std::chrono::duration_cast<std::chrono::microseconds>(stop_range_query_time - start_range_query_time);
    range_query_time = duration_knn_or_range_query_time.count() + duration_range_query_time.count();
    range_query_time /= 1000;
    range_query_time *= (range / 1000) * (d_value / 2000);
}

void write_terrain_off_file(std::vector<double> terrain_vertex,
                            std::vector<unsigned> terrain_face)
{
    std::string terrain_write_path = "../build/simplified_terrain.off";
    std::ofstream ofs(&terrain_write_path[0], std::ofstream::trunc);

    ofs << "OFF\n";

    int total_vertex_num = terrain_vertex.size() / 3;
    int total_face_num = terrain_face.size() / 3;
    int total_edge_num = 0;

    ofs << total_vertex_num << " " << total_face_num << " " << total_edge_num << " \n";

    for (unsigned i = 0; i < terrain_vertex.size() / 3; ++i)
    {
        ofs << std::fixed << int(terrain_vertex[i * 3]) << " " << int(terrain_vertex[i * 3 + 1])
            << " " << int(terrain_vertex[i * 3 + 2]) << " \n";
    }

    for (unsigned i = 0; i < terrain_face.size() / 3; ++i)
    {
        ofs << "3 " << terrain_face[i * 3] << " " << terrain_face[i * 3 + 1]
            << " " << terrain_face[i * 3 + 2] << " \n";
    }
}