#include "algorithms.h"
#include <filesystem>

struct input_struct
{
    std::string height_map_or_point_cloud_or_terrain;
    int source_index;
    int destination_index;

    input_struct() {}

    input_struct(std::string _height_map_or_point_cloud_or_terrain,
                 int _source_index, int _destination_index)
    {
        height_map_or_point_cloud_or_terrain = _height_map_or_point_cloud_or_terrain;
        source_index = _source_index;
        destination_index = _destination_index;
    }
};

int main(int argc, char *argv[])
{
    int input_file_index = std::stoi(argv[1]);
    double epsilon = std::stod(argv[2]);
    int run_knn_and_range_query_index = std::stod(argv[3]);

    std::string input_folder = "../input/";

    std::vector<input_struct> input_file;
    input_file.push_back(input_struct("BH_1024", 0, 1023));
    input_file.push_back(input_struct("EP_1024", 0, 1023));
    input_file.push_back(input_struct("EP_10000", 0, 9999));
    input_file.push_back(input_struct("EP_20164", 0, 20163));
    input_file.push_back(input_struct("EP_30276", 0, 30275));
    input_file.push_back(input_struct("EP_40000", 0, 39999));
    input_file.push_back(input_struct("EP_50176", 0, 50175));
    input_file.push_back(input_struct("GF_1024", 0, 1023));
    input_file.push_back(input_struct("LM_1024", 0, 1023));
    input_file.push_back(input_struct("RM_1024", 0, 1023));
    input_file.push_back(input_struct("BH_501264", 0, 501263));
    input_file.push_back(input_struct("BH_5004169", 0, 5004168));
    input_file.push_back(input_struct("BH_10004569", 0, 10004568));
    input_file.push_back(input_struct("BH_15000129", 0, 15000128));
    input_file.push_back(input_struct("BH_20007729", 0, 20007728));
    input_file.push_back(input_struct("BH_25000000", 0, 24999999));
    input_file.push_back(input_struct("EP_501264", 0, 501263));
    input_file.push_back(input_struct("EP_5004169", 0, 5004168));
    input_file.push_back(input_struct("EP_10004569", 0, 10004568));
    input_file.push_back(input_struct("EP_15000129", 0, 15000128));
    input_file.push_back(input_struct("EP_20007729", 0, 20007728));
    input_file.push_back(input_struct("EP_25000000", 0, 24999999));
    input_file.push_back(input_struct("GF_501264", 0, 501263));
    input_file.push_back(input_struct("GF_5004169", 0, 5004168));
    input_file.push_back(input_struct("GF_10004569", 0, 10004568));
    input_file.push_back(input_struct("GF_15000129", 0, 15000128));
    input_file.push_back(input_struct("GF_20007729", 0, 20007728));
    input_file.push_back(input_struct("GF_25000000", 0, 24999999));
    input_file.push_back(input_struct("LM_501264", 0, 501263));
    input_file.push_back(input_struct("LM_5004169", 0, 5004168));
    input_file.push_back(input_struct("LM_10004569", 0, 10004568));
    input_file.push_back(input_struct("LM_15000129", 0, 15000128));
    input_file.push_back(input_struct("LM_20007729", 0, 20007728));
    input_file.push_back(input_struct("LM_25000000", 0, 24999999));
    input_file.push_back(input_struct("RM_501264", 0, 501263));
    input_file.push_back(input_struct("RM_5004169", 0, 5004168));
    input_file.push_back(input_struct("RM_10004569", 0, 10004568));
    input_file.push_back(input_struct("RM_15000129", 0, 15000128));
    input_file.push_back(input_struct("RM_20007729", 0, 20007728));
    input_file.push_back(input_struct("RM_25000000", 0, 24999999));

    std::vector<std::string> input_height_map_file;

    int knn_and_range_query_obj_num = 500;
    int k_value = 5;
    double range = 1000;
    double d_value = 2000;
    double epsilon_prime_ds = 0.1;
    bool run_knn_and_range_query = false;
    if (run_knn_and_range_query_index == 1)
    {
        run_knn_and_range_query = true;
    }

    std::string input_height_map = input_folder + input_file[input_file_index].height_map_or_point_cloud_or_terrain + ".png";
    std::string input_point_cloud = input_folder + input_file[input_file_index].height_map_or_point_cloud_or_terrain + ".xyz";
    std::string input_terrain = input_folder + input_file[input_file_index].height_map_or_point_cloud_or_terrain + ".off";
    int source_index = input_file[input_file_index].source_index;
    int destination_index = input_file[input_file_index].destination_index;

    std::cout.precision(10);

    std::vector<int> pixels;
    int width;
    int height;
    height_map_geodesic::read_org_height_map_from_file(input_height_map, pixels, width, height);
    height_map_geodesic::HeightMap org_height_map;
    org_height_map.initialize_height_map_data(pixels);

    std::vector<double> points;
    point_cloud_geodesic::read_point_cloud_from_file(&input_point_cloud[0], points);
    point_cloud_geodesic::PointCloud point_cloud;
    point_cloud.initialize_point_cloud_data(points);

    std::vector<double> terrain_points;
    std::vector<unsigned> terrain_faces;
    geodesic::read_mesh_from_file(&input_terrain[0], terrain_points, terrain_faces);
    geodesic::Mesh mesh;
    mesh.initialize_mesh_data(terrain_points, terrain_faces);

    std::string write_file_header = input_file[input_file_index].height_map_or_point_cloud_or_terrain + "\t" +
                                    std::to_string(org_height_map.hm_points().size()) + "\t" +
                                    std::to_string(epsilon);

    assert(org_height_map.hm_points().size() > knn_and_range_query_obj_num + 1);

    std::cout << "dataset: " << input_file[input_file_index].height_map_or_point_cloud_or_terrain << "\tdataset_size: " << org_height_map.hm_points().size() << "\tepsilon: " << epsilon << std::endl;
    std::cout << std::endl;

    std::string output_file = "../output/output.txt";
    std::ofstream ofs(output_file, std::ofstream::app);
    ofs << "# dataset\tdataset_size\tepsilon\theight_map_to_point_cloud_or_terrain_time\theight_map_to_point_cloud_or_terrain_memory_usage\tpreprocessing_time\tmemory_usage\toutput_size\tnum_of_cell_point_vertex\tquery_time\tdistance_error_height_map_or_point_cloud\tdistance_error_terrain\tknn_query_time\tknn_error_height_map_or_point_cloud\tknn_error_terrain\trange_query_time\trange_error_height_map_or_point_cloud\trange_error_terrain\n\n";
    ofs.close();

    double height_map_or_point_cloud_exact_distance = 0;
    double terrain_exact_distance = 0;
    std::vector<int> height_map_or_point_cloud_exact_knn_list;
    std::vector<int> terrain_exact_knn_list;
    std::vector<int> height_map_or_point_cloud_exact_range_list;
    std::vector<int> terrain_exact_range_list;
    height_map_or_point_cloud_exact_knn_list.clear();
    terrain_exact_knn_list.clear();
    height_map_or_point_cloud_exact_range_list.clear();
    terrain_exact_range_list.clear();

    calculate_height_map_or_point_cloud_exact_distance(&org_height_map, source_index,
                                                       destination_index, height_map_or_point_cloud_exact_distance);
    calculate_terrain_exact_distance(&org_height_map, source_index,
                                     destination_index, terrain_exact_distance);
    if (run_knn_and_range_query)
    {
        calculate_height_map_or_point_cloud_exact_knn_and_range_query(&org_height_map, source_index,
                                                                      knn_and_range_query_obj_num,
                                                                      k_value, range,
                                                                      height_map_or_point_cloud_exact_knn_list,
                                                                      height_map_or_point_cloud_exact_range_list);
        calculate_terrain_exact_knn_and_range_query(&org_height_map, source_index,
                                                    knn_and_range_query_obj_num, k_value, range,
                                                    terrain_exact_knn_list, terrain_exact_range_list);
    }

    if (epsilon > 0)
    {
        if (input_file_index >= 0 && input_file_index <= 9)
        {
            std::cout << "== TIN_SSimplify_Adapt(HM) and TIN_ESSP_Adapt(HM) on the simplified TIN ==" << std::endl;
            simplified_terrain_face_exact_and_face_appr_and_vertex_with_output(
                output_file, &org_height_map, epsilon, source_index, destination_index,
                1, height_map_or_point_cloud_exact_distance, terrain_exact_distance,
                run_knn_and_range_query, knn_and_range_query_obj_num,
                k_value, range, d_value,
                height_map_or_point_cloud_exact_knn_list, terrain_exact_knn_list,
                height_map_or_point_cloud_exact_range_list, terrain_exact_range_list,
                write_file_header, 1);
            std::cout << std::endl;

            std::cout << "== TIN_NSimplify_Adapt(HM) and TIN_NSP_Adapt(HM) on the simplified TIN ==" << std::endl;
            simplified_terrain_face_exact_and_face_appr_and_vertex_with_output(
                output_file, &org_height_map, epsilon, source_index, destination_index,
                3, height_map_or_point_cloud_exact_distance, terrain_exact_distance,
                run_knn_and_range_query, knn_and_range_query_obj_num,
                k_value, range, d_value,
                height_map_or_point_cloud_exact_knn_list, terrain_exact_knn_list,
                height_map_or_point_cloud_exact_range_list, terrain_exact_range_list,
                write_file_header, 1);
            std::cout << std::endl;

            std::cout << "== PC_Simplify_Adapt(HM) and PC_SP_Adapt(HM) on the simplified point cloud ==" << std::endl;
            simplified_height_map_or_point_cloud_with_output(output_file, &org_height_map, epsilon, epsilon_prime_ds, source_index,
                                                             destination_index, 6,
                                                             height_map_or_point_cloud_exact_distance, terrain_exact_distance,
                                                             run_knn_and_range_query, knn_and_range_query_obj_num,
                                                             k_value, range, d_value,
                                                             height_map_or_point_cloud_exact_knn_list, terrain_exact_knn_list,
                                                             height_map_or_point_cloud_exact_range_list, terrain_exact_range_list,
                                                             write_file_header, 1);
            std::cout << std::endl;

            std::cout << "== HM_Simplify_LS and HM_SP_LS on the simplified height map ==" << std::endl;
            simplified_height_map_or_point_cloud_with_output(output_file, &org_height_map, epsilon, epsilon_prime_ds, source_index,
                                                             destination_index, 4,
                                                             height_map_or_point_cloud_exact_distance, terrain_exact_distance,
                                                             run_knn_and_range_query, knn_and_range_query_obj_num,
                                                             k_value, range, d_value,
                                                             height_map_or_point_cloud_exact_knn_list, terrain_exact_knn_list,
                                                             height_map_or_point_cloud_exact_range_list, terrain_exact_range_list,
                                                             write_file_header, 1);
            std::cout << std::endl;

            std::cout << "== HM_Simplify_LST and HM_SP_LST on the simplified height map ==" << std::endl;
            simplified_height_map_or_point_cloud_with_output(output_file, &org_height_map, epsilon, epsilon_prime_ds, source_index,
                                                             destination_index, 5,
                                                             height_map_or_point_cloud_exact_distance, terrain_exact_distance,
                                                             run_knn_and_range_query, knn_and_range_query_obj_num,
                                                             k_value, range, d_value,
                                                             height_map_or_point_cloud_exact_knn_list, terrain_exact_knn_list,
                                                             height_map_or_point_cloud_exact_range_list, terrain_exact_range_list,
                                                             write_file_header, 1);
            std::cout << std::endl;

            std::cout << "== HM_Simplify_DS and HM_SP_DS on the simplified height map ==" << std::endl;
            simplified_height_map_or_point_cloud_with_output(output_file, &org_height_map, epsilon, epsilon_prime_ds, source_index,
                                                             destination_index, 6,
                                                             height_map_or_point_cloud_exact_distance, terrain_exact_distance,
                                                             run_knn_and_range_query, knn_and_range_query_obj_num,
                                                             k_value, range, d_value,
                                                             height_map_or_point_cloud_exact_knn_list, terrain_exact_knn_list,
                                                             height_map_or_point_cloud_exact_range_list, terrain_exact_range_list,
                                                             write_file_header, 1);
            std::cout << std::endl;

            std::cout << "== TIN_SSimplify_Adapt(PC) and TIN_ESSP_Adapt(PC) on the simplified TIN ==" << std::endl;
            simplified_terrain_face_exact_and_face_appr_and_vertex_with_output(
                output_file, &org_height_map, epsilon, source_index, destination_index,
                1, height_map_or_point_cloud_exact_distance, terrain_exact_distance,
                run_knn_and_range_query, knn_and_range_query_obj_num,
                k_value, range, d_value,
                height_map_or_point_cloud_exact_knn_list, terrain_exact_knn_list,
                height_map_or_point_cloud_exact_range_list, terrain_exact_range_list,
                write_file_header, 2);
            std::cout << std::endl;

            std::cout << "== TIN_NSimplify_Adapt(PC) and TIN_NSP_Adapt(PC) on the simplified TIN ==" << std::endl;
            simplified_terrain_face_exact_and_face_appr_and_vertex_with_output(
                output_file, &org_height_map, epsilon, source_index, destination_index,
                3, height_map_or_point_cloud_exact_distance, terrain_exact_distance,
                run_knn_and_range_query, knn_and_range_query_obj_num,
                k_value, range, d_value,
                height_map_or_point_cloud_exact_knn_list, terrain_exact_knn_list,
                height_map_or_point_cloud_exact_range_list, terrain_exact_range_list,
                write_file_header, 2);
            std::cout << std::endl;

            std::cout << "== PC_Simplify and PC_SP on the simplified point cloud ==" << std::endl;
            simplified_height_map_or_point_cloud_with_output(output_file, &org_height_map, epsilon, epsilon_prime_ds, source_index,
                                                             destination_index, 6,
                                                             height_map_or_point_cloud_exact_distance, terrain_exact_distance,
                                                             run_knn_and_range_query, knn_and_range_query_obj_num,
                                                             k_value, range, d_value,
                                                             height_map_or_point_cloud_exact_knn_list, terrain_exact_knn_list,
                                                             height_map_or_point_cloud_exact_range_list, terrain_exact_range_list,
                                                             write_file_header, 2);
            std::cout << std::endl;

            std::cout << "== TIN_SSimplify and TIN_ESSP on the simplified TIN ==" << std::endl;
            simplified_terrain_face_exact_and_face_appr_and_vertex_with_output(
                output_file, &org_height_map, epsilon, source_index, destination_index,
                1, height_map_or_point_cloud_exact_distance, terrain_exact_distance,
                run_knn_and_range_query, knn_and_range_query_obj_num,
                k_value, range, d_value,
                height_map_or_point_cloud_exact_knn_list, terrain_exact_knn_list,
                height_map_or_point_cloud_exact_range_list, terrain_exact_range_list,
                write_file_header, 3);
            std::cout << std::endl;

            std::cout << "== TIN_NSimplify and TIN_NSP on the simplified TIN ==" << std::endl;
            simplified_terrain_face_exact_and_face_appr_and_vertex_with_output(
                output_file, &org_height_map, epsilon, source_index, destination_index,
                3, height_map_or_point_cloud_exact_distance, terrain_exact_distance,
                run_knn_and_range_query, knn_and_range_query_obj_num,
                k_value, range, d_value,
                height_map_or_point_cloud_exact_knn_list, terrain_exact_knn_list,
                height_map_or_point_cloud_exact_range_list, terrain_exact_range_list,
                write_file_header, 3);
            std::cout << std::endl;

            std::cout << "== PC_Simplify_Adapt(TIN) and PC_SP_Adapt(TIN) on the simplified point cloud ==" << std::endl;
            simplified_height_map_or_point_cloud_with_output(output_file, &org_height_map, epsilon, epsilon_prime_ds, source_index,
                                                             destination_index, 6,
                                                             height_map_or_point_cloud_exact_distance, terrain_exact_distance,
                                                             run_knn_and_range_query, knn_and_range_query_obj_num,
                                                             k_value, range, d_value,
                                                             height_map_or_point_cloud_exact_knn_list, terrain_exact_knn_list,
                                                             height_map_or_point_cloud_exact_range_list, terrain_exact_range_list,
                                                             write_file_header, 3);
            std::cout << std::endl;
        }

        std::cout << "== HM_Simplify and HM_SP on the simplified height map ==" << std::endl;
        simplified_height_map_or_point_cloud_with_output(output_file, &org_height_map, epsilon, epsilon_prime_ds, source_index,
                                                         destination_index, 1,
                                                         height_map_or_point_cloud_exact_distance, terrain_exact_distance,
                                                         run_knn_and_range_query, knn_and_range_query_obj_num,
                                                         k_value, range, d_value,
                                                         height_map_or_point_cloud_exact_knn_list, terrain_exact_knn_list,
                                                         height_map_or_point_cloud_exact_range_list, terrain_exact_range_list,
                                                         write_file_header, 1);
        std::cout << std::endl;

        std::cout << "== HM_Simplify_LQT1 and HM_SP_LQT1 on the simplified height map ==" << std::endl;
        simplified_height_map_or_point_cloud_with_output(output_file, &org_height_map, epsilon, epsilon_prime_ds, source_index,
                                                         destination_index, 2,
                                                         height_map_or_point_cloud_exact_distance, terrain_exact_distance,
                                                         run_knn_and_range_query, knn_and_range_query_obj_num,
                                                         k_value, range, d_value,
                                                         height_map_or_point_cloud_exact_knn_list, terrain_exact_knn_list,
                                                         height_map_or_point_cloud_exact_range_list, terrain_exact_range_list,
                                                         write_file_header, 1);
        std::cout << std::endl;

        std::cout << "== HM_Simplify_LQT2 and HM_SP_LQT2 on the simplified height map ==" << std::endl;
        simplified_height_map_or_point_cloud_with_output(output_file, &org_height_map, epsilon, epsilon_prime_ds, source_index,
                                                         destination_index, 3,
                                                         height_map_or_point_cloud_exact_distance, terrain_exact_distance,
                                                         run_knn_and_range_query, knn_and_range_query_obj_num,
                                                         k_value, range, d_value,
                                                         height_map_or_point_cloud_exact_knn_list, terrain_exact_knn_list,
                                                         height_map_or_point_cloud_exact_range_list, terrain_exact_range_list,
                                                         write_file_header, 1);
        std::cout << std::endl;

        std::cout << "== TIN_ASSP_Adapt(HM) ==" << std::endl;
        terrain_face_exact_and_face_appr_and_vertex_with_output(
            output_file, &org_height_map, epsilon, source_index,
            destination_index, 2, height_map_or_point_cloud_exact_distance,
            terrain_exact_distance, run_knn_and_range_query,
            knn_and_range_query_obj_num, k_value, range, d_value,
            height_map_or_point_cloud_exact_knn_list, terrain_exact_knn_list,
            height_map_or_point_cloud_exact_range_list, terrain_exact_range_list,
            write_file_header, 1);
        std::cout << std::endl;

        std::cout << "== HM_Simplify_Adapt(PC) and HM_SP_Adapt(PC) on the simplified height map ==" << std::endl;
        simplified_height_map_or_point_cloud_with_output(output_file, &org_height_map, epsilon, epsilon_prime_ds, source_index,
                                                         destination_index, 1,
                                                         height_map_or_point_cloud_exact_distance, terrain_exact_distance,
                                                         run_knn_and_range_query, knn_and_range_query_obj_num,
                                                         k_value, range, d_value,
                                                         height_map_or_point_cloud_exact_knn_list, terrain_exact_knn_list,
                                                         height_map_or_point_cloud_exact_range_list, terrain_exact_range_list,
                                                         write_file_header, 2);
        std::cout << std::endl;

        std::cout << "== TIN_ASSP_Adapt(PC) ==" << std::endl;
        terrain_face_exact_and_face_appr_and_vertex_with_output(
            output_file, &org_height_map, epsilon, source_index,
            destination_index, 2, height_map_or_point_cloud_exact_distance,
            terrain_exact_distance, run_knn_and_range_query,
            knn_and_range_query_obj_num, k_value, range, d_value,
            height_map_or_point_cloud_exact_knn_list, terrain_exact_knn_list,
            height_map_or_point_cloud_exact_range_list, terrain_exact_range_list,
            write_file_header, 2);
        std::cout << std::endl;

        std::cout << "== HM_Simplify_Adapt(TIN) and HM_SP_Adapt(TIN) on the simplified height map ==" << std::endl;
        simplified_height_map_or_point_cloud_with_output(output_file, &org_height_map, epsilon, epsilon_prime_ds, source_index,
                                                         destination_index, 1,
                                                         height_map_or_point_cloud_exact_distance, terrain_exact_distance,
                                                         run_knn_and_range_query, knn_and_range_query_obj_num,
                                                         k_value, range, d_value,
                                                         height_map_or_point_cloud_exact_knn_list, terrain_exact_knn_list,
                                                         height_map_or_point_cloud_exact_range_list, terrain_exact_range_list,
                                                         write_file_header, 3);
        std::cout << std::endl;

        std::cout << "== TIN_ASSP ==" << std::endl;
        terrain_face_exact_and_face_appr_and_vertex_with_output(
            output_file, &org_height_map, epsilon, source_index,
            destination_index, 2, height_map_or_point_cloud_exact_distance,
            terrain_exact_distance, run_knn_and_range_query,
            knn_and_range_query_obj_num, k_value, range, d_value,
            height_map_or_point_cloud_exact_knn_list, terrain_exact_knn_list,
            height_map_or_point_cloud_exact_range_list, terrain_exact_range_list,
            write_file_header, 3);
        std::cout << std::endl;
    }
    else if (epsilon == 0)
    {
        std::cout << "== TIN_ESSP_Adapt(HM) on the original height map ==" << std::endl;
        terrain_face_exact_and_face_appr_and_vertex_with_output(
            output_file, &org_height_map, epsilon, source_index,
            destination_index, 1, height_map_or_point_cloud_exact_distance,
            terrain_exact_distance, run_knn_and_range_query,
            knn_and_range_query_obj_num, k_value, range, d_value,
            height_map_or_point_cloud_exact_knn_list, terrain_exact_knn_list,
            height_map_or_point_cloud_exact_range_list, terrain_exact_range_list,
            write_file_header, 1);
        std::cout << std::endl;

        std::cout << "== TIN_NSP_Adapt(HM) on the original height map ==" << std::endl;
        terrain_face_exact_and_face_appr_and_vertex_with_output(
            output_file, &org_height_map, epsilon, source_index,
            destination_index, 3, height_map_or_point_cloud_exact_distance,
            terrain_exact_distance, run_knn_and_range_query,
            knn_and_range_query_obj_num, k_value, range, d_value,
            height_map_or_point_cloud_exact_knn_list, terrain_exact_knn_list,
            height_map_or_point_cloud_exact_range_list, terrain_exact_range_list,
            write_file_header, 1);
        std::cout << std::endl;

        std::cout << "== PC_SP_Adapt(HM) on the original height map ==" << std::endl;
        height_map_or_point_cloud_with_output(output_file, &org_height_map, source_index, destination_index,
                                              height_map_or_point_cloud_exact_distance, terrain_exact_distance,
                                              run_knn_and_range_query, knn_and_range_query_obj_num,
                                              k_value, range, d_value,
                                              height_map_or_point_cloud_exact_knn_list, terrain_exact_knn_list,
                                              height_map_or_point_cloud_exact_range_list, terrain_exact_range_list, 2,
                                              write_file_header, 1);
        std::cout << std::endl;

        std::cout << "== HM_SP on the original height map ==" << std::endl;
        height_map_or_point_cloud_with_output(output_file, &org_height_map, source_index, destination_index,
                                              height_map_or_point_cloud_exact_distance, terrain_exact_distance,
                                              run_knn_and_range_query, knn_and_range_query_obj_num,
                                              k_value, range, d_value,
                                              height_map_or_point_cloud_exact_knn_list, terrain_exact_knn_list,
                                              height_map_or_point_cloud_exact_range_list, terrain_exact_range_list, 1,
                                              write_file_header, 1);
        std::cout << std::endl;

        std::cout << "== TIN_ESSP_Adapt(PC) on the original point cloud ==" << std::endl;
        terrain_face_exact_and_face_appr_and_vertex_with_output(
            output_file, &org_height_map, epsilon, source_index,
            destination_index, 1, height_map_or_point_cloud_exact_distance,
            terrain_exact_distance, run_knn_and_range_query,
            knn_and_range_query_obj_num, k_value, range, d_value,
            height_map_or_point_cloud_exact_knn_list, terrain_exact_knn_list,
            height_map_or_point_cloud_exact_range_list, terrain_exact_range_list,
            write_file_header, 2);
        std::cout << std::endl;

        std::cout << "== TIN_NSP_Adapt(PC) on the original point cloud ==" << std::endl;
        terrain_face_exact_and_face_appr_and_vertex_with_output(
            output_file, &org_height_map, epsilon, source_index,
            destination_index, 3, height_map_or_point_cloud_exact_distance,
            terrain_exact_distance, run_knn_and_range_query,
            knn_and_range_query_obj_num, k_value, range, d_value,
            height_map_or_point_cloud_exact_knn_list, terrain_exact_knn_list,
            height_map_or_point_cloud_exact_range_list, terrain_exact_range_list,
            write_file_header, 2);
        std::cout << std::endl;

        std::cout << "== PC_SP on the original point cloud ==" << std::endl;
        height_map_or_point_cloud_with_output(output_file, &org_height_map, source_index, destination_index,
                                              height_map_or_point_cloud_exact_distance, terrain_exact_distance,
                                              run_knn_and_range_query, knn_and_range_query_obj_num,
                                              k_value, range, d_value,
                                              height_map_or_point_cloud_exact_knn_list, terrain_exact_knn_list,
                                              height_map_or_point_cloud_exact_range_list, terrain_exact_range_list, 2,
                                              write_file_header, 2);
        std::cout << std::endl;

        std::cout << "== HM_SP_Adapt(PC) on the original point cloud ==" << std::endl;
        height_map_or_point_cloud_with_output(output_file, &org_height_map, source_index, destination_index,
                                              height_map_or_point_cloud_exact_distance, terrain_exact_distance,
                                              run_knn_and_range_query, knn_and_range_query_obj_num,
                                              k_value, range, d_value,
                                              height_map_or_point_cloud_exact_knn_list, terrain_exact_knn_list,
                                              height_map_or_point_cloud_exact_range_list, terrain_exact_range_list, 1,
                                              write_file_header, 2);
        std::cout << std::endl;

        std::cout << "== TIN_ESSP on the original TIN ==" << std::endl;
        terrain_face_exact_and_face_appr_and_vertex_with_output(
            output_file, &org_height_map, epsilon, source_index,
            destination_index, 1, height_map_or_point_cloud_exact_distance,
            terrain_exact_distance, run_knn_and_range_query,
            knn_and_range_query_obj_num, k_value, range, d_value,
            height_map_or_point_cloud_exact_knn_list, terrain_exact_knn_list,
            height_map_or_point_cloud_exact_range_list, terrain_exact_range_list,
            write_file_header, 3);
        std::cout << std::endl;

        std::cout << "== TIN_NSP on the original TIN ==" << std::endl;
        terrain_face_exact_and_face_appr_and_vertex_with_output(
            output_file, &org_height_map, epsilon, source_index,
            destination_index, 3, height_map_or_point_cloud_exact_distance,
            terrain_exact_distance, run_knn_and_range_query,
            knn_and_range_query_obj_num, k_value, range, d_value,
            height_map_or_point_cloud_exact_knn_list, terrain_exact_knn_list,
            height_map_or_point_cloud_exact_range_list, terrain_exact_range_list,
            write_file_header, 3);
        std::cout << std::endl;

        std::cout << "== PC_SP_Adapt(TIN) on the original TIN ==" << std::endl;
        height_map_or_point_cloud_with_output(output_file, &org_height_map, source_index, destination_index,
                                              height_map_or_point_cloud_exact_distance, terrain_exact_distance,
                                              run_knn_and_range_query, knn_and_range_query_obj_num,
                                              k_value, range, d_value,
                                              height_map_or_point_cloud_exact_knn_list, terrain_exact_knn_list,
                                              height_map_or_point_cloud_exact_range_list, terrain_exact_range_list, 2,
                                              write_file_header, 3);
        std::cout << std::endl;

        std::cout << "== HM_SP_Adapt(TIN) on the original TIN ==" << std::endl;
        height_map_or_point_cloud_with_output(output_file, &org_height_map, source_index, destination_index,
                                              height_map_or_point_cloud_exact_distance, terrain_exact_distance,
                                              run_knn_and_range_query, knn_and_range_query_obj_num,
                                              k_value, range, d_value,
                                              height_map_or_point_cloud_exact_knn_list, terrain_exact_knn_list,
                                              height_map_or_point_cloud_exact_range_list, terrain_exact_range_list, 1,
                                              write_file_header, 3);
        std::cout << std::endl;
    }

    std::__fs::filesystem::remove("../build/temp_terrain.off");
    std::__fs::filesystem::remove("../build/simplified_terrain.off");
}