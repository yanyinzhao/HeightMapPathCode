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
    input_file.push_back(input_struct("BH_1014", 0, 1013));
    input_file.push_back(input_struct("BH_10086", 0, 10085));
    input_file.push_back(input_struct("EP_10062", 0, 10061));
    input_file.push_back(input_struct("EP_20130", 0, 20129));
    input_file.push_back(input_struct("EP_30098", 0, 30097));
    input_file.push_back(input_struct("EP_40076", 0, 40075));
    input_file.push_back(input_struct("EP_50373", 0, 50372));
    input_file.push_back(input_struct("GF_10092", 0, 10091));
    input_file.push_back(input_struct("LM_10092", 0, 10091));
    input_file.push_back(input_struct("RM_10092", 0, 10091));
    input_file.push_back(input_struct("BH_500835", 0, 500834));
    input_file.push_back(input_struct("BH_1000414", 0, 1000413));
    input_file.push_back(input_struct("BH_1500996", 0, 1500995));
    input_file.push_back(input_struct("BH_2001610", 0, 2001609));
    input_file.push_back(input_struct("BH_2502596", 0, 2502595));
    input_file.push_back(input_struct("EP_500384", 0, 500383));
    input_file.push_back(input_struct("EP_1001040", 0, 1001039));
    input_file.push_back(input_struct("EP_1501578", 0, 1501577));
    input_file.push_back(input_struct("EP_2001536", 0, 2001535));
    input_file.push_back(input_struct("EP_2500560", 0, 2500559));
    input_file.push_back(input_struct("GF_500208", 0, 50027));
    input_file.push_back(input_struct("GF_1000518", 0, 1000517));
    input_file.push_back(input_struct("GF_1501668", 0, 1501667));
    input_file.push_back(input_struct("GF_2000832", 0, 2000831));
    input_file.push_back(input_struct("GF_2502075", 0, 2502074));
    input_file.push_back(input_struct("LM_500208", 0, 50027));
    input_file.push_back(input_struct("LM_1000518", 0, 1000517));
    input_file.push_back(input_struct("LM_1501668", 0, 1501667));
    input_file.push_back(input_struct("LM_2000832", 0, 2000831));
    input_file.push_back(input_struct("LM_2502075", 0, 2502074));
    input_file.push_back(input_struct("RM_500208", 0, 50027));
    input_file.push_back(input_struct("RM_1000518", 0, 1000517));
    input_file.push_back(input_struct("RM_1501668", 0, 1501667));
    input_file.push_back(input_struct("RM_2000832", 0, 2000831));
    input_file.push_back(input_struct("RM_2502075", 0, 2502074));

    std::vector<std::string> input_height_map_file;

    int knn_and_range_query_obj_num = 1000;
    int k_value = 5;
    double range = 1000;
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
    ofs << "# dataset\tdataset_size\tepsilon\theight_map_to_point_cloud_or_terrain_time\theight_map_to_point_cloud_or_terrain_memory_usage\tsimplification_time\tmemory_usage\toutput_size\tquery_time\tdistance_error_height_map_or_point_cloud\tdistance_error_terrain\tknn_query_time\tknn_error_height_map_or_point_cloud\tknn_error_terrain\trange_query_time\trange_error_height_map_or_point_cloud\trange_error_terrain\n\n";
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

    if (input_file_index >= 0 && input_file_index <= 9)
    {
        std::cout << "== TIN_SurSimQ_Adapt(HM) ==" << std::endl;
        simplified_terrain_face_exact_and_face_appr_and_vertex_with_output(
            output_file, &org_height_map, epsilon, source_index, destination_index,
            1, height_map_or_point_cloud_exact_distance, terrain_exact_distance,
            run_knn_and_range_query, knn_and_range_query_obj_num,
            k_value, range,
            height_map_or_point_cloud_exact_knn_list, terrain_exact_knn_list,
            height_map_or_point_cloud_exact_range_list, terrain_exact_range_list,
            write_file_header, 1);
        std::cout << std::endl;

        std::cout << "== TIN_NetSimQ_Adapt(HM) ==" << std::endl;
        simplified_terrain_face_exact_and_face_appr_and_vertex_with_output(
            output_file, &org_height_map, epsilon, source_index, destination_index,
            3, height_map_or_point_cloud_exact_distance, terrain_exact_distance,
            run_knn_and_range_query, knn_and_range_query_obj_num,
            k_value, range,
            height_map_or_point_cloud_exact_knn_list, terrain_exact_knn_list,
            height_map_or_point_cloud_exact_range_list, terrain_exact_range_list,
            write_file_header, 1);
        std::cout << std::endl;

        std::cout << "== HM_MesSimQ_LS ==" << std::endl;
        simplified_height_map_or_point_cloud_with_output(output_file, &org_height_map, epsilon, source_index,
                                                         destination_index, 4,
                                                         height_map_or_point_cloud_exact_distance, terrain_exact_distance,
                                                         run_knn_and_range_query, knn_and_range_query_obj_num,
                                                         k_value, range,
                                                         height_map_or_point_cloud_exact_knn_list, terrain_exact_knn_list,
                                                         height_map_or_point_cloud_exact_range_list, terrain_exact_range_list,
                                                         write_file_header, 1);
        std::cout << std::endl;

        std::cout << "== HM_MesSimQ_LST ==" << std::endl;
        simplified_height_map_or_point_cloud_with_output(output_file, &org_height_map, epsilon, source_index,
                                                         destination_index, 5,
                                                         height_map_or_point_cloud_exact_distance, terrain_exact_distance,
                                                         run_knn_and_range_query, knn_and_range_query_obj_num,
                                                         k_value, range,
                                                         height_map_or_point_cloud_exact_knn_list, terrain_exact_knn_list,
                                                         height_map_or_point_cloud_exact_range_list, terrain_exact_range_list,
                                                         write_file_header, 1);
        std::cout << std::endl;

        std::cout << "== TIN_SurSimQ_Adapt(PC) ==" << std::endl;
        simplified_terrain_face_exact_and_face_appr_and_vertex_with_output(
            output_file, &org_height_map, epsilon, source_index, destination_index,
            1, height_map_or_point_cloud_exact_distance, terrain_exact_distance,
            run_knn_and_range_query, knn_and_range_query_obj_num,
            k_value, range,
            height_map_or_point_cloud_exact_knn_list, terrain_exact_knn_list,
            height_map_or_point_cloud_exact_range_list, terrain_exact_range_list,
            write_file_header, 2);
        std::cout << std::endl;

        std::cout << "== TIN_NetSimQ_Adapt(PC) ==" << std::endl;
        simplified_terrain_face_exact_and_face_appr_and_vertex_with_output(
            output_file, &org_height_map, epsilon, source_index, destination_index,
            3, height_map_or_point_cloud_exact_distance, terrain_exact_distance,
            run_knn_and_range_query, knn_and_range_query_obj_num,
            k_value, range,
            height_map_or_point_cloud_exact_knn_list, terrain_exact_knn_list,
            height_map_or_point_cloud_exact_range_list, terrain_exact_range_list,
            write_file_header, 2);
        std::cout << std::endl;

        std::cout << "== TIN_SurSimQ ==" << std::endl;
        simplified_terrain_face_exact_and_face_appr_and_vertex_with_output(
            output_file, &org_height_map, epsilon, source_index, destination_index,
            1, height_map_or_point_cloud_exact_distance, terrain_exact_distance,
            run_knn_and_range_query, knn_and_range_query_obj_num,
            k_value, range,
            height_map_or_point_cloud_exact_knn_list, terrain_exact_knn_list,
            height_map_or_point_cloud_exact_range_list, terrain_exact_range_list,
            write_file_header, 3);
        std::cout << std::endl;

        std::cout << "== TIN_NetSimQ ==" << std::endl;
        simplified_terrain_face_exact_and_face_appr_and_vertex_with_output(
            output_file, &org_height_map, epsilon, source_index, destination_index,
            3, height_map_or_point_cloud_exact_distance, terrain_exact_distance,
            run_knn_and_range_query, knn_and_range_query_obj_num,
            k_value, range,
            height_map_or_point_cloud_exact_knn_list, terrain_exact_knn_list,
            height_map_or_point_cloud_exact_range_list, terrain_exact_range_list,
            write_file_header, 3);
        std::cout << std::endl;
    }

    std::cout << "== PC_MesSimQ_Adapt(HM) ==" << std::endl;
    simplified_height_map_or_point_cloud_with_output(output_file, &org_height_map, epsilon, source_index,
                                                     destination_index, 6,
                                                     height_map_or_point_cloud_exact_distance, terrain_exact_distance,
                                                     run_knn_and_range_query, knn_and_range_query_obj_num,
                                                     k_value, range,
                                                     height_map_or_point_cloud_exact_knn_list, terrain_exact_knn_list,
                                                     height_map_or_point_cloud_exact_range_list, terrain_exact_range_list,
                                                     write_file_header, 1);
    std::cout << std::endl;

    std::cout << "== HM_MesSimQ ==" << std::endl;
    simplified_height_map_or_point_cloud_with_output(output_file, &org_height_map, epsilon, source_index,
                                                     destination_index, 1,
                                                     height_map_or_point_cloud_exact_distance, terrain_exact_distance,
                                                     run_knn_and_range_query, knn_and_range_query_obj_num,
                                                     k_value, range,
                                                     height_map_or_point_cloud_exact_knn_list, terrain_exact_knn_list,
                                                     height_map_or_point_cloud_exact_range_list, terrain_exact_range_list,
                                                     write_file_header, 1);
    std::cout << std::endl;

    std::cout << "== HM_MesSimQ_LQT1 ==" << std::endl;
    simplified_height_map_or_point_cloud_with_output(output_file, &org_height_map, epsilon, source_index,
                                                     destination_index, 2,
                                                     height_map_or_point_cloud_exact_distance, terrain_exact_distance,
                                                     run_knn_and_range_query, knn_and_range_query_obj_num,
                                                     k_value, range,
                                                     height_map_or_point_cloud_exact_knn_list, terrain_exact_knn_list,
                                                     height_map_or_point_cloud_exact_range_list, terrain_exact_range_list,
                                                     write_file_header, 1);
    std::cout << std::endl;

    std::cout << "== HM_MesSimQ_LQT2 ==" << std::endl;
    simplified_height_map_or_point_cloud_with_output(output_file, &org_height_map, epsilon, source_index,
                                                     destination_index, 3,
                                                     height_map_or_point_cloud_exact_distance, terrain_exact_distance,
                                                     run_knn_and_range_query, knn_and_range_query_obj_num,
                                                     k_value, range,
                                                     height_map_or_point_cloud_exact_knn_list, terrain_exact_knn_list,
                                                     height_map_or_point_cloud_exact_range_list, terrain_exact_range_list,
                                                     write_file_header, 1);
    std::cout << std::endl;

    std::cout << "== TIN_UnfQ_Adapt(HM) ==" << std::endl;
    terrain_face_exact_and_face_appr_and_vertex_with_output(
        output_file, &org_height_map, epsilon, source_index,
        destination_index, 1, height_map_or_point_cloud_exact_distance,
        terrain_exact_distance, run_knn_and_range_query,
        knn_and_range_query_obj_num, k_value, range,
        height_map_or_point_cloud_exact_knn_list, terrain_exact_knn_list,
        height_map_or_point_cloud_exact_range_list, terrain_exact_range_list,
        write_file_header, 1);
    std::cout << std::endl;

    std::cout << "== TIN_SteQ_Adapt(HM) ==" << std::endl;
    terrain_face_exact_and_face_appr_and_vertex_with_output(
        output_file, &org_height_map, epsilon, source_index,
        destination_index, 2, height_map_or_point_cloud_exact_distance,
        terrain_exact_distance, run_knn_and_range_query,
        knn_and_range_query_obj_num, k_value, range,
        height_map_or_point_cloud_exact_knn_list, terrain_exact_knn_list,
        height_map_or_point_cloud_exact_range_list, terrain_exact_range_list,
        write_file_header, 1);
    std::cout << std::endl;

    std::cout << "== TIN_DijQ_Adapt(HM) ==" << std::endl;
    terrain_face_exact_and_face_appr_and_vertex_with_output(
        output_file, &org_height_map, epsilon, source_index,
        destination_index, 3, height_map_or_point_cloud_exact_distance,
        terrain_exact_distance, run_knn_and_range_query,
        knn_and_range_query_obj_num, k_value, range,
        height_map_or_point_cloud_exact_knn_list, terrain_exact_knn_list,
        height_map_or_point_cloud_exact_range_list, terrain_exact_range_list,
        write_file_header, 1);
    std::cout << std::endl;

    std::cout << "== PC_ConQ_Adapt(HM) ==" << std::endl;
    height_map_or_point_cloud_with_output(output_file, &org_height_map, source_index, destination_index,
                                          height_map_or_point_cloud_exact_distance, terrain_exact_distance,
                                          run_knn_and_range_query, knn_and_range_query_obj_num,
                                          k_value, range,
                                          height_map_or_point_cloud_exact_knn_list, terrain_exact_knn_list,
                                          height_map_or_point_cloud_exact_range_list, terrain_exact_range_list, 2,
                                          write_file_header, 1);
    std::cout << std::endl;

    std::cout << "== HM_EffQ ==" << std::endl;
    height_map_or_point_cloud_with_output(output_file, &org_height_map, source_index, destination_index,
                                          height_map_or_point_cloud_exact_distance, terrain_exact_distance,
                                          run_knn_and_range_query, knn_and_range_query_obj_num,
                                          k_value, range,
                                          height_map_or_point_cloud_exact_knn_list, terrain_exact_knn_list,
                                          height_map_or_point_cloud_exact_range_list, terrain_exact_range_list, 1,
                                          write_file_header, 1);
    std::cout << std::endl;

    std::cout << "== PC_MesSimQ ==" << std::endl;
    simplified_height_map_or_point_cloud_with_output(output_file, &org_height_map, epsilon, source_index,
                                                     destination_index, 6,
                                                     height_map_or_point_cloud_exact_distance, terrain_exact_distance,
                                                     run_knn_and_range_query, knn_and_range_query_obj_num,
                                                     k_value, range,
                                                     height_map_or_point_cloud_exact_knn_list, terrain_exact_knn_list,
                                                     height_map_or_point_cloud_exact_range_list, terrain_exact_range_list,
                                                     write_file_header, 2);
    std::cout << std::endl;

    std::cout << "== HM_MesSimQ_Adapt(PC) ==" << std::endl;
    simplified_height_map_or_point_cloud_with_output(output_file, &org_height_map, epsilon, source_index,
                                                     destination_index, 1,
                                                     height_map_or_point_cloud_exact_distance, terrain_exact_distance,
                                                     run_knn_and_range_query, knn_and_range_query_obj_num,
                                                     k_value, range,
                                                     height_map_or_point_cloud_exact_knn_list, terrain_exact_knn_list,
                                                     height_map_or_point_cloud_exact_range_list, terrain_exact_range_list,
                                                     write_file_header, 2);
    std::cout << std::endl;

    std::cout << "== TIN_UnfQ_Adapt(PC) ==" << std::endl;
    terrain_face_exact_and_face_appr_and_vertex_with_output(
        output_file, &org_height_map, epsilon, source_index,
        destination_index, 1, height_map_or_point_cloud_exact_distance,
        terrain_exact_distance, run_knn_and_range_query,
        knn_and_range_query_obj_num, k_value, range,
        height_map_or_point_cloud_exact_knn_list, terrain_exact_knn_list,
        height_map_or_point_cloud_exact_range_list, terrain_exact_range_list,
        write_file_header, 2);
    std::cout << std::endl;

    std::cout << "== TIN_SteQ_Adapt(PC) ==" << std::endl;
    terrain_face_exact_and_face_appr_and_vertex_with_output(
        output_file, &org_height_map, epsilon, source_index,
        destination_index, 2, height_map_or_point_cloud_exact_distance,
        terrain_exact_distance, run_knn_and_range_query,
        knn_and_range_query_obj_num, k_value, range,
        height_map_or_point_cloud_exact_knn_list, terrain_exact_knn_list,
        height_map_or_point_cloud_exact_range_list, terrain_exact_range_list,
        write_file_header, 2);
    std::cout << std::endl;

    std::cout << "== TIN_DijQ_Adapt(PC) ==" << std::endl;
    terrain_face_exact_and_face_appr_and_vertex_with_output(
        output_file, &org_height_map, epsilon, source_index,
        destination_index, 3, height_map_or_point_cloud_exact_distance,
        terrain_exact_distance, run_knn_and_range_query,
        knn_and_range_query_obj_num, k_value, range,
        height_map_or_point_cloud_exact_knn_list, terrain_exact_knn_list,
        height_map_or_point_cloud_exact_range_list, terrain_exact_range_list,
        write_file_header, 2);
    std::cout << std::endl;

    std::cout << "== PC_ConQ ==" << std::endl;
    height_map_or_point_cloud_with_output(output_file, &org_height_map, source_index, destination_index,
                                          height_map_or_point_cloud_exact_distance, terrain_exact_distance,
                                          run_knn_and_range_query, knn_and_range_query_obj_num,
                                          k_value, range,
                                          height_map_or_point_cloud_exact_knn_list, terrain_exact_knn_list,
                                          height_map_or_point_cloud_exact_range_list, terrain_exact_range_list, 2,
                                          write_file_header, 2);
    std::cout << std::endl;

    std::cout << "== HM_EffQ_Adapt(PC) ==" << std::endl;
    height_map_or_point_cloud_with_output(output_file, &org_height_map, source_index, destination_index,
                                          height_map_or_point_cloud_exact_distance, terrain_exact_distance,
                                          run_knn_and_range_query, knn_and_range_query_obj_num,
                                          k_value, range,
                                          height_map_or_point_cloud_exact_knn_list, terrain_exact_knn_list,
                                          height_map_or_point_cloud_exact_range_list, terrain_exact_range_list, 1,
                                          write_file_header, 2);
    std::cout << std::endl;

    std::cout << "== PC_MesSimQ_Adapt(TIN) ==" << std::endl;
    simplified_height_map_or_point_cloud_with_output(output_file, &org_height_map, epsilon, source_index,
                                                     destination_index, 6,
                                                     height_map_or_point_cloud_exact_distance, terrain_exact_distance,
                                                     run_knn_and_range_query, knn_and_range_query_obj_num,
                                                     k_value, range,
                                                     height_map_or_point_cloud_exact_knn_list, terrain_exact_knn_list,
                                                     height_map_or_point_cloud_exact_range_list, terrain_exact_range_list,
                                                     write_file_header, 3);
    std::cout << std::endl;

    std::cout << "== HM_MesSimQ_Adapt(TIN) ==" << std::endl;
    simplified_height_map_or_point_cloud_with_output(output_file, &org_height_map, epsilon, source_index,
                                                     destination_index, 1,
                                                     height_map_or_point_cloud_exact_distance, terrain_exact_distance,
                                                     run_knn_and_range_query, knn_and_range_query_obj_num,
                                                     k_value, range,
                                                     height_map_or_point_cloud_exact_knn_list, terrain_exact_knn_list,
                                                     height_map_or_point_cloud_exact_range_list, terrain_exact_range_list,
                                                     write_file_header, 3);
    std::cout << std::endl;

    std::cout << "== TIN_UnfQ ==" << std::endl;
    terrain_face_exact_and_face_appr_and_vertex_with_output(
        output_file, &org_height_map, epsilon, source_index,
        destination_index, 1, height_map_or_point_cloud_exact_distance,
        terrain_exact_distance, run_knn_and_range_query,
        knn_and_range_query_obj_num, k_value, range,
        height_map_or_point_cloud_exact_knn_list, terrain_exact_knn_list,
        height_map_or_point_cloud_exact_range_list, terrain_exact_range_list,
        write_file_header, 3);
    std::cout << std::endl;

    std::cout << "== TIN_SteQ ==" << std::endl;
    terrain_face_exact_and_face_appr_and_vertex_with_output(
        output_file, &org_height_map, epsilon, source_index,
        destination_index, 2, height_map_or_point_cloud_exact_distance,
        terrain_exact_distance, run_knn_and_range_query,
        knn_and_range_query_obj_num, k_value, range,
        height_map_or_point_cloud_exact_knn_list, terrain_exact_knn_list,
        height_map_or_point_cloud_exact_range_list, terrain_exact_range_list,
        write_file_header, 3);
    std::cout << std::endl;

    std::cout << "== TIN_DijQ ==" << std::endl;
    terrain_face_exact_and_face_appr_and_vertex_with_output(
        output_file, &org_height_map, epsilon, source_index,
        destination_index, 3, height_map_or_point_cloud_exact_distance,
        terrain_exact_distance, run_knn_and_range_query,
        knn_and_range_query_obj_num, k_value, range,
        height_map_or_point_cloud_exact_knn_list, terrain_exact_knn_list,
        height_map_or_point_cloud_exact_range_list, terrain_exact_range_list,
        write_file_header, 3);
    std::cout << std::endl;

    std::cout << "== PC_ConQ_Adapt(TIN) ==" << std::endl;
    height_map_or_point_cloud_with_output(output_file, &org_height_map, source_index, destination_index,
                                          height_map_or_point_cloud_exact_distance, terrain_exact_distance,
                                          run_knn_and_range_query, knn_and_range_query_obj_num,
                                          k_value, range,
                                          height_map_or_point_cloud_exact_knn_list, terrain_exact_knn_list,
                                          height_map_or_point_cloud_exact_range_list, terrain_exact_range_list, 2,
                                          write_file_header, 3);
    std::cout << std::endl;

    std::cout << "== HM_EffQ_Adapt(TIN) ==" << std::endl;
    height_map_or_point_cloud_with_output(output_file, &org_height_map, source_index, destination_index,
                                          height_map_or_point_cloud_exact_distance, terrain_exact_distance,
                                          run_knn_and_range_query, knn_and_range_query_obj_num,
                                          k_value, range,
                                          height_map_or_point_cloud_exact_knn_list, terrain_exact_knn_list,
                                          height_map_or_point_cloud_exact_range_list, terrain_exact_range_list, 1,
                                          write_file_header, 3);
    std::cout << std::endl;

    std::__fs::filesystem::remove("../build/temp_terrain.off");
    std::__fs::filesystem::remove("../build/simplified_terrain.off");
}