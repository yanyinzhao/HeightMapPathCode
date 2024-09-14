#include <assert.h>
#include <cstddef>
#include <vector>
#include <math.h>
#include <memory>
#include <limits>
#include <fstream>
#include <cmath>
#include <iostream>
#include <algorithm>
#include <ctime>
#include <set>
#include <unordered_map>

extern "C"
{
#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"
}

namespace height_map_geodesic
{
    // ========= memory =========
    template <class T> // quickly allocates multiple elements of a given type; no deallocation
    class HeightMapSimlpeMemoryAllocator
    {
    public:
        typedef T *pointer;

        HeightMapSimlpeMemoryAllocator(unsigned block_size = 0,
                                       unsigned max_number_of_blocks = 0)
        {
            reset(block_size,
                  max_number_of_blocks);
        };

        ~HeightMapSimlpeMemoryAllocator() {};

        void reset(unsigned block_size,
                   unsigned max_number_of_blocks)
        {
            m_block_size = block_size;
            m_max_number_of_blocks = max_number_of_blocks;

            m_current_position = 0;

            m_storage.reserve(max_number_of_blocks);
            m_storage.resize(1);
            m_storage[0].resize(block_size);
        };

        pointer allocate(unsigned const n) // allocate n units
        {
            assert(n < m_block_size);

            if (m_current_position + n >= m_block_size)
            {
                m_storage.push_back(std::vector<T>());
                m_storage.back().resize(m_block_size);
                m_current_position = 0;
            }
            pointer result = &m_storage.back()[m_current_position];
            m_current_position += n;

            return result;
        };

    private:
        std::vector<std::vector<T>> m_storage;
        unsigned m_block_size;           // size of a single block
        unsigned m_max_number_of_blocks; // maximum allowed number of blocks
        unsigned m_current_position;     // first unused element inside the current block
    };

    template <class T> // quickly allocates and deallocates single elements of a given type
    class HeightMapMemoryAllocator
    {
    public:
        typedef T *pointer;

        HeightMapMemoryAllocator(unsigned block_size = 1024,
                                 unsigned max_number_of_blocks = 1024)
        {
            reset(block_size,
                  max_number_of_blocks);
        };

        ~HeightMapMemoryAllocator() {};

        void clear()
        {
            reset(m_block_size,
                  m_max_number_of_blocks);
        }

        void reset(unsigned block_size,
                   unsigned max_number_of_blocks)
        {
            m_block_size = block_size;
            m_max_number_of_blocks = max_number_of_blocks;

            assert(m_block_size > 0);
            assert(m_max_number_of_blocks > 0);

            m_current_position = 0;

            m_storage.reserve(max_number_of_blocks);
            m_storage.resize(1);
            m_storage[0].resize(block_size);

            m_deleted.clear();
            m_deleted.reserve(2 * block_size);
        };

        pointer allocate() // allocates single unit of memory
        {
            pointer result;
            if (m_deleted.empty())
            {
                if (m_current_position + 1 >= m_block_size)
                {
                    m_storage.push_back(std::vector<T>());
                    m_storage.back().resize(m_block_size);
                    m_current_position = 0;
                }
                result = &m_storage.back()[m_current_position];
                ++m_current_position;
            }
            else
            {
                result = m_deleted.back();
                m_deleted.pop_back();
            }

            return result;
        };

        void deallocate(pointer p) // allocate n units
        {
            if (m_deleted.size() < m_deleted.capacity())
            {
                m_deleted.push_back(p);
            }
        };

    private:
        std::vector<std::vector<T>> m_storage;
        unsigned m_block_size;           // size of a single block
        unsigned m_max_number_of_blocks; // maximum allowed number of blocks
        unsigned m_current_position;     // first unused element inside the current block

        std::vector<pointer> m_deleted; // pointers to deleted elemets
    };

    class HeightMapOutputBuffer
    {
    public:
        HeightMapOutputBuffer() : m_num_bytes(0)
        {
        }

        void clear()
        {
            m_num_bytes = 0;
            m_buffer = std::unique_ptr<double>();
        }

        template <class T>
        T *allocate(unsigned n)
        {
            double wanted = n * sizeof(T);
            if (wanted > m_num_bytes)
            {
                unsigned new_size = (unsigned)ceil(wanted / (double)sizeof(double));
                m_buffer = std::unique_ptr<double>(new double[new_size]);
                m_num_bytes = new_size * sizeof(double);
            }

            return (T *)m_buffer.get();
        }

        template <class T>
        T *get()
        {
            return (T *)m_buffer.get();
        }

        template <class T>
        unsigned capacity()
        {
            return (unsigned)floor((double)m_num_bytes / (double)sizeof(T));
        };

    private:
        std::unique_ptr<double> m_buffer;
        unsigned m_num_bytes;
    };
    // ========= memory =========

    // ========= element =========
    class HM_Point;
    class HeightMap;
    class HeightMapElementBase;

    typedef HM_Point *hm_point_pointer;
    typedef HeightMap *height_map_pointer;
    typedef HeightMapElementBase *base_pointer;

    template <class Data>       // simple vector that stores info about height map references
    class HeightMapSimpleVector // for efficiency, it uses an outside memory allocator
    {
    public:
        HeightMapSimpleVector() : m_size(0),
                                  m_begin(NULL) {};

        typedef Data *iterator;

        unsigned size() { return m_size; };
        iterator begin() { return m_begin; };
        iterator end() { return m_begin + m_size; };

        template <class DataPointer>
        void set_allocation(DataPointer begin, unsigned size)
        {
            assert(begin != NULL || size == 0);
            m_size = size;
            m_begin = (iterator)begin;
        }

        Data &operator[](unsigned i)
        {
            assert(i < m_size);
            return *(m_begin + i);
        }

        void clear()
        {
            m_size = 0;
            m_begin = NULL;
        }

    private:
        unsigned m_size;
        Data *m_begin;
    };

    enum PointType
    {
        HM_POINT,
        UNDEFINED_POINT
    };

    class HeightMapElementBase // prototype of hm points
    {
    public:
        typedef std::unordered_map<int, double> double_map;

        HeightMapElementBase() : m_id(0),
                                 m_type(UNDEFINED_POINT) {};

        double_map &adj_hm_points_dist() { return m_adjacent_hm_points_distance; };

        unsigned &id() { return m_id; };
        PointType type() { return m_type; };

    protected:
        double_map m_adjacent_hm_points_distance; // list of the adjacent hm points distance

        unsigned m_id;    // unique id
        PointType m_type; // hm point, edge or face
    };

    class HeightMapPoint3D // point in 3D and corresponding operations
    {
    public:
        HeightMapPoint3D() {};
        HeightMapPoint3D(HeightMapPoint3D *p)
        {
            x() = p->x();
            y() = p->y();
            z() = p->z();
        };

        double *xyz() { return m_coordinates; };
        double &x() { return *m_coordinates; };
        double &y() { return *(m_coordinates + 1); };
        double &z() { return *(m_coordinates + 2); };

        double getx() { return *m_coordinates; };
        double gety() { return *(m_coordinates + 1); };
        double getz() { return *(m_coordinates + 2); };

        void set(double new_x, double new_y, double new_z)
        {
            x() = new_x;
            y() = new_y;
            z() = new_z;
        }

        void set(double *data)
        {
            x() = *data;
            y() = *(data + 1);
            z() = *(data + 2);
        }

        double distance(double *v)
        {
            double dx = m_coordinates[0] - v[0];
            double dy = m_coordinates[1] - v[1];
            double dz = m_coordinates[2] - v[2];

            return sqrt(dx * dx + dy * dy + dz * dz);
        };

        double distance(HeightMapPoint3D *v)
        {
            return distance(v->xyz());
        };

        void add(HeightMapPoint3D *v)
        {
            x() += v->x();
            y() += v->y();
            z() += v->z();
        };

        void multiply(double v)
        {
            x() *= v;
            y() *= v;
            z() *= v;
        };

        double m_coordinates[3]; // xyz
    };

    class HM_Point : public HeightMapElementBase, public HeightMapPoint3D
    {
    public:
        HM_Point()
        {
            m_type = HM_POINT;
        };

        ~HM_Point() {};

        bool &boundary() { return m_boundary; };
        bool &deleted() { return m_deleted; };

    private:
        bool m_boundary;
        bool m_deleted;
    };

    class PathPoint : public HeightMapPoint3D // point on the surface of the height map
    {
    public:
        PathPoint() : m_p(NULL) {};

        PathPoint(hm_point_pointer v) : // set the path point in the hm point
                                        PathPoint::HeightMapPoint3D(v),
                                        m_p(v) {};

        PathPoint(base_pointer g,
                  double x,
                  double y,
                  double z,
                  PointType t = UNDEFINED_POINT) : m_p(g)
        {
            set(x, y, z);
        };

        void initialize(PathPoint const &p)
        {
            *this = p;
        }

        ~PathPoint() {};

        PointType type() { return m_p ? m_p->type() : UNDEFINED_POINT; };
        base_pointer &base_element() { return m_p; };

    protected:
        base_pointer m_p; // could be face, hm point or edge pointer
    };
    // ========= element =========

    // ========= simple functions =========
    double const INFIN = 1e100;

    // https://cplusplus.com/forum/beginner/267364/
    // Loads as RGBA... even if file is only RGB
    // Feel free to adjust this if you so please, by changing the 4 to a 0.
    bool load_image(std::vector<unsigned char> &image, const std::string &filename, int &x, int &y)
    {
        int n;
        unsigned char *data = stbi_load(filename.c_str(), &x, &y, &n, 4);
        if (data != nullptr)
        {
            image = std::vector<unsigned char>(data, data + x * y * 4);
        }
        stbi_image_free(data);
        return (data != nullptr);
    }

    void read_org_height_map_from_file(std::string filename, std::vector<int> &points,
                                       int width, int height)
    {
        std::vector<unsigned char> image;
        bool success = load_image(image, filename, width, height);
        if (!success)
        {
            std::cout << "Error loading image\n";
            assert(false);
        }

        // std::cout << "Image width = " << width << '\n';
        // std::cout << "Image height = " << height << '\n';

        points.clear();

        const size_t RGBA = 4;

        for (int j = 0; j < height; j++)
        {
            for (int i = 0; i < width; i++)
            {
                size_t index = RGBA * (j * width + i);
                points.push_back(2 * i);
                points.push_back(2 * j);
                points.push_back(static_cast<int>(0.3 * 2 * height * image[index + 0] / 255));
            }
        }
    }

    // ========= simple functions =========

    // ========= height map =========
    class HeightMap
    {
    public:
        HeightMap() {};

        ~HeightMap() {};

        template <class Points>
        void initialize_height_map_data(Points &p); // build height map

        void height_map_to_terrain(double &memory_usage); // convert height map to terrain

        void copy_height_map(HeightMap *height_map);

        void same_height_map(HeightMap *height_map);

        void update_height_map_merge_four_point(double added_center_point_x, double added_center_point_y,
                                                double added_center_point_z, int added_center_point_index,
                                                std::unordered_map<int, std::unordered_map<int, int>> dominate_table_map,
                                                std::unordered_map<int, int> del_p_dom_by_map,
                                                std::unordered_map<int, bool> &prev_deleted,
                                                std::unordered_map<int, std::unordered_map<int, double>> &prev_adjacent_hm_points_and_distance);

        void update_height_map_exp_four_direct(double added_point_updated_x, double added_point_updated_y,
                                               double added_point_updated_z, int added_center_point_index,
                                               std::unordered_map<int, int> del_p_dom_by_map,
                                               std::vector<int> &exp_four_direct_v_index_list,
                                               double &added_point_prev_x, double &added_point_prev_y, double &added_point_prev_z,
                                               std::unordered_map<int, bool> &prev_deleted,
                                               std::unordered_map<int, std::unordered_map<int, double>> &prev_adjacent_hm_points_and_distance);

        void update_height_map_exp_three_two_one_direct(double new_merged_bottom_left_x, double new_merged_bottom_left_y,
                                                        double new_merged_top_right_x, double new_merged_top_right_y,
                                                        double added_point_updated_x, double added_point_updated_y,
                                                        double added_point_updated_z, int added_center_point_index,
                                                        std::unordered_map<int, int> del_p_dom_by_map,
                                                        std::vector<int> &exp_three_two_one_direct_v_index_list,
                                                        double &added_point_prev_x, double &added_point_prev_y, double &added_point_prev_z,
                                                        std::unordered_map<int, bool> &prev_deleted,
                                                        std::unordered_map<int, std::unordered_map<int, double>> &prev_adjacent_hm_points_and_distance);

        void restore_height_map_merge_four_point(
            std::unordered_map<int, bool> &prev_deleted,
            std::unordered_map<int, std::unordered_map<int, double>> &prev_adjacent_hm_points_and_distance);

        void restore_height_map_exp_four_three_two_one_direct(
            int added_center_point_index, double added_point_prev_x,
            double added_point_prev_y, double added_point_prev_z,
            std::unordered_map<int, bool> &prev_deleted,
            std::unordered_map<int, std::unordered_map<int, double>> &prev_adjacent_hm_points_and_distance);

        std::vector<HM_Point> &hm_points() { return m_hm_points; };

        unsigned closest_hm_points(PathPoint *p,
                                   std::vector<hm_point_pointer> *storage = NULL); // list hm points closest to the point
        double m_xmin, m_xlength, m_xincrement, m_ymin, m_ylength, m_yincrement;
        int m_xpointnum, m_ypointnum;

    private:
        typedef void *void_pointer;
        void_pointer allocate_pointers(unsigned n)
        {
            return m_pointer_allocator.allocate(n);
        }

        std::vector<HM_Point> m_hm_points;

        HeightMapSimlpeMemoryAllocator<void_pointer> m_pointer_allocator; // fast memory allocating for Face/HM_Point/Edge cross-references
    };

    inline unsigned HeightMap::closest_hm_points(PathPoint *p,
                                                 std::vector<hm_point_pointer> *storage)
    {
        if (p->type() == HM_POINT)
        {
            if (storage)
            {
                storage->push_back(static_cast<hm_point_pointer>(p->base_element()));
            }
            return 1;
        }

        assert(0);
        return 0;
    }

    template <class Points>
    void HeightMap::initialize_height_map_data(Points &p)
    {
        assert(p.size() % 3 == 0);
        unsigned const num_hm_points = p.size() / 3;

        unsigned const approximate_number_of_internal_pointers = (num_hm_points) * 3;
        unsigned const max_number_of_pointer_blocks = 100;
        m_pointer_allocator.reset(approximate_number_of_internal_pointers,
                                  max_number_of_pointer_blocks);

        m_hm_points.clear();
        m_hm_points.resize(num_hm_points);

        for (unsigned i = 0; i < num_hm_points; ++i) // copy coordinates to hm points
        {
            HM_Point &v = m_hm_points[i];
            v.id() = i;

            unsigned shift = 3 * i;
            v.x() = p[shift];
            v.y() = p[shift + 1];
            v.z() = p[shift + 2];
            v.boundary() = false;
            v.deleted() = false;
        }

        double x_min = 1e100;
        double x_max = -1e100;
        double y_min = 1e100;
        double y_max = -1e100;
        double z_min = 1e100;
        double z_max = -1e100;
        for (unsigned i = 0; i < m_hm_points.size(); ++i)
        {
            HM_Point &v = m_hm_points[i];
            x_min = std::min(x_min, v.x());
            x_max = std::max(x_max, v.x());
            y_min = std::min(y_min, v.y());
            y_max = std::max(y_max, v.y());
            z_min = std::min(z_min, v.z());
            z_max = std::max(z_max, v.z());
        }

        m_xmin = x_min;
        m_xlength = x_max - x_min;
        m_ymin = y_min;
        m_ylength = y_max - y_min;

        double x_length = x_max - x_min;
        double y_length = y_max - y_min;

        int x_point_num = 0;
        int y_point_num = 0;

        for (unsigned i = 0; i < m_hm_points.size(); ++i)
        {
            HM_Point &v = m_hm_points[i];
            if (v.y() == y_min)
            {
                x_point_num++;
            }
            else
            {
                break;
            }
        }
        y_point_num = m_hm_points.size() / x_point_num;

        m_xpointnum = x_point_num;
        m_ypointnum = y_point_num;

        m_xincrement = m_xlength / (x_point_num - 1);
        m_yincrement = m_ylength / (y_point_num - 1);

        for (int j = 0; j < y_point_num; ++j)
        {
            for (int i = 0; i < x_point_num; ++i)
            {
                int center_index = i + j * x_point_num;

                int bottom_left_index = (i - 1) + (j - 1) * x_point_num;
                int bottom_index = i + (j - 1) * x_point_num;
                int bottom_right_index = (i + 1) + (j - 1) * x_point_num;

                int left_index = (i - 1) + j * x_point_num;
                int right_index = (i + 1) + j * x_point_num;

                int top_left_index = (i - 1) + (j + 1) * x_point_num;
                int top_index = i + (j + 1) * x_point_num;
                int top_right_index = (i + 1) + (j + 1) * x_point_num;

                if (i == 0 && j == 0)
                {
                    m_hm_points[center_index].boundary() = true;
                    m_hm_points[center_index].adj_hm_points_dist()[right_index] = m_hm_points[center_index].distance(&m_hm_points[right_index]);
                    m_hm_points[center_index].adj_hm_points_dist()[top_index] = m_hm_points[center_index].distance(&m_hm_points[top_index]);
                    m_hm_points[center_index].adj_hm_points_dist()[top_right_index] = m_hm_points[center_index].distance(&m_hm_points[top_right_index]);
                }
                else if (i == x_point_num - 1 && j == 0)
                {
                    m_hm_points[center_index].boundary() = true;
                    m_hm_points[center_index].adj_hm_points_dist()[left_index] = m_hm_points[center_index].distance(&m_hm_points[left_index]);
                    m_hm_points[center_index].adj_hm_points_dist()[top_left_index] = m_hm_points[center_index].distance(&m_hm_points[top_left_index]);
                    m_hm_points[center_index].adj_hm_points_dist()[top_index] = m_hm_points[center_index].distance(&m_hm_points[top_index]);
                }
                else if (i == 0 && j == y_point_num - 1)
                {
                    m_hm_points[center_index].boundary() = true;
                    m_hm_points[center_index].adj_hm_points_dist()[bottom_index] = m_hm_points[center_index].distance(&m_hm_points[bottom_index]);
                    m_hm_points[center_index].adj_hm_points_dist()[bottom_right_index] = m_hm_points[center_index].distance(&m_hm_points[bottom_right_index]);
                    m_hm_points[center_index].adj_hm_points_dist()[right_index] = m_hm_points[center_index].distance(&m_hm_points[right_index]);
                }
                else if (i == x_point_num - 1 && j == y_point_num - 1)
                {
                    m_hm_points[center_index].boundary() = true;
                    m_hm_points[center_index].adj_hm_points_dist()[bottom_left_index] = m_hm_points[center_index].distance(&m_hm_points[bottom_left_index]);
                    m_hm_points[center_index].adj_hm_points_dist()[bottom_index] = m_hm_points[center_index].distance(&m_hm_points[bottom_index]);
                    m_hm_points[center_index].adj_hm_points_dist()[left_index] = m_hm_points[center_index].distance(&m_hm_points[left_index]);
                }
                else if (i != 0 && i != x_point_num - 1 && j == 0)
                {
                    m_hm_points[center_index].boundary() = true;
                    m_hm_points[center_index].adj_hm_points_dist()[left_index] = m_hm_points[center_index].distance(&m_hm_points[left_index]);
                    m_hm_points[center_index].adj_hm_points_dist()[right_index] = m_hm_points[center_index].distance(&m_hm_points[right_index]);
                    m_hm_points[center_index].adj_hm_points_dist()[top_left_index] = m_hm_points[center_index].distance(&m_hm_points[top_left_index]);
                    m_hm_points[center_index].adj_hm_points_dist()[top_index] = m_hm_points[center_index].distance(&m_hm_points[top_index]);
                    m_hm_points[center_index].adj_hm_points_dist()[top_right_index] = m_hm_points[center_index].distance(&m_hm_points[top_right_index]);
                }
                else if (i == 0 && j != 0 && j != y_point_num - 1)
                {
                    m_hm_points[center_index].boundary() = true;
                    m_hm_points[center_index].adj_hm_points_dist()[bottom_index] = m_hm_points[center_index].distance(&m_hm_points[bottom_index]);
                    m_hm_points[center_index].adj_hm_points_dist()[bottom_right_index] = m_hm_points[center_index].distance(&m_hm_points[bottom_right_index]);
                    m_hm_points[center_index].adj_hm_points_dist()[right_index] = m_hm_points[center_index].distance(&m_hm_points[right_index]);
                    m_hm_points[center_index].adj_hm_points_dist()[top_index] = m_hm_points[center_index].distance(&m_hm_points[top_index]);
                    m_hm_points[center_index].adj_hm_points_dist()[top_right_index] = m_hm_points[center_index].distance(&m_hm_points[top_right_index]);
                }
                else if (i == x_point_num - 1 && j != 0 && j != y_point_num - 1)
                {
                    m_hm_points[center_index].boundary() = true;
                    m_hm_points[center_index].adj_hm_points_dist()[bottom_left_index] = m_hm_points[center_index].distance(&m_hm_points[bottom_left_index]);
                    m_hm_points[center_index].adj_hm_points_dist()[bottom_index] = m_hm_points[center_index].distance(&m_hm_points[bottom_index]);
                    m_hm_points[center_index].adj_hm_points_dist()[left_index] = m_hm_points[center_index].distance(&m_hm_points[left_index]);
                    m_hm_points[center_index].adj_hm_points_dist()[top_left_index] = m_hm_points[center_index].distance(&m_hm_points[top_left_index]);
                    m_hm_points[center_index].adj_hm_points_dist()[top_index] = m_hm_points[center_index].distance(&m_hm_points[top_index]);
                }
                else if (i != 0 && i != x_point_num - 1 && j == y_point_num - 1)
                {
                    m_hm_points[center_index].boundary() = true;
                    m_hm_points[center_index].adj_hm_points_dist()[bottom_left_index] = m_hm_points[center_index].distance(&m_hm_points[bottom_left_index]);
                    m_hm_points[center_index].adj_hm_points_dist()[bottom_index] = m_hm_points[center_index].distance(&m_hm_points[bottom_index]);
                    m_hm_points[center_index].adj_hm_points_dist()[bottom_right_index] = m_hm_points[center_index].distance(&m_hm_points[bottom_right_index]);
                    m_hm_points[center_index].adj_hm_points_dist()[left_index] = m_hm_points[center_index].distance(&m_hm_points[left_index]);
                    m_hm_points[center_index].adj_hm_points_dist()[right_index] = m_hm_points[center_index].distance(&m_hm_points[right_index]);
                }
                else if (i != 0 && i != x_point_num - 1 && j != 0 && j != y_point_num - 1)
                {
                    m_hm_points[center_index].adj_hm_points_dist()[bottom_left_index] = m_hm_points[center_index].distance(&m_hm_points[bottom_left_index]);
                    m_hm_points[center_index].adj_hm_points_dist()[bottom_index] = m_hm_points[center_index].distance(&m_hm_points[bottom_index]);
                    m_hm_points[center_index].adj_hm_points_dist()[bottom_right_index] = m_hm_points[center_index].distance(&m_hm_points[bottom_right_index]);
                    m_hm_points[center_index].adj_hm_points_dist()[left_index] = m_hm_points[center_index].distance(&m_hm_points[left_index]);
                    m_hm_points[center_index].adj_hm_points_dist()[right_index] = m_hm_points[center_index].distance(&m_hm_points[right_index]);
                    m_hm_points[center_index].adj_hm_points_dist()[top_left_index] = m_hm_points[center_index].distance(&m_hm_points[top_left_index]);
                    m_hm_points[center_index].adj_hm_points_dist()[top_index] = m_hm_points[center_index].distance(&m_hm_points[top_index]);
                    m_hm_points[center_index].adj_hm_points_dist()[top_right_index] = m_hm_points[center_index].distance(&m_hm_points[top_right_index]);
                }
                else
                {
                    std::cout << "adjacent hm points error!" << std::endl;
                    exit(0);
                }
            }
        }

        // std::cout << std::endl;
        // std::cout << "Height map has " << m_hm_points.size() << " pixels" << std::endl;
        // std::cout << "enclosing XYZ box:"
        //           << " X[" << x_min << "," << x_max << "]"
        //           << " Y[" << y_min << "," << y_max << "]"
        //           << " Z[" << z_min << "," << z_max << "]"
        //           << std::endl;
        // std::cout << std::endl;
    }

    void HeightMap::height_map_to_terrain(double &memory_usage)
    {
        double x_min = 1e100;
        double x_max = -1e100;
        double y_min = 1e100;
        double y_max = -1e100;
        double z_min = 1e100;
        double z_max = -1e100;
        for (unsigned i = 0; i < m_hm_points.size(); ++i)
        {
            HM_Point &v = m_hm_points[i];
            x_min = std::min(x_min, v.x());
            x_max = std::max(x_max, v.x());
            y_min = std::min(y_min, v.y());
            y_max = std::max(y_max, v.y());
            z_min = std::min(z_min, v.z());
            z_max = std::max(z_max, v.z());
        }

        double x_length = x_max - x_min;
        double y_length = y_max - y_min;

        int x_point_num = 0;
        int y_point_num = 0;

        for (unsigned i = 0; i < m_hm_points.size(); ++i)
        {
            HM_Point &v = m_hm_points[i];
            if (v.y() == y_min)
            {
                x_point_num++;
            }
            else
            {
                break;
            }
        }
        y_point_num = m_hm_points.size() / x_point_num;

        std::string terrain_write_path = "temp_terrain.off";
        std::ofstream ofs(&terrain_write_path[0], std::ofstream::trunc);

        ofs << "OFF\n";

        int total_vertex_num = x_point_num * y_point_num;
        int total_face_num = (x_point_num - 1) * (y_point_num - 1) * 2;
        int total_edge_num = (x_point_num - 1) * y_point_num + x_point_num * (y_point_num - 1) + (x_point_num - 1) * (y_point_num - 1);

        ofs << total_vertex_num << " " << total_face_num << " " << total_edge_num << " \n";

        for (unsigned i = 0; i < m_hm_points.size(); ++i)
        {
            HM_Point &v = m_hm_points[i];
            ofs << std::fixed << int(v.x()) << " " << int(v.y()) << " " << int(v.z()) << " \n";
        }

        for (int j = 0; j < y_point_num - 1; ++j)
        {
            for (int i = 0; i < x_point_num - 1; ++i)
            {
                int bottom_left_index = i + j * x_point_num;
                int bottom_right_index = (i + 1) + j * x_point_num;
                int top_left_index = i + (j + 1) * x_point_num;
                int top_right_index = (i + 1) + (j + 1) * x_point_num;

                ofs << "3 " << bottom_left_index << " " << bottom_right_index << " " << top_left_index << " \n";
                ofs << "3 " << bottom_right_index << " " << top_right_index << " " << top_left_index << " \n";
            }
        }

        ofs.close();

        std::ifstream in_file(&terrain_write_path[0], std::ios::binary);
        in_file.seekg(0, std::ios::end);
        memory_usage = in_file.tellg();
    }

    void HeightMap::copy_height_map(HeightMap *height_map)
    {
        unsigned const num_hm_points = height_map->hm_points().size();

        unsigned const approximate_number_of_internal_pointers = (num_hm_points) * 3;
        unsigned const max_number_of_pointer_blocks = 100;
        m_pointer_allocator.reset(approximate_number_of_internal_pointers,
                                  max_number_of_pointer_blocks);

        m_hm_points.clear();
        m_hm_points.resize(num_hm_points);

        for (unsigned i = 0; i < num_hm_points; i++)
        {
            HM_Point &v = m_hm_points[i];

            HM_Point u = height_map->hm_points()[i];
            v.id() = i;
            v.x() = u.x();
            v.y() = u.y();
            v.z() = u.z();
            v.boundary() = u.boundary();
            v.deleted() = u.deleted();
        }

        for (unsigned i = 0; i < num_hm_points; i++)
        {
            HM_Point &v = m_hm_points[i];
            HM_Point u = height_map->hm_points()[i];
            for (auto ite : u.adj_hm_points_dist())
            {
                v.adj_hm_points_dist()[ite.first] = u.adj_hm_points_dist()[ite.first];
            }
        }

        m_xmin = height_map->m_xmin;
        m_xlength = height_map->m_xlength;
        m_xincrement = height_map->m_xincrement;
        m_ymin = height_map->m_ymin;
        m_ylength = height_map->m_ylength;
        m_yincrement = height_map->m_yincrement;
        m_xpointnum = height_map->m_xpointnum;
        m_ypointnum = height_map->m_ypointnum;
    }

    void HeightMap::same_height_map(HeightMap *height_map)
    {
        assert(m_hm_points.size() == height_map->hm_points().size());
        for (unsigned i = 0; i < m_hm_points.size(); i++)
        {
            HM_Point &v = m_hm_points[i];
            HM_Point u = height_map->hm_points()[i];
            assert(v.x() == u.x());
            assert(v.y() == u.y());
            assert(v.z() == u.z());
            assert(v.boundary() == u.boundary());
            assert(v.deleted() == u.deleted());
            for (auto ite : u.adj_hm_points_dist())
            {
                assert(v.adj_hm_points_dist()[ite.first] == u.adj_hm_points_dist()[ite.first]);
            }
        }
    }

    void HeightMap::update_height_map_merge_four_point(
        double added_center_point_x, double added_center_point_y,
        double added_center_point_z, int added_center_point_index,
        std::unordered_map<int, std::unordered_map<int, int>> dominate_table_map,
        std::unordered_map<int, int> del_p_dom_by_map,
        std::unordered_map<int, bool> &prev_deleted,
        std::unordered_map<int, std::unordered_map<int, double>> &prev_adjacent_hm_points_and_distance)
    {
        prev_deleted.clear();
        prev_adjacent_hm_points_and_distance.clear();

        std::unordered_map<int, int> neig_of_deleted_p;
        neig_of_deleted_p.clear();

        // calculate the neighbour of deleted pixel
        for (auto ite : dominate_table_map[added_center_point_index])
        {
            int neig_index = ite.second;
            for (auto ite2 : m_hm_points[neig_index].adj_hm_points_dist())
            {
                int adj_p_of_neig_index = ite2.first;
                if (del_p_dom_by_map.count(adj_p_of_neig_index) == 0)
                {
                    if (neig_of_deleted_p.count(adj_p_of_neig_index) == 0)
                    {
                        neig_of_deleted_p[adj_p_of_neig_index] = adj_p_of_neig_index;
                    }
                }
            }
        }

        // store the adjacent information of the original neighbour of the deleted point
        for (auto ite : neig_of_deleted_p)
        {
            int neig_of_deleted_index = ite.second;
            for (auto ite2 : m_hm_points[neig_of_deleted_index].adj_hm_points_dist())
            {
                prev_adjacent_hm_points_and_distance[neig_of_deleted_index][ite2.first] = ite2.second;
            }
        }

        // store the information of the four deleted point
        for (auto ite : dominate_table_map[added_center_point_index])
        {
            int deleted_index = ite.second;
            prev_deleted[deleted_index] = m_hm_points[deleted_index].deleted();
            for (auto ite2 : m_hm_points[deleted_index].adj_hm_points_dist())
            {
                prev_adjacent_hm_points_and_distance[deleted_index][ite2.first] = ite2.second;
            }
        }

        HM_Point h;
        m_hm_points.push_back(h);
        HM_Point &v = m_hm_points[m_hm_points.size() - 1];
        v.id() = (unsigned)m_hm_points.size() - 1;
        assert((int)v.id() == added_center_point_index);
        v.x() = added_center_point_x;
        v.y() = added_center_point_y;
        v.z() = added_center_point_z;
        v.boundary() = false;
        v.deleted() = false;

        // update the added point
        for (auto ite : neig_of_deleted_p)
        {
            assert(m_hm_points[ite.first].id() == ite.first);
            v.adj_hm_points_dist()[ite.first] = v.distance(&m_hm_points[ite.first]);
        }

        // update for the original neighbour of the deleted point
        for (auto ite : neig_of_deleted_p)
        {
            int neig_of_deleted_index = ite.second;
            for (auto ite2 : dominate_table_map[added_center_point_index])
            {
                int neig_index = ite2.second;
                if (m_hm_points[neig_of_deleted_index].adj_hm_points_dist().count(neig_index) != 0)
                {
                    assert(m_hm_points[neig_of_deleted_index].adj_hm_points_dist().count(neig_index) != 0);
                    // remove the deleted point from the adjacent point of these neighbours
                    m_hm_points[neig_of_deleted_index].adj_hm_points_dist().erase(neig_index);
                }
            }
            // add the newly addded point as adjacent point of these neighbours
            assert(m_hm_points[v.id()].id() == v.id());
            m_hm_points[neig_of_deleted_index].adj_hm_points_dist()[v.id()] = m_hm_points[neig_of_deleted_index].distance(&m_hm_points[v.id()]);
        }

        // update the four deleted point
        for (auto ite : dominate_table_map[added_center_point_index])
        {
            int deleted_index = ite.second;
            m_hm_points[deleted_index].deleted() = true;
            m_hm_points[deleted_index].adj_hm_points_dist().clear();
        }
    }

    void HeightMap::update_height_map_exp_four_direct(
        double added_point_updated_x, double added_point_updated_y,
        double added_point_updated_z, int added_center_point_index,
        std::unordered_map<int, int> del_p_dom_by_map,
        std::vector<int> &exp_four_direct_v_index_list,
        double &added_point_prev_x, double &added_point_prev_y, double &added_point_prev_z,
        std::unordered_map<int, bool> &prev_deleted,
        std::unordered_map<int, std::unordered_map<int, double>> &prev_adjacent_hm_points_and_distance)
    {
        HM_Point &v = m_hm_points[added_center_point_index];
        added_point_prev_x = v.x();
        added_point_prev_y = v.y();
        added_point_prev_z = v.z();

        prev_deleted.clear();
        prev_adjacent_hm_points_and_distance.clear();

        std::unordered_map<int, int> neig_of_deleted_p;
        neig_of_deleted_p.clear();

        // calculate the neighbour of the previous added center point
        for (int i = 0; i < exp_four_direct_v_index_list.size(); i++)
        {
            int neig_index = exp_four_direct_v_index_list[i];
            for (auto ite : m_hm_points[neig_index].adj_hm_points_dist())
            {
                int adj_p_of_neig_index = ite.first;
                if (adj_p_of_neig_index == added_center_point_index)
                {
                    continue;
                }
                if (del_p_dom_by_map.count(adj_p_of_neig_index) == 0)
                {
                    if (neig_of_deleted_p.count(adj_p_of_neig_index) == 0)
                    {
                        neig_of_deleted_p[adj_p_of_neig_index] = adj_p_of_neig_index;
                    }
                }
            }
        }

        // store the adjacent information of the previous added center point
        for (auto ite : m_hm_points[added_center_point_index].adj_hm_points_dist())
        {
            prev_adjacent_hm_points_and_distance[added_center_point_index][ite.first] = ite.second;
        }

        // store the adjacent information of the original neighbour of the deleted point
        for (auto ite : neig_of_deleted_p)
        {
            int neig_of_deleted_index = ite.second;
            for (auto ite2 : m_hm_points[neig_of_deleted_index].adj_hm_points_dist())
            {
                prev_adjacent_hm_points_and_distance[neig_of_deleted_index][ite2.first] = ite2.second;
            }
        }

        // store the information of the deleted point
        for (int i = 0; i < exp_four_direct_v_index_list.size(); i++)
        {
            int deleted_index = exp_four_direct_v_index_list[i];
            prev_deleted[deleted_index] = m_hm_points[deleted_index].deleted();
            for (auto ite2 : m_hm_points[deleted_index].adj_hm_points_dist())
            {
                prev_adjacent_hm_points_and_distance[deleted_index][ite2.first] = ite2.second;
            }
        }

        v.x() = added_point_updated_x;
        v.y() = added_point_updated_y;
        v.z() = added_point_updated_z;

        // update the previous added center point
        v.adj_hm_points_dist().clear();
        for (auto ite : neig_of_deleted_p)
        {
            assert(m_hm_points[ite.first].id() == ite.first);
            v.adj_hm_points_dist()[ite.first] = v.distance(&m_hm_points[ite.first]);
        }

        // update for the original neighbour of the deleted point
        for (auto ite : neig_of_deleted_p)
        {
            int neig_of_deleted_index = ite.second;
            for (int i = 0; i < exp_four_direct_v_index_list.size(); i++)
            {
                int neig_index = exp_four_direct_v_index_list[i];
                if (m_hm_points[neig_of_deleted_index].adj_hm_points_dist().count(neig_index) != 0)
                {
                    assert(m_hm_points[neig_of_deleted_index].adj_hm_points_dist().count(neig_index) != 0);
                    // remove the deleted point from the adjacent point of these neighbours
                    m_hm_points[neig_of_deleted_index].adj_hm_points_dist().erase(neig_index);
                }
            }
            // add the previous addded center point as adjacent point of these neighbours
            assert(m_hm_points[ite.first].id() == ite.first);
            m_hm_points[neig_of_deleted_index].adj_hm_points_dist()[v.id()] = m_hm_points[neig_of_deleted_index].distance(&m_hm_points[v.id()]);
        }

        // update the deleted point
        for (int i = 0; i < exp_four_direct_v_index_list.size(); i++)
        {
            int deleted_index = exp_four_direct_v_index_list[i];
            m_hm_points[deleted_index].deleted() = true;
            m_hm_points[deleted_index].adj_hm_points_dist().clear();
        }
    }

    void HeightMap::update_height_map_exp_three_two_one_direct(
        double new_merged_bottom_left_x, double new_merged_bottom_left_y,
        double new_merged_top_right_x, double new_merged_top_right_y,
        double added_point_updated_x, double added_point_updated_y,
        double added_point_updated_z, int added_center_point_index,
        std::unordered_map<int, int> del_p_dom_by_map,
        std::vector<int> &exp_three_two_one_direct_v_index_list,
        double &added_point_prev_x, double &added_point_prev_y, double &added_point_prev_z,
        std::unordered_map<int, bool> &prev_deleted,
        std::unordered_map<int, std::unordered_map<int, double>> &prev_adjacent_hm_points_and_distance)
    {
        HM_Point &v = m_hm_points[added_center_point_index];
        added_point_prev_x = v.x();
        added_point_prev_y = v.y();
        added_point_prev_z = v.z();

        prev_deleted.clear();
        prev_adjacent_hm_points_and_distance.clear();

        std::unordered_map<int, int> neig_of_deleted_p_and_not_adj_center;
        neig_of_deleted_p_and_not_adj_center.clear();
        std::unordered_map<int, int> neig_of_deleted_p_and_adj_center;
        neig_of_deleted_p_and_adj_center.clear();

        // calculate the original neighbour of the deleted point which are both adjacent and not adjacent to the previous centere point
        for (int i = 0; i < exp_three_two_one_direct_v_index_list.size(); i++)
        {
            int neig_index = exp_three_two_one_direct_v_index_list[i];
            for (auto ite : m_hm_points[neig_index].adj_hm_points_dist())
            {
                int adj_p_of_neig_index = ite.first;
                if (adj_p_of_neig_index == added_center_point_index)
                {
                    continue;
                }
                if (del_p_dom_by_map.count(adj_p_of_neig_index) == 0)
                {
                    if (v.adj_hm_points_dist().count(adj_p_of_neig_index) == 0)
                    {
                        if (neig_of_deleted_p_and_not_adj_center.count(adj_p_of_neig_index) == 0)
                        {
                            neig_of_deleted_p_and_not_adj_center[adj_p_of_neig_index] = adj_p_of_neig_index;
                        }
                    }
                    else
                    {
                        if (neig_of_deleted_p_and_adj_center.count(adj_p_of_neig_index) == 0)
                        {
                            neig_of_deleted_p_and_adj_center[adj_p_of_neig_index] = adj_p_of_neig_index;
                        }
                    }
                }
            }
        }

        // store the adjacent information of the previous added center point
        for (auto ite : m_hm_points[added_center_point_index].adj_hm_points_dist())
        {
            prev_adjacent_hm_points_and_distance[added_center_point_index][ite.first] = ite.second;
        }

        // store the adjacent information of the original neighbour of the deleted point which are not adjacent to the previous centere point
        for (auto ite : neig_of_deleted_p_and_not_adj_center)
        {
            int neig_of_deleted_index = ite.second;
            for (auto ite2 : m_hm_points[neig_of_deleted_index].adj_hm_points_dist())
            {
                prev_adjacent_hm_points_and_distance[neig_of_deleted_index][ite2.first] = ite2.second;
            }
        }

        // store the adjacent information of the original neighbour of the deleted point which are adjacent to the previous centere point
        for (auto ite : neig_of_deleted_p_and_adj_center)
        {
            int neig_of_deleted_index = ite.second;
            for (auto ite2 : m_hm_points[neig_of_deleted_index].adj_hm_points_dist())
            {
                prev_adjacent_hm_points_and_distance[neig_of_deleted_index][ite2.first] = ite2.second;
            }
        }

        // store the information of the deleted point
        for (int i = 0; i < exp_three_two_one_direct_v_index_list.size(); i++)
        {
            int deleted_index = exp_three_two_one_direct_v_index_list[i];
            prev_deleted[deleted_index] = m_hm_points[deleted_index].deleted();
            for (auto ite2 : m_hm_points[deleted_index].adj_hm_points_dist())
            {
                prev_adjacent_hm_points_and_distance[deleted_index][ite2.first] = ite2.second;
            }
        }

        // update the previous added center point
        std::vector<int> not_del_adj_of_v_index_list;
        not_del_adj_of_v_index_list.clear();

        for (auto ite : v.adj_hm_points_dist())
        {
            int adj_p_of_added_p_index = ite.first;
            double adj_p_of_added_p_x = m_hm_points[adj_p_of_added_p_index].getx();
            double adj_p_of_added_p_y = m_hm_points[adj_p_of_added_p_index].gety();

            if (adj_p_of_added_p_x >= new_merged_bottom_left_x &&
                adj_p_of_added_p_x <= new_merged_top_right_x &&
                adj_p_of_added_p_y >= new_merged_bottom_left_y &&
                adj_p_of_added_p_y <= new_merged_top_right_y)
            {
                continue;
            }
            not_del_adj_of_v_index_list.push_back(adj_p_of_added_p_index);
        }

        v.x() = added_point_updated_x;
        v.y() = added_point_updated_y;
        v.z() = added_point_updated_z;

        v.adj_hm_points_dist().clear();

        for (int i = 0; i < not_del_adj_of_v_index_list.size(); i++)
        {
            assert(m_hm_points[not_del_adj_of_v_index_list[i]].id() == not_del_adj_of_v_index_list[i]);
            v.adj_hm_points_dist()[not_del_adj_of_v_index_list[i]] = v.distance(&m_hm_points[not_del_adj_of_v_index_list[i]]);
        }

        for (auto ite : neig_of_deleted_p_and_not_adj_center)
        {
            assert(m_hm_points[ite.first].id() == ite.first);
            v.adj_hm_points_dist()[ite.first] = v.distance(&m_hm_points[ite.first]);
        }

        // update for the original neighbour of the deleted point which are not adjacent to the previous centere point
        for (auto ite : neig_of_deleted_p_and_not_adj_center)
        {
            int neig_of_deleted_index = ite.second;
            for (int i = 0; i < exp_three_two_one_direct_v_index_list.size(); i++)
            {
                int neig_index = exp_three_two_one_direct_v_index_list[i];
                if (m_hm_points[neig_of_deleted_index].adj_hm_points_dist().count(neig_index) != 0)
                {
                    assert(m_hm_points[neig_of_deleted_index].adj_hm_points_dist().count(neig_index) != 0);
                    // remove the deleted point from the adjacent point of these neighbours
                    m_hm_points[neig_of_deleted_index].adj_hm_points_dist().erase(neig_index);
                }
            }
            // add the previous addded center point as adjacent point of these neighbours
            assert(m_hm_points[v.id()].id() == v.id());
            m_hm_points[neig_of_deleted_index].adj_hm_points_dist()[v.id()] = m_hm_points[neig_of_deleted_index].distance(&m_hm_points[v.id()]);
        }

        // update for the original neighbour of the deleted point which are adjacent to the previous centere point
        for (auto ite : neig_of_deleted_p_and_adj_center)
        {
            int neig_of_deleted_index = ite.second;
            for (int i = 0; i < exp_three_two_one_direct_v_index_list.size(); i++)
            {
                int neig_index = exp_three_two_one_direct_v_index_list[i];
                if (m_hm_points[neig_of_deleted_index].adj_hm_points_dist().count(neig_index) != 0)
                {
                    assert(m_hm_points[neig_of_deleted_index].adj_hm_points_dist().count(neig_index) != 0);
                    // remove the deleted point from the adjacent point of these neighbours
                    m_hm_points[neig_of_deleted_index].adj_hm_points_dist().erase(neig_index);
                }
            }
        }

        // update the deleted point
        for (int i = 0; i < exp_three_two_one_direct_v_index_list.size(); i++)
        {
            int neig_index = exp_three_two_one_direct_v_index_list[i];
            m_hm_points[neig_index].deleted() = true;
            m_hm_points[neig_index].adj_hm_points_dist().clear();
        }
    }

    void HeightMap::restore_height_map_merge_four_point(
        std::unordered_map<int, bool> &prev_deleted,
        std::unordered_map<int, std::unordered_map<int, double>> &prev_adjacent_hm_points_and_distance)
    {
        m_hm_points.pop_back();

        for (auto ite : prev_deleted)
        {
            m_hm_points[ite.first].deleted() = ite.second;
        }
        for (auto ite : prev_adjacent_hm_points_and_distance)
        {
            int neig_of_deleted_index = ite.first;
            m_hm_points[neig_of_deleted_index].adj_hm_points_dist().clear();
            for (auto ite2 : prev_adjacent_hm_points_and_distance[neig_of_deleted_index])
            {
                assert(m_hm_points[ite2.first].id() == ite2.first);
                m_hm_points[neig_of_deleted_index].adj_hm_points_dist()[ite2.first] = ite2.second;
            }
        }
    }

    void HeightMap::restore_height_map_exp_four_three_two_one_direct(
        int added_center_point_index, double added_point_prev_x,
        double added_point_prev_y, double added_point_prev_z,
        std::unordered_map<int, bool> &prev_deleted,
        std::unordered_map<int, std::unordered_map<int, double>> &prev_adjacent_hm_points_and_distance)
    {
        HM_Point &v = m_hm_points[added_center_point_index];
        v.x() = added_point_prev_x;
        v.y() = added_point_prev_y;
        v.z() = added_point_prev_z;

        for (auto ite : prev_deleted)
        {
            m_hm_points[ite.first].deleted() = ite.second;
        }
        for (auto ite : prev_adjacent_hm_points_and_distance)
        {
            int neig_of_deleted_index = ite.first;
            m_hm_points[neig_of_deleted_index].adj_hm_points_dist().clear();
            for (auto ite2 : prev_adjacent_hm_points_and_distance[neig_of_deleted_index])
            {
                assert(m_hm_points[ite2.first].id() == ite2.first);
                m_hm_points[neig_of_deleted_index].adj_hm_points_dist()[ite2.first] = ite2.second;
            }
        }
    }

    // ========= height map =========

    // ========= algorithm base =========
    class HeightMapGeodesicAlgorithmBase
    {
    public:
        enum AlgorithmType
        {
            DIJKSTRA,
            UNDEFINED_ALGORITHM
        };

        HeightMapGeodesicAlgorithmBase(height_map_geodesic::HeightMap *height_map) : m_type(UNDEFINED_ALGORITHM),
                                                                                     m_max_propagation_distance(1e100),
                                                                                     m_height_map(height_map) {};

        virtual ~HeightMapGeodesicAlgorithmBase() {};

        virtual unsigned best_source(PathPoint &point, // after propagation step is done, quickly find what source this point belongs to and what is the distance to this source
                                     double &best_source_distance) = 0;

        AlgorithmType type() { return m_type; };

        virtual std::string name();

        height_map_geodesic::HeightMap *height_map() { return m_height_map; };

    protected:
        void set_stop_conditions(std::vector<PathPoint> *stop_points);

        void set_stop_conditions(double stop_distance);

        double stop_distance()
        {
            return m_max_propagation_distance;
        }

        AlgorithmType m_type; // type of the algorithm

        typedef std::pair<hm_point_pointer, double> stop_hm_point_with_distace_type;
        std::vector<stop_hm_point_with_distace_type> m_stop_hm_points; // algorithm stops propagation after covering certain hm points
        double m_max_propagation_distance;                             // or reaching the certain distance

        height_map_geodesic::HeightMap *m_height_map;

        double m_time_consumed;                // how much time does the propagation step takes
        double m_propagation_distance_stopped; // at what distance (if any) the propagation algorithm stopped
    };

    inline double length(std::vector<PathPoint> &path)
    {
        double length = 0;
        if (!path.empty())
        {
            for (unsigned i = 0; i < path.size() - 1; ++i)
            {
                length += path[i].distance(&path[i + 1]);
            }
        }
        return length;
    }

    inline void print_info_about_path(std::vector<PathPoint> &path)
    {
        std::cout << "number of the points in the path = " << path.size()
                  << ", length of the path = " << length(path)
                  << std::endl;
    }

    inline std::string HeightMapGeodesicAlgorithmBase::name()
    {
        switch (m_type)
        {
        case DIJKSTRA:
            return "dijkstra";
        default:
        case UNDEFINED_ALGORITHM:
            return "undefined";
        }
    }

    inline void HeightMapGeodesicAlgorithmBase::set_stop_conditions(std::vector<PathPoint> *stop_points)
    {
        if (!stop_points)
        {
            m_stop_hm_points.clear();
            return;
        }

        m_stop_hm_points.resize(stop_points->size());

        std::vector<hm_point_pointer> possible_hm_points;
        for (unsigned i = 0; i < stop_points->size(); ++i)
        {
            PathPoint *point = &(*stop_points)[i];

            possible_hm_points.clear();
            m_height_map->closest_hm_points(point, &possible_hm_points);

            hm_point_pointer closest_vertex = NULL;
            double min_distance = 1e100;
            for (unsigned j = 0; j < possible_hm_points.size(); ++j)
            {
                double distance = point->distance(possible_hm_points[j]);
                if (distance < min_distance)
                {
                    min_distance = distance;
                    closest_vertex = possible_hm_points[j];
                }
            }
            assert(closest_vertex);

            m_stop_hm_points[i].first = closest_vertex;
            m_stop_hm_points[i].second = min_distance;
        }
    }

    inline void HeightMapGeodesicAlgorithmBase::set_stop_conditions(double stop_distance)
    {
        m_max_propagation_distance = stop_distance;
    }
    // ========= algorithm base =========

    // ========= algorithm graph base =========
    template <class Node>
    class HeightMapGeodesicAlgorithmGraphBase : public HeightMapGeodesicAlgorithmBase
    {
    public:
        typedef Node *node_pointer;

        HeightMapGeodesicAlgorithmGraphBase(height_map_geodesic::HeightMap *height_map) : HeightMapGeodesicAlgorithmBase(height_map) {};

        ~HeightMapGeodesicAlgorithmGraphBase() {};

        void propagate(std::vector<PathPoint> &sources,
                       double max_propagation_distance, int iter_num);

        void propagate(std::vector<PathPoint> &sources,
                       std::vector<PathPoint> *stop_points);

        void trace_back(PathPoint &destination,
                        std::vector<PathPoint> &path);

        unsigned best_source(PathPoint &point,
                             double &best_source_distance);

        void get_memory(int &total_node_size, size_t &node_size);

    protected:
        unsigned node_index(hm_point_pointer v) // gives index of the node that corresponds to this hm point
        {
            return v->id();
        };

        void set_sources(std::vector<PathPoint> &sources)
        {
            m_sources = sources;
        }

        node_pointer best_first_node(PathPoint &point, double &best_total_distance)
        {
            node_pointer best_node = NULL;
            if (point.type() == HM_POINT)
            {
                hm_point_pointer v = (hm_point_pointer)point.base_element();
                best_node = &m_nodes[node_index(v)];
                best_total_distance = best_node->distance_from_source();
            }

            if (best_total_distance > m_propagation_distance_stopped) // result is unreliable
            {
                best_total_distance = INFIN;
                return NULL;
            }
            else
            {
                return best_node;
            }
        }; // quickly find what node will be the next one in geodesic path

        bool check_stop_conditions_distance();                    // check when propagation should stop
        bool check_stop_conditions_cover_points(unsigned &index); // check when propagation should stop

        virtual void list_nodes_visible_from_source(HeightMapElementBase *p,
                                                    std::vector<node_pointer> &storage) = 0; // list all nodes that belong to this height map element

        virtual void list_nodes_visible_from_node(node_pointer node, // list all nodes that belong to this height map element
                                                  std::vector<node_pointer> &storage,
                                                  std::vector<double> &distances,
                                                  double threshold_distance) = 0; // list only the nodes whose current distance is larger than the threshold

        std::vector<Node> m_nodes; // nodes of the graph

        typedef std::set<node_pointer, Node> queue_type;
        queue_type m_queue;

        std::vector<PathPoint> m_sources; // for simplicity, we keep sources as they are
    };

    template <class Node>
    void HeightMapGeodesicAlgorithmGraphBase<Node>::propagate(std::vector<PathPoint> &sources,
                                                              double max_propagation_distance,
                                                              int iter_num)
    {
        for (int m = 0; m < iter_num; m++)
        {
            set_stop_conditions(max_propagation_distance);
            set_sources(sources);

            m_queue.clear();
            m_propagation_distance_stopped = INFIN;
            for (unsigned i = 0; i < m_nodes.size(); ++i)
            {
                m_nodes[i].clear();
            }

            clock_t start = clock();

            std::vector<node_pointer> visible_nodes; // initialize hm points directly visible from sources
            for (unsigned i = 0; i < m_sources.size(); ++i)
            {
                PathPoint *source = &m_sources[i];
                list_nodes_visible_from_source(source->base_element(),
                                               visible_nodes);

                for (unsigned j = 0; j < visible_nodes.size(); ++j)
                {
                    node_pointer node = visible_nodes[j];
                    double distance = node->distance(source);
                    if (distance < node->distance_from_source())
                    {
                        node->distance_from_source() = distance;
                        node->source_index() = i;
                        node->previous() = NULL;
                    }
                }
                visible_nodes.clear();
            }

            for (unsigned i = 0; i < m_nodes.size(); ++i) // initialize the queue
            {
                if (m_nodes[i].distance_from_source() < INFIN)
                {
                    m_queue.insert(&m_nodes[i]);
                }
            }

            unsigned counter = 0;

            std::vector<double> distances_between_nodes;
            while (!m_queue.empty()) // main cycle
            {
                if (counter++ % 10 == 0) // check if we covered all required hm points
                {
                    if (check_stop_conditions_distance())
                    {
                        std::cout << "break" << std::endl;
                        break;
                    }
                }

                node_pointer min_node = *m_queue.begin();
                m_queue.erase(m_queue.begin());
                assert(min_node->distance_from_source() < INFIN);

                visible_nodes.clear();
                distances_between_nodes.clear();
                list_nodes_visible_from_node(min_node,
                                             visible_nodes,
                                             distances_between_nodes,
                                             min_node->distance_from_source());

                for (unsigned i = 0; i < visible_nodes.size(); ++i) // update all the adgecent hm points
                {
                    node_pointer next_node = visible_nodes[i];

                    if (next_node->distance_from_source() > min_node->distance_from_source() +
                                                                distances_between_nodes[i])
                    {
                        if (next_node->distance_from_source() < INFIN) // remove it from the queue
                        {
                            typename queue_type::iterator iter = m_queue.find(next_node);
                            assert(iter != m_queue.end());
                            m_queue.erase(iter);
                        }
                        next_node->distance_from_source() = min_node->distance_from_source() +
                                                            distances_between_nodes[i];
                        next_node->source_index() = min_node->source_index();
                        next_node->previous() = min_node;
                        m_queue.insert(next_node);
                    }
                }
            }

            m_propagation_distance_stopped = m_queue.empty() ? INFIN : (*m_queue.begin())->distance_from_source();
            clock_t finish = clock();
            m_time_consumed = (static_cast<double>(finish) - static_cast<double>(start)) / CLOCKS_PER_SEC;
        }
    }

    template <class Node>
    void HeightMapGeodesicAlgorithmGraphBase<Node>::propagate(std::vector<PathPoint> &sources,
                                                              std::vector<PathPoint> *stop_points)
    {
        set_stop_conditions(stop_points);
        set_sources(sources);

        m_queue.clear();
        m_propagation_distance_stopped = INFIN;
        for (unsigned i = 0; i < m_nodes.size(); ++i)
        {
            m_nodes[i].clear();
        }

        clock_t start = clock();

        std::vector<node_pointer> visible_nodes; // initialize hm points directly visible from sources
        for (unsigned i = 0; i < m_sources.size(); ++i)
        {
            PathPoint *source = &m_sources[i];
            list_nodes_visible_from_source(source->base_element(),
                                           visible_nodes);

            for (unsigned j = 0; j < visible_nodes.size(); ++j)
            {
                node_pointer node = visible_nodes[j];
                double distance = node->distance(source);
                if (distance < node->distance_from_source())
                {
                    node->distance_from_source() = distance;
                    node->source_index() = i;
                    node->previous() = NULL;
                }
            }
            visible_nodes.clear();
        }

        for (unsigned i = 0; i < m_nodes.size(); ++i) // initialize the queue
        {
            if (m_nodes[i].distance_from_source() < INFIN)
            {
                m_queue.insert(&m_nodes[i]);
            }
        }

        unsigned counter = 0;
        unsigned satisfied_index = 0;

        std::vector<double> distances_between_nodes;
        while (!m_queue.empty()) // main cycle
        {
            if (counter++ % 10 == 0) // check if we covered all required hm points
            {
                if (check_stop_conditions_cover_points(satisfied_index))
                {
                    break;
                }
            }

            node_pointer min_node = *m_queue.begin();
            m_queue.erase(m_queue.begin());
            assert(min_node->distance_from_source() < INFIN);

            visible_nodes.clear();
            distances_between_nodes.clear();
            list_nodes_visible_from_node(min_node,
                                         visible_nodes,
                                         distances_between_nodes,
                                         min_node->distance_from_source());

            for (unsigned i = 0; i < visible_nodes.size(); ++i) // update all the adgecent hm points
            {
                node_pointer next_node = visible_nodes[i];

                if (next_node->distance_from_source() > min_node->distance_from_source() +
                                                            distances_between_nodes[i])
                {
                    if (next_node->distance_from_source() < INFIN) // remove it from the queue
                    {
                        typename queue_type::iterator iter = m_queue.find(next_node);
                        assert(iter != m_queue.end());
                        m_queue.erase(iter);
                    }
                    next_node->distance_from_source() = min_node->distance_from_source() +
                                                        distances_between_nodes[i];
                    next_node->source_index() = min_node->source_index();
                    next_node->previous() = min_node;
                    m_queue.insert(next_node);
                }
            }
        }

        m_propagation_distance_stopped = m_queue.empty() ? INFIN : (*m_queue.begin())->distance_from_source();
        clock_t finish = clock();
        m_time_consumed = (static_cast<double>(finish) - static_cast<double>(start)) / CLOCKS_PER_SEC;
    }

    template <class Node>
    inline bool HeightMapGeodesicAlgorithmGraphBase<Node>::check_stop_conditions_distance()
    {
        double queue_min_distance = (*m_queue.begin())->distance_from_source();
        if (queue_min_distance < m_max_propagation_distance)
        {
            return false;
        }
        return true;
    }

    template <class Node>
    inline bool HeightMapGeodesicAlgorithmGraphBase<Node>::check_stop_conditions_cover_points(unsigned &index)
    {
        double queue_min_distance = (*m_queue.begin())->distance_from_source();
        while (index < m_stop_hm_points.size())
        {
            hm_point_pointer v = m_stop_hm_points[index].first;
            Node &node = m_nodes[node_index(v)];
            if (queue_min_distance < node.distance_from_source() + m_stop_hm_points[index].second)
            {
                return false;
            }
            ++index;
        }
        return true;
    }

    template <class Node>
    inline void HeightMapGeodesicAlgorithmGraphBase<Node>::trace_back(PathPoint &destination, // trace back piecewise-linear path
                                                                      std::vector<PathPoint> &path)
    {
        path.clear();

        double total_path_length;
        node_pointer node = best_first_node(destination, total_path_length);

        if (total_path_length > INFIN / 2.0) // unable to find the path
        {
            return;
        }

        path.push_back(destination);

        if (node->distance(&destination) > 1e-50)
        {
            path.push_back(node->path_point());
        }

        while (node->previous()) // follow the path
        {
            node = node->previous();
            path.push_back(node->path_point());
        }

        PathPoint &source = m_sources[node->source_index()]; // add source to the path if it is not already there
        if (node->distance(&source) > 1e-50)
        {
            path.push_back(source);
        }
    }

    template <class Node>
    inline void HeightMapGeodesicAlgorithmGraphBase<Node>::get_memory(int &total_node_size, size_t &node_size)
    {
        total_node_size = m_nodes.size();
        node_size = sizeof(Node);
    }

    template <class Node>
    inline unsigned HeightMapGeodesicAlgorithmGraphBase<Node>::best_source(PathPoint &point, // quickly find what source this point belongs to and what is the distance to this source
                                                                           double &best_source_distance)
    {
        node_pointer node = best_first_node(point, best_source_distance);
        return node ? node->source_index() : 0;
    };

    // ========= algorithm graph base =========

    // ========= Dijkstra shortest path algorithm =========
    class HeightMapDijkstraNode
    {
        typedef HeightMapDijkstraNode *node_pointer;

    public:
        HeightMapDijkstraNode() {};
        ~HeightMapDijkstraNode() {};

        double &distance_from_source() { return m_distance; };
        node_pointer &previous() { return m_previous; };
        unsigned &source_index() { return m_source_index; };
        hm_point_pointer &hm_point() { return m_hm_point; };

        void clear()
        {
            m_distance = INFIN;
            m_previous = NULL;
        }

        bool operator()(node_pointer const s1, node_pointer const s2) const
        {
            return s1->distance_from_source() != s2->distance_from_source() ? s1->distance_from_source() < s2->distance_from_source() : s1->hm_point()->id() < s2->hm_point()->id();
        };

        double distance(PathPoint *p)
        {
            return m_hm_point->distance(p);
        }

        PathPoint path_point()
        {
            return PathPoint(m_hm_point);
        }

    private:
        double m_distance;           // distance to the closest source
        unsigned m_source_index;     // closest source index
        node_pointer m_previous;     // previous node in the geodesic path
        hm_point_pointer m_hm_point; // correspoding hm point
    };

    class HeightMapGeodesicAlgorithmDijkstra : public HeightMapGeodesicAlgorithmGraphBase<HeightMapDijkstraNode>
    {
    public:
        typedef HeightMapDijkstraNode Node;
        typedef Node *node_pointer;

        HeightMapGeodesicAlgorithmDijkstra(height_map_geodesic::HeightMap *height_map) : HeightMapGeodesicAlgorithmGraphBase<Node>(height_map)
        {
            m_type = DIJKSTRA;

            m_nodes.resize(height_map->hm_points().size());
            for (unsigned i = 0; i < m_nodes.size(); ++i)
            {
                m_nodes[i].hm_point() = &m_height_map->hm_points()[i];
            }
        };

        ~HeightMapGeodesicAlgorithmDijkstra() {};

    protected:
        void list_nodes_visible_from_source(HeightMapElementBase *p,
                                            std::vector<node_pointer> &storage); // list all nodes that belong to this height_map element

        void list_nodes_visible_from_node(node_pointer node, // list all nodes that belong to this height_map element
                                          std::vector<node_pointer> &storage,
                                          std::vector<double> &distances,
                                          double threshold_distance); // list only the nodes whose current distance is larger than the threshold
    };

    void HeightMapGeodesicAlgorithmDijkstra::list_nodes_visible_from_source(HeightMapElementBase *p,
                                                                            std::vector<node_pointer> &storage)
    {
        assert(p->type() != UNDEFINED_POINT);

        if (p->type() == HM_POINT)
        {
            hm_point_pointer v = static_cast<hm_point_pointer>(p);
            storage.push_back(&m_nodes[node_index(v)]);
        }
    }

    inline void HeightMapGeodesicAlgorithmDijkstra::list_nodes_visible_from_node(node_pointer node, // list all nodes that belong to this height map element
                                                                                 std::vector<node_pointer> &storage,
                                                                                 std::vector<double> &distances,
                                                                                 double threshold_distance)
    {

        hm_point_pointer v = node->hm_point();
        assert(storage.size() == distances.size());
        for (auto ite : v->adj_hm_points_dist())
        {
            if (!m_height_map->hm_points()[ite.first].deleted())
            {
                node_pointer new_node = &m_nodes[ite.first];
                if (new_node->distance_from_source() > threshold_distance + v->adj_hm_points_dist()[ite.first])
                {
                    storage.push_back(new_node);
                    distances.push_back(v->adj_hm_points_dist()[ite.first]);
                }
            }
        }
    }
    // ========= Dijkstra shortest path algorithm =========

} // geodesic