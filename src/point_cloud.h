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

namespace point_cloud_geodesic
{
    // ========= memory =========
    template <class T> // quickly allocates multiple elements of a given type; no deallocation
    class PointCloudSimlpeMemoryAllocator
    {
    public:
        typedef T *pointer;

        PointCloudSimlpeMemoryAllocator(unsigned block_size = 0,
                                        unsigned max_number_of_blocks = 0)
        {
            reset(block_size,
                  max_number_of_blocks);
        };

        ~PointCloudSimlpeMemoryAllocator() {};

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
    class PointCloudMemoryAllocator
    {
    public:
        typedef T *pointer;

        PointCloudMemoryAllocator(unsigned block_size = 1024,
                                  unsigned max_number_of_blocks = 1024)
        {
            reset(block_size,
                  max_number_of_blocks);
        };

        ~PointCloudMemoryAllocator() {};

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

    class PointCloudOutputBuffer
    {
    public:
        PointCloudOutputBuffer() : m_num_bytes(0)
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
    class PC_Point;
    class PointCloud;
    class PointCloudElementBase;

    typedef PC_Point *pc_point_pointer;
    typedef PointCloud *point_cloud_pointer;
    typedef PointCloudElementBase *base_pointer;

    template <class Data>        // simple vector that stores info about point cloud references
    class PointCloudSimpleVector // for efficiency, it uses an outside memory allocator
    {
    public:
        PointCloudSimpleVector() : m_size(0),
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
        PC_POINT,
        UNDEFINED_POINT
    };

    class PointCloudElementBase // prototype of pc points, edges and faces
    {
    public:
        typedef PointCloudSimpleVector<pc_point_pointer> pc_point_pointer_vector;
        typedef PointCloudSimpleVector<double> double_vector;

        PointCloudElementBase() : m_id(0),
                                  m_type(UNDEFINED_POINT) {};

        pc_point_pointer_vector &adjacent_pc_points() { return m_adjacent_pc_points; };
        double_vector &adjacent_pc_points_distance() { return m_adjacent_pc_points_distance; };

        unsigned &id() { return m_id; };
        PointType type() { return m_type; };

    protected:
        pc_point_pointer_vector m_adjacent_pc_points; // list of the adjacent pc points
        double_vector m_adjacent_pc_points_distance;  // list of the adjacent pc points distance

        unsigned m_id;    // unique id
        PointType m_type; // pc point, edge or face
    };

    class PointCloudPoint3D // point in 3D and corresponding operations
    {
    public:
        PointCloudPoint3D() {};
        PointCloudPoint3D(PointCloudPoint3D *p)
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

        double distance(PointCloudPoint3D *v)
        {
            return distance(v->xyz());
        };

        void add(PointCloudPoint3D *v)
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

    class PC_Point : public PointCloudElementBase, public PointCloudPoint3D
    {
    public:
        PC_Point()
        {
            m_type = PC_POINT;
        };

        ~PC_Point() {};
    };

    class PathPoint : public PointCloudPoint3D // point on the surface of the point cloud
    {
    public:
        PathPoint() : m_p(NULL) {};

        PathPoint(pc_point_pointer v) : // set the path point in the pc point
                                        PathPoint::PointCloudPoint3D(v),
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
        base_pointer m_p; // could be face, pc point or edge pointer
    };
    // ========= element =========

    // ========= simple functions =========
    double const INFIN = 1e100;

    template <class Points>
    inline bool read_point_cloud_from_file(char *filename,
                                           Points &points)
    {
        std::ifstream file(filename);
        assert(file.is_open());
        if (!file.is_open())
            return false;

        std::string dum;
        unsigned num_points;
        file >> num_points;
        assert(num_points >= 1);

        points.resize(num_points * 3);
        for (typename Points::iterator i = points.begin(); i != points.end(); ++i)
        {
            file >> *i;
        }

        file.close();

        return true;
    }
    // ========= simple functions =========

    // ========= point cloud =========
    class PointCloud
    {
    public:
        PointCloud() {};

        ~PointCloud() {};

        template <class Points>
        void initialize_point_cloud_data(Points &p); // build point cloud

        void point_cloud_to_terrain(double &memory_usage); // convert point cloud to terrain

        std::vector<PC_Point> &pc_points() { return m_pc_points; };

        unsigned closest_pc_points(PathPoint *p,
                                   std::vector<pc_point_pointer> *storage = NULL); // list pc points closest to the point
        double m_xmin, m_xlength, m_xincrement, m_ymin, m_ylength, m_yincrement;
        int m_xpointnum, m_ypointnum;

    private:
        typedef void *void_pointer;
        void_pointer allocate_pointers(unsigned n)
        {
            return m_pointer_allocator.allocate(n);
        }

        std::vector<PC_Point> m_pc_points;

        PointCloudSimlpeMemoryAllocator<void_pointer> m_pointer_allocator; // fast memory allocating for Face/PC_Point/Edge cross-references
    };

    inline unsigned PointCloud::closest_pc_points(PathPoint *p,
                                                  std::vector<pc_point_pointer> *storage)
    {
        //	assert(p->type() != UNDEFINED_POINT);

        if (p->type() == PC_POINT)
        {
            if (storage)
            {
                storage->push_back(static_cast<pc_point_pointer>(p->base_element()));
            }
            return 1;
        }

        assert(0);
        return 0;
    }

    template <class Points>
    void PointCloud::initialize_point_cloud_data(Points &p)
    {
        assert(p.size() % 3 == 0);
        unsigned const num_pc_points = p.size() / 3;

        unsigned const approximate_number_of_internal_pointers = (num_pc_points) * 3;
        unsigned const max_number_of_pointer_blocks = 100;
        m_pointer_allocator.reset(approximate_number_of_internal_pointers,
                                  max_number_of_pointer_blocks);

        m_pc_points.resize(num_pc_points);
        for (unsigned i = 0; i < num_pc_points; ++i) // copy coordinates to pc points
        {
            PC_Point &v = m_pc_points[i];
            v.id() = i;

            unsigned shift = 3 * i;
            v.x() = p[shift];
            v.y() = p[shift + 1];
            v.z() = p[shift + 2];
        }

        double x_min = 1e100;
        double x_max = -1e100;
        double y_min = 1e100;
        double y_max = -1e100;
        double z_min = 1e100;
        double z_max = -1e100;
        for (unsigned i = 0; i < m_pc_points.size(); ++i)
        {
            PC_Point &v = m_pc_points[i];
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

        for (unsigned i = 0; i < m_pc_points.size(); ++i)
        {
            PC_Point &v = m_pc_points[i];
            if (v.y() == y_min)
            {
                x_point_num++;
            }
            else
            {
                break;
            }
        }
        y_point_num = m_pc_points.size() / x_point_num;

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
                    m_pc_points[center_index].adjacent_pc_points().set_allocation(allocate_pointers(3), 3);
                    m_pc_points[center_index].adjacent_pc_points_distance().set_allocation(allocate_pointers(3), 3);
                    m_pc_points[center_index].adjacent_pc_points()[0] = &m_pc_points[right_index];
                    m_pc_points[center_index].adjacent_pc_points()[1] = &m_pc_points[top_index];
                    m_pc_points[center_index].adjacent_pc_points()[2] = &m_pc_points[top_right_index];
                    m_pc_points[center_index].adjacent_pc_points_distance()[0] = m_pc_points[center_index].distance(&m_pc_points[right_index]);
                    m_pc_points[center_index].adjacent_pc_points_distance()[1] = m_pc_points[center_index].distance(&m_pc_points[top_index]);
                    m_pc_points[center_index].adjacent_pc_points_distance()[2] = m_pc_points[center_index].distance(&m_pc_points[top_right_index]);
                }
                else if (i == x_point_num - 1 && j == 0)
                {
                    m_pc_points[center_index].adjacent_pc_points().set_allocation(allocate_pointers(3), 3);
                    m_pc_points[center_index].adjacent_pc_points_distance().set_allocation(allocate_pointers(3), 3);
                    m_pc_points[center_index].adjacent_pc_points()[0] = &m_pc_points[left_index];
                    m_pc_points[center_index].adjacent_pc_points()[1] = &m_pc_points[top_left_index];
                    m_pc_points[center_index].adjacent_pc_points()[2] = &m_pc_points[top_index];
                    m_pc_points[center_index].adjacent_pc_points_distance()[0] = m_pc_points[center_index].distance(&m_pc_points[left_index]);
                    m_pc_points[center_index].adjacent_pc_points_distance()[1] = m_pc_points[center_index].distance(&m_pc_points[top_left_index]);
                    m_pc_points[center_index].adjacent_pc_points_distance()[2] = m_pc_points[center_index].distance(&m_pc_points[top_index]);
                }
                else if (i == 0 && j == y_point_num - 1)
                {
                    m_pc_points[center_index].adjacent_pc_points().set_allocation(allocate_pointers(3), 3);
                    m_pc_points[center_index].adjacent_pc_points_distance().set_allocation(allocate_pointers(3), 3);
                    m_pc_points[center_index].adjacent_pc_points()[0] = &m_pc_points[bottom_index];
                    m_pc_points[center_index].adjacent_pc_points()[1] = &m_pc_points[bottom_right_index];
                    m_pc_points[center_index].adjacent_pc_points()[2] = &m_pc_points[right_index];
                    m_pc_points[center_index].adjacent_pc_points_distance()[0] = m_pc_points[center_index].distance(&m_pc_points[bottom_index]);
                    m_pc_points[center_index].adjacent_pc_points_distance()[1] = m_pc_points[center_index].distance(&m_pc_points[bottom_right_index]);
                    m_pc_points[center_index].adjacent_pc_points_distance()[2] = m_pc_points[center_index].distance(&m_pc_points[right_index]);
                }
                else if (i == x_point_num - 1 && j == y_point_num - 1)
                {
                    m_pc_points[center_index].adjacent_pc_points().set_allocation(allocate_pointers(3), 3);
                    m_pc_points[center_index].adjacent_pc_points_distance().set_allocation(allocate_pointers(3), 3);
                    m_pc_points[center_index].adjacent_pc_points()[0] = &m_pc_points[bottom_left_index];
                    m_pc_points[center_index].adjacent_pc_points()[1] = &m_pc_points[bottom_index];
                    m_pc_points[center_index].adjacent_pc_points()[2] = &m_pc_points[left_index];
                    m_pc_points[center_index].adjacent_pc_points_distance()[0] = m_pc_points[center_index].distance(&m_pc_points[bottom_left_index]);
                    m_pc_points[center_index].adjacent_pc_points_distance()[1] = m_pc_points[center_index].distance(&m_pc_points[bottom_index]);
                    m_pc_points[center_index].adjacent_pc_points_distance()[2] = m_pc_points[center_index].distance(&m_pc_points[left_index]);
                }
                else if (i != 0 && i != x_point_num - 1 && j == 0)
                {
                    m_pc_points[center_index].adjacent_pc_points().set_allocation(allocate_pointers(5), 5);
                    m_pc_points[center_index].adjacent_pc_points_distance().set_allocation(allocate_pointers(5), 5);
                    m_pc_points[center_index].adjacent_pc_points()[0] = &m_pc_points[left_index];
                    m_pc_points[center_index].adjacent_pc_points()[1] = &m_pc_points[right_index];
                    m_pc_points[center_index].adjacent_pc_points()[2] = &m_pc_points[top_left_index];
                    m_pc_points[center_index].adjacent_pc_points()[3] = &m_pc_points[top_index];
                    m_pc_points[center_index].adjacent_pc_points()[4] = &m_pc_points[top_right_index];
                    m_pc_points[center_index].adjacent_pc_points_distance()[0] = m_pc_points[center_index].distance(&m_pc_points[left_index]);
                    m_pc_points[center_index].adjacent_pc_points_distance()[1] = m_pc_points[center_index].distance(&m_pc_points[right_index]);
                    m_pc_points[center_index].adjacent_pc_points_distance()[2] = m_pc_points[center_index].distance(&m_pc_points[top_left_index]);
                    m_pc_points[center_index].adjacent_pc_points_distance()[3] = m_pc_points[center_index].distance(&m_pc_points[top_index]);
                    m_pc_points[center_index].adjacent_pc_points_distance()[4] = m_pc_points[center_index].distance(&m_pc_points[top_right_index]);
                }
                else if (i == 0 && j != 0 && j != y_point_num - 1)
                {
                    m_pc_points[center_index].adjacent_pc_points().set_allocation(allocate_pointers(5), 5);
                    m_pc_points[center_index].adjacent_pc_points_distance().set_allocation(allocate_pointers(5), 5);
                    m_pc_points[center_index].adjacent_pc_points()[0] = &m_pc_points[bottom_index];
                    m_pc_points[center_index].adjacent_pc_points()[1] = &m_pc_points[bottom_right_index];
                    m_pc_points[center_index].adjacent_pc_points()[2] = &m_pc_points[right_index];
                    m_pc_points[center_index].adjacent_pc_points()[3] = &m_pc_points[top_index];
                    m_pc_points[center_index].adjacent_pc_points()[4] = &m_pc_points[top_right_index];
                    m_pc_points[center_index].adjacent_pc_points_distance()[0] = m_pc_points[center_index].distance(&m_pc_points[bottom_index]);
                    m_pc_points[center_index].adjacent_pc_points_distance()[1] = m_pc_points[center_index].distance(&m_pc_points[bottom_right_index]);
                    m_pc_points[center_index].adjacent_pc_points_distance()[2] = m_pc_points[center_index].distance(&m_pc_points[right_index]);
                    m_pc_points[center_index].adjacent_pc_points_distance()[3] = m_pc_points[center_index].distance(&m_pc_points[top_index]);
                    m_pc_points[center_index].adjacent_pc_points_distance()[4] = m_pc_points[center_index].distance(&m_pc_points[top_right_index]);
                }
                else if (i == x_point_num - 1 && j != 0 && j != y_point_num - 1)
                {
                    m_pc_points[center_index].adjacent_pc_points().set_allocation(allocate_pointers(5), 5);
                    m_pc_points[center_index].adjacent_pc_points_distance().set_allocation(allocate_pointers(5), 5);
                    m_pc_points[center_index].adjacent_pc_points()[0] = &m_pc_points[bottom_left_index];
                    m_pc_points[center_index].adjacent_pc_points()[1] = &m_pc_points[bottom_index];
                    m_pc_points[center_index].adjacent_pc_points()[2] = &m_pc_points[left_index];
                    m_pc_points[center_index].adjacent_pc_points()[3] = &m_pc_points[top_left_index];
                    m_pc_points[center_index].adjacent_pc_points()[4] = &m_pc_points[top_index];
                    m_pc_points[center_index].adjacent_pc_points_distance()[0] = m_pc_points[center_index].distance(&m_pc_points[bottom_left_index]);
                    m_pc_points[center_index].adjacent_pc_points_distance()[1] = m_pc_points[center_index].distance(&m_pc_points[bottom_index]);
                    m_pc_points[center_index].adjacent_pc_points_distance()[2] = m_pc_points[center_index].distance(&m_pc_points[left_index]);
                    m_pc_points[center_index].adjacent_pc_points_distance()[3] = m_pc_points[center_index].distance(&m_pc_points[top_left_index]);
                    m_pc_points[center_index].adjacent_pc_points_distance()[4] = m_pc_points[center_index].distance(&m_pc_points[top_index]);
                }
                else if (i != 0 && i != x_point_num - 1 && j == y_point_num - 1)
                {
                    m_pc_points[center_index].adjacent_pc_points().set_allocation(allocate_pointers(5), 5);
                    m_pc_points[center_index].adjacent_pc_points_distance().set_allocation(allocate_pointers(5), 5);
                    m_pc_points[center_index].adjacent_pc_points()[0] = &m_pc_points[bottom_left_index];
                    m_pc_points[center_index].adjacent_pc_points()[1] = &m_pc_points[bottom_index];
                    m_pc_points[center_index].adjacent_pc_points()[2] = &m_pc_points[bottom_right_index];
                    m_pc_points[center_index].adjacent_pc_points()[3] = &m_pc_points[left_index];
                    m_pc_points[center_index].adjacent_pc_points()[4] = &m_pc_points[right_index];
                    m_pc_points[center_index].adjacent_pc_points_distance()[0] = m_pc_points[center_index].distance(&m_pc_points[bottom_left_index]);
                    m_pc_points[center_index].adjacent_pc_points_distance()[1] = m_pc_points[center_index].distance(&m_pc_points[bottom_index]);
                    m_pc_points[center_index].adjacent_pc_points_distance()[2] = m_pc_points[center_index].distance(&m_pc_points[bottom_right_index]);
                    m_pc_points[center_index].adjacent_pc_points_distance()[3] = m_pc_points[center_index].distance(&m_pc_points[left_index]);
                    m_pc_points[center_index].adjacent_pc_points_distance()[4] = m_pc_points[center_index].distance(&m_pc_points[right_index]);
                }
                else if (i != 0 && i != x_point_num - 1 && j != 0 && j != y_point_num - 1)
                {
                    m_pc_points[center_index].adjacent_pc_points().set_allocation(allocate_pointers(8), 8);
                    m_pc_points[center_index].adjacent_pc_points_distance().set_allocation(allocate_pointers(8), 8);
                    m_pc_points[center_index].adjacent_pc_points()[0] = &m_pc_points[bottom_left_index];
                    m_pc_points[center_index].adjacent_pc_points()[1] = &m_pc_points[bottom_index];
                    m_pc_points[center_index].adjacent_pc_points()[2] = &m_pc_points[bottom_right_index];
                    m_pc_points[center_index].adjacent_pc_points()[3] = &m_pc_points[left_index];
                    m_pc_points[center_index].adjacent_pc_points()[4] = &m_pc_points[right_index];
                    m_pc_points[center_index].adjacent_pc_points()[5] = &m_pc_points[top_left_index];
                    m_pc_points[center_index].adjacent_pc_points()[6] = &m_pc_points[top_index];
                    m_pc_points[center_index].adjacent_pc_points()[7] = &m_pc_points[top_right_index];
                    m_pc_points[center_index].adjacent_pc_points_distance()[0] = m_pc_points[center_index].distance(&m_pc_points[bottom_left_index]);
                    m_pc_points[center_index].adjacent_pc_points_distance()[1] = m_pc_points[center_index].distance(&m_pc_points[bottom_index]);
                    m_pc_points[center_index].adjacent_pc_points_distance()[2] = m_pc_points[center_index].distance(&m_pc_points[bottom_right_index]);
                    m_pc_points[center_index].adjacent_pc_points_distance()[3] = m_pc_points[center_index].distance(&m_pc_points[left_index]);
                    m_pc_points[center_index].adjacent_pc_points_distance()[4] = m_pc_points[center_index].distance(&m_pc_points[right_index]);
                    m_pc_points[center_index].adjacent_pc_points_distance()[5] = m_pc_points[center_index].distance(&m_pc_points[top_left_index]);
                    m_pc_points[center_index].adjacent_pc_points_distance()[6] = m_pc_points[center_index].distance(&m_pc_points[top_index]);
                    m_pc_points[center_index].adjacent_pc_points_distance()[7] = m_pc_points[center_index].distance(&m_pc_points[top_right_index]);
                }
                else
                {
                    std::cout << "adjacent pc points error!" << std::endl;
                    exit(0);
                }
            }
        }

        // std::cout << std::endl;
        // std::cout << "Point cloud has " << m_pc_points.size() << " points" << std::endl;
        // std::cout << "enclosing XYZ box:"
        //           << " X[" << x_min << "," << x_max << "]"
        //           << " Y[" << y_min << "," << y_max << "]"
        //           << " Z[" << z_min << "," << z_max << "]"
        //           << std::endl;
        // std::cout << std::endl;
    }

    void PointCloud::point_cloud_to_terrain(double &memory_usage)
    {
        double x_min = 1e100;
        double x_max = -1e100;
        double y_min = 1e100;
        double y_max = -1e100;
        double z_min = 1e100;
        double z_max = -1e100;
        for (unsigned i = 0; i < m_pc_points.size(); ++i)
        {
            PC_Point &v = m_pc_points[i];
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

        for (unsigned i = 0; i < m_pc_points.size(); ++i)
        {
            PC_Point &v = m_pc_points[i];
            if (v.y() == y_min)
            {
                x_point_num++;
            }
            else
            {
                break;
            }
        }
        y_point_num = m_pc_points.size() / x_point_num;

        std::string terrain_write_path = "temp_terrain.off";
        std::ofstream ofs(&terrain_write_path[0], std::ofstream::trunc);

        ofs << "OFF\n";

        int total_vertex_num = x_point_num * y_point_num;
        int total_face_num = (x_point_num - 1) * (y_point_num - 1) * 2;
        int total_edge_num = (x_point_num - 1) * y_point_num + x_point_num * (y_point_num - 1) + (x_point_num - 1) * (y_point_num - 1);

        ofs << total_vertex_num << " " << total_face_num << " " << total_edge_num << " \n";

        for (unsigned i = 0; i < m_pc_points.size(); ++i)
        {
            PC_Point &v = m_pc_points[i];
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

    // ========= point cloud =========

    // ========= algorithm base =========
    class PointCloudGeodesicAlgorithmBase
    {
    public:
        enum AlgorithmType
        {
            DIJKSTRA,
            UNDEFINED_ALGORITHM
        };

        PointCloudGeodesicAlgorithmBase(point_cloud_geodesic::PointCloud *point_cloud) : m_type(UNDEFINED_ALGORITHM),
                                                                                         m_max_propagation_distance(1e100),
                                                                                         m_point_cloud(point_cloud) {};

        virtual ~PointCloudGeodesicAlgorithmBase() {};

        // virtual void propagate(std::vector<PathPoint> &sources,
        //                        double max_propagation_distance = INFIN,           // propagation algorithm stops after reaching the certain distance from the source
        //                        std::vector<PathPoint> *stop_points = NULL) = 0; // or after ensuring that all the stop_points are covered

        // virtual void trace_back(PathPoint &destination, // trace back piecewise-linear path
        //                         std::vector<PathPoint> &path) = 0;

        virtual unsigned best_source(PathPoint &point, // after propagation step is done, quickly find what source this point belongs to and what is the distance to this source
                                     double &best_source_distance) = 0;

        virtual void print_statistics() // print info about timing and memory usage in the propagation step of the algorithm
        {
            std::cout << "propagation step took " << m_time_consumed << " seconds " << std::endl;
        };

        AlgorithmType type() { return m_type; };

        virtual std::string name();

        point_cloud_geodesic::PointCloud *point_cloud() { return m_point_cloud; };

    protected:
        void set_stop_conditions(std::vector<PathPoint> *stop_points);

        void set_stop_conditions(double stop_distance);

        double stop_distance()
        {
            return m_max_propagation_distance;
        }

        AlgorithmType m_type; // type of the algorithm

        typedef std::pair<pc_point_pointer, double> stop_pc_point_with_distace_type;
        std::vector<stop_pc_point_with_distace_type> m_stop_pc_points; // algorithm stops propagation after covering certain pc points
        double m_max_propagation_distance;                             // or reaching the certain distance

        point_cloud_geodesic::PointCloud *m_point_cloud;

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

    inline std::string PointCloudGeodesicAlgorithmBase::name()
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

    inline void PointCloudGeodesicAlgorithmBase::set_stop_conditions(std::vector<PathPoint> *stop_points)
    {
        if (!stop_points)
        {
            m_stop_pc_points.clear();
            return;
        }

        m_stop_pc_points.resize(stop_points->size());

        std::vector<pc_point_pointer> possible_pc_points;
        for (unsigned i = 0; i < stop_points->size(); ++i)
        {
            PathPoint *point = &(*stop_points)[i];

            possible_pc_points.clear();
            m_point_cloud->closest_pc_points(point, &possible_pc_points);

            pc_point_pointer closest_vertex = NULL;
            double min_distance = 1e100;
            for (unsigned j = 0; j < possible_pc_points.size(); ++j)
            {
                double distance = point->distance(possible_pc_points[j]);
                if (distance < min_distance)
                {
                    min_distance = distance;
                    closest_vertex = possible_pc_points[j];
                }
            }
            assert(closest_vertex);

            m_stop_pc_points[i].first = closest_vertex;
            m_stop_pc_points[i].second = min_distance;
        }
    }

    inline void PointCloudGeodesicAlgorithmBase::set_stop_conditions(double stop_distance)
    {
        m_max_propagation_distance = stop_distance;
    }
    // ========= algorithm base =========

    // ========= algorithm graph base =========
    template <class Node>
    class PointCloudGeodesicAlgorithmGraphBase : public PointCloudGeodesicAlgorithmBase
    {
    public:
        typedef Node *node_pointer;

        PointCloudGeodesicAlgorithmGraphBase(point_cloud_geodesic::PointCloud *point_cloud) : PointCloudGeodesicAlgorithmBase(point_cloud) {};

        ~PointCloudGeodesicAlgorithmGraphBase() {};

        void propagate(std::vector<PathPoint> &sources,
                       double max_propagation_distance);

        void propagate(std::vector<PathPoint> &sources,
                       std::vector<PathPoint> *stop_points);

        void trace_back(PathPoint &destination,
                        std::vector<PathPoint> &path);

        unsigned best_source(PathPoint &point,
                             double &best_source_distance);

        void print_statistics()
        {
            PointCloudGeodesicAlgorithmBase::print_statistics();

            double memory = m_nodes.size() * sizeof(Node);
            std::cout << "uses about " << memory / 1e6 << "Mb of memory" << std::endl;
        }

        double get_memory();

    protected:
        unsigned node_index(pc_point_pointer v) // gives index of the node that corresponds to this pc point
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
            if (point.type() == PC_POINT)
            {
                pc_point_pointer v = (pc_point_pointer)point.base_element();
                best_node = &m_nodes[node_index(v)];
                best_total_distance = best_node->distance_from_source();
            }

            // assert(best_node);
            // assert(best_total_distance<INFIN);
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

        virtual void list_nodes_visible_from_source(PointCloudElementBase *p,
                                                    std::vector<node_pointer> &storage) = 0; // list all nodes that belong to this point cloud element

        virtual void list_nodes_visible_from_node(node_pointer node, // list all nodes that belong to this point cloud element
                                                  std::vector<node_pointer> &storage,
                                                  std::vector<double> &distances,
                                                  double threshold_distance) = 0; // list only the nodes whose current distance is larger than the threshold

        std::vector<Node> m_nodes; // nodes of the graph

        typedef std::set<node_pointer, Node> queue_type;
        queue_type m_queue;

        std::vector<PathPoint> m_sources; // for simplicity, we keep sources as they are
    };

    template <class Node>
    void PointCloudGeodesicAlgorithmGraphBase<Node>::propagate(std::vector<PathPoint> &sources,
                                                               double max_propagation_distance)
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

        std::vector<node_pointer> visible_nodes; // initialize pc points directly visible from sources
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
            if (counter++ % 10 == 0) // check if we covered all required pc points
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

            for (unsigned i = 0; i < visible_nodes.size(); ++i) // update all the adgecent pc points
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
        // std::cout << std::endl;
    }

    template <class Node>
    void PointCloudGeodesicAlgorithmGraphBase<Node>::propagate(std::vector<PathPoint> &sources,
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

        std::vector<node_pointer> visible_nodes; // initialize pc points directly visible from sources
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
            if (counter++ % 10 == 0) // check if we covered all required pc points
            {
                if (check_stop_conditions_cover_points(satisfied_index))
                {
                    // std::cout << "break" << std::endl;
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

            for (unsigned i = 0; i < visible_nodes.size(); ++i) // update all the adgecent pc points
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
        // std::cout << std::endl;
    }

    template <class Node>
    inline bool PointCloudGeodesicAlgorithmGraphBase<Node>::check_stop_conditions_distance()
    {
        double queue_min_distance = (*m_queue.begin())->distance_from_source();
        if (queue_min_distance < m_max_propagation_distance)
        {
            return false;
        }
        return true;
    }

    template <class Node>
    inline bool PointCloudGeodesicAlgorithmGraphBase<Node>::check_stop_conditions_cover_points(unsigned &index)
    {
        double queue_min_distance = (*m_queue.begin())->distance_from_source();
        while (index < m_stop_pc_points.size())
        {
            pc_point_pointer v = m_stop_pc_points[index].first;
            Node &node = m_nodes[node_index(v)];
            if (queue_min_distance < node.distance_from_source() + m_stop_pc_points[index].second)
            {
                return false;
            }
            ++index;
        }
        return true;
    }

    template <class Node>
    inline void PointCloudGeodesicAlgorithmGraphBase<Node>::trace_back(PathPoint &destination, // trace back piecewise-linear path
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
    inline double PointCloudGeodesicAlgorithmGraphBase<Node>::get_memory()
    {
        double memory_usage = m_nodes.size() * sizeof(Node);
        return memory_usage;
    }

    template <class Node>
    inline unsigned PointCloudGeodesicAlgorithmGraphBase<Node>::best_source(PathPoint &point, // quickly find what source this point belongs to and what is the distance to this source
                                                                            double &best_source_distance)
    {
        node_pointer node = best_first_node(point, best_source_distance);
        return node ? node->source_index() : 0;
    };

    // ========= algorithm graph base =========

    // ========= Dijkstra shortest path algorithm =========
    class PointCloudDijkstraNode
    {
        typedef PointCloudDijkstraNode *node_pointer;

    public:
        PointCloudDijkstraNode() {};
        ~PointCloudDijkstraNode() {};

        double &distance_from_source() { return m_distance; };
        node_pointer &previous() { return m_previous; };
        unsigned &source_index() { return m_source_index; };
        pc_point_pointer &pc_point() { return m_pc_point; };

        void clear()
        {
            m_distance = INFIN;
            m_previous = NULL;
        }

        bool operator()(node_pointer const s1, node_pointer const s2) const
        {
            return s1->distance_from_source() != s2->distance_from_source() ? s1->distance_from_source() < s2->distance_from_source() : s1->pc_point()->id() < s2->pc_point()->id();
        };

        double distance(PathPoint *p)
        {
            return m_pc_point->distance(p);
        }

        PathPoint path_point()
        {
            return PathPoint(m_pc_point);
        }

    private:
        double m_distance;           // distance to the closest source
        unsigned m_source_index;     // closest source index
        node_pointer m_previous;     // previous node in the geodesic path
        pc_point_pointer m_pc_point; // correspoding pc point
    };

    class PointCloudGeodesicAlgorithmDijkstra : public PointCloudGeodesicAlgorithmGraphBase<PointCloudDijkstraNode>
    {
    public:
        typedef PointCloudDijkstraNode Node;
        typedef Node *node_pointer;

        PointCloudGeodesicAlgorithmDijkstra(point_cloud_geodesic::PointCloud *point_cloud) : PointCloudGeodesicAlgorithmGraphBase<Node>(point_cloud)
        {
            m_type = DIJKSTRA;

            m_nodes.resize(point_cloud->pc_points().size());
            for (unsigned i = 0; i < m_nodes.size(); ++i)
            {
                m_nodes[i].pc_point() = &m_point_cloud->pc_points()[i];
            }
        };

        ~PointCloudGeodesicAlgorithmDijkstra() {};

    protected:
        void list_nodes_visible_from_source(PointCloudElementBase *p,
                                            std::vector<node_pointer> &storage); // list all nodes that belong to this point_cloud element

        void list_nodes_visible_from_node(node_pointer node, // list all nodes that belong to this point_cloud element
                                          std::vector<node_pointer> &storage,
                                          std::vector<double> &distances,
                                          double threshold_distance); // list only the nodes whose current distance is larger than the threshold
    };

    void PointCloudGeodesicAlgorithmDijkstra::list_nodes_visible_from_source(PointCloudElementBase *p,
                                                                             std::vector<node_pointer> &storage)
    {
        assert(p->type() != UNDEFINED_POINT);

        if (p->type() == PC_POINT)
        {
            pc_point_pointer v = static_cast<pc_point_pointer>(p);
            storage.push_back(&m_nodes[node_index(v)]);
        }
    }

    inline void PointCloudGeodesicAlgorithmDijkstra::list_nodes_visible_from_node(node_pointer node, // list all nodes that belong to this point cloud element
                                                                                  std::vector<node_pointer> &storage,
                                                                                  std::vector<double> &distances,
                                                                                  double threshold_distance)
    {
        pc_point_pointer v = node->pc_point();
        assert(storage.size() == distances.size());

        for (unsigned i = 0; i < v->adjacent_pc_points().size(); ++i)
        {
            pc_point_pointer new_p = v->adjacent_pc_points()[i];
            node_pointer new_node = &m_nodes[node_index(new_p)];

            if (new_node->distance_from_source() > threshold_distance + v->adjacent_pc_points_distance()[i])
            {
                storage.push_back(new_node);
                distances.push_back(v->adjacent_pc_points_distance()[i]);
            }
        }
    }
    // ========= Dijkstra shortest path algorithm =========

} // geodesic