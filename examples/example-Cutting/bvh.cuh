#include "aabb.cuh"
#include <stdint.h>
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>



struct node
{
    unsigned int parent_idx; // parent node
    unsigned int left_idx;   // index of left  child node
    unsigned int right_idx;  // index of right child node
    unsigned int object_idx; // == 0xFFFFFFFF if internal node.
};

struct buildNode
{
    __device__
    node operator()(int idx)
    {
        node n;
        n.parent_idx = 0xFFFFFFFF;
        n.left_idx   = 0xFFFFFFFF;
        n.right_idx  = 0xFFFFFFFF;
        n.object_idx = idx;
        return n;
    }
};


template<typename Real, typename Object>
struct default_morton_code_calculator
{
    default_morton_code_calculator(aabb w): whole(w) {}
    default_morton_code_calculator()  = default;
    ~default_morton_code_calculator() = default;
    default_morton_code_calculator(default_morton_code_calculator const&) = default;
    default_morton_code_calculator(default_morton_code_calculator&&)      = default;
    default_morton_code_calculator& operator=(default_morton_code_calculator const&) = default;
    default_morton_code_calculator& operator=(default_morton_code_calculator&&)      = default;

    __device__ __host__
    inline unsigned int operator()(const Object&, const aabb& box) noexcept
    {
        auto p = centroid(box);
        p.x -= whole.lower.x;
        p.y -= whole.lower.y;
        p.z -= whole.lower.z;
        p.x /= (whole.upper.x - whole.lower.x);
        p.y /= (whole.upper.y - whole.lower.y);
        p.z /= (whole.upper.z - whole.lower.z);
        return morton_code(p);
    }
    aabb whole;
};

template<typename Real, typename Object, typename AABBGetter,
        typename MortonCodeCalculator = default_morton_code_calculator<Real, Object>>
class bvh
{

public:

    template<typename InputIterator>
    bvh(InputIterator first, InputIterator last)
            : objects_d_(first, last)
    {
        this->construct();
    }

    void clear()
    {
        this->objects_h_.clear();
        this->objects_d_.clear();
        this->aabbs_h_.clear();
        this->aabbs_.clear();
        this->nodes_h_.clear();
        this->nodes_.clear();
        return ;
    }

    void construct()
    {
        objects_h_ = objects_d_;

        assert(objects_h_.size() == objects_d_.size());
        if(objects_h_.size() == 0u) {
            return;}

        const unsigned int num_objects        = objects_h_.size();
        const unsigned int num_internal_nodes = num_objects - 1;
        const unsigned int num_nodes          = num_objects * 2 - 1;

        // --------------------------------------------------------------------
        // calculate morton code of each points

        const auto inf = std::numeric_limits<Real>::infinity();
        aabb default_aabb;
        default_aabb.upper.x = -inf; default_aabb.lower.x = inf;
        default_aabb.upper.y = -inf; default_aabb.lower.y = inf;
        default_aabb.upper.z = -inf; default_aabb.lower.z = inf;

        this->aabbs_.resize(num_nodes, default_aabb);

        thrust::transform(this->objects_d_.begin(), this->objects_d_.end(),
                          aabbs_.begin() + num_internal_nodes, AABBGetter());

        const auto aabb_whole = thrust::reduce(
                aabbs_.begin() + num_internal_nodes, aabbs_.end(), default_aabb,
                ::merge());

        thrust::device_vector<::uint32_t> morton(num_objects);
        thrust::transform(this->objects_d_.begin(), this->objects_d_.end(),
                          aabbs_.begin() + num_internal_nodes, morton.begin(),
                          MortonCodeCalculator(aabb_whole));


        // --------------------------------------------------------------------
        // sort object-indices by morton code

        thrust::device_vector<::uint32_t> indices(num_objects);
        thrust::copy(thrust::make_counting_iterator<int>(0),
                     thrust::make_counting_iterator<int>(num_objects),
                     indices.begin());
        // keep indices ascending order
        thrust::stable_sort_by_key(morton.begin(), morton.end(),
                                   thrust::make_zip_iterator(
                                           thrust::make_tuple(aabbs_.begin() + num_internal_nodes,
                                                              indices.begin())));

        // --------------------------------------------------------------------
        // check morton codes are unique

        thrust::device_vector<unsigned long long int> morton64(num_objects);
        const auto uniqued = thrust::unique_copy(morton.begin(), morton.end(),
                                                 morton64.begin());

        const bool morton_code_is_unique = (morton64.end() == uniqued);
        if(!morton_code_is_unique)
        {
            thrust::transform(morton.begin(), morton.end(), indices.begin(),
                              morton64.begin(),m64());
        }

        // --------------------------------------------------------------------
        // construct leaf nodes and aabbs

        node default_node;
        default_node.parent_idx = 0xFFFFFFFF;
        default_node.left_idx   = 0xFFFFFFFF;
        default_node.right_idx  = 0xFFFFFFFF;
        default_node.object_idx = 0xFFFFFFFF;
        this->nodes_.resize(num_nodes, default_node);

        thrust::transform(indices.begin(), indices.end(),
                          this->nodes_.begin() + num_internal_nodes, buildNode());


        // --------------------------------------------------------------------
        // construct internal nodes

        const auto self = this->get_device_repr();
        if(morton_code_is_unique)
        {
            const unsigned int* node_code = morton.data().get();
            construct_internal_nodes(self, node_code, num_objects);
        }
        else // 64bit version
        {
            const unsigned long long int* node_code = morton64.data().get();
            construct_internal_nodes(self, node_code, num_objects);
        }

        return;
    }

private:

    // check to see if host vector is needed
    thrust::host_vector  <object_type>   objects_h_;
    thrust::device_vector<object_type>   objects_d_;
    thrust::host_vector  <aabb>     aabbs_h_;
    thrust::device_vector<aabb>     aabbs_;
    thrust::host_vector  <node>     nodes_h_;
    thrust::device_vector<node>     nodes_;
};