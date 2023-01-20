
struct aabb
{
    float3 upper;
    float3 lower;
};

__device__ __host__
inline bool intersects(const aabb& lhs, const aabb& rhs) noexcept
{
if(lhs.upper.x < rhs.lower.x || rhs.upper.x < lhs.lower.x) {return false;}
if(lhs.upper.y < rhs.lower.y || rhs.upper.y < lhs.lower.y) {return false;}
if(lhs.upper.z < rhs.lower.z || rhs.upper.z < lhs.lower.z) {return false;}
return true;
}

struct merge {
    __device__
     aabb operator()(const aabb &lhs, const aabb &rhs) const noexcept {
        aabb merged;
        merged.upper.x = ::fmax(lhs.upper.x, rhs.upper.x);
        merged.upper.y = ::fmax(lhs.upper.y, rhs.upper.y);
        merged.upper.z = ::fmax(lhs.upper.z, rhs.upper.z);
        merged.lower.x = ::fmin(lhs.lower.x, rhs.lower.x);
        merged.lower.y = ::fmin(lhs.lower.y, rhs.lower.y);
        merged.lower.z = ::fmin(lhs.lower.z, rhs.lower.z);
        return merged;
    }
};

__device__ __host__
inline float3 centroid(const aabb& box) noexcept
{
float3 c;
c.x = (box.upper.x + box.lower.x) * 0.5;
c.y = (box.upper.y + box.lower.y) * 0.5;
c.z = (box.upper.z + box.lower.z) * 0.5;
return c;
}

__device__ __host__
inline std::uint32_t expand_bits(std::uint32_t v) noexcept
{
    v = (v * 0x00010001u) & 0xFF0000FFu;
    v = (v * 0x00000101u) & 0x0F00F00Fu;
    v = (v * 0x00000011u) & 0xC30C30C3u;
    v = (v * 0x00000005u) & 0x49249249u;
    return v;
}

struct m64
{
__device__
    unsigned long long int operator()(const unsigned int m, const unsigned int idx)
    {
        unsigned long long int m64 = m;
        m64 <<= 32;
        m64 |= idx;
        return m64;
    }
};

// Calculates a 30-bit Morton code for the
// given 3D point located within the unit cube [0,1].
__device__ __host__
inline std::uint32_t morton_code(float3 xyz, float resolution = 1024.0f) noexcept
{
    xyz.x = ::fminf(::fmaxf(xyz.x * resolution, 0.0f), resolution - 1.0f);
    xyz.y = ::fminf(::fmaxf(xyz.y * resolution, 0.0f), resolution - 1.0f);
    xyz.z = ::fminf(::fmaxf(xyz.z * resolution, 0.0f), resolution - 1.0f);
    const std::uint32_t xx = expand_bits(static_cast<std::uint32_t>(xyz.x));
    const std::uint32_t yy = expand_bits(static_cast<std::uint32_t>(xyz.y));
    const std::uint32_t zz = expand_bits(static_cast<std::uint32_t>(xyz.z));
    return xx * 4 + yy * 2 + zz;
}

__device__ __host__
inline std::uint32_t morton_code(double3 xyz, double resolution = 1024.0) noexcept
{
    xyz.x = ::fmin(::fmax(xyz.x * resolution, 0.0), resolution - 1.0);
    xyz.y = ::fmin(::fmax(xyz.y * resolution, 0.0), resolution - 1.0);
    xyz.z = ::fmin(::fmax(xyz.z * resolution, 0.0), resolution - 1.0);
    const std::uint32_t xx = expand_bits(static_cast<std::uint32_t>(xyz.x));
    const std::uint32_t yy = expand_bits(static_cast<std::uint32_t>(xyz.y));
    const std::uint32_t zz = expand_bits(static_cast<std::uint32_t>(xyz.z));
    return xx * 4 + yy * 2 + zz;
}