#include <Eigen/Dense>
#include <thrust/functional.h>

#define _gamma 5.828427124 // FOUR_GAMMA_SQUARED = sqrt(8)+3;
#define _cstar 0.923879532 // cos(pi/8)
#define _sstar 0.3826834323 // sin(p/8)
#define EPSILON 1e-6

struct object_type
{
    // end triangle
    double end[3][3];

    // start triangle
    double start[3][3];

};

__host__ __device__ __forceinline__
double accurateSqrt(double x)
{
    return x * rsqrt(x);
}

__host__ __device__ __forceinline__
void condSwap(bool c, double &X, double &Y)
{
    // used in step 2
    double Z = X;
    X = c ? Y : X;
    Y = c ? Z : Y;
}

__host__ __device__ __forceinline__
void condNegSwap(bool c, double &X, double &Y)
{
    // used in step 2 and 3
    double Z = -X;
    X = c ? Y : X;
    Y = c ? Z : Y;
}

// matrix multiplication M = A * B
__host__ __device__ __forceinline__
void multAB(double a11, double a12, double a13,
            double a21, double a22, double a23,
            double a31, double a32, double a33,
        //
            double b11, double b12, double b13,
            double b21, double b22, double b23,
            double b31, double b32, double b33,
        //
            double &m11, double &m12, double &m13,
            double &m21, double &m22, double &m23,
            double &m31, double &m32, double &m33)
{

    m11=a11*b11 + a12*b21 + a13*b31; m12=a11*b12 + a12*b22 + a13*b32; m13=a11*b13 + a12*b23 + a13*b33;
    m21=a21*b11 + a22*b21 + a23*b31; m22=a21*b12 + a22*b22 + a23*b32; m23=a21*b13 + a22*b23 + a23*b33;
    m31=a31*b11 + a32*b21 + a33*b31; m32=a31*b12 + a32*b22 + a33*b32; m33=a31*b13 + a32*b23 + a33*b33;
}

// matrix multiplication M = Transpose[A] * B
__host__ __device__ __forceinline__
void multAtB(double a11, double a12, double a13,
             double a21, double a22, double a23,
             double a31, double a32, double a33,
        //
             double b11, double b12, double b13,
             double b21, double b22, double b23,
             double b31, double b32, double b33,
        //
             double &m11, double &m12, double &m13,
             double &m21, double &m22, double &m23,
             double &m31, double &m32, double &m33)
{
    m11=a11*b11 + a21*b21 + a31*b31; m12=a11*b12 + a21*b22 + a31*b32; m13=a11*b13 + a21*b23 + a31*b33;
    m21=a12*b11 + a22*b21 + a32*b31; m22=a12*b12 + a22*b22 + a32*b32; m23=a12*b13 + a22*b23 + a32*b33;
    m31=a13*b11 + a23*b21 + a33*b31; m32=a13*b12 + a23*b22 + a33*b32; m33=a13*b13 + a23*b23 + a33*b33;
}

__host__ __device__ __forceinline__
void quatToMat3(const double * qV,
                double &m11, double &m12, double &m13,
                double &m21, double &m22, double &m23,
                double &m31, double &m32, double &m33
)
{
    double w = qV[3];
    double x = qV[0];
    double y = qV[1];
    double z = qV[2];

    double qxx = x*x;
    double qyy = y*y;
    double qzz = z*z;
    double qxz = x*z;
    double qxy = x*y;
    double qyz = y*z;
    double qwx = w*x;
    double qwy = w*y;
    double qwz = w*z;

    m11=1 - 2*(qyy + qzz); m12=2*(qxy - qwz); m13=2*(qxz + qwy);
    m21=2*(qxy + qwz); m22=1 - 2*(qxx + qzz); m23=2*(qyz - qwx);
    m31=2*(qxz - qwy); m32=2*(qyz + qwx); m33=1 - 2*(qxx + qyy);
}

__host__ __device__ __forceinline__
void approximateGivensQuaternion(double a11, double a12, double a22, double &ch, double &sh)
{
/*
     * Given givens angle computed by approximateGivensAngles,
     * compute the corresponding rotation quaternion.
     */
    ch = 2*(a11-a22);
    sh = a12;
    bool b = _gamma*sh*sh < ch*ch;
    double w = rsqrt(ch*ch+sh*sh);
    ch=b?w*ch:_cstar;
    sh=b?w*sh:_sstar;
}

__host__ __device__ __forceinline__
void jacobiConjugation( const int x, const int y, const int z,
                        double &s11,
                        double &s21, double &s22,
                        double &s31, double &s32, double &s33,
                        double * qV)
{
    double ch,sh;
    approximateGivensQuaternion(s11,s21,s22,ch,sh);

    double scale = ch*ch+sh*sh;
    double a = (ch*ch-sh*sh)/scale;
    double b = (2*sh*ch)/scale;

    // make temp copy of S
    double _s11 = s11;
    double _s21 = s21; double _s22 = s22;
    double _s31 = s31; double _s32 = s32; double _s33 = s33;

    // perform conjugation S = Q'*S*Q
    // Q already implicitly solved from a, b
    s11 =a*(a*_s11 + b*_s21) + b*(a*_s21 + b*_s22);
    s21 =a*(-b*_s11 + a*_s21) + b*(-b*_s21 + a*_s22);	s22=-b*(-b*_s11 + a*_s21) + a*(-b*_s21 + a*_s22);
    s31 =a*_s31 + b*_s32;								s32=-b*_s31 + a*_s32; s33=_s33;

    // update cumulative rotation qV
    double tmp[3];
    tmp[0]=qV[0]*sh;
    tmp[1]=qV[1]*sh;
    tmp[2]=qV[2]*sh;
    sh *= qV[3];

    qV[0] *= ch;
    qV[1] *= ch;
    qV[2] *= ch;
    qV[3] *= ch;

    // (x,y,z) corresponds to ((0,1,2),(1,2,0),(2,0,1))
    // for (p,q) = ((0,1),(1,2),(0,2))
    qV[z] += sh;
    qV[3] -= tmp[z]; // w
    qV[x] += tmp[y];
    qV[y] -= tmp[x];

    // re-arrange matrix for next iteration
    _s11 = s22;
    _s21 = s32; _s22 = s33;
    _s31 = s21; _s32 = s31; _s33 = s11;
    s11 = _s11;
    s21 = _s21; s22 = _s22;
    s31 = _s31; s32 = _s32; s33 = _s33;

}

__host__ __device__ __forceinline__
double dist2(double x, double y, double z)
{
    return x*x+y*y+z*z;
}

// finds transformation that diagonalizes a symmetric matrix
__host__ __device__ __forceinline__
void jacobiEigenanlysis( // symmetric matrix
        double &s11,
        double &s21, double &s22,
        double &s31, double &s32, double &s33,
        // quaternion representation of V
        double * qV)
{
    qV[3]=1; qV[0]=0;qV[1]=0;qV[2]=0; // follow same indexing convention as GLM
    for (int i=0;i<4;i++)
    {
        // we wish to eliminate the maximum off-diagonal element
        // on every iteration, but cycling over all 3 possible rotations
        // in fixed order (p,q) = (1,2) , (2,3), (1,3) still retains
        //  asymptotic convergence
        jacobiConjugation(0,1,2,s11,s21,s22,s31,s32,s33,qV); // p,q = 0,1
        jacobiConjugation(1,2,0,s11,s21,s22,s31,s32,s33,qV); // p,q = 1,2
        jacobiConjugation(2,0,1,s11,s21,s22,s31,s32,s33,qV); // p,q = 0,2
    }
}

__host__ __device__ __forceinline__
void sortSingularValues(// matrix that we want to decompose
        double &b11, double &b12, double &b13,
        double &b21, double &b22, double &b23,
        double &b31, double &b32, double &b33,
        // sort V simultaneously
        double &v11, double &v12, double &v13,
        double &v21, double &v22, double &v23,
        double &v31, double &v32, double &v33)
{
    double rho1 = dist2(b11,b21,b31);
    double rho2 = dist2(b12,b22,b32);
    double rho3 = dist2(b13,b23,b33);
    bool c;
    c = rho1 < rho2;
    condNegSwap(c,b11,b12); condNegSwap(c,v11,v12);
    condNegSwap(c,b21,b22); condNegSwap(c,v21,v22);
    condNegSwap(c,b31,b32); condNegSwap(c,v31,v32);
    condSwap(c,rho1,rho2);
    c = rho1 < rho3;
    condNegSwap(c,b11,b13); condNegSwap(c,v11,v13);
    condNegSwap(c,b21,b23); condNegSwap(c,v21,v23);
    condNegSwap(c,b31,b33); condNegSwap(c,v31,v33);
    condSwap(c,rho1,rho3);
    c = rho2 < rho3;
    condNegSwap(c,b12,b13); condNegSwap(c,v12,v13);
    condNegSwap(c,b22,b23); condNegSwap(c,v22,v23);
    condNegSwap(c,b32,b33); condNegSwap(c,v32,v33);
}

__host__ __device__ __forceinline__
void QRGivensQuaternion(double a1, double a2, double &ch, double &sh)
{
    // a1 = pivot point on diagonal
    // a2 = lower triangular entry we want to annihilate
    double epsilon = EPSILON;
    double rho = accurateSqrt(a1*a1 + a2*a2);

    sh = rho > epsilon ? a2 : 0;
    ch = fabs(a1) + fmax(rho,epsilon);
    bool b = a1 < 0;
    condSwap(b,sh,ch);
    double w = rsqrt(ch*ch+sh*sh);
    ch *= w;
    sh *= w;
}

__host__ __device__ __forceinline__
void QRDecomposition(// matrix that we want to decompose
        double b11, double b12, double b13,
        double b21, double b22, double b23,
        double b31, double b32, double b33,
        // output Q
        double &q11, double &q12, double &q13,
        double &q21, double &q22, double &q23,
        double &q31, double &q32, double &q33,
        // output R
        double &r11, double &r12, double &r13,
        double &r21, double &r22, double &r23,
        double &r31, double &r32, double &r33)
{
    double ch1,sh1,ch2,sh2,ch3,sh3;
    double a,b;

    // first givens rotation (ch,0,0,sh)
    QRGivensQuaternion(b11,b21,ch1,sh1);
    a=1-2*sh1*sh1;
    b=2*ch1*sh1;
    // apply B = Q' * B
    r11=a*b11+b*b21;  r12=a*b12+b*b22;  r13=a*b13+b*b23;
    r21=-b*b11+a*b21; r22=-b*b12+a*b22; r23=-b*b13+a*b23;
    r31=b31;          r32=b32;          r33=b33;

    // second givens rotation (ch,0,-sh,0)
    QRGivensQuaternion(r11,r31,ch2,sh2);
    a=1-2*sh2*sh2;
    b=2*ch2*sh2;
    // apply B = Q' * B;
    b11=a*r11+b*r31;  b12=a*r12+b*r32;  b13=a*r13+b*r33;
    b21=r21;           b22=r22;           b23=r23;
    b31=-b*r11+a*r31; b32=-b*r12+a*r32; b33=-b*r13+a*r33;

    // third givens rotation (ch,sh,0,0)
    QRGivensQuaternion(b22,b32,ch3,sh3);
    a=1-2*sh3*sh3;
    b=2*ch3*sh3;
    // R is now set to desired value
    r11=b11;             r12=b12;           r13=b13;
    r21=a*b21+b*b31;     r22=a*b22+b*b32;   r23=a*b23+b*b33;
    r31=-b*b21+a*b31;    r32=-b*b22+a*b32;  r33=-b*b23+a*b33;

    // construct the cumulative rotation Q=Q1 * Q2 * Q3
    // the number of doubleing point operations for three quaternion multiplications
    // is more or less comparable to the explicit form of the joined matrix.
    // certainly more memory-efficient!
    double sh12=sh1*sh1;
    double sh22=sh2*sh2;
    double sh32=sh3*sh3;

    q11=(-1+2*sh12)*(-1+2*sh22);
    q12=4*ch2*ch3*(-1+2*sh12)*sh2*sh3+2*ch1*sh1*(-1+2*sh32);
    q13=4*ch1*ch3*sh1*sh3-2*ch2*(-1+2*sh12)*sh2*(-1+2*sh32);

    q21=2*ch1*sh1*(1-2*sh22);
    q22=-8*ch1*ch2*ch3*sh1*sh2*sh3+(-1+2*sh12)*(-1+2*sh32);
    q23=-2*ch3*sh3+4*sh1*(ch3*sh1*sh3+ch1*ch2*sh2*(-1+2*sh32));

    q31=2*ch2*sh2;
    q32=2*ch3*(1-2*sh22)*sh3;
    q33=(-1+2*sh22)*(-1+2*sh32);
}

__host__ __device__ double det(Eigen::Matrix3d A)
{
    double r1 = A(0,0) * ((A(1,1) * A(2,2))
                    - (A(2,1) * A(1,2)));

    double r2 = A(0,1) * ((A(1,0) * A(2,2))
                    - (A(2,0) * A(1,2)));

    double r3 = A(0,2) * ((A(1,0) * A(2,1))
                    - (A(2,0) * A(1,1)));

    return (r1 - r2 + r3);
}

__host__ __device__ Eigen::Matrix3d inv(Eigen::Matrix3d A)
{
    double d = det(A);
    Eigen::Matrix3d A_inv;

    A(0,0) = A(1,1)*A(2,2) - A(1,2)*A(2,1);
    A(0,1) = A(1,0)*A(2,2) - A(1,2)*A(2,0);
    A(0,2) = A(1,0)*A(2,1) - A(1,1)*A(2,0);

    A(1,0) = A(0,1)*A(2,2) - A(0,2)*A(2,1);
    A(1,1) = A(0,0)*A(2,2) - A(0,2)*A(2,0);
    A(1,2) = A(0,0)*A(2,1) - A(0,1)*A(2,0);

    A(2,0) = A(0,1)*A(1,2) - A(0,2)*A(1,1);
    A(2,1) = A(0,0)*A(1,2) - A(0,2)*A(1,0);
    A(2,2) = A(0,0)*A(1,1) - A(0,1)*A(1,0);

    return A;

}

__host__ __device__ __forceinline__ Eigen::Matrix3d mul(Eigen::Matrix3d A, Eigen::Matrix3d B)
{
    Eigen::Matrix3d ret;

    for (int i = 0; i < 3; i++) {
        for (int k = 0; k < 3; k++) {
            ret(i,k) = A(i,0) * B(0,k) + A(i,1) * B(1,k) + A(i,2) * B(2,k);
        }
    }

    return ret;
}

__host__ __device__ __forceinline__
void svd(// input A
        double a11, double a12, double a13,
        double a21, double a22, double a23,
        double a31, double a32, double a33,
        // output U
        double &u11, double &u12, double &u13,
        double &u21, double &u22, double &u23,
        double &u31, double &u32, double &u33,
        // output S
        double &s11, double &s12, double &s13,
        double &s21, double &s22, double &s23,
        double &s31, double &s32, double &s33,
        // output V
        double &v11, double &v12, double &v13,
        double &v21, double &v22, double &v23,
        double &v31, double &v32, double &v33)
{
    // normal equations matrix
    double ATA11, ATA12, ATA13;
    double ATA21, ATA22, ATA23;
    double ATA31, ATA32, ATA33;

    multAtB(a11,a12,a13,a21,a22,a23,a31,a32,a33,
            a11,a12,a13,a21,a22,a23,a31,a32,a33,
            ATA11,ATA12,ATA13,ATA21,ATA22,ATA23,ATA31,ATA32,ATA33);

    // symmetric eigenalysis
    double qV[4];
    jacobiEigenanlysis( ATA11,ATA21,ATA22, ATA31,ATA32,ATA33,qV);
    quatToMat3(qV,v11,v12,v13,v21,v22,v23,v31,v32,v33);

    double b11, b12, b13;
    double b21, b22, b23;
    double b31, b32, b33;
    multAB(a11,a12,a13,a21,a22,a23,a31,a32,a33,
           v11,v12,v13,v21,v22,v23,v31,v32,v33,
           b11, b12, b13, b21, b22, b23, b31, b32, b33);

    // sort singular values and find V
    sortSingularValues(b11, b12, b13, b21, b22, b23, b31, b32, b33,
                       v11,v12,v13,v21,v22,v23,v31,v32,v33);

    // QR decomposition
    QRDecomposition(b11, b12, b13, b21, b22, b23, b31, b32, b33,
                    u11, u12, u13, u21, u22, u23, u31, u32, u33,
                    s11, s12, s13, s21, s22, s23, s31, s32, s33
    );
}

/// polar decomposition can be reconstructed trivially from SVD result
/// A = UP
__host__ __device__ __forceinline__
void pd(double a11, double a12, double a13,
        double a21, double a22, double a23,
        double a31, double a32, double a33,
        // output U
        double &u11, double &u12, double &u13,
        double &u21, double &u22, double &u23,
        double &u31, double &u32, double &u33,
        // output P
        double &p11, double &p12, double &p13,
        double &p21, double &p22, double &p23,
        double &p31, double &p32, double &p33)
{
    double w11, w12, w13, w21, w22, w23, w31, w32, w33;
    double s11, s12, s13, s21, s22, s23, s31, s32, s33;
    double v11, v12, v13, v21, v22, v23, v31, v32, v33;

    svd(a11, a12, a13, a21, a22, a23, a31, a32, a33,
        w11, w12, w13, w21, w22, w23, w31, w32, w33,
        s11, s12, s13, s21, s22, s23, s31, s32, s33,
        v11, v12, v13, v21, v22, v23, v31, v32, v33);

    // P = VSV'
    double t11, t12, t13, t21, t22, t23, t31, t32, t33;
    multAB(v11, v12, v13, v21, v22, v23, v31, v32, v33,
           s11, s12, s13, s21, s22, s23, s31, s32, s33,
           t11, t12, t13, t21, t22, t23, t31, t32, t33);

    multAB(t11, t12, t13, t21, t22, t23, t31, t32, t33,
           v11, v21, v31, v12, v22, v32, v13, v23, v33,
           p11, p12, p13, p21, p22, p23, p31, p32, p33);

    // U = WV'
    multAB(w11, w12, w13, w21, w22, w23, w31, w32, w33,
           v11, v21, v31, v12, v22, v32, v13, v23, v33,
           u11, u12, u13, u21, u22, u23, u31, u32, u33);
}

__global__ void timestepSimulation(double d_dt, int d_n, double* d_v, double* d_a, double* d_x, double* d_p)
{
    for (int i = 0; i < d_n; i++)
    {
        double vexplicit[3] = {d_v[3*i + 0] + d_dt*d_a[3*i + 0],
                           d_v[3*i + 1] + d_dt*d_a[3*i + 1],
                           d_v[3*i + 2] + d_dt*d_a[3*i + 2]};

        d_p[3*i+0] = d_x[3*i+0] + d_dt * vexplicit[0];
        d_p[3*i+1] = d_x[3*i+1] + d_dt * vexplicit[1];
        d_p[3*i+2] = d_x[3*i+2] + d_dt * vexplicit[2];
    }
}

__global__ void fixedPointConstraint(int d_n, int* d_i, double* d_p, double* d_p0)
{


    int N = blockIdx.x * blockDim.x + threadIdx.x;

    if (N < d_n)
    {
        int v1 = d_i[N + 0];
        d_p[3 * v1 + 0] = d_p0[3 * v1 + 0];
        d_p[3 * v1 + 1] = d_p0[3 * v1 + 1];
        d_p[3 * v1 + 2] = d_p0[3 * v1 + 2];
    }

}

__global__ void neohookeanConstraint(double dt, int d_n, double d_mu, double d_lambda, double d_alpha, double d_beta,
                                     int* d_i, double* d_p, double* d_p0, double* d_DmInv, double* d_v0,
                                     double* m, double* lagrange)
{


    for (int i = 0; i < d_n; i++)
    {
        //tested
        int v1 = d_i[4*i + 0];
        int v2 = d_i[4*i + 1];
        int v3 = d_i[4*i + 2];
        int v4 = d_i[4*i + 3];

        //printf("%i\n",i);
        //printf("%i,%i,%i,%i\n",v1,v2,v3,v4);

        // tested
        double V0 = d_v0[i];

        //tested
        Eigen::Matrix3d DmInv;
        DmInv << d_DmInv[9*i+0],d_DmInv[9*i+1],d_DmInv[9*i+2],
                d_DmInv[9*i+3], d_DmInv[9*i+4], d_DmInv[9*i+5],
                d_DmInv[9*i+6], d_DmInv[9*i+7], d_DmInv[9*i+8];

        //printf("%f,%f,%f,%f,%f,%f,%f,%f,%f\n",d_DmInv[9*i+0],d_DmInv[9*i+1],d_DmInv[9*i+2],
        //       d_DmInv[9*i+3], d_DmInv[9*i+4], d_DmInv[9*i+5],
        //       d_DmInv[9*i+6], d_DmInv[9*i+7], d_DmInv[9*i+8]);

        Eigen::Vector3d p1 = {d_p[3*v1+0],d_p[3*v1+1], d_p[3*v1+2]};
        Eigen::Vector3d p2 = {d_p[3*v2+0],d_p[3*v2+1], d_p[3*v2+2]};
        Eigen::Vector3d p3 = {d_p[3*v3+0],d_p[3*v3+1], d_p[3*v3+2]};
        Eigen::Vector3d p4 = {d_p[3*v4+0],d_p[3*v4+1], d_p[3*v4+2]};

        //printf("%i,%i,%i\n",3*v1+0,3*v1+1, 3*v1+2);
        //printf("%i,%i,%i\n",3*v2+0,3*v2+1, 3*v2+2);
        //printf("%i,%i,%i\n",3*v3+0,3*v3+1, 3*v3+2);
        //printf("%i,%i,%i\n",3*v4+0,3*v4+1, 3*v4+2);

        //printf("%f,%f,%f\n",d_p[3*v1+0],d_p[3*v1+1], d_p[3*v1+2]);
        //printf("%f,%f,%f\n",d_p[3*v2+0],d_p[3*v2+1], d_p[3*v2+2]);
        //printf("%f,%f,%f\n",d_p[3*v3+0],d_p[3*v3+1], d_p[3*v3+2]);
        //printf("%f,%f,%f\n",d_p[3*v4+0],d_p[3*v4+1], d_p[3*v4+2]);


        Eigen::Vector3d p1_0 = {d_p0[3*v1+0],d_p0[3*v1+1], d_p0[3*v1+2]};
        Eigen::Vector3d p2_0 = {d_p0[3*v2+0],d_p0[3*v2+1], d_p0[3*v2+2]};
        Eigen::Vector3d p3_0 = {d_p0[3*v3+0],d_p0[3*v3+1], d_p0[3*v3+2]};
        Eigen::Vector3d p4_0 = {d_p0[3*v4+0],d_p0[3*v4+1], d_p0[3*v4+2]};

        p1 = p1_0;
        p2 = p2_0;
        p3 = p3_0;
        p4 = p4_0;

        //printf("%f,%f,%f\n",p1(0),p1(1),p1(2));
        //printf("%f,%f,%f\n",p2(0),p2(1),p2(2));
        //printf("%f,%f,%f\n",p3(0),p3(1),p3(2));
        //printf("%f,%f,%f\n",p4(0),p4(1),p4(2));

        //printf("%f,%f,%f\n",d_p0[3*v1+0],d_p0[3*v1+1], d_p0[3*v1+2]);
        //printf("%f,%f,%f\n",d_p0[3*v2+0],d_p0[3*v2+1], d_p0[3*v2+2]);
        //printf("%f,%f,%f\n",d_p0[3*v3+0],d_p0[3*v3+1], d_p0[3*v3+2]);
        //printf("%f,%f,%f\n",d_p0[3*v4+0],d_p0[3*v4+1], d_p0[3*v4+2]);

        double w1 = 1. / m[v1];
        double w2 = 1. / m[v2];
        double w3 = 1. / m[v3];
        double w4 = 1. / m[v4];

        Eigen::Matrix3d Ds;
        Ds.col(0) = (p1 - p4).transpose();
        Ds.col(1) = (p2 - p4).transpose();
        Ds.col(2) = (p3 - p4).transpose();

        double Vsigned = (1. / 6.) * det(Ds);

        bool is_V_positive = Vsigned >= 0 ;

        bool is_V0_positive = V0 >= 0;
        bool is_tet_inverted = (is_V_positive && !is_V0_positive) || (!is_V_positive && is_V0_positive);

        Eigen::Matrix3d F = mul(Ds, DmInv);
        // printf("%f,%f,%f,",F(0,0),F(1,1),F(2,2));
        Eigen::Matrix3d I = Eigen::Matrix3d::Identity();

        double epsilon = 1e-20;
        Eigen::Matrix3d Piola;
        double psi{};

        // TODO: Implement correct inversion handling described in
        // Irving, Geoffrey, Joseph Teran, and Ronald Fedkiw. "Invertible finite elements for robust
        // simulation of large deformation." Proceedings of the 2004 ACM SIGGRAPH/Eurographics symposium
        // on Computer animation. 2004.
        if (is_tet_inverted)
        {

            // TODO: Get singular values
        double U[3][3]; double S[3][3]; double V[3][3];
        svd(F(0,0),F(0,1),F(0,2),
             F(1,0),F(1,1),F(1,2),
             F(2,0),F(2,1),F(2,2),
             U[0][0], U[0][1], U[0][2],
             U[1][0], U[1][1], U[1][2],
             U[2][0], U[2][1], U[2][2],
             S[0][0], S[0][1], S[0][2],
             S[1][0], S[1][1], S[1][2],
             S[2][0], S[2][1], S[2][2],
             V[0][0], V[0][1], V[0][2],
             V[1][0], V[1][1], V[1][2],
             V[2][0], V[2][1], V[2][2]);

        Eigen::Vector3d Fsigma = {S[0][0],S[1][1],S[2][2]};
        Eigen::Matrix3d Fhat;
        Fhat.setZero();
        Fhat(0, 0) = Fsigma(0);
        Fhat(1, 1) = Fsigma(1);
        Fhat(2, 2) = Fsigma(2);

        Eigen::Matrix3d Umat;
        Umat << U[0][0],U[0][1],U[0][2],
                U[1][0],U[1][1],U[1][2],
                U[2][0],U[2][1],U[2][2];

        Eigen::Matrix3d Vmat;
        Vmat << V[0][0],V[0][1],V[0][2],
                V[1][0],V[1][1],V[1][2],
                V[2][0],V[2][1],V[2][2];

        auto smallest_element_idx = 0;
        if (Fsigma(0) < Fsigma(1) && Fsigma(0) < Fsigma(2))
            smallest_element_idx = 0;
        if (Fsigma(1) < Fsigma(0) && Fsigma(1) < Fsigma(2))
            smallest_element_idx = 1;
        if (Fsigma(2) < Fsigma(0) && Fsigma(2) < Fsigma(1))
            smallest_element_idx = 2;

        Fhat(smallest_element_idx, smallest_element_idx) =
                -Fhat(smallest_element_idx, smallest_element_idx);
        Umat.col(smallest_element_idx) = -Umat.col(smallest_element_idx);

        // stress reaches maximum at 58% compression
        double constexpr min_singular_value = 0.577;
        Fhat(0, 0)                               = min(Fhat(0, 0), min_singular_value);
        Fhat(1, 1)                               = min(Fhat(1, 1), min_singular_value);
        Fhat(2, 2)                               = min(Fhat(2, 2), min_singular_value);

        Eigen::Matrix3d const Fprime = Umat * Fhat * Vmat.transpose();
        Eigen::Matrix3d const F2     = Fprime.transpose() * Fprime;
        Eigen::Matrix3d const Finv   = inv(Fprime);
        Eigen::Matrix3d const FinvT  = Finv.transpose();
        double const I1         = F2.trace();
        double const J          = det(Fprime);

        double const logJ = std::log(J);
        // psi(I1, J) = (mu/2)*(I1 - 3) - mu*log(J) + (lambda/2)*log^2(J)
        psi = static_cast<double>(0.5) * d_mu * (I1 - static_cast<double>(3.)) -
              d_mu * logJ + static_cast<double>(0.5) * d_lambda * logJ * logJ;
        // P(F) = mu*(F - mu*F^-T) + lambda*log(J)*F^-T
        Piola = d_mu * (Fprime - d_mu * FinvT) + d_lambda * logJ * FinvT;
    }
    else
    {
        Eigen::Matrix3d const F2    = mul(F.transpose() ,F);
        Eigen::Matrix3d const Finv  = inv(F);
        Eigen::Matrix3d const FinvT = Finv.transpose();
        double const I1        = F2.trace();

        double const J         = det(F);

        double const logJ = std::log(J);
        // psi(I1, J) = (mu/2)*(I1 - 3) - mu*log(J) + (lambda/2)*log^2(J)
        psi = static_cast<double>(0.5) * d_mu * (I1 - static_cast<double>(3.)) -
              d_mu * logJ + static_cast<double>(0.5) * d_lambda * logJ * logJ;
        // P(F) = mu*(F - mu*F^-T) + lambda*log(J)*F^-T
        Piola = d_mu * (F - d_mu * FinvT) + d_lambda * logJ * FinvT;
        // printf("%f,%f,%f,",Piola(0,0),Piola(1,1),Piola(2,2));
    }


        // H is the negative gradient of the elastic potential
        V0 = std::abs(V0);
        Eigen::Matrix3d H = -V0 * mul(Piola , DmInv.transpose());
        Eigen::Vector3d f1 = H.col(0);
        Eigen::Vector3d f2 = H.col(1);
        Eigen::Vector3d f3 = H.col(2);
        Eigen::Vector3d f4 = -(f1 + f2 + f3);

        double weighted_sum_of_gradients =
                w1 * f1.squaredNorm() +
                w2 * f2.squaredNorm() +
                w3 * f3.squaredNorm() +
                w4 * f4.squaredNorm();

        if (weighted_sum_of_gradients < epsilon)
            continue;

        double C = V0 * psi;
        double alpha_tilde = d_alpha / (dt * dt);
        double beta_tilde  = d_beta * dt * dt;
        double gamma       = alpha_tilde * beta_tilde * dt;

        double delta_lagrange =
                -(C + alpha_tilde * lagrange[i]) / (weighted_sum_of_gradients + alpha_tilde);

        lagrange[i] += delta_lagrange;

        // because f = - grad(potential), then grad(potential) = -f and thus grad(C) = -f
        d_p[3*v1 + 0] += w1 * -f1(0) * delta_lagrange;
        d_p[3*v1 + 1] += w1 * -f1(1) * delta_lagrange;
        d_p[3*v1 + 2] += w1 * -f1(2) * delta_lagrange;

        d_p[3*v2 + 0] += w2 * -f2(0) * delta_lagrange;
        d_p[3*v2 + 1] += w2 * -f2(1) * delta_lagrange;
        d_p[3*v2 + 2] += w2 * -f2(2) * delta_lagrange;

        d_p[3*v3 + 0] += w3 * -f3(0) * delta_lagrange;
        d_p[3*v3 + 1] += w3 * -f3(1) * delta_lagrange;
        d_p[3*v3 + 2] += w3 * -f3(2) * delta_lagrange;

        d_p[3*v4 + 0] += w4 * -f4(0) * delta_lagrange;
        d_p[3*v4 + 1] += w4 * -f4(1) * delta_lagrange;
        d_p[3*v4 + 2] += w4 * -f4(2) * delta_lagrange;

    }
}

__global__ void volumeConstraint(double dt, int d_n, double d_alpha, double d_beta,
                                 int* d_i, double* d_p, double* d_p0, double* d_v0,
                                 double* m, double* lagrange, int* d_graph, int color, int maxelem)
{

    int globalIdx = blockIdx.x * blockDim.x + threadIdx.x;

    if (globalIdx < maxelem) {

        int i = d_graph[color * maxelem + globalIdx];

        if (i != -1) {

            int v1 = d_i[4 * i + 0];
            int v2 = d_i[4 * i + 1];
            int v3 = d_i[4 * i + 2];
            int v4 = d_i[4 * i + 3];

            double V0 = d_v0[i];

            Eigen::Vector3d p1 = {d_p[3 * v1 + 0], d_p[3 * v1 + 1], d_p[3 * v1 + 2]};
            Eigen::Vector3d p2 = {d_p[3 * v2 + 0], d_p[3 * v2 + 1], d_p[3 * v2 + 2]};
            Eigen::Vector3d p3 = {d_p[3 * v3 + 0], d_p[3 * v3 + 1], d_p[3 * v3 + 2]};
            Eigen::Vector3d p4 = {d_p[3 * v4 + 0], d_p[3 * v4 + 1], d_p[3 * v4 + 2]};


            double w1 = 1. / m[v1];
            double w2 = 1. / m[v2];
            double w3 = 1. / m[v3];
            double w4 = 1. / m[v4];

            double v = (1. / 6.) * abs((p2 - p1).cross(p3 - p1).dot(p4 - p1));

            auto const C = v - V0;

            Eigen::RowVector3d const grad1 = (1. / 6.) * (p2 - p3).cross(p4 - p3);
            Eigen::RowVector3d const grad2 = (1. / 6.) * (p3 - p1).cross(p4 - p1);
            Eigen::RowVector3d const grad3 = (1. / 6.) * (p1 - p2).cross(p4 - p2);
            Eigen::RowVector3d const grad4 = (1. / 6.) * (p2 - p1).cross(p3 - p1);

            auto const weighted_sum_of_gradients = w1 * grad1.squaredNorm() + w2 * grad2.squaredNorm() +
                                                   w3 * grad3.squaredNorm() + w4 * grad4.squaredNorm();


            double const alpha_tilde = d_alpha / (dt * dt);
            double const beta_tilde = d_beta * dt * dt;
            double const gamma = (alpha_tilde * beta_tilde) / dt;

            if (abs(alpha_tilde + weighted_sum_of_gradients) < 1e-6)
                return;

            double const delta_lagrange =
                    -(C + alpha_tilde * lagrange[i]) / (weighted_sum_of_gradients + alpha_tilde);

            lagrange[i] += delta_lagrange;

            d_p[3 * v1 + 0] += w1 * grad1(0) * delta_lagrange;
            d_p[3 * v1 + 1] += w1 * grad1(1) * delta_lagrange;
            d_p[3 * v1 + 2] += w1 * grad1(2) * delta_lagrange;

            d_p[3 * v2 + 0] += w2 * grad2(0) * delta_lagrange;
            d_p[3 * v2 + 1] += w2 * grad2(1) * delta_lagrange;
            d_p[3 * v2 + 2] += w2 * grad2(2) * delta_lagrange;

            d_p[3 * v3 + 0] += w3 * grad3(0) * delta_lagrange;
            d_p[3 * v3 + 1] += w3 * grad3(1) * delta_lagrange;
            d_p[3 * v3 + 2] += w3 * grad3(2) * delta_lagrange;

            d_p[3 * v4 + 0] += w4 * grad4(0) * delta_lagrange;
            d_p[3 * v4 + 1] += w4 * grad4(1) * delta_lagrange;
            d_p[3 * v4 + 2] += w4 * grad4(2) * delta_lagrange;

            lagrange[i] = 0;
        }
    }

}

__global__ void edgeConstraint(double dt, int d_n, double d_alpha, double d_beta,
                                  int* d_i, double* d_p, double* d_p0, double* d_e0,
                                  double* m, double* lagrange, int* d_graph, int color, int maxelem)
{

    int globalIdx = blockIdx.x * blockDim.x + threadIdx.x;

    if (globalIdx < maxelem) {

        int i = d_graph[color * maxelem + globalIdx];
        //printf("%i \n",i);
        if (i != -1) {

            int v1 = d_i[2 * i + 0];
            int v2 = d_i[2 * i + 1];

            double e0 = d_e0[i];


            Eigen::Vector3d p1 = {d_p[3 * v1 + 0], d_p[3 * v1 + 1], d_p[3 * v1 + 2]};
            Eigen::Vector3d p2 = {d_p[3 * v2 + 0], d_p[3 * v2 + 1], d_p[3 * v2 + 2]};

            double w1 = 1. / m[v1];
            double w2 = 1. / m[v2];

            double e = (p2 - p1).norm();

            auto const C = e - e0;

            // <n,n> = 1 and <-n,-n> = 1

            auto const grad1 = p1 - p2;
            auto const grad2 = p2 - p1;

            auto const weighted_sum_of_gradients = w1 + w2;
            auto const alpha_tilde = d_alpha / (dt * dt);

            if (abs(alpha_tilde + weighted_sum_of_gradients) < 1e-6)
                return;

            double const delta_lagrange =
                    -(C + alpha_tilde * lagrange[i]) / (weighted_sum_of_gradients + alpha_tilde);

            lagrange[i] += delta_lagrange;

            d_p[3 * v1 + 0] += w1 * grad1(0) * delta_lagrange;
            d_p[3 * v1 + 1] += w1 * grad1(1) * delta_lagrange;
            d_p[3 * v1 + 2] += w1 * grad1(2) * delta_lagrange;

            d_p[3 * v2 + 0] += w2 * grad2(0) * delta_lagrange;
            d_p[3 * v2 + 1] += w2 * grad2(1) * delta_lagrange;
            d_p[3 * v2 + 2] += w2 * grad2(2) * delta_lagrange;

            lagrange[i] = 0;
        }

    }

}

__global__ void explicitEuler(Eigen::Vector3i d_tri, Eigen::Vector3d d_f, int d_n, double d_dt,
                              double* d_p, double* d_x, double* d_pdot, double* d_m )
{
    int globalIdx = blockIdx.x * blockDim.x + threadIdx.x;


    if (globalIdx < d_n)
    {
        double a[3] = {0.,0.,0.};
        int d_i0 = d_tri[0];
        int d_i1 = d_tri[1];
        int d_i2 = d_tri[2];

        // apply force if in index
        if (globalIdx == d_i0)
        {
            a[0] += d_f[0] / d_m[d_i0];
            a[1] += d_f[1] / d_m[d_i0];
            a[2] += d_f[2] / d_m[d_i0];
        }
        if (globalIdx == d_i1)
        {
            a[0] += d_f[0] / d_m[d_i1];
            a[1] += d_f[1] / d_m[d_i1];
            a[2] += d_f[2] / d_m[d_i1];
        }
        if ( globalIdx == d_i2)
        {
            a[0] += d_f[0] / d_m[d_i2];
            a[1] += d_f[1] / d_m[d_i2];
            a[2] += d_f[2] / d_m[d_i2];
        }

        // apply gravity
        a[2] -= 9.81;

        // explicit euler
        double pdot_exp[3] = {0.,0.,0.};
        pdot_exp[0] = d_pdot[3*globalIdx + 0] + d_dt * a[0];
        pdot_exp[1] = d_pdot[3*globalIdx + 1] + d_dt * a[1];
        pdot_exp[2] = d_pdot[3*globalIdx + 2] + d_dt * a[2];

        // trivial damping times 0.99
        d_p[3*globalIdx + 0] = d_x[3*globalIdx + 0] + 0.99* d_dt * pdot_exp[0];
        d_p[3*globalIdx + 1] = d_x[3*globalIdx + 1] + 0.99* d_dt * pdot_exp[1];
        d_p[3*globalIdx + 2] = d_x[3*globalIdx + 2] + 0.99* d_dt * pdot_exp[2];
    }

}

__global__ void updateSolution(double d_dt, int d_n, double* d_xdot, double* d_x, double* d_xlast, double* d_p)
{
    int globalIdx = blockIdx.x * blockDim.x + threadIdx.x;

    if (globalIdx < d_n)
    {
        // set the last
        d_xlast[3*globalIdx + 0] = d_x[3*globalIdx + 0];
        d_xlast[3*globalIdx + 1] = d_x[3*globalIdx + 1];
        d_xlast[3*globalIdx + 2] = d_x[3*globalIdx + 2];

        d_xdot[3 * globalIdx + 0] = (d_p[3 * globalIdx + 0] - d_x[3 * globalIdx + 0]) / d_dt;
        d_xdot[3 * globalIdx + 1] = (d_p[3 * globalIdx + 1] - d_x[3 * globalIdx + 1]) / d_dt;
        d_xdot[3 * globalIdx + 2] = (d_p[3 * globalIdx + 2] - d_x[3 * globalIdx + 2]) / d_dt;

        d_x[3 * globalIdx + 0] = d_p[3 * globalIdx + 0];
        d_x[3 * globalIdx + 1] = d_p[3 * globalIdx + 1];
        d_x[3 * globalIdx + 2] = d_p[3 * globalIdx + 2];
    }

}

__global__ void computeNormals(const double* d_x, double* d_N0, double* d_N , int* d_F, int nfaces)
{

    int globalIdx = blockIdx.x * blockDim.x + threadIdx.x;

    if (globalIdx < nfaces)
    {
        d_N0[3*globalIdx + 0] = d_N[3*globalIdx + 0];
        d_N0[3*globalIdx + 1] = d_N[3*globalIdx + 1];
        d_N0[3*globalIdx + 2] = d_N[3*globalIdx + 2];

        int i1 = d_F[3*globalIdx + 0];
        int i2 = d_F[3*globalIdx + 1];
        int i3 = d_F[3*globalIdx + 2];

        Eigen::Vector3d A(d_x[3*i1 + 0], d_x[3*i1 + 1], d_x[3*i1 + 2]);
        Eigen::Vector3d B(d_x[3*i2 + 0], d_x[3*i2 + 1], d_x[3*i2 + 2]);
        Eigen::Vector3d C(d_x[3*i3 + 0], d_x[3*i3 + 1], d_x[3*i3 + 2]);

        Eigen::Vector3d N = (A - B).cross(B - C);

        d_N[3*globalIdx + 0] = N(0);
        d_N[3*globalIdx + 1] = N(1);
        d_N[3*globalIdx + 2] = N(2);
    }

}


/*
__global__ void transform2objecttype(object_type* tris, int* f, double* start, double* end, int ntris)
{
    //int globalIdx = blockIdx.x * blockDim.x + threadIdx.x;

    for (int globalIdx = 0; globalIdx < ntris; globalIdx++)//if (globalIdx < ntris)
    {
        for (int i = 0; i < 3; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                int idx = f[3 * globalIdx + i];

                tris[globalIdx].start[i][j] = start[3 * idx + j];
                tris[globalIdx].end[i][j] = end[3 * idx + j];
            }
        }
    }
}
 */