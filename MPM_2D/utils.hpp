//
//  utils.hpp
//  MPM_2D
//
//  Created by LeeSangheon on 2024/02/17.
//

#ifndef utils_hpp
#define utils_hpp

#include <stdio.h>
#include <iostream>
#include <random>
#include <algorithm>

///
/// Simple
///

template <int I, typename T>
constexpr inline T pow(T a) noexcept {
  T ret(1);
  for (int i = 0; i < I; i++) ret *= a;
  return ret;
};

template <typename T>
inline T clamp(const T &a, const T &min, const T &max) noexcept {
  if (a < min) return min;
  if (a > max) return max;
  return a;
}

///
/// Vector
///
template <int dim, typename  T>
struct VectorNDBase {
    T d[dim];
    
    VectorNDBase(T t) {
        for (int i=0; i<dim; i++) { d[i] = t; }
    }
};

template <typename T>
struct VectorNDBase<1, T> {
    union {
        T d[1];
        struct {
            T x;
        };
    };
};

template <typename T>
struct VectorNDBase<2, T> {
    union {
        T d[2];
        struct {
            T x, y;
        };
    };
};

template <typename T>
struct VectorNDBase<3, T> {
    union {
        T d[3];
        struct {
            T x, y, z;
        };
    };
};

template <typename T>
struct VectorNDBase<4, T> {
    union {
        T d[4];
        struct {
            T x, y, z, w;
        };
    };
};

template <int dim, typename T>
struct VectorND : public VectorNDBase<dim, T> {
    using VectorNDBase<dim,T>::d;
    
    inline VectorND() {
        for (int i=0; i<dim; i++) {
            this->d[i] = T(0);
        }
    }
    
    static inline VectorND from_array(const T new_val[dim]) {
        VectorND ret;
        for (int i=0; i<dim; i++) {
            ret.d[i] = new_val[i];
        }
        return ret;
    }
    
    template <int dim_, typename T_>
    explicit inline VectorND(const VectorND<dim_, T_> &o) : VectorND() {
        for (int i=0; i<std::min(dim_, dim); i++) { d[i] = o[i]; }
    }
    
    explicit inline VectorND(const std::array<T, dim> &o) {
        for (int i=0; i<dim; i++) { d[i] = o[i]; }
    }
    
    inline VectorND(T x) {
        for (int i=0; i<dim; i++) { this->d[i] = x;
        }
    }
    explicit inline VectorND(T x, T y) {
        static_assert(dim == 2, "Vector dim must be 2");
        this->d[0]=x; this->d[1]=y;
    }
    explicit inline VectorND(T x, T y, T z) {
        static_assert(dim == 3, "Vector dim must be 3");
        this->d[0]=x; this->d[1]=y; this->d[2]=z;
    }
    explicit inline VectorND(T x, T y, T z, T w) {
        static_assert(dim >= 4, "Vector dim must be 4");
        this->d[0]=x; this->d[1]=y; this->d[2]=z; this->d[3]=w;
    }
    explicit inline VectorND(VectorND<dim-1, T> v, T x) {
        for (int i=0; i<dim-1; i++) this->d[i] = v.d[i];
        this->d[dim-1] = x;
    }
    
    template<typename F,
        std::enable_if_t<std::is_convertible<F, std::function<T(int)>>::value, int> = 0>
    explicit inline VectorND(const F &f) {
        for (int i=0; i<dim; i++) this->d[i] = f(i);
    }
    
    inline T &operator[](int i) { return this->d[i]; }
    inline const T &operator[](int i) const { return this->d[i]; }
    
    inline T &operator()(int i) { return this->d[i]; }
    inline const T &operator()(int i) const { return this->d[i]; }
    
    inline T dot(VectorND<dim, T> o) const  {
        T ret = T(0);
        for (int i = 0; i < dim; i++) ret += this->d[i] * o[i];
        return ret;
    }
    
    inline VectorND &operator=(const VectorND &o) {
        memcpy(this, &o, sizeof(*this));
        return *this;
    }
    
    inline auto map(T(f)(T)) const
        -> VectorND<dim, decltype(f(T(0)))> {
      VectorND<dim, decltype(f(T(0)))> ret;
      for (int i = 0; i < dim; i++)
        ret[i] = f(this->d[i]);
      return ret;
    }
    
    inline VectorND operator+(const VectorND &o) const {
        return VectorND([=](int i) { return this->d[i] + o[i]; });
    }
    inline VectorND operator-(const VectorND &o) const {
        return VectorND([=](int i) { return this->d[i] - o[i]; });
    }
    inline VectorND operator*(const VectorND &o) const {
        return VectorND([=](int i) { return this->d[i] * o[i]; });
    }
    inline VectorND operator/(const VectorND &o) const {
        return VectorND([=](int i) { return this->d[i] / o[i]; });
    }
    inline VectorND operator*(const T &o) const {
        return VectorND([=](int i) { return this->d[i] * o; });
    }
    inline VectorND operator/(const T &o) const {
        return VectorND([=](int i) { return this->d[i] / o; });
    }
    
    
    VectorND &operator+=(const VectorND &o) {
        (*this) = (*this) + o; return *this;
    }
    VectorND &operator-=(const VectorND &o) {
        (*this) = (*this) - o; return *this;
    }
    VectorND &operator*=(const VectorND &o) {
        (*this) = (*this) * o; return *this;
    }
    VectorND &operator/=(const VectorND &o) {
        (*this) = (*this) / o; return *this;
    }
    VectorND &operator*=(const T &o) {
        (*this) = (*this) * o; return *this;
    }
    VectorND &operator/=(const T &o) {
        (*this) = (*this) / o; return *this;
    }
    inline VectorND operator-() const {
        return VectorND([=](int i) { return -this->d[i]; });
    }
    
    inline bool operator==(const VectorND &o) const {
      for (int i = 0; i < dim; i++)
        if (this->d[i] != o[i]) return false;
      return true;
    }
    inline bool operator!=(const VectorND &o) const {
      for (int i = 0; i < dim; i++)
        if (this->d[i] != o[i]) return true;
      return false;
    }
    inline bool operator<(const VectorND &o) const {
      for (int i = 0; i < dim; i++)
        if (this->d[i] >= o[i]) return false;
      return true;
    }
    inline bool operator<=(const VectorND &o) const {
      for (int i = 0; i < dim; i++)
        if (this->d[i] > o[i]) return false;
      return true;
    }
    inline bool operator>(const VectorND &o) const {
      for (int i = 0; i < dim; i++)
        if (this->d[i] <= o[i]) return false;
      return true;
    }
    inline bool operator>=(const VectorND &o) const {
      for (int i = 0; i < dim; i++)
        if (this->d[i] < o[i]) return false;
      return true;
    }
    
    inline VectorND sqr() const {
        return VectorND([&](int i) { return d[i]*d[i]; });
    }
    
    template <typename G>
    inline VectorND<dim, G> cast() const {
      return VectorND<dim, G>(
          [this](int i) { return static_cast<G>(this->d[i]); });
    }
    
    void print() const {
      for (int i = 0; i < dim; i++) { std::cout << this->d[i] << " "; }
      std::cout << std::endl;
    }
    
    static VectorND randomVector() {
        VectorND ret;
        for (int i=0; i<dim; i++) { ret[i] = rand() / static_cast<float>(RAND_MAX); }
        return ret;
    }
};

template <int dim, typename T>
inline VectorND<dim, T> operator*(T a, const VectorND<dim, T> &v) {
    return v * a;
}

///
/// Matrix
///

template <int dim, typename T>
struct MatrixND {
    using Vec = VectorND<dim, T>;
    Vec d[dim];
    
    inline MatrixND() {
        for (int i=0; i<dim; i++) { d[i] = Vec{}; }
    }
    
    inline MatrixND(T v) : MatrixND() {
        for (int i=0; i<dim; i++) {
            d[i][i] = v;
        }
    }
    
    inline explicit MatrixND(Vec v0, Vec v1) {
      static_assert(dim == 2, "Matrix dim must be 2");
      this->d[0] = v0; this->d[1] = v1;
    }
    inline explicit MatrixND(Vec v0, Vec v1, Vec v2) {
      static_assert(dim == 3, "Matrix dim must be 3");
      this->d[0] = v0; this->d[1] = v1; this->d[2] = v2;
    }
    inline explicit MatrixND(Vec v0, Vec v1, Vec v2, Vec v3) {
      static_assert(dim == 4, "Matrix dim must be 4");
      this->d[0] = v0; this->d[1] = v1; this->d[2] = v2; this->d[3] = v3;
    }
    
    template <typename F,
        std::enable_if_t<std::is_convertible<F, std::function<VectorND<dim, T>(int)>>::value,
        int> = 0>
    inline explicit MatrixND(const F &f) {
        for (int i = 0; i < dim; i++) this->d[i] = f(i);
    }
    template <typename F,
        std::enable_if_t<std::is_convertible<F, std::function<VectorND<dim, T>(int)>>::value,
        int> = 0>
    inline MatrixND &set(const F &f) {
        for (int i = 0; i < dim; i++) this->d[i] = f(i);
        return *this;
    }
    
    inline MatrixND &operator=(const MatrixND &o) {
        for (int i=0; i<dim; i++) this->d[i] = o[i];
        return *this;
    }
    inline Vec &operator[](int i) {
      return d[i];
    }
    inline const Vec &operator[](int i) const {
      return d[i];
    }
    inline T &operator()(int i, int j) {
      return d[j][i];
    }
    inline const T &operator()(int i, int j) const {
      return d[j][i];
    }
    
    inline MatrixND operator*(const T &o) const {
        return MatrixND([=](int i) { return this->d[i]*o; });
    }
    inline Vec operator*(const Vec &o) const {
        Vec ret = d[0] * o[0];
        for (int i=1; i<dim; i++) { ret += d[i]*o[i]; }
        return ret;
    }
    inline MatrixND operator*(const MatrixND &o) const {
      return MatrixND([&](int i) { return (*this) * o[i]; });
    }
    inline MatrixND operator+(const MatrixND &o) const {
      return MatrixND([=](int i) { return this->d[i] + o[i]; });
    }
    inline MatrixND operator-(const MatrixND &o) const {
      return MatrixND([=](int i) { return this->d[i] - o[i]; });
    }
    inline MatrixND &operator+=(const MatrixND &o) {
      return this->set([&](int i) { return this->d[i] + o[i]; });
    }
    inline MatrixND &operator-=(const MatrixND &o) {
      return this->set([&](int i) { return this->d[i] - o[i]; });
    }
    inline MatrixND operator-() const {
      return MatrixND([=](int i) { return -this->d[i]; });
    }
    inline bool operator==(const MatrixND &o) const {
      for (int i = 0; i < dim; i++)
        for (int j = 0; j < dim; j++)
          if (d[i][j] != o[i][j]) return false;
      return true;
    }
    inline bool operator!=(const MatrixND &o) const {
      for (int i = 0; i < dim; i++)
        for (int j = 0; j < dim; j++)
          if (d[i][j] != o[i][j]) return true;
      return false;
    }
    
    inline MatrixND transposed() const {
      MatrixND ret;
      for (int i = 0; i < dim; i++) {
        for (int j = 0; j < dim; j++) {
          ret[i][j] = d[j][i];
        }
      }
      return ret;
    }
    
    template <typename G>
    inline MatrixND<dim, G> cast() const {
      return MatrixND<dim, G>(
          [=](int i) { return d[i].template cast<G>(); });
    }
    
    inline static MatrixND outer_product(Vec column, Vec row) {
        return MatrixND([&](int i) { return column * row[i]; });
    }
};

template <typename T>
inline T determinant(const MatrixND<2, T> &mat) {
    return mat[0][0] * mat[1][1] - mat[0][1] * mat[1][0];
}

template <typename T>
inline T determinant(const MatrixND<3, T> &mat) {
    return mat[0][0] * (mat[1][1] * mat[2][2] - mat[2][1] * mat[1][2]) -
         mat[1][0] * (mat[0][1] * mat[2][2] - mat[2][1] * mat[0][2]) +
         mat[2][0] * (mat[0][1] * mat[1][2] - mat[1][1] * mat[0][2]);
}

template <typename T>
T determinant(const MatrixND<4, T> &m) {
    T Coef00 = m[2][2] * m[3][3] - m[3][2] * m[2][3];
    T Coef02 = m[1][2] * m[3][3] - m[3][2] * m[1][3];
    T Coef03 = m[1][2] * m[2][3] - m[2][2] * m[1][3];

    T Coef04 = m[2][1] * m[3][3] - m[3][1] * m[2][3];
    T Coef06 = m[1][1] * m[3][3] - m[3][1] * m[1][3];
    T Coef07 = m[1][1] * m[2][3] - m[2][1] * m[1][3];

    T Coef08 = m[2][1] * m[3][2] - m[3][1] * m[2][2];
    T Coef10 = m[1][1] * m[3][2] - m[3][1] * m[1][2];
    T Coef11 = m[1][1] * m[2][2] - m[2][1] * m[1][2];

    T Coef12 = m[2][0] * m[3][3] - m[3][0] * m[2][3];
    T Coef14 = m[1][0] * m[3][3] - m[3][0] * m[1][3];
    T Coef15 = m[1][0] * m[2][3] - m[2][0] * m[1][3];

    T Coef16 = m[2][0] * m[3][2] - m[3][0] * m[2][2];
    T Coef18 = m[1][0] * m[3][2] - m[3][0] * m[1][2];
    T Coef19 = m[1][0] * m[2][2] - m[2][0] * m[1][2];

    T Coef20 = m[2][0] * m[3][1] - m[3][0] * m[2][1];
    T Coef22 = m[1][0] * m[3][1] - m[3][0] * m[1][1];
    T Coef23 = m[1][0] * m[2][1] - m[2][0] * m[1][1];

    using Vec = VectorND<4, T>;

    Vec Fac0(Coef00, Coef00, Coef02, Coef03);
    Vec Fac1(Coef04, Coef04, Coef06, Coef07);
    Vec Fac2(Coef08, Coef08, Coef10, Coef11);
    Vec Fac3(Coef12, Coef12, Coef14, Coef15);
    Vec Fac4(Coef16, Coef16, Coef18, Coef19);
    Vec Fac5(Coef20, Coef20, Coef22, Coef23);

    Vec Vec0(m[1][0], m[0][0], m[0][0], m[0][0]);
    Vec Vec1(m[1][1], m[0][1], m[0][1], m[0][1]);
    Vec Vec2(m[1][2], m[0][2], m[0][2], m[0][2]);
    Vec Vec3(m[1][3], m[0][3], m[0][3], m[0][3]);

    Vec Inv0(Vec1 * Fac0 - Vec2 * Fac1 + Vec3 * Fac2);
    Vec Inv1(Vec0 * Fac0 - Vec2 * Fac3 + Vec3 * Fac4);
    Vec Inv2(Vec0 * Fac1 - Vec1 * Fac3 + Vec3 * Fac5);
    Vec Inv3(Vec0 * Fac2 - Vec1 * Fac4 + Vec2 * Fac5);

    Vec SignA(+1, -1, +1, -1);
    Vec SignB(-1, +1, -1, +1);
    MatrixND<4, T> Inverse(Inv0 * SignA, Inv1 * SignB, Inv2 * SignA, Inv3 * SignB);

    Vec Row0(Inverse[0][0], Inverse[1][0], Inverse[2][0], Inverse[3][0]);

    Vec Dot0(m[0] * Row0);
    T Dot1 = (Dot0.x + Dot0.y) + (Dot0.z + Dot0.w);

    return Dot1;
}

template <typename T>
inline MatrixND<2, T> inversed(const MatrixND<2, T> &mat) {
    T det = determinant(mat);
    return static_cast<T>(1) / det *
         MatrixND<2, T>(VectorND<2, T>(mat[1][1], -mat[0][1]),
                             VectorND<2, T>(-mat[1][0], mat[0][0]));
}

template <typename T>
MatrixND<3, T> inversed(const MatrixND<3, T> &mat) {
    T det = determinant(mat);
    return T(1.0) / det *
         MatrixND<3, T>(
             VectorND<3, T>(mat[1][1] * mat[2][2] - mat[2][1] * mat[1][2],
                            mat[2][1] * mat[0][2] - mat[0][1] * mat[2][2],
                            mat[0][1] * mat[1][2] - mat[1][1] * mat[0][2]),
             VectorND<3, T>(mat[2][0] * mat[1][2] - mat[1][0] * mat[2][2],
                            mat[0][0] * mat[2][2] - mat[2][0] * mat[0][2],
                            mat[1][0] * mat[0][2] - mat[0][0] * mat[1][2]),
             VectorND<3, T>(mat[1][0] * mat[2][1] - mat[2][0] * mat[1][1],
                            mat[2][0] * mat[0][1] - mat[0][0] * mat[2][1],
                            mat[0][0] * mat[1][1] - mat[1][0] * mat[0][1]));
}

template <typename T>
MatrixND<4, T> inversed(const MatrixND<4, T> &m) {
    T Coef00 = m[2][2] * m[3][3] - m[3][2] * m[2][3];
    T Coef02 = m[1][2] * m[3][3] - m[3][2] * m[1][3];
    T Coef03 = m[1][2] * m[2][3] - m[2][2] * m[1][3];

    T Coef04 = m[2][1] * m[3][3] - m[3][1] * m[2][3];
    T Coef06 = m[1][1] * m[3][3] - m[3][1] * m[1][3];
    T Coef07 = m[1][1] * m[2][3] - m[2][1] * m[1][3];

    T Coef08 = m[2][1] * m[3][2] - m[3][1] * m[2][2];
    T Coef10 = m[1][1] * m[3][2] - m[3][1] * m[1][2];
    T Coef11 = m[1][1] * m[2][2] - m[2][1] * m[1][2];

    T Coef12 = m[2][0] * m[3][3] - m[3][0] * m[2][3];
    T Coef14 = m[1][0] * m[3][3] - m[3][0] * m[1][3];
    T Coef15 = m[1][0] * m[2][3] - m[2][0] * m[1][3];

    T Coef16 = m[2][0] * m[3][2] - m[3][0] * m[2][2];
    T Coef18 = m[1][0] * m[3][2] - m[3][0] * m[1][2];
    T Coef19 = m[1][0] * m[2][2] - m[2][0] * m[1][2];

    T Coef20 = m[2][0] * m[3][1] - m[3][0] * m[2][1];
    T Coef22 = m[1][0] * m[3][1] - m[3][0] * m[1][1];
    T Coef23 = m[1][0] * m[2][1] - m[2][0] * m[1][1];

    using Vec = VectorND<4, T>;

    Vec Fac0(Coef00, Coef00, Coef02, Coef03);
    Vec Fac1(Coef04, Coef04, Coef06, Coef07);
    Vec Fac2(Coef08, Coef08, Coef10, Coef11);
    Vec Fac3(Coef12, Coef12, Coef14, Coef15);
    Vec Fac4(Coef16, Coef16, Coef18, Coef19);
    Vec Fac5(Coef20, Coef20, Coef22, Coef23);

    Vec Vec0(m[1][0], m[0][0], m[0][0], m[0][0]);
    Vec Vec1(m[1][1], m[0][1], m[0][1], m[0][1]);
    Vec Vec2(m[1][2], m[0][2], m[0][2], m[0][2]);
    Vec Vec3(m[1][3], m[0][3], m[0][3], m[0][3]);

    Vec Inv0(Vec1 * Fac0 - Vec2 * Fac1 + Vec3 * Fac2);
    Vec Inv1(Vec0 * Fac0 - Vec2 * Fac3 + Vec3 * Fac4);
    Vec Inv2(Vec0 * Fac1 - Vec1 * Fac3 + Vec3 * Fac5);
    Vec Inv3(Vec0 * Fac2 - Vec1 * Fac4 + Vec2 * Fac5);

    Vec SignA(+1, -1, +1, -1);
    Vec SignB(-1, +1, -1, +1);
    MatrixND<4, T> Inverse(Inv0*SignA, Inv1*SignB, Inv2*SignA, Inv3*SignB);

    Vec Row0(Inverse[0][0], Inverse[1][0], Inverse[2][0], Inverse[3][0]);

    Vec Dot0(m[0] * Row0);
    T Dot1 = (Dot0.x + Dot0.y) + (Dot0.z + Dot0.w);

    T OneOverDeterminant = static_cast<T>(1) / Dot1;

    return Inverse * OneOverDeterminant;
}

template <int dim, typename T>
inline MatrixND<dim, T> operator*(const T a, const MatrixND<dim, T> &M) {
    return M * a;
}

template <typename T>
inline void polar_decomp(MatrixND<2, T> m, MatrixND<2,T> &R, MatrixND<2,T> &S) {
    auto x = m(0, 0) + m(1, 1);
    auto y = m(1, 0) - m(0, 1);
    auto scale = 1.0 / std::sqrt(x * x + y * y);
    auto c = x * scale, s = y * scale;
    R(0, 0) = c; R(0, 1) = -s; R(1, 0) = s; R(1, 1) = c;
    S = R.transposed() * m;
}

inline void svd(MatrixND<2, float> m, MatrixND<2, float> &U,
                MatrixND<2, float> &sig, MatrixND<2, float> &V) {
  MatrixND<2, float> S;
  polar_decomp(m, U, S);
  float c, s;
  if (std::abs(S(0, 1)) < 1e-6) {
    sig = S;
    c = 1;
    s = 0;
  } else {
    auto tao = 0.5 * (S(0, 0) - S(1, 1));
    auto w = std::sqrt(tao * tao + S(0, 1) * S(0, 1));
    auto t = tao > 0 ? S(0, 1) / (tao + w) : S(0, 1) / (tao - w);
    c = 1.0 / std::sqrt(t * t + 1);
    s = -t * c;
    sig(0, 0) = pow<2>(c) * S(0, 0) - 2 * c * s * S(0, 1) + pow<2>(s) * S(1, 1);
    sig(1, 1) = pow<2>(s) * S(0, 0) + 2 * c * s * S(0, 1) + pow<2>(c) * S(1, 1);
  }
  if (sig(0, 0) < sig(1, 1)) {
    std::swap(sig(0, 0), sig(1, 1));
    V(0, 0) = -s;
    V(0, 1) = -c;
    V(1, 0) = c;
    V(1, 1) = -s;
  } else {
    V(0, 0) = c;
    V(0, 1) = -s;
    V(1, 0) = s;
    V(1, 1) = c;
  }
  V = V.transposed();
  U = U * V;
}

#endif /* utils_hpp */
