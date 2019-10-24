#include <cmath>
#include "Vector.h"

using namespace mat_vec;

Vector::Vector(size_t size, double value) : vecSize(size), vec(new double[size]) {
    for (size_t i = 0; i < size; i++) {
        vec[i] = value;
    }
}

Vector::~Vector() {
    delete [] vec;
    vecSize = 0;
}

Vector::Vector(const Vector &src) : vecSize(src.vecSize), vec(new double[src.vecSize]) {
    for (size_t i = 0; i < src.vecSize; i++) {
        vec[i] = src.vec[i];
    }
}

size_t Vector::size() const {
    return vecSize;
}

bool Vector::operator==(const Vector &rhs) const {
    if (vecSize != rhs.vecSize) {
        return false;
    }
    else {
        for (size_t i = 0; i < vecSize; i++) {
            if (vec[i] != rhs.vec[i]) {
                return false;
            }
        }
        return true;
    }
}

bool Vector::operator!=(const Vector &rhs) const {
    if (vecSize != rhs.vecSize) {
        return true;
    }
    else {
        for (size_t i = 0; i < vecSize; i++) {
            if (vec[i] == rhs.vec[i]) {
                return false;
            }
        }
        return true;
    }
}

double Vector::norm() const {
    double res = 0;
    for (size_t i = 0; i < vecSize; i++) {
        res += vec[i] * vec[i];
    }
    return sqrt(res);
}

Vector& Vector::operator=(const Vector &rhs) {
    delete [] vec;
    vecSize = rhs.vecSize;
    vec = new double[rhs.vecSize];
    for (size_t i = 0; i < rhs.vecSize; i++) {
        vec[i] = rhs.vec[i];
    }
    return *this;
}

Vector Vector::operator*(double k) const {
    Vector res(vecSize);
    for (size_t i = 0; i < vecSize; i++) {
        res.vec[i] = vec[i] * k;
    }
    return res;
}

Vector& Vector::operator*=(double k) {
    for (size_t i = 0; i < vecSize; i++) {
        vec[i] *= k;
    }
    return *this;
}

void Vector::normalize() {
    double countNorm = norm();
    for (size_t i = 0; i < vecSize; i++) {
        vec[i] /= countNorm;
    }
}

Vector Vector::normalized() const {
    Vector res(*this);
    double countNorm = norm();
    for (size_t i = 0; i < vecSize; i++) {
        res.vec[i] /= countNorm;
    }
    return res;
}

double Vector::operator[](size_t n) const {
    return vec[n];
}

double& Vector::operator[](size_t n) {
    return vec[n];
}

Vector Vector::operator/(double k) const {
    Vector res(vecSize);
    for (size_t i = 0; i < vecSize; i++) {
        res.vec[i] = vec[i] / k;
    }
    return res;
}

Vector& Vector::operator/=(double k) {
    for (size_t i = 0; i < vecSize; i++) {
        vec[i] /= k;
    }
    return *this;
}

Vector Vector::operator+(const Vector &rhs) const {
    Vector res(vecSize);
    for (size_t i = 0; i < vecSize; i++) {
        res.vec[i] = vec[i] + rhs.vec[i];
    }
    return res;
}

Vector& Vector::operator+=(const Vector &rhs) {
    for (size_t i = 0; i < vecSize; i++) {
        vec[i] = vec[i] + rhs.vec[i];
    }
    return *this;
}

Vector Vector::operator-(const Vector &rhs) const {
    Vector res(vecSize);
    for (size_t i = 0; i < vecSize; i++) {
        res.vec[i] = vec[i] - rhs.vec[i];
    }
    return res;
}

Vector& Vector::operator-=(const Vector &rhs) {
    for (size_t i = 0; i < vecSize; i++) {
        vec[i] = vec[i] - rhs.vec[i];
    }
    return *this;
}

Vector Vector::operator^(const mat_vec::Vector &rhs) const {
    Vector res(vecSize);
    for (size_t i = 0; i < vecSize; i++) {
        res.vec[i] = vec[i] * rhs.vec[i];
    }
    return res;
}

Vector& Vector::operator^=(const mat_vec::Vector &rhs) {
    for (size_t i = 0; i < vecSize; i++) {
        vec[i] = vec[i] * rhs.vec[i];
    }
    return *this;
}

double Vector::operator*(const Vector &rhs) const {
    double res = 0;
    for (size_t i = 0; i < vecSize; i++) {
        res += vec[i] * rhs.vec[i];
    }
    return res;
}

Vector mat_vec::operator*(double k, const Vector &v) {
    Vector res(v.size());
    for (size_t i = 0; i < v.size(); i++) {
        res[i] = v[i] * k;
    }
    return res;
}

Vector Vector::operator*(const Matrix &mat) const {
    std::pair<size_t, size_t> p = mat.shape();
    Vector res(p.second);
    for (size_t i = 0; i < p.second; ++i) {
        for(size_t j = 0; j < vecSize; j++) {
            res.vec[i] += vec[j] * mat.get(j, i);
        }
    }
    return res;
}

Vector& Vector::operator*=(const Matrix &mat) {
    std::pair<size_t, size_t> p = mat.shape();
    Vector res(p.second);
    for (size_t i = 0; i < p.second; ++i) {
        for(size_t j = 0; j < vecSize; j++) {
            res.vec[i] += vec[j] * mat.get(j, i);
        }
    }
    *this = res;
    return *this;
}