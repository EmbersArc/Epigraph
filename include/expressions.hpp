/**
 * @file expressions.hpp
 * 
 */

#pragma once

#include "parameter.hpp"
#include "variable.hpp"

#include <Eigen/Sparse>

namespace Eigen
{

    template <>
    struct NumTraits<cvx::Scalar>
        : NumTraits<double>
    {
        using Real = cvx::Scalar;
        using NonInteger = cvx::Scalar;
        using Nested = cvx::Scalar;

        enum
        {
            IsComplex = 0,
            IsInteger = 0,
            IsSigned = 1,
            RequireInitialization = 1,
            ReadCost = 10,
            AddCost = 200,
            MulCost = 200,
        };
    };

} // namespace Eigen

namespace cvx
{
    class Constraint;
    class OptimizationProblem;
    class Scalar;

    namespace internal
    {
        class Affine;
        class SOCPWrapperBase;
        class QPWrapperBase;

        class Term
        {
        public:
            Term();

            Parameter parameter;
            Variable variable;

            bool operator==(const Term &other) const;
            Term &operator*=(const Parameter &param);
            Term &operator/=(const Parameter &param);

            operator Affine() const;

            friend std::ostream &operator<<(std::ostream &os, const Term &term);
            double evaluate() const;
        };

        class Affine
        {
        public:
            bool operator==(const Affine &other) const;

            Parameter constant = Parameter(0.);
            std::vector<Term> terms;

            friend std::ostream &operator<<(std::ostream &os, const Affine &affine);
            double evaluate() const;
            Affine &operator+=(const Affine &other);
            Affine &operator-=(const Affine &other);
            Affine &operator*=(const Parameter &param);
            Affine &operator/=(const Parameter &param);
            // Affine operator+(const Affine &other) const;
            Affine operator-(const Affine &other) const;
            Affine operator-() const;

            void cleanUp();

            bool isZero() const;
            bool isConstant() const;
            bool isFirstOrder() const;
        };

        class Product
        {
        public:
            explicit Product(const Affine &term);
            Product(const Affine &lhs, const Affine &rhs);
            Affine &firstTerm();
            Affine &secondTerm();
            const Affine &firstTerm() const;
            const Affine &secondTerm() const;
            void toSquaredTerm();
            double evaluate() const;
            bool isSquare() const;

            bool operator==(const Product &other) const;

            friend std::ostream &operator<<(std::ostream &os, const Product &product);

        private:
            std::vector<Affine> factors;
        };

    } // namespace internal

    class Scalar
    {
    public:
        Scalar() = default;
        explicit Scalar(int x);
        Scalar(double x);
        explicit Scalar(double *x);

        Scalar &operator+=(const Scalar &other);
        Scalar &operator-=(const Scalar &other);
        Scalar &operator*=(const Scalar &other);
        Scalar &operator/=(const Scalar &other);
        Scalar operator-() const;
        friend Scalar operator+(const Scalar &lhs, const Scalar &rhs);
        friend Scalar operator-(const Scalar &lhs, const Scalar &rhs);
        friend Scalar operator*(const Scalar &lhs, const Scalar &rhs);
        friend Scalar operator/(const Scalar &lhs, const Scalar &rhs);

        bool operator==(const cvx::Scalar &other) const;

        double evaluate() const;
        size_t getOrder() const;
        bool isNorm() const;

        friend OptimizationProblem;
        friend internal::SOCPWrapperBase;
        friend internal::QPWrapperBase;

        explicit operator double() const;

        friend Scalar sqrt(const Scalar &scalar);
        friend Scalar square(const Scalar &scalar);
        friend Scalar abs2(const Scalar &scalar);

        // friend internal::Parameter::operator Scalar() const;
        friend internal::Variable::operator Scalar() const;

        friend Constraint equalTo(const Scalar &lhs, const Scalar &rhs);
        friend Constraint lessThan(const Scalar &lhs, const Scalar &rhs);
        friend Constraint greaterThan(const Scalar &lhs, const Scalar &rhs);
        friend Constraint box(const Scalar &lower, const Scalar &middle, const Scalar &upper);

    private:
        internal::Affine affine;
        std::vector<internal::Product> products;
        bool norm = false;

        friend std::ostream &operator<<(std::ostream &os, const Scalar &scalar);
    };

    using MatrixX = Eigen::Matrix<cvx::Scalar, Eigen::Dynamic, Eigen::Dynamic>;
    using VectorX = Eigen::Matrix<cvx::Scalar, Eigen::Dynamic, 1>;

    /**
     * @brief Creates a constant parameter.
     * 
     * @param p The value of the parameter
     * @return Scalar The constant parameter
     */
    Scalar par(double p);

    /**
     * @brief Creates a dynamic parameter.
     * 
     * @details Internally stores pointers to the original values.
     * 
     * @warning Do not delete the source before the parameter is no longer required.
     * 
     * @param p The value of the parameter
     * @return Scalar The constant parameter
     */
    Scalar dynpar(double &p);

    /**
     * @brief Creates a constant parameter from a dense Eigen type.
     * 
     * @tparam Derived 
     * @param m A dense Eigen type containing problem parameters
     * @return auto A dense Eigen type with cvx::Scalar as scalar type
     */
    template <typename Derived>
    inline auto par(const Eigen::MatrixBase<Derived> &m)
    {
        return m.template cast<Scalar>().eval();
    }

    /**
     * @brief Creates a constant parameter from a sparse Eigen type.
     * 
     * @tparam Derived 
     * @param m A sparse Eigen type containing problem parameters
     * @return auto A sparse Eigen type with cvx::Scalar as scalar type
     */
    template <typename Derived>
    inline auto par(const Eigen::SparseMatrixBase<Derived> &m)
    {
        return m.template cast<Scalar>().eval();
    }

    /**
     * @brief Creates a dynamic parameter from a dense Eigen type.
     * 
     * @details Internally stores pointers to the original values.
     * 
     * @warning Do not delete the source before the parameter is no longer required.
     * 
     * @tparam Derived 
     * @param m A dense Eigen type containing problem parameters
     * @return auto A dense Eigen type with cvx::Scalar as scalar type
     */
    template <typename Derived>
    auto dynpar(Eigen::MatrixBase<Derived> &m)
    {
        auto result = m.template cast<Scalar>().eval();

        for (int row = 0; row < m.rows(); row++)
        {
            for (int col = 0; col < m.cols(); col++)
            {
                result.coeffRef(row, col) = dynpar(m.coeffRef(row, col));
            }
        }

        return result;
    }

    /**
     * @brief Creates a dynamic parameter from a sparse Eigen type.
     * 
     * @details Internally stores pointers to the original values.
     * 
     * @warning Do not delete the source before the parameter is no longer required.
     * 
     * @tparam T 
     * @param m A sparse Eigen type containing problem parameters
     * @return auto A sparse Eigen type with cvx::Scalar as scalar type
     */
    template <typename T>
    auto dynpar(Eigen::SparseMatrix<T> &m)
    {
        auto result = m.template cast<Scalar>().eval();

        for (int k = 0; k < result.nonZeros(); k++)
        {
            result.valuePtr()[k] = dynpar(m.valuePtr()[k]);
        }

        return result;
    }

    /**
     * @brief Evaluates the scalar.
     * 
     * @param s The scalar to be evaluated
     * @return double The value of the scalar
     */
    double eval(const Scalar &s);

    /**
     * @brief Evaluates a dense Eigen type
     * 
     * @tparam Derived Has to be cvx::Scalar
     * @param m A dense Eigen type to be evaluated
     * @return auto The evaluated dense Eigen type
     */
    template <typename Derived>
    inline auto eval(const Eigen::MatrixBase<Derived> &m)
    {
        return m.template cast<double>();
    }

    /**
     * @brief Evaluates a sparse Eigen type
     * 
     * @tparam Derived Has to be cvx::Scalar
     * @param m A sparse Eigen type to be evaluated
     * @return auto The evaluated sparse Eigen type
     */
    template <typename Derived>
    inline auto eval(const Eigen::SparseMatrixBase<Derived> &m)
    {
        return m.template cast<double>();
    }

    inline const Scalar &conj(const Scalar &x) { return x; }
    inline const Scalar &real(const Scalar &x) { return x; }
    inline Scalar imag(const Scalar &) { return Scalar(0.); }
    inline Scalar square(const Scalar &x)
    {
        Scalar new_scalar;
        new_scalar.products = {internal::Product(x.affine)};
        return new_scalar;
    }
    inline Scalar abs2(const Scalar &x)
    {
        return square(x);
    }

} // namespace cvx
