#pragma once

#include "parameter.hpp"
#include "variable.hpp"

#include <Eigen/Core>

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

    class Affine;
    class Constraint;
    class Scalar;
    class OptimizationProblem;
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
        friend SOCPWrapperBase;
        friend QPWrapperBase;

        explicit operator double() const;

    private:
        Affine affine;
        std::vector<Product> products;
        bool norm = false;

        friend std::ostream &operator<<(std::ostream &os, const Scalar &scalar);

        /**
         * Useful to call .norm() on a matrix.
         * 
         * Possible when only squared expressions are present.
         * 
         */
        friend Scalar sqrt(const Scalar &scalar);
        friend Scalar square(const Scalar &scalar);
        friend Scalar abs2(const Scalar &scalar);

        // friend Parameter::operator Scalar() const;
        friend Variable::operator Scalar() const;

        friend Constraint equalTo(const Scalar &lhs, const Scalar &rhs);
        friend Constraint lessThan(const Scalar &lhs, const Scalar &rhs);
        friend Constraint box(const Scalar &lower, const Scalar &middle, const Scalar &upper);
    };

    using MatrixX = Eigen::Matrix<cvx::Scalar, Eigen::Dynamic, Eigen::Dynamic>;
    using VectorX = Eigen::Matrix<cvx::Scalar, Eigen::Dynamic, 1>;

    /**
     * @brief Creates a single variable.
     * 
     * @warning Do not share variables between different problems
     * 
     * @param name The name of the variable
     * @return Scalar The variable
     */
    Scalar var(const std::string &name);

    /**
     * @brief Creates a vector of variables.
     * 
     * @warning Do not share variables between different problems
     * 
     * @param name The name of the variable
     * @param rows The number of elements in the vector
     * @return VectorX The vector of variables
     */
    VectorX var(const std::string &name,
                size_t rows);

    /**
     * @brief Creates a matrix of variables.
     * 
     * @warning Do not share variables between different problems
     * 
     * @param name The name of the variable
     * @param rows The number of rows of the matrix
     * @param cols The number of columns in the matrix
     * @return MatrixX The matrix of variables
     */
    MatrixX var(const std::string &name,
                size_t rows,
                size_t cols);

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
     * @details Internally store pointers to the original values.
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
     * @brief Creates a dynamic parameter from a dense Eigen type.
     * 
     * @details Internally store pointers to the original values.
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
        new_scalar.products = {Product(x.affine)};
        return new_scalar;
    }
    inline Scalar abs2(const Scalar &x)
    {
        return square(x);
    }

} // namespace cvx
