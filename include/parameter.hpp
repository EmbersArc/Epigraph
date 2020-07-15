#pragma once

#include <memory>

namespace cvx
{

    enum class ParamOpcode
    {
        Add,
        Mul,
        Div,
        Sqrt,
    };

    enum class ParameterType
    {
        Constant,
        Pointer,
        Operation,
    };

    class ParameterSource
    {
    public:
        virtual double getValue() const = 0;
        virtual ParameterType getType() const = 0;
    };

    class ConstantSource final : public ParameterSource
    {
    public:
        explicit ConstantSource(double const_value);
        double getValue() const override;
        ParameterType getType() const override;
        bool operator==(const ConstantSource &other) const;

    private:
        double value;
    };

    class PointerSource final : public ParameterSource
    {
    public:
        explicit PointerSource(const double *value_ptr);
        double getValue() const override;
        ParameterType getType() const override;
        bool operator==(const PointerSource &other) const;

    private:
        const double *ptr;
    };

    class OperationSource final : public ParameterSource
    {
    public:
        OperationSource(ParamOpcode op,
                        std::shared_ptr<ParameterSource> p1,
                        std::shared_ptr<ParameterSource> p2);
        double getValue() const override;
        ParameterType getType() const override;
        bool operator==(const OperationSource &other) const;

    private:
        ParamOpcode op;
        std::shared_ptr<ParameterSource> p1;
        std::shared_ptr<ParameterSource> p2;
    };

    class Scalar;
    class Affine;

    class Parameter
    {
    public:
        Parameter();
        explicit Parameter(int const_value);
        explicit Parameter(double const_value);
        explicit Parameter(double *value_ptr);

        bool isZero() const;
        bool isOne() const;
        double getValue() const;

        bool operator==(const Parameter &other) const;

        Parameter operator+(const Parameter &other) const;
        Parameter operator-(const Parameter &other) const;
        Parameter operator-() const;
        Parameter operator*(const Parameter &other) const;
        Parameter operator/(const Parameter &other) const;
        Parameter &operator+=(const Parameter &other);
        Parameter &operator*=(const Parameter &other);
        Parameter &operator/=(const Parameter &other);
        explicit operator Affine() const;
        explicit operator Scalar() const;
        explicit operator double() const;

        friend Parameter sqrt(const Parameter &param);
        friend std::ostream &operator<<(std::ostream &os, const Parameter &parameter);

    private:
        std::shared_ptr<ParameterSource> source;
    };

} // namespace cvx
