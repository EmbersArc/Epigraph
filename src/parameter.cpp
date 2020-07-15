#include "parameter.hpp"

#include <sstream>
#include <cassert>
#include <cmath>

namespace cvx
{

    ConstantSource::ConstantSource(double const_value)
        : value(const_value) {}

    PointerSource::PointerSource(const double *value_ptr)
        : ptr(value_ptr) {}

    OperationSource::OperationSource(ParamOpcode op,
                                     std::shared_ptr<ParameterSource> p1,
                                     std::shared_ptr<ParameterSource> p2)
        : op(op), p1(p1), p2(p2) {}

    ParameterType ConstantSource::getType() const
    {
        return ParameterType::Constant;
    }
    ParameterType PointerSource::getType() const
    {
        return ParameterType::Pointer;
    }
    ParameterType OperationSource::getType() const
    {
        return ParameterType::Operation;
    }

    double ConstantSource::getValue() const
    {
        return value;
    }

    double PointerSource::getValue() const
    {
        return *ptr;
    }

    double OperationSource::getValue() const
    {
        switch (op)
        {
        case ParamOpcode::Add:
            return p1->getValue() +
                   p2->getValue();
        case ParamOpcode::Mul:
            return p1->getValue() *
                   p2->getValue();
        case ParamOpcode::Div:
            assert(p2->getValue() != 0.);
            return p1->getValue() /
                   p2->getValue();
        default: // ParamOpcode::Sqrt:
            assert(p1->getValue() >= 0.);
            return std::sqrt(p1->getValue());
        }
    }

    Parameter::Parameter()
        : source(std::make_shared<ConstantSource>(0.))
    {
    }

    Parameter::Parameter(int const_value)
        : source(std::make_shared<ConstantSource>(double(const_value)))
    {
    }

    Parameter::Parameter(double const_value)
        : source(std::make_shared<ConstantSource>(const_value))
    {
    }

    Parameter::Parameter(double *value_ptr)
        : source(std::make_shared<PointerSource>(value_ptr))
    {
    }

    double Parameter::getValue() const
    {
        return source->getValue();
    }

    bool compare_sources(std::shared_ptr<ParameterSource> p1, std::shared_ptr<ParameterSource> p2)
    {
        if (p1 == p2)
        {
            return true;
        }
        else if (p1->getType() == p2->getType())
        {
            if (p1->getType() == ParameterType::Constant)
            {
                return *std::dynamic_pointer_cast<ConstantSource>(p1) ==
                       *std::dynamic_pointer_cast<ConstantSource>(p2);
            }
            else if (p1->getType() == ParameterType::Pointer)
            {
                return *std::dynamic_pointer_cast<PointerSource>(p1) ==
                       *std::dynamic_pointer_cast<PointerSource>(p2);
            }
            else if (p1->getType() == ParameterType::Operation)
            {
                return *std::dynamic_pointer_cast<OperationSource>(p1) ==
                       *std::dynamic_pointer_cast<OperationSource>(p2);
            }
        }

        return false;
    }

    bool ConstantSource::operator==(const ConstantSource &other) const
    {
        return this->value == other.value;
    }
    bool PointerSource::operator==(const PointerSource &other) const
    {
        return this->ptr == other.ptr;
    }
    bool OperationSource::operator==(const OperationSource &other) const
    {
        // All available operations are commutative
        if (this->op == other.op)
        {
            switch (this->op)
            {
            case ParamOpcode::Add:
            case ParamOpcode::Mul:
                return ((compare_sources(this->p1, other.p1) and compare_sources(this->p2, other.p2)) or
                        (compare_sources(this->p1, other.p2) and compare_sources(this->p2, other.p1)));
            case ParamOpcode::Div:
                return compare_sources(this->p1, other.p1) and compare_sources(this->p2, other.p2);
            case ParamOpcode::Sqrt:
                return compare_sources(this->p1, other.p1);
            }
        }
        return false;
    }

    bool Parameter::operator==(const Parameter &other) const
    {
        return compare_sources(this->source, other.source);
    }

    Parameter::operator double() const
    {
        return getValue();
    }

    std::ostream &operator<<(std::ostream &os, const Parameter &parameter)
    {
        os << parameter.getValue();

        return os;
    }

    bool Parameter::isZero() const
    {
        return source->getType() == ParameterType::Constant and getValue() == 0.;
    }

    bool Parameter::isOne() const
    {
        return source->getType() == ParameterType::Constant and getValue() == 1.;
    }

    Parameter Parameter::operator+(const Parameter &other) const
    {
        Parameter result = *this;
        result += other;
        return result;
    }

    Parameter &Parameter::operator+=(const Parameter &other)
    {
        if (other.isZero())
        {
            return *this;
        }
        else if (this->isZero())
        {
            *this = other;
            return *this;
        }
        else if (source->getType() == ParameterType::Constant and
                 other.source->getType() == ParameterType::Constant)
        {
            this->source = std::make_shared<ConstantSource>(getValue() + other.getValue());
            return *this;
        }
        else
        {
            this->source = std::make_shared<OperationSource>(ParamOpcode::Add,
                                                             this->source, other.source);
            return *this;
        }
    }

    Parameter Parameter::operator-() const
    {
        return Parameter(-1.) * *this;
    }

    Parameter Parameter::operator-(const Parameter &other) const
    {
        return *this + -other;
    }

    Parameter Parameter::operator*(const Parameter &other) const
    {
        Parameter result = *this;
        result *= other;
        return result;
    }

    Parameter Parameter::operator/(const Parameter &other) const
    {
        Parameter result = *this;
        result /= other;
        return result;
    }

    Parameter &Parameter::operator*=(const Parameter &other)
    {
        if (this->isZero())
        {
            return *this;
        }
        else if (other.isZero())
        {
            this->source = other.source;
            return *this;
        }
        else if (source->getType() == ParameterType::Constant and
                 other.source->getType() == ParameterType::Constant)
        {
            this->source = std::make_shared<ConstantSource>(getValue() * other.getValue());
            return *this;
        }
        else
        {
            this->source = std::make_shared<OperationSource>(ParamOpcode::Mul,
                                                             this->source, other.source);
            return *this;
        }
    }

    Parameter &Parameter::operator/=(const Parameter &other)
    {
        if (other.isZero())
        {
            throw std::runtime_error("Found a division by zero.");
        }

        if (this->isZero() or other.isOne())
        {
            return *this;
        }
        else if (source->getType() == ParameterType::Constant and
                 other.source->getType() == ParameterType::Constant)
        {
            this->source = std::make_shared<ConstantSource>(getValue() / other.getValue());
            return *this;
        }
        else
        {
            this->source = std::make_shared<OperationSource>(ParamOpcode::Div,
                                                             this->source, other.source);
            return *this;
        }
    }

    Parameter sqrt(const Parameter &param)
    {
        if (param.isZero())
        {
            return Parameter(0.);
        }
        else if (param.isOne())
        {
            return Parameter(1.);
        }
        else if (param.source->getType() == ParameterType::Constant)
        {
            assert(param.getValue() >= 0.);
            return Parameter(std::sqrt(param.getValue()));
        }
        else
        {
            Parameter param_sqrt;
            param_sqrt.source = std::make_shared<OperationSource>(ParamOpcode::Sqrt, param.source, nullptr);
            return param_sqrt;
        }
    }

} // namespace cvx