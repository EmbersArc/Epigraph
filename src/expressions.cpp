#include "expressions.hpp"

namespace cvx
{
    namespace internal
    {
        std::ostream &operator<<(std::ostream &os, const Term &term)
        {
            if (not term.parameter.isOne())
            {
                os << term.parameter.getValue() << " * ";
            }
            os << term.variable;
            return os;
        }

        std::ostream &operator<<(std::ostream &os, const Affine &affine)
        {
            for (size_t i = 0; i < affine.terms.size(); i++)
            {
                os << affine.terms[i];
                if (i != affine.terms.size() - 1)
                    os << " + ";
            }

            if (not affine.terms.empty() and not affine.constant.isZero())
            {
                os << " + ";
            }

            if (affine.terms.empty() or not affine.constant.isZero())
            {
                os << affine.constant;
            }

            return os;
        }

        std::ostream &operator<<(std::ostream &os, const Product &product)
        {
            if (product.factors.size() == 1)
            {
                os << "(" << product.factors[0] << ")^2";
            }
            else if (product.factors.size() == 2)
            {
                os << "(" << product.factors[0] << ") * (" << product.factors[1] << ")";
            }

            return os;
        }
    } // namespace internal

    using namespace internal;

    std::ostream &operator<<(std::ostream &os, const Scalar &exp)
    {
        if (not exp.products.empty())
        {
            if (exp.isNorm())
            {
                os << "(";
            }

            for (size_t i = 0; i < exp.products.size(); i++)
            {
                os << exp.products[i];

                if (i < exp.products.size() - 1)
                {
                    os << " + ";
                }
            }

            if (exp.isNorm())
            {
                os << ")^(1/2)";
            }
        }
        if (not exp.affine.isZero() and not exp.products.empty())
        {
            os << " + ";
        }
        if (not exp.affine.isZero())
        {
            os << exp.affine;
        }
        return os;
    }

    // Term

    Term::Term()
        : parameter(0.) {}

    bool Term::operator==(const Term &other) const
    {
        return this->parameter == other.parameter and
               this->variable == other.variable;
    }

    double Term::evaluate() const
    {
        return parameter.getValue() * variable.getSolution();
    }

    Term &Term::operator*=(const Parameter &param)
    {
        this->parameter *= param;
        return *this;
    }

    Term &Term::operator/=(const Parameter &param)
    {
        this->parameter /= param;
        return *this;
    }

    Term operator*(const Parameter &parameter, const Variable &variable)
    {
        Term term;

        term.parameter = parameter;
        term.variable = variable;

        return term;
    }

    Variable::operator Term() const
    {
        return Parameter(1.) * *this;
    }

    // Affine

    Affine operator*(const Parameter &param, const Affine &affine)
    {
        Affine result = affine;
        for (Term &term : result.terms)
        {
            term *= param;
        }
        result.constant *= param;
        return result;
    }

    double Affine::evaluate() const
    {
        double sum = this->constant.getValue();
        for (const Term &term : terms)
        {
            sum += term.evaluate();
        }
        return sum;
    }

    bool Affine::operator==(const Affine &other) const
    {
        return this->terms == other.terms and this->constant == other.constant;
    }

    Affine &Affine::operator+=(const Affine &other)
    {
        this->terms.insert(this->terms.end(), other.terms.cbegin(), other.terms.cend());

        this->constant += other.constant;

        return *this;
    }

    Affine &Affine::operator-=(const Affine &other)
    {
        *this += -other;

        return *this;
    }

    Affine &Affine::operator*=(const Parameter &param)
    {
        this->constant *= param;
        for (Term &term : this->terms)
        {
            term *= param;
        }
        return *this;
    }

    Affine &Affine::operator/=(const Parameter &param)
    {
        this->constant /= param;
        for (Term &term : this->terms)
        {
            term /= param;
        }
        return *this;
    }

    // Affine Affine::operator+(const Affine &other) const
    // {
    //     Affine result = *this;

    //     result += other;

    //     return result;
    // }

    Affine Affine::operator-() const
    {
        return Parameter(-1.) * *this;
    }

    Affine Affine::operator-(const Affine &other) const
    {
        Affine result = *this;
        result += Parameter(-1.) * other;
        return result;
    }

    bool Affine::isZero() const
    {
        return this->terms.empty() and this->constant.isZero();
    }

    bool Affine::isConstant() const
    {
        return this->terms.empty();
    }

    bool Affine::isFirstOrder() const
    {
        return not this->terms.empty();
    }

    Term::operator Affine() const
    {
        Affine affine;
        affine.terms = {*this};
        return affine;
    }

    Parameter::operator Affine() const
    {
        Affine affine;
        affine.constant = *this;
        return affine;
    }

    void Affine::cleanUp()
    {
        std::vector<Term> new_terms;
        for (const Term &term : this->terms)
        {
            auto same_variable = [&term](const Term &t) { return t.variable == term.variable; };

            auto found_term = std::find_if(new_terms.begin(), new_terms.end(), same_variable);

            if (found_term != new_terms.end())
            {
                found_term->parameter += term.parameter;
            }
            else
            {
                new_terms.push_back(term);
            }
        }

        auto parameter_is_zero = [](const Term &t) { return t.parameter.isZero(); };
        new_terms.erase(std::remove_if(new_terms.begin(), new_terms.end(), parameter_is_zero), new_terms.end());

        this->terms = new_terms;
    }

    // Product

    Product::Product(const Affine &term)
    {
        factors.push_back(term);
    }

    Product::Product(const Affine &lhs, const Affine &rhs)
    {
        factors.push_back(lhs);

        if (not(lhs == rhs))
        {
            factors.push_back(rhs);
        }
    }

    void Product::toSquaredTerm()
    {
        if (not isSquare())
        {
            // This is an edge case but necessary for cases like (p1 * x1) * (p2 * x1)
            if (firstTerm().terms.size() == 1 and
                secondTerm().terms.size() == 1 and
                firstTerm().constant.isZero() and
                secondTerm().constant.isZero() and
                firstTerm().terms[0].variable == secondTerm().terms[0].variable)
            {
                Parameter new_param = sqrt(firstTerm().terms[0].parameter * secondTerm().terms[0].parameter);
                Variable variable = firstTerm().terms[0].variable;

                factors = {new_param * variable};
            }
            else
            {
                throw std::runtime_error("Could not convert product expression into a squared expression.");
            }
        }
    }

    Affine &Product::firstTerm()
    {
        return factors[0];
    }

    Affine &Product::secondTerm()
    {
        if (isSquare())
        {
            return factors[0];
        }
        else
        {
            return factors[1];
        }
    }

    const Affine &Product::firstTerm() const
    {
        return factors[0];
    }

    const Affine &Product::secondTerm() const
    {
        if (isSquare())
        {
            return factors[0];
        }
        else
        {
            return factors[1];
        }
    }

    double Product::evaluate() const
    {
        if (isSquare())
        {
            return std::pow(factors[0].evaluate(), 2);
        }
        else
        {
            return factors[0].evaluate() * factors[1].evaluate();
        }
    }

    bool Product::operator==(const Product &other) const
    {
        return (this->firstTerm() == other.firstTerm() and
                this->secondTerm() == other.secondTerm()) or
               (this->firstTerm() == other.secondTerm() and
                this->secondTerm() == other.firstTerm());
    }

    bool Product::isSquare() const
    {
        return factors.size() == 1;
    }

    // Scalar

    Scalar::Scalar(int x)
    {
        this->affine.constant = Parameter(double(x));
    }

    Scalar::Scalar(double x)
    {
        this->affine.constant = Parameter(x);
    }

    Scalar::Scalar(double *x)
    {
        this->affine.constant = Parameter(x);
    }

    bool Scalar::operator==(const cvx::Scalar &other) const
    {
        bool equal = true;

        equal &= this->affine == other.affine;
        equal &= this->products == other.products;
        equal &= this->norm == other.norm;

        return equal;
    }

    double Scalar::evaluate() const
    {
        double sum = 0.;

        for (const Product &product : this->products)
        {
            sum += product.evaluate();
        }

        if (this->norm)
        {
            sum = std::sqrt(sum);
        }

        sum += this->affine.evaluate();

        return sum;
    }

    Scalar::operator double() const
    {
        return this->evaluate();
    }

    Scalar &Scalar::operator+=(const Scalar &other)
    {
        if ((this->isNorm() and other.getOrder() == 2) or
            (this->getOrder() == 2 and other.isNorm()) or
            (this->isNorm() and other.isNorm()))
        {
            throw std::runtime_error("Incompatible addition.");
        }

        this->affine += other.affine;

        this->products.insert(this->products.end(),
                              other.products.cbegin(),
                              other.products.cend());

        return *this;
    }

    Scalar &Scalar::operator-=(const Scalar &other)
    {
        if (other.getOrder() > 1)
        {
            throw std::runtime_error("Subtraction is not supported for higher-order terms.");
        }

        this->affine -= other.affine;

        return *this;
    }

    Scalar &Scalar::operator*=(const Scalar &other)
    {
        if (this->getOrder() == 2 or other.getOrder() == 2)
        {
            // TODO: Maybe allow scaling sums of squares, but not norms for performance reasons
            throw std::runtime_error("Factors in a multiplication have to be constant or linear.");
        }

        if (this->affine.isFirstOrder() and other.affine.isFirstOrder())
        {
            this->products.emplace_back(this->affine, other.affine);
            this->affine = Affine();
        }
        else if (this->affine.isConstant())
        {
            this->affine = this->affine.constant * other.affine;
        }
        else if (other.affine.isConstant())
        {
            this->affine *= other.affine.constant;
        }

        return *this;
    }

    Scalar &Scalar::operator/=(const Scalar &other)
    {
        if (this->getOrder() == 2)
        {
            throw std::runtime_error("The dividend has to be constant or linear.");
        }
        if (other.getOrder() > 0)
        {
            throw std::runtime_error("The divisor has to be constant.");
        }

        this->affine /= other.affine.constant;

        return *this;
    }

    Scalar Scalar::operator-() const
    {
        return Scalar(-1.) * *this;
    }

    Scalar operator+(const Scalar &lhs, const Scalar &rhs)
    {
        Scalar result = lhs;

        result += rhs;

        return result;
    }

    Scalar operator-(const Scalar &lhs, const Scalar &rhs)
    {
        Scalar result = lhs;

        result -= rhs;

        return result;
    }

    Scalar operator*(const Scalar &lhs, const Scalar &rhs)
    {
        Scalar result = lhs;

        result *= rhs;

        return result;
    }

    Scalar operator/(const Scalar &lhs, const Scalar &rhs)
    {
        Scalar result = lhs;

        result /= rhs;

        return result;
    }

    size_t Scalar::getOrder() const
    {
        if (not this->products.empty())
        {
            return 2;
        }
        else if (this->affine.isFirstOrder())
        {
            return 1;
        }
        else // if (not this->affine.isConstant())
        {
            return 0;
        }
    }

    bool Scalar::isNorm() const
    {
        return this->norm;
    }

    Scalar sqrt(const Scalar &scalar)
    {
        // Turn all higher order terms into squared terms
        Scalar e = scalar;

        for (Product &product : e.products)
        {
            product.toSquaredTerm();
        }

        if (e.affine.isConstant())
        {
            if (not e.affine.constant.isZero())
            {
                e.products.emplace_back(Affine(sqrt(e.affine.constant)));
                e.affine = Affine();
            }
            e.norm = true;
        }
        else
        {
            throw std::runtime_error("Can only take the square root when no linear terms are present.");
        }

        return e;
    }

    // Parameter::operator Scalar() const
    // {
    //     Scalar scalar;
    //     scalar.affine.constant = *this;
    //     return scalar;
    // }

    Variable::operator Scalar() const
    {
        Scalar scalar;
        scalar.affine.terms = {*this};
        return scalar;
    }

    Scalar par(double p)
    {
        return Scalar(p);
    }

    Scalar dynpar(double &p)
    {
        return Scalar(&p);
    }

    double eval(const Scalar &s)
    {
        return double(s);
    }

} // namespace cvx