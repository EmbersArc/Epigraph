#include "variable.hpp"

#include <sstream>

namespace cvx::internal
{
    Variable::Variable(const std::string &name)
    {
        this->source = std::make_shared<VariableSource>();
        this->source->name = name;
        this->source->type = VariableSource::Type::Scalar;
    }

    Variable::Variable(const std::string &name, size_t row)
    {
        this->source = std::make_shared<VariableSource>();
        this->source->index = {row, 0};
        this->source->name = name;
        this->source->type = VariableSource::Type::Vector;
    }

    Variable::Variable(const std::string &name, size_t row, size_t col)
    {
        this->source = std::make_shared<VariableSource>();
        this->source->index = {row, col};
        this->source->name = name;
        this->source->type = VariableSource::Type::Matrix;
    }

    void Variable::unlink()
    {
        this->source->solution_ptr = nullptr;
        this->source->solution_idx = 0;
    }

    bool Variable::operator==(const Variable &other) const
    {
        return this->source == other.source;
    }

    bool Variable::isLinkedToSolver() const
    {
        return source->solution_ptr != nullptr;
    }

    bool Variable::linkToSolver(std::vector<double> *solution_ptr, size_t solution_idx)
    {
        if (isLinkedToSolver())
        {
            if (this->source->solution_ptr != solution_ptr)
            {
                throw std::runtime_error("Linking variables to multiple solvers is not supported.");
            }
            return false;
        }
        else
        {
            source->solution_ptr = solution_ptr;
            source->solution_idx = solution_idx;
            return true;
        }
    }

    double Variable::getSolution() const
    {
        if (not this->isLinkedToSolver())
        {
            // Don't throw here since variables might indeed be unused.
            return 0.;
        }
        else
        {
            return source->solution_ptr->at(source->solution_idx);
        }
    }

    size_t Variable::getProblemIndex() const
    {
        if (not this->isLinkedToSolver())
        {
            throw std::runtime_error("Variable must be linked to a problem first!");
        }

        return source->solution_idx;
    }

    std::ostream &operator<<(std::ostream &os, const Variable &variable)
    {
        os << variable.source->name;

        if (not(variable.source->type == VariableSource::Type::Scalar))
        {
            os << "[" << variable.source->index.first;

            if (variable.source->type == VariableSource::Type::Matrix)
            {
                os << ", " << variable.source->index.second;
            }

            os << "]";
        }
        if (variable.isLinkedToSolver())
            os << "@(" << variable.source->solution_idx << ")";

        return os;
    }

} // namespace cvx::internal
