#include "variable.hpp"

#include <sstream>

namespace cvx
{
    Variable::Variable(const std::string &name, size_t row, size_t col, VariableSource::Type type)
    {
        this->source = std::make_shared<VariableSource>();
        this->source->index = {row, col};
        this->source->name = name;
        this->source->type = type;
    }

    bool Variable::operator==(const Variable &other) const
    {
        return this->source == other.source;
    }

    bool Variable::isLinkedToProblem() const
    {
        return source->solution_ptr != nullptr;
    }

    void Variable::linkToProblem(std::vector<double> *solution_ptr, size_t solution_idx)
    {
        source->solution_ptr = solution_ptr;
        source->solution_idx = solution_idx;
    }

    double Variable::getSolution() const
    {
        if (not this->isLinkedToProblem())
        {
            // Don't throw error here since variables might indeed be unused.
            return 0.;
        }
        else
        {
            return source->solution_ptr->at(source->solution_idx);
        }
    }

    size_t Variable::getProblemIndex() const
    {
        if (not this->isLinkedToProblem())
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
        if (variable.isLinkedToProblem())
            os << "@(" << variable.source->solution_idx << ")";

        return os;
    }

} // namespace cvx
