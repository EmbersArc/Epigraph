#pragma once

#include <utility>
#include <ostream>
#include <vector>
#include <string>
#include <memory>

namespace cvx
{

    class Term;
    class Scalar;

    struct VariableSource
    {
        enum class Type
        {
            Scalar,
            Vector,
            Matrix,
        };
        std::vector<double> *solution_ptr = nullptr;
        size_t solution_idx = 0;
        std::string name;
        std::pair<size_t, size_t> index;
        Type type;
    };

    class Variable
    {
    public:
        Variable() = default;
        Variable(const std::string &name);
        Variable(const std::string &name, size_t row);
        Variable(const std::string &name, size_t row, size_t col);

        bool operator==(const Variable &other) const;

        bool isLinkedToProblem() const;
        void linkToProblem(std::vector<double> *solution_ptr, size_t solution_idx);
        double getSolution() const;
        size_t getProblemIndex() const;

        friend std::ostream &operator<<(std::ostream &os,
                                        const Variable &variable);
        operator Term() const;
        operator Scalar() const;

    private:
        std::shared_ptr<VariableSource> source;
    };

} // namespace cvx