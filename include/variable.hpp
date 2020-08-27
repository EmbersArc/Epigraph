#pragma once

#include <utility>
#include <ostream>
#include <vector>
#include <string>
#include <memory>

namespace cvx
{
    class Scalar;

    namespace internal
    {

        class Term;

        struct VariableSource
        {
            enum class Type
            {
                Scalar,
                Vector,
                Matrix,
            };

            // The pointer gets deleted with the problem.
            std::vector<double> *solution_ptr = nullptr;
            size_t solution_idx = 0;
            std::string name;
            std::pair<size_t, size_t> index = {0, 0};
            Type type;
        };

        class Variable
        {
        public:
            Variable() = default;
            explicit Variable(const std::string &name);
            Variable(const std::string &name, size_t row);
            Variable(const std::string &name, size_t row, size_t col);

            bool operator==(const Variable &other) const;

            bool isLinkedToSolver() const;
            bool linkToSolver(std::vector<double> *solution_ptr, size_t solution_idx);
            double getSolution() const;
            size_t getProblemIndex() const;
            void unlink();

            friend std::ostream &operator<<(std::ostream &os,
                                            const Variable &variable);
            operator Term() const;
            operator Scalar() const;

        private:
            std::shared_ptr<VariableSource> source;
        };

    } // namespace internal

} // namespace cvx