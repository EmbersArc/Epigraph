#include "wrapperBase.hpp"

namespace cvx
{

    size_t WrapperBase::getNumVariables() const
    {
        return variables.size();
    }

    WrapperBase::~WrapperBase()
    {
        for (Variable &var : variables)
        {
            var.unlink();
        }
    }

} // namespace cvx
