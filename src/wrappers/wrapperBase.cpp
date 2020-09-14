#include "wrappers/wrapperBase.hpp"

namespace cvx::internal
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
