/**
 * @file epigraph.hpp
 * @author Sven Niederberger (s-niederberger@outlook.com)
 * @brief The Epigraph library. Enabled solvers will be included here.
 * @date 2020-08-11
 * 
 */
#pragma once

#ifdef ENABLE_ECOS
#include "ecosWrapper.hpp"
#endif

#ifdef ENABLE_EICOS
#include "eicosWrapper.hpp"
#endif

#ifdef ENABLE_OSQP
#include "osqpWrapper.hpp"
#endif