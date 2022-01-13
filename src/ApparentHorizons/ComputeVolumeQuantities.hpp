// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstdint>
#include <type_traits>

#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tags/TempTensor.hpp"
#include "DataStructures/TempBuffer.hpp"
#include "DataStructures/Tensor/EagerMath/DeterminantAndInverse.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "DataStructures/Tensor/TypeAliases.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/Tags.hpp"
#include "PointwiseFunctions/GeneralRelativity/Tags.hpp"
#include "Utilities/ContainerHelpers.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/TMPL.hpp"

namespace detail {

template <typename Tag, typename DestTags, typename TempTags>
typename Tag::type& get_from_target_or_temp(
    const gsl::not_null<Variables<DestTags>*> target_vars,
    const gsl::not_null<TempBuffer<TempTags>*> temp_vars) {
  if constexpr (tmpl::list_contains<DestTags, Tag>::value) {
    (void) temp_vars; // silence gcc 'unused variable' warning
    return get<Tag>(*target_vars);
  } else {
    (void) target_vars; // silence gcc 'unused variable' warning
    return get<Tag>(*temp_vars);
  }
}
}  // namespace detail
