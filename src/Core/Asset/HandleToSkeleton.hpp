#pragma once

#include <Core/RaCore.hpp>
#include <map>

namespace Ra {
namespace Core {
namespace Animation {
class Skeleton;
}
} // namespace Core
} // namespace Ra
namespace Ra {
namespace Core {
namespace Asset {

class HandleData;

/// Create a skeleton from a Handle Data extracted from a file.
/// Outputs the skeleton.
void RA_CORE_API createSkeleton( const Ra::Core::Asset::HandleData& data,
                                 Core::Animation::Skeleton& skelOut );

} // namespace Asset
} // namespace Core
} // namespace Ra
