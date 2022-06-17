#include <Engine/Scene/ComponentMessenger.hpp>

#include <functional>

namespace Ra {
namespace Engine {
namespace Scene {

RA_SINGLETON_IMPLEMENTATION( ComponentMessenger );

void ComponentMessenger::unregisterAll( const Entity* entity, Component* comp ) {
    for ( auto& entityList : std::array<std::reference_wrapper<EntityMap>, 3> {
              m_entitySetLists, m_entityGetLists, m_entityRwLists } ) {
        auto listItr = entityList.get().find( entity );
        if ( listItr != entityList.get().end() ) {
            CallbackMap& list = listItr->second;
            auto pred = [comp]( const Key& k ) -> bool { return k.first == comp->getName(); };
            for ( auto it = list.begin(); it != list.end(); ) {
                if ( pred( it->first ) )
                    it = list.erase( it );
                else
                    ++it;
            }
        }
    }
}

} // namespace Scene
} // namespace Engine
} // namespace Ra
