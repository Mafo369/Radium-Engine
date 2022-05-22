#pragma once
#include <Core/RaCore.hpp>

namespace Ra {
namespace Core {
namespace Math {
class Quadric;
class DualQuaternion;
} // namespace Math

namespace Asset {
class AbstractVolume;
class AnimationData;
class AnimationTime;
class AssetData;
class BlinnPhongMaterialData;
class Camera;
template <class DATA>
class DataLoader;
class FileData;
class FileLoaderInterface;
class GeometryData;
class HandleData;
class LightData;
class MaterialData;
class Skeleton;
} // namespace Asset

namespace Containers {
class AdjacencyList;
class Grid;
class Tex;
} // namespace Containers

namespace Utils {
template <typename T>
class Attrib;
class AttribBase;
template <typename T>
class AttribHandle;
class CircularIndex;
class ContainerIntrospectionInterface;
class FunctionTask;
class Index;
template <typename T>
class IndexMap;
class IndexedObject;
template <typename T>
class Log;
class ObjectWithSemantic;
class ObservableVoid;
class Task;
class TaskQueue;
} // namespace Utils

namespace Geometry {
class AbstractDiscreteVolume;
class AbstractVolume;
class AttribArrayGeometry;
class CatmullClarkSubdivider;
class CubicBezier;
class Curve2D;
class GeometryIndexLayerBase;
template <typename T>
class IndexedGeometry;
class IndexedPointCloud;
class Line;
class LineMesh;
class LineStrip;
class LoopSubdivider;
class MultiIndexedGeometry;
class Obb;
class PointCloud;
class PolyLine;
class PolyMesh;
class QuadraSpline;
template <uint D, uint K>
class Spline;
class SplineCurve;
class TopologicalMesh;
class TriangleMesh;
class VolumeGrid;
class VolumeSparse;
} // namespace Geometry

namespace Animation {
class Cage;
class HandleArray;
template <typename VALUE_TYPE>
class KeyFramedValue;
class KeyFramedValueBase;
class KeyFramedValueController;
class Sequence;
class Skeleton;
} // namespace Animation
} // namespace Core
} // namespace Ra
