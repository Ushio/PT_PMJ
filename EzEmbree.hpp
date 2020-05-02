#pragma once

/*
	Premake Example Setting

	-- embree
	includedirs { "libs/embree/include" }
	libdirs { "libs/embree/lib" }
	filter {"Debug"}
		links { "embree3", "tbb_debug", "tbbmalloc_debug" }
	filter {"Release"}
		links { "embree3",  "tbb", "tbbmalloc" }
	filter{}
	postbuildcommands {
		"{COPY} %{prj.location}../libs/embree/bin/embree3.dll %{cfg.targetdir}",
		"{COPY} %{prj.location}../libs/embree/bin/tbb.dll %{cfg.targetdir}",
		"{COPY} %{prj.location}../libs/embree/bin/tbb_debug.dll %{cfg.targetdir}",
		"{COPY} %{prj.location}../libs/embree/bin/tbbmalloc.dll %{cfg.targetdir}",
		"{COPY} %{prj.location}../libs/embree/bin/tbbmalloc_debug.dll %{cfg.targetdir}",
	}
*/

#include <embree3/rtcore.h>
#include "pr.hpp"

namespace EzEmbree {
	struct TriangleBufferDescriptor
	{
		uint32_t geometryID = 0;

		const glm::vec3 *vertices = nullptr;
		uint32_t vertexStrideInBytes = sizeof(glm::vec3);
		uint32_t vertexCount = 0;

		const uint32_t *indices = nullptr;
		uint32_t primitiveCount = 0;

		glm::mat4 xform = glm::identity<glm::mat4>();
	};

	template <class T>
	class RTCPtr {
	public:
		typedef void(*Functor)(T);
		RTCPtr(T p, Functor retain, Functor release)
			: _ptr(p)
			, _retain(retain)
			, _release(release)
		{
		}
		RTCPtr() {

		}
		~RTCPtr() {
			if (_ptr) {
				_release(_ptr);
				_ptr = nullptr;
			}
		}
		RTCPtr(const RTCPtr& rhs) :
			: _ptr(rhs._ptr)
			, _retain(rhs._retain)
			, _release(rhs._release)
		{
			if (_ptr) {
				_retain(_ptr);
			}
		}

		RTCPtr& operator=(const RTCPtr& rhs)
		{
			if (rhs._ptr) {
				rhs._retain(rhs._ptr);
			}
			if (_ptr) {
				_release(_ptr);
			}

			_ptr = rhs._ptr;
			_retain = rhs._retain;
			_release = rhs._release;

			return *this;
		}
		T get() const {
			return _ptr;
		}
	private:
		T _ptr = nullptr;
		Functor _retain = [](T) {};
		Functor _release = [](T) {};
	};
	struct EmbreeHit {
		EmbreeHit() {

		}
		bool hasHit = false;

		float t = 0.0f;

		uint32_t geometryID;
		uint32_t primitiveID;

		// If you provide vertices as ClockWise (e.g. Houdini), you might want to flip the normal.
		glm::vec3 Ng;

		/*
		The parametrization of a triangle uses the first vertex p0 as base point, the
		vector p1 - p0 as u-direction and the vector p2 - p0 as v-direction. Thus vertex
		attributes t0,t1,t2 can be linearly interpolated over the triangle the following
		way:
			t_uv = (1-u-v)*t0 + u*t1 + v*t2
			     = t0 + u*(t1-t0) + v*(t2-t0)
		*/
		float u = 0.0f;
		float v = 0.0f;
	};
	class Embree {
	public:
		Embree() {
			_device = RTCPtr<RTCDevice>(rtcNewDevice("set_affinity=1"), rtcRetainDevice, rtcReleaseDevice);
			_scene = RTCPtr<RTCScene>(rtcNewScene(_device.get()), rtcRetainScene, rtcReleaseScene);
			rtcSetSceneBuildQuality(_scene.get(), RTC_BUILD_QUALITY_HIGH);
		}
		void add(TriangleBufferDescriptor descriptor)
		{
			RTCPtr<RTCScene> instanceSrcScene(rtcNewScene(_device.get()), rtcRetainScene, rtcReleaseScene);
			rtcSetSceneBuildQuality(instanceSrcScene.get(), RTC_BUILD_QUALITY_HIGH);

			RTCPtr<RTCGeometry> g(rtcNewGeometry(_device.get(), RTC_GEOMETRY_TYPE_TRIANGLE), rtcRetainGeometry, rtcReleaseGeometry);
			rtcSetSharedGeometryBuffer(g.get(), RTC_BUFFER_TYPE_VERTEX, 0, RTC_FORMAT_FLOAT3, descriptor.vertices, 0 /*byteoffset*/, descriptor.vertexStrideInBytes, descriptor.vertexCount);
			rtcSetSharedGeometryBuffer(g.get(), RTC_BUFFER_TYPE_INDEX, 0, RTC_FORMAT_UINT3, descriptor.indices, 0 /*byteoffset*/, sizeof(uint32_t) * 3, descriptor.primitiveCount);
			rtcCommitGeometry(g.get());
			rtcAttachGeometryByID(instanceSrcScene.get(), g.get(), descriptor.geometryID);
			rtcCommitScene(instanceSrcScene.get());

			RTCPtr<RTCGeometry> instance(rtcNewGeometry(_device.get(), RTC_GEOMETRY_TYPE_INSTANCE), rtcRetainGeometry, rtcReleaseGeometry);
			rtcSetGeometryInstancedScene(instance.get(), instanceSrcScene.get());
			rtcSetGeometryTransform(instance.get(), 0, RTC_FORMAT_FLOAT4X4_COLUMN_MAJOR, glm::value_ptr(descriptor.xform));
			rtcSetGeometryTimeStepCount(instance.get(), 1);
			rtcCommitGeometry(instance.get());
			rtcAttachGeometry(_scene.get(), instance.get());
		}
		void build() {
			rtcCommitScene(_scene.get());
		}

		EmbreeHit intersect( const glm::vec3& ro, const glm::vec3& rd ) {
			RTCIntersectContext context;
			rtcInitIntersectContext(&context);

			RTCRayHit rayhit;
			rayhit.ray.org_x = ro.x;
			rayhit.ray.org_y = ro.y;
			rayhit.ray.org_z = ro.z;
			rayhit.ray.dir_x = rd.x;
			rayhit.ray.dir_y = rd.y;
			rayhit.ray.dir_z = rd.z;
			rayhit.ray.time = 0.0f;

			rayhit.ray.tfar = FLT_MAX;
			rayhit.ray.tnear = 0.0f;

			rayhit.ray.mask = 0;
			rayhit.ray.id = 0;
			rayhit.ray.flags = 0;
			rayhit.hit.geomID = RTC_INVALID_GEOMETRY_ID;
			rayhit.hit.instID[0] = RTC_INVALID_GEOMETRY_ID;
			rtcIntersect1(_scene.get(), &context, &rayhit);

			if (rayhit.hit.geomID == RTC_INVALID_GEOMETRY_ID) {
				return EmbreeHit();
			}

			EmbreeHit h;
			h.hasHit = true;
			h.t = rayhit.ray.tfar;
			h.geometryID = rayhit.hit.geomID;
			h.primitiveID = rayhit.hit.primID;
			h.Ng = glm::vec3(rayhit.hit.Ng_x, rayhit.hit.Ng_y, rayhit.hit.Ng_z);

			h.u = rayhit.hit.u;
			h.v = rayhit.hit.v;

			return h;
		}
	private:
		RTCPtr<RTCDevice> _device;
		RTCPtr<RTCScene> _scene;
	};
}
