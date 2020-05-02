#include "pr.hpp"
#include "rapidjson/document.h"
#include "rapidjson/rapidjson.h"
#include <iostream>
#include <memory>

#include "EzEmbree.hpp"

struct Polygon
{
public:
	uint32_t geometryID = 0;
	glm::mat4 xform = glm::identity<glm::mat4>();
	std::vector<glm::vec3> P;	   // Points
	std::vector<uint32_t> indices; // Vertices

	EzEmbree::TriangleBufferDescriptor descriptor() const {
		EzEmbree::TriangleBufferDescriptor d;
		d.geometryID = geometryID;
		d.vertices = P.data();
		d.vertexCount = P.size();
		d.vertexStrideInBytes = sizeof(glm::vec3);
		d.indices = indices.data();
		d.primitiveCount = indices.size() / 3;
		d.xform = xform;
		return d;
	}
};

class Scene {
public:
	void loadGeometry(const char *name) {
		std::shared_ptr<Polygon> polygon(new Polygon());

		pr::BinaryLoader loader;
		loader.load(name);
		loader.push_back('\0');

		rapidjson::Document d;
		d.ParseInsitu((char*)loader.data());

		PR_ASSERT(d.HasParseError() == false);
		PR_ASSERT(d.HasMember("type"));

		const rapidjson::Value& objectType = d["type"];
		PR_ASSERT(objectType.IsString());

		if (objectType == "Polygon")
		{
			PR_ASSERT(d.HasMember("xform"));
			{
				const rapidjson::Value& xform = d["xform"];
				PR_ASSERT(xform.IsArray());
				PR_ASSERT(xform.Size() == 16);

				glm::mat4 m;
				for (rapidjson::SizeType i = 0; i < xform.Size(); i++)
				{
					PR_ASSERT(xform[i].IsNumber());
					glm::value_ptr(m)[i] = xform[i].GetFloat();
				}
				polygon->xform = m;
			}

			PR_ASSERT(d.HasMember("Points"));
			const rapidjson::Value& Points = d["Points"];
			PR_ASSERT(Points.IsObject());
			{
				PR_ASSERT(Points.HasMember("P"));
				const rapidjson::Value& P = Points["P"];
				PR_ASSERT(P.IsArray());
				PR_ASSERT(P.Size() % 3 == 0); // xyz
				polygon->P.reserve(P.Size() / 3);	// xyz
				for (rapidjson::SizeType i = 0; i < P.Size(); i += 3)
				{
					PR_ASSERT(P[i].IsNumber());
					PR_ASSERT(P[i + 1].IsNumber());
					PR_ASSERT(P[i + 2].IsNumber());
					float x = P[i].GetFloat();
					float y = P[i + 1].GetFloat();
					float z = P[i + 2].GetFloat();
					polygon->P.push_back(glm::vec3(x, y, z));
				}
			}

			const rapidjson::Value& Vertices = d["Vertices"];
			PR_ASSERT(Vertices.IsObject());
			{
				PR_ASSERT(Vertices.HasMember("Point Num"));
				const rapidjson::Value& PointNum = Vertices["Point Num"];
				PR_ASSERT(PointNum.IsArray());
				PR_ASSERT(PointNum.Size() % 3 == 0); // Triangle

				polygon->indices.reserve(PointNum.Size());
				for (rapidjson::SizeType i = 0; i < PointNum.Size(); i++)
				{
					PR_ASSERT(PointNum[i].IsNumber());
					polygon->indices.push_back(PointNum[i].GetInt());
				}
			}
		}

		polygon->geometryID = polygons.size();
		polygons.push_back(polygon);
	}
	void build() {
		for (auto polygon : polygons)
		{
			embree.add(polygon->descriptor());
		}
		embree.build();
	}
	std::vector<std::shared_ptr<Polygon>> polygons;
	EzEmbree::Embree embree;
};


int main()
{
	using namespace pr;

	SetDataDir( JoinPath( ExecutableDir(), "..", "data" ) );

	Config config;
	config.ScreenWidth = 1920;
	config.ScreenHeight = 1080;
	config.SwapInterval = 1;
	Initialize( config );

	Camera3D camera;
	camera.origin = {4, 4, 4};
	camera.lookat = {0, 0, 0};
	camera.zUp = false;

	Scene scene;

	Stopwatch sw;
	scene.loadGeometry("out/box.json");
	scene.loadGeometry("out/sphere1.json");
	printf( "load: %f ms", sw.elapsed() * 1000.0f );
	scene.build();

	pr::ITexture *rtTexture = CreateTexture();
	int stride = 4;

	double e = GetElapsedTime();

	while ( pr::NextFrame() == false )
	{
		if (IsImGuiUsingMouse() == false)
		{
			UpdateCameraBlenderLike(&camera);
		}

		BeginCamera(camera);
		PushGraphicState();

		Image2DRGBA8 image;
		image.allocate(GetScreenWidth() / stride, GetScreenHeight() / stride);

		CameraRayGenerator rayGenerator(GetCurrentViewMatrix(), GetCurrentProjMatrix(), image.width(), image.height());

		ParallelFor(image.height(), [&](int j) {
			for (int i = 0; i < image.width(); ++i)
			{
				glm::vec3 ro, rd;
				rayGenerator.shoot(&ro, &rd, i, j);
				EzEmbree::EmbreeHit hit = scene.embree.intersect(ro, rd);

				if (hit.hasHit) {
					glm::vec3 n = glm::normalize(hit.Ng);
					//if (hit.geometryID == 1) {
					//	n = glm::vec3(1, 0, 0);
					//}
					glm::vec3 color = (n + glm::vec3(1.0f)) * 0.5f;
					image(i, j) = { 255 * color.r, 255 * color.g, 255 * color.b, 255 };
				}
				else {
					image(i, j) = { 0, 0, 0, 255 };
				}

				//if (focus)
				//{
				//	image(i, j) = { 255, 0, 0, 255 };
				//}
			}
		});

		rtTexture->upload(image);

		// ClearBackground(0, 0, 0, 1);
		ClearBackground(rtTexture);


		DrawGrid( GridAxis::XZ, 1.0f, 10, {128, 128, 128} );
		DrawXYZAxis( 1.0f );

		//for (auto polygon : scene.polygons)
		//{
		//	BeginCameraWithObjectTransform( camera, polygon->xform );
		//	PrimBegin( PrimitiveMode::Lines );
		//	for ( int i = 0; i < polygon->indices.size(); i += 3 )
		//	{
		//		int a = polygon->indices[i];
		//		int b = polygon->indices[i + 1];
		//		int c = polygon->indices[i + 2];
		//		glm::u8vec3 color = {255, 255, 255};
		//		PrimVertex( polygon->P[a], color );
		//		PrimVertex( polygon->P[b], color );
		//		PrimVertex( polygon->P[b], color );
		//		PrimVertex( polygon->P[c], color );
		//		PrimVertex( polygon->P[c], color );
		//		PrimVertex( polygon->P[a], color );
		//	}
		//	PrimEnd();
		//	EndCamera();
		//}

		PopGraphicState();
		EndCamera();

		BeginImGui();

		ImGui::SetNextWindowSize( {500, 800}, ImGuiCond_Once );
		ImGui::Begin( "Panel" );
		ImGui::Text( "fps = %f", GetFrameRate() );
		ImGui::Image(rtTexture, ImVec2((float)rtTexture->width(), (float)rtTexture->height()));
		ImGui::End();

		EndImGui();
	}

	pr::CleanUp();
}
