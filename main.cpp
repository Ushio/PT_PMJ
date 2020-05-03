#include "pr.hpp"
#include "rapidjson/document.h"
#include "rapidjson/rapidjson.h"
#include <iostream>
#include <memory>

#include "EzEmbree.hpp"

template <int C>
int nextPowerOf(int M)
{
	int a = 1;
	while (a < M)
	{
		a *= C;
	}
	return a;
}
class PMJSequence
{
public:
	// Highly reccommend to use integer for sample. 
	// We have to avoid numerical error because the stratum check must be strict.
	enum {
		RANDOM_MAX = 0x7FFFFF,
		RANDOM_LENGTH,
	};

	// this method should be called before first extend. 
	void setSeed(uint32_t s)
	{
		_seed = s;
	}
	void clear()
	{
		_samples.clear();
	}
	const glm::ivec2* samples() const
	{
		return _samples.data();
	}
	int size() const
	{
		return _samples.size();
	}

	void extend(int M)
	{
		int N = _samples.size();
		if (N == 0)
		{
			_random = decltype(_random)(_seed);
			_samples.emplace_back(
				_random.uniformi() % RANDOM_LENGTH,
				_random.uniformi() % RANDOM_LENGTH
			);
			N = 1;
		}
		_samples.resize(nextPowerOf<4>(M));

		// number of cells
		// N = 1, 4, 16, 64, 256...
		while (N < M)
		{
			std::vector<bool> xstratum, ystratum;
			buildOccupied(N, _samples, xstratum, ystratum);
			extendSequenceDiagonal(N, _samples, _random, xstratum, ystratum);
			// check whether stratums are filled.
			// PR_ASSERT(std::all_of(xstratum.begin(), xstratum.end(), [](bool b) { return b; }), "");
			// PR_ASSERT(std::all_of(ystratum.begin(), ystratum.end(), [](bool b) { return b; }), "");

			buildOccupied(N * 2, _samples, xstratum, ystratum);
			extendSequenceNonDiagonal(N * 2, _samples, _random, xstratum, ystratum);
			// check whether stratums are filled.
			// PR_ASSERT(std::all_of(xstratum.begin(), xstratum.end(), [](bool b) { return b; }), "");
			// PR_ASSERT(std::all_of(ystratum.begin(), ystratum.end(), [](bool b) { return b; }), "");

			// printf("generated: %d -> %d\n", N, N * 4);
			N = N * 4;
		}
	}
	glm::vec2 to01(glm::ivec2 s) const
	{
		return glm::vec2(s) / glm::vec2(RANDOM_LENGTH);
	}
private:
	/*
		i, j : cell index
		xhalf, yhalf: 0 or 1, these indicate the sub-cell
		n : rows, cols count

		(i,j) cell and sub cell
		+------------------+------------------+
		|(xhalf=0, yhalf=0)|(xhalf=1, yhalf=0)|
		+------------------+------------------+
		|(xhalf=0, yhalf=0)|(xhalf=1, yhalf=1)|
		+------------------+------------------+
	*/
	static glm::ivec2 generateSamplePoint(int N, int i, int j, int xhalf, int yhalf, int n, pr::IRandomNumberGenerator& random, std::vector<bool>& xstratum, std::vector<bool>& ystratum)
	{
		/*
		 This is the stratum count.
		 We'll generate [N, Nx2) samples on the current step.
		 So N stratums are already generated and filled so stratums should be Nx2 and these will be filled on the current step.
		*/
		int Nx2 = N * 2;

		int squareLength = (RANDOM_LENGTH / n);
		int halfSquareLength = squareLength / 2;
		int stratumLength = (RANDOM_LENGTH / Nx2);
		int x;
		for (;;)
		{
			x = i * squareLength + xhalf * halfSquareLength + (random.uniformi() % halfSquareLength);
			int xstratum_index = x / stratumLength;
			if (xstratum[xstratum_index] == false)
			{
				xstratum[xstratum_index] = true;
				break;
			}
		}

		int y;
		for (;;)
		{
			y = j * squareLength + yhalf * halfSquareLength + (random.uniformi() % halfSquareLength);
			int ystratum_index = y / stratumLength;
			if (ystratum[ystratum_index] == false)
			{
				ystratum[ystratum_index] = true;
				break;
			}
		}

		return {
			x,
			y
		};
	}

	static void buildOccupied(int N, std::vector<glm::ivec2>& samples, std::vector<bool>& xstratum, std::vector<bool>& ystratum)
	{
		/*
		 This is the stratum count.
		 We'll generate [N, Nx2) samples on the current step.
		 So N stratums are already generated and filled so stratums should be Nx2 and these will be filled on the current step.
		*/
		int Nx2 = N * 2;

		xstratum.clear();
		ystratum.clear();
		xstratum.resize(Nx2);
		ystratum.resize(Nx2);
		std::fill(xstratum.begin(), xstratum.end(), false);
		std::fill(ystratum.begin(), ystratum.end(), false);

		for (int i = 0; i < N; ++i)
		{
			int xstratum_index = samples[i].x / (RANDOM_LENGTH / Nx2);
			int ystratum_index = samples[i].y / (RANDOM_LENGTH / Nx2);
			xstratum[xstratum_index] = true;
			ystratum[ystratum_index] = true;
		}
	}

	// Generate [N, Nx2) sequence.
	// samples: sample sequence
	static void extendSequenceDiagonal(int N, std::vector<glm::ivec2>& samples, pr::IRandomNumberGenerator& random, std::vector<bool>& xstratum, std::vector<bool>& ystratum)
	{
		// number of rows, cols
		// n = 1, 2, 4, 8, 16
		int n = sqrt(N);
		for (int s = 0; s < N; ++s)
		{
			glm::ivec2 oldpt = samples[s];
			int squareLength = (RANDOM_LENGTH / n);
			int i = oldpt.x / squareLength;
			int j = oldpt.y / squareLength;
			int i_mod = oldpt.x % squareLength;
			int j_mod = oldpt.y % squareLength;

			// local sub-square index
			int xhalf = i_mod < (squareLength / 2) ? 0 : 1;
			int yhalf = j_mod < (squareLength / 2) ? 0 : 1;

			/* choose a diagonal child cell
			+-+-+
			|o| |
			+-+-+
			| |x|
			+-+-+

			o: first cell
			x: diagonal cell
			*/
			xhalf = 1 - xhalf;
			yhalf = 1 - yhalf;
			samples[N + s] = generateSamplePoint(N, i, j, xhalf, yhalf, n, random, xstratum, ystratum);
		}
	}
	// Generate [Nx2, Nx3) sequence.
	static void extendSequenceNonDiagonal(int Nx2, std::vector<glm::ivec2>& samples, pr::IRandomNumberGenerator& random, std::vector<bool>& xstratum, std::vector<bool>& ystratum)
	{
		int N = Nx2 / 2;
		// number of rows, cols
		// n = 1, 2, 4, 8, 16
		int n = sqrt(N);
		for (int s = 0; s < N; ++s)
		{
			glm::ivec2 oldpt = samples[s];
			int squareLength = (RANDOM_LENGTH / n);
			int i = oldpt.x / squareLength;
			int j = oldpt.y / squareLength;
			int i_mod = oldpt.x % squareLength;
			int j_mod = oldpt.y % squareLength;

			// local sub-square index
			int xhalf = i_mod < (squareLength / 2) ? 0 : 1;
			int yhalf = j_mod < (squareLength / 2) ? 0 : 1;

			/* choose a or b
			+-+-+
			|o|a|
			+-+-+
			|b|o|
			+-+-+
			*/
			if (random.uniformf() < 0.5f)
			{
				xhalf = 1 - xhalf;
			}
			else
			{
				yhalf = 1 - yhalf;
			}

			samples[N * 2 + s] = generateSamplePoint(Nx2, i, j, xhalf, yhalf, n, random, xstratum, ystratum);

			// b -> a or a -> b
			xhalf = 1 - xhalf;
			yhalf = 1 - yhalf;
			samples[N * 3 + s] = generateSamplePoint(Nx2, i, j, xhalf, yhalf, n, random, xstratum, ystratum);
		}
	}
private:
	uint32_t _seed = 1;
	pr::Xoshiro128StarStar _random;
	std::vector<glm::ivec2> _samples;
};

struct Polygon
{
public:
	uint32_t geometryID = 0;
	glm::mat4 xform = glm::identity<glm::mat4>();
	std::vector<glm::vec3> P;	   // Points
	std::vector<uint32_t> indices; // Vertices

	std::vector<glm::vec3> primitiveColors;
	std::vector<glm::vec3> primitiveLes;

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

			PR_ASSERT(d.HasMember("Vertices"));
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

			auto loadV3 = [](const rapidjson::Value& attrib, std::vector<glm::vec3> &values) {
				PR_ASSERT(attrib.IsArray());
				PR_ASSERT(attrib.Size() % 3 == 0);

				values.clear();
				values.reserve(attrib.Size() / 3);
				for (rapidjson::SizeType i = 0; i < attrib.Size(); i += 3)
				{
					PR_ASSERT(attrib[i].IsNumber());
					values.emplace_back(attrib[i].GetFloat(), attrib[i + 1].GetFloat(), attrib[i + 2].GetFloat());
				}
			};

			PR_ASSERT(d.HasMember("Primitives"));
			const rapidjson::Value& Primitives = d["Primitives"];
			PR_ASSERT(Primitives.IsObject());
			{
				const rapidjson::Value& Cd = Primitives["Cd"];
				PR_ASSERT(Cd.IsArray());
				loadV3(Cd, polygon->primitiveColors);

				const rapidjson::Value& Le = Primitives["Le"];
				PR_ASSERT(Le.IsArray());
				loadV3(Le, polygon->primitiveLes);
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

class Sampler {
public:
	Sampler() :_seed(0) {
		_dimOffsetter = pr::Xoshiro128StarStar(_seed);
		_baseOffset = { _dimOffsetter.uniformf(), _dimOffsetter.uniformf() };
	}
	Sampler(int seed) :_seed(seed) {
		_dimOffsetter = pr::Xoshiro128StarStar(_seed);
		_baseOffset = { _dimOffsetter.uniformf(), _dimOffsetter.uniformf() };
	}
	void clearDimension() {
		_dimOffsetter = pr::Xoshiro128StarStar(_seed);
		_baseOffset = { _dimOffsetter.uniformf(), _dimOffsetter.uniformf() };
	}
	glm::vec2 sample(const PMJSequence *pmj, int i) {
		float x = _dimOffsetter.uniformf();
		float y = _dimOffsetter.uniformf();
		auto sample = pmj->samples()[i % pmj->size()];
		return glm::fract(pmj->to01(sample) + _baseOffset + glm::vec2(x, y));
	}
	int _seed = 0;
	glm::vec2 _baseOffset;
	pr::Xoshiro128StarStar _dimOffsetter;
};


// z up
static glm::vec3 sampleLambertian(float a, float b, const glm::vec3 &Ng) {
	float r = std::sqrt(a);
	float theta = b * glm::pi<float>() * 2.0f;

	// uniform in xy circle, a = r * r
	float x = r * cos(theta);
	float y = r * sin(theta);

	// unproject to hemisphere
	float z = std::sqrt(std::max(1.0f - a, 0.0f));

	// local to global
	glm::vec3 xaxis;
	glm::vec3 yaxis;
	pr::GetOrthonormalBasis(Ng, &xaxis, &yaxis);
	return xaxis * x + yaxis * y + Ng * z;
}

int main()
{
	using namespace pr;

	SetDataDir( JoinPath( ExecutableDir(), "..", "data" ) );

	Config config;
	config.ScreenWidth = 1280;
	config.ScreenHeight = 720;
	config.SwapInterval = 1;
	Initialize( config );

	Camera3D camera;
	camera.origin = {0, 2.5f, 10};
	camera.lookat = {0, 2.5f, 0};
	camera.zUp = false;

	int focusX = 0;
	int focusY = 0;

	Scene scene;

	Stopwatch sw;
	scene.loadGeometry("out/box.json");
	 //scene.loadGeometry("out/sphere1.json");
	printf( "load: %f ms", sw.elapsed() * 1000.0f );
	scene.build();

	PMJSequence pmj;
	pmj.extend(8096);

	pr::ITexture *rtTexture = CreateTexture();
	int stride = 1;

	int iteration = 0;
	struct RayPayload {
		glm::dvec3 color;
		int samples = 0;
		Xoshiro128StarStar random;
		Sampler sampler;
	};
	std::vector<RayPayload> rayPayloads;

	double e = GetElapsedTime();

	while ( pr::NextFrame() == false )
	{
		bool restart = false;

		if (IsImGuiUsingMouse() == false)
		{
			bool e = UpdateCameraBlenderLike(&camera);
			if (e) {
				restart = true;
			}
		}

		if (IsKeyPressed(KEY_C))
		{
			focusX = GetMousePosition().x / stride;
			focusY = GetMousePosition().y / stride;
		}

		BeginCamera(camera);
		PushGraphicState();

		Image2DRGBA8 image;
		image.allocate(GetScreenWidth() / stride, GetScreenHeight() / stride);

		if (rayPayloads.size() != image.width() * image.height()) {
			rayPayloads.resize(image.width() * image.height());
			restart = true;
		}

		// Restart
		if (restart) {
			for (int i = 0; i < rayPayloads.size(); ++i) {
				rayPayloads[i].color = glm::vec3(0.0f);
				rayPayloads[i].samples = 0;
				rayPayloads[i].random = Xoshiro128StarStar(i + 1);
				rayPayloads[i].sampler = Sampler(i + 1);
			}
		}

		CameraRayGenerator rayGenerator(GetCurrentViewMatrix(), GetCurrentProjMatrix(), image.width(), image.height());

		ParallelFor(image.height(), [&](int j) {
			for (int i = 0; i < image.width(); ++i)
			{
				int index = j * image.width() + i;

				bool focus = focusX == i && focusY == j;
				if (focus) {
					printf("");
				}

				rayPayloads[index].sampler.clearDimension();
				// glm::vec2 rForShoot = { rayPayloads[index].random.uniformf(), rayPayloads[index].random.uniformf() };
				glm::vec2 rForShoot = rayPayloads[index].sampler.sample(&pmj, rayPayloads[index].samples);
				glm::vec3 ro, rd;
				rayGenerator.shoot(&ro, &rd, i, j, rForShoot.x, rForShoot.y);

				glm::vec3 accumColor = glm::vec3(0.0f);

				glm::vec3 T = glm::vec3(1.0f);
				for (int depth = 0; depth < 8; ++depth)
				{
					EzEmbree::EmbreeHit hit = scene.embree.intersect(ro, rd);
					if (hit.hasHit == false) {
						break;
					}
					auto hitGeom = scene.polygons[hit.geometryID];
					glm::vec3 Ng = glm::normalize(hit.Ng);
					if (0.0f < glm::dot(rd, Ng)) {
						Ng = -Ng;
					}
					glm::vec3 Cd = hitGeom->primitiveColors[hit.primitiveID];
					glm::vec3 Le = hitGeom->primitiveLes[hit.primitiveID];

					accumColor += T * Le;
					T *= Cd;

					glm::vec3 pHit = ro + rd * hit.t;

					// glm::vec2 rForSample = rayPayloads[index].sampler.sample(&pmj, rayPayloads[index].samples);
					glm::vec2 rForSample = { rayPayloads[index].random.uniformf(), rayPayloads[index].random.uniformf() };
					rd = sampleLambertian(
						rForSample.x,
						rForSample.y,
						Ng
					);
					ro = pHit + Ng * 0.001f;

					//glm::vec3 n = glm::normalize(hit.Ng);
					//accumColor = (n + glm::vec3(1.0f)) * 0.5f;
				}

				rayPayloads[index].color += accumColor;
				rayPayloads[index].samples += 1;

				// update pixels
				glm::vec3 color = glm::vec3(rayPayloads[index].color) / (float)rayPayloads[index].samples;
				color.x = std::pow(color.x, 1.0f / 2.2f);
				color.y = std::pow(color.y, 1.0f / 2.2f);
				color.z = std::pow(color.z, 1.0f / 2.2f);
				glm::ivec3 icolor = glm::ivec3( color * 255.0f + glm::vec3(0.5f) );
				icolor.x = std::min(icolor.x, 255);
				icolor.y = std::min(icolor.y, 255);
				icolor.z = std::min(icolor.z, 255);
				image(i, j) = { icolor.r, icolor.g, icolor.b, 255 };

				if (focus)
				{
					image(i, j) = { 255, 0, 0, 255 };
				}
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
		ImGui::Text("%d spp", rayPayloads[0].samples);
		ImGui::End();

		int nSample = rayPayloads[0].samples;
		if (32 <= nSample && __popcnt(rayPayloads[0].samples) == 1)
		{
			char buffer[256];
			// sprintf(buffer, "rand_%06d.png", nSample);
			sprintf(buffer, "pmj_%06d.png", nSample);
			image.save(buffer);
		}

		EndImGui();
	}

	pr::CleanUp();
}
