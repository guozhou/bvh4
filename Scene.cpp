#include "stdafx.h"
#include "Scene.h"
#include "BVH4.h"
#include "Camera.h"
#include "Image.h"
#include "Ray.h"
#include "Sampler.h"
#include "Util.h"
#include "Color.h"
#include <tbb/parallel_for.h>
#include <tbb/tbb_thread.h>
#include <tbb/enumerable_thread_specific.h>

Scene::Scene(const char* _fileName)
	:fileName(_fileName)
{
	mesh = new TriMesh;
	std::ifstream is(fileName, std::ios::in);
	OpenMesh::IO::Options ropt;
	if (!OpenMesh::IO::read_mesh(*mesh, is, ".obj", ropt))
	{
		throw "error loading mesh from file " + fileName;
	}
	printf("loaded mesh %s\n", _fileName);
	std::cout << "# Vertices: " << mesh->n_vertices() << std::endl;
	std::cout << "# Edges   : " << mesh->n_edges() << std::endl;
	std::cout << "# Faces   : " << mesh->n_faces() << std::endl;

	bvh = new BVH4(mesh);
	auto build = [this] { bvh->build(); };
	std::cout << "elapsed milliseconds " << Util::measure<>::execution(build) << std::endl;
	
	camera = new Camera(width, height);
	camera->lookAt({ 0.f, 0.795f, 4.f }, { 0.f, 0.795f, 0.f }, { 0.f, 1.f, 0.f });
	camera->perspective();

	image = new Image(width, height);

	sampler = new Sampler;
}

Scene::~Scene()
{
	delete sampler;
	delete image;
	delete camera;
	delete bvh;
	delete mesh;
}

void Scene::render() const
{
	union Id
	{
		Id() { j = tbb::this_tbb_thread::get_id(); }
		pcg32::state_type i;
		tbb::tbb_thread::id j;
	};
	tbb::enumerable_thread_specific<pcg32> rng;

	auto normal = [this] {
		for (int y = 0; y < height; y++)
		{
			tbb::parallel_for<int>(0, width, [&](int x) {
				Ray ray;
				camera->sampleRay(ray, { x + 0.5f, height - (y + 0.5f) });
				if (bvh->intersect(ray))
				{
					image->write(y, x, ray.Ng.normalized().array().abs());
				}
			});
			printf("row:%4d\t", y);
		}
	};

	auto directPoint = [this, &rng] {
		const int nSamples = 4; // no. of samples per pixel
		const Vector3f p(0.f, 1.f, 0.f); // point emitter
		for (int y = 0; y < height; y++)
		{
			tbb::parallel_for<int>(0, width, [&](int x) {
				auto& _rng = rng.local();
				_rng.set_stream(Id().i);
				auto roll = [&_rng]() { return unif(_rng); };

				float Li = 0.f;
				for (int z = 0; z < nSamples; z++)
				{
					Ray ray;
					camera->sampleRay(ray, { x + roll(), height - (y + roll()) });
					if (!bvh->intersect(ray)) { continue; }

					Vector3f o = ray.org + ray.tfar * ray.dir, d = p - o;
					float r = d.norm();
					Ray shadowRay(o, d / r);
					shadowRay.tfar = r;
					if (!bvh->intersect(shadowRay, true)) {
						Li += 0.8f / M_PI * ray.Ng.normalized().dot(d) / r;
					}
				}
				image->write(y, x, Array3f::Constant(Li / float(nSamples)));
			});
			printf("row:%4d\t", y);
		}
	};

	auto ao = [this, &rng](const float t) {
		const int nSamples = 4, nShadowRays = 16;
		for (int y = 0; y < height; y++)
		{
			tbb::parallel_for<int>(0, width, [&](int x) {
				auto& _rng = rng.local();
				_rng.set_stream(Id().i);
				auto roll = [&_rng]() { return unif(_rng); };
				
				float Li = 0.f;
				for (int z = 0; z < nSamples; z++)
				{
					Ray ray;
					camera->sampleRay(ray, { x + roll(), height - (y + roll()) });
					if (!bvh->intersect(ray)) { continue; }

					const Vector3f o = ray.org + ray.tfar * ray.dir;
					for (int w = 0; w < nShadowRays; w++)
					{
						Vector3f d = Ray::frame(ray.Ng.normalized()) *
							Sampler::square2CosHemisphere({ roll(), roll() });
						Ray shadowRay(o, d);
						shadowRay.tfar = t;
						if (!bvh->intersect(shadowRay, true))
						{
							Li += 1.f / nShadowRays;
						}
					}
				}
				image->write(y, x, Array3f::Constant(Li / float(nSamples)));
			});
			printf("row:%4d\t", y);
		}
	};

	auto spectral = [this, &rng] {
		const int nSamples = 8; // no. of samples per pixel
		const int nWavelengths = 4; // no. of wavelengths per sample
		const Vector3f p(0.f, 1.f, 0.f); // point emitter
		for (int y = 0; y < height; y++)
		{
			tbb::parallel_for<int>(0, width, [&](int x) {
				auto& _rng = rng.local();
				_rng.set_stream(Id().i);
				auto roll = [&_rng]() { return unif(_rng); };

				Array3f frag(Array3f::Zero());
				for (int z = 0; z < nSamples; z++)
				{
					Ray ray;
					camera->sampleRay(ray, { x + roll(), height - (y + roll()) });
					if (!bvh->intersect(ray)) { continue; }

					Vector3f o = ray.org + ray.tfar * ray.dir, d = p - o;
					float r = d.norm();
					Ray shadowRay(o, d / r);
					shadowRay.tfar = r;
					if (bvh->intersect(shadowRay, true)) { continue; }

					const float throughput = 0.8f / M_PI *
						ray.Ng.normalized().dot(d) / r;
					for (int i = 0; i < nWavelengths; i++)
					{
						const float lambda = Color::sampleWavelength(roll());
						// the pdf of sampling the visible spectrum is 1 / Color::WAVELENGTH_LEN
						frag += Color::SPD2XYZ(lambda,
							Color::WAVELENGTH_LEN / float(nWavelengths * nSamples) * throughput * Color::illuminantA(lambda));
					}
				}
				image->write(y, x, Color::XYZ2RGB(frag));
			});
			printf("row:%4d\t", y);
		}
	};

	std::cout << "elapsed milliseconds " << Util::measure<>::execution(spectral) << std::endl;
	//std::cout << "elapsed milliseconds " << Util::measure<>::execution(ao, 0.2f) << std::endl;
	//std::cout << "elapsed milliseconds " << Util::measure<>::execution(normal) << std::endl;
	//std::cout << "elapsed milliseconds " << Util::measure<>::execution(directPoint) << std::endl;
	image->write((fileName + ".exr").c_str());

#if 0
	Image* image0 = new Image(width, height);
	for (int y = 0; y < height; y++)
	{
#pragma omp parallel for
		for (int x = 0; x < width; x++)
		{
			bool hit = false;
			Ray ray;
			camera->sampleRay(ray, { x + 0.5f, height - (y + 0.5f) });
			for (auto o : bvh->outers)
			{
				hit |= o->intersect(ray);
			}

			if (hit) {
				image0->write(y, x, ray.Ng.array());
			}
		}
		printf("row:%4d\t", y);
	}
	image0->write((fileName + "direct.exr").c_str());
	delete image0;
#endif
}
