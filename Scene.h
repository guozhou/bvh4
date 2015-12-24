#pragma once
using namespace Eigen;

class Camera;
class BVH4;
class Image;
class Sampler;
class Scene
{
public:
	explicit Scene(const char* fileName);
	~Scene();

	void render() const;

	static const int width = 640, height = 360;
	const std::string fileName;
private:
	TriMesh* mesh;
	BVH4* bvh;
	Camera* camera;
	Image* image;
	Sampler* sampler;
};

