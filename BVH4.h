#pragma once
using namespace Eigen;

typedef Array<float, 8, 3> Array8x3f;
typedef Array<float, 16, 3> Array16x3f;
typedef Array<float, Dynamic, 3> ArrayXx3f;
typedef Array<bool, 4, 1> Array4b;
typedef Array<float, 4, 3> Array4x3f;

struct Ray;

struct Node
{
	virtual ~Node() {}
	virtual bool isInner() const = 0;
};

struct Inner : public Node
{
	Inner();
	~Inner() {
		for (int i = 0; i < 4; i++) { 
			if (child[i]) { delete child[i]; }
		} 
	}
	bool isInner() const { return true; }
	void store(int i, const AlignedBox3f& ab);
	void load(int i, AlignedBox3f& box) const;
	int intersect(const Ray& ray, Array4f& tnear) const;
	
	Node* child[4];
private:
	Array8x3f data; // (4 lower, 4 upper) * XYZ
};

struct Outer : public Node
{
	Outer(int offset);
	bool isInner() const { return false; }
	void store(int i, const Array3f& p0, const Array3f& p1, const Array3f& p2);
	void store(int offset); // padding the last node
	bool intersect(Ray& ray, const bool shadow = false) const;
	int size() const { 
		assert(begin >= 0 && end >= 0 && begin < end);
		return end - begin; 
	}

private:
	bool intersect(Ray& ray, int i, const bool shadow = false) const;

	int begin, end;
	static ArrayXx3f data[4]; // Ng, p0, e1 (p1-p0), e2 (p2-p0)
};

typedef std::pair<const Node*, float> NodeRecord; // node pointer and signed distance
typedef std::vector<NodeRecord> NodeVector; // BVH traversal stack

typedef std::tuple<Vector3f, AlignedBox3f, int> Reference; // center, box, ID
typedef std::vector<Reference> ReferenceVector;
struct ReferenceRecord
{
	int size() const {
		assert(begin >= 0 && end >= 0 && begin < end);
		return end - begin;
	}

	int begin, end;
	
	union {
		int pivot;
		float plane;
	};
	enum Dim { X, Y, Z, N } dim;
	
	AlignedBox3f box, lBox, rBox;
	// pointer to pointer to corresponding node of this record
	Node** pptr; 
};
typedef std::deque<std::pair<std::array<ReferenceRecord, 4>, int>> ReferenceRecordQueue;

class BVH4
{
public:
	BVH4(const TriMesh* mesh);
	~BVH4();

	bool intersect(Ray& ray, const bool shadow = false) const;
	void build();
	// On Quality Metrics of Bounding Volume Hierarchies
	// the cost of 2 ray-box tests and 1 ray-triangle test respectively
	static constexpr float Ct = 1.2f;
	static constexpr float Ci = 1.f;

#if _DEBUG
	std::vector<const Outer*> outers;
#endif
private:
	float halfArea(const AlignedBox3f& box) const;
	float cost() const;
	// attempt to split into 4 children at the most
	int split(Reference* refVec, int* intVec, AlignedBox3f* boxVec, ReferenceRecord* refRec4) const;
	float split(const Reference* refVec, int* intVec, AlignedBox3f* boxVec, ReferenceRecord& refRec) const;

	Node* root;
	AlignedBox3f rootBox;
	const TriMesh* mesh;
};