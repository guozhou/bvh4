#include "stdafx.h"
#include "BVH4.h"
#include "Ray.h"
#include <tbb/parallel_sort.h>
#include <tbb/parallel_reduce.h>
#include <tbb/parallel_scan.h>

BVH4::BVH4(const TriMesh* _mesh) : root(nullptr), mesh(_mesh)
{
}

BVH4::~BVH4()
{
	mesh = nullptr;
	if (root) { delete root; }
}

// Combining Single and Packet-Ray Tracing for Arbitrary Ray Distributions on the Intel MIC Architecture
bool BVH4::intersect(Ray& ray, const bool shadow) const
{
	NodeVector nodeVec;
	nodeVec.reserve(16);
	nodeVec.push_back({ root, -finf });
	bool hit = false;

	auto bsf = [](int& v) {
		unsigned long i = 0;
		_BitScanForward(&i, v);
		v &= v - 1;
		return i;
	};
	// http://stackoverflow.com/a/2789530/695128
	auto sort3 = [&nodeVec]() {
		assert(nodeVec.size() >= 3);
		auto a = nodeVec.end() - 3, b = nodeVec.end() - 2;
		if (std::get<1>(*a) < std::get<1>(*b)) { std::iter_swap(a, b); }
		b = nodeVec.end() - 1;
		if (std::get<1>(*a) < std::get<1>(*b)) { std::iter_swap(a, b); }
		a = nodeVec.end() - 2;
		if (std::get<1>(*a) < std::get<1>(*b)) { std::iter_swap(a, b); }
	};
	auto sort4 = [&nodeVec]() {
		assert(nodeVec.size() >= 4);
		auto a = nodeVec.end() - 4, b = nodeVec.end() - 3;
		if (std::get<1>(*a) < std::get<1>(*b)) { std::iter_swap(a, b); }
		a = nodeVec.end() - 2, b = nodeVec.end() - 1;
		if (std::get<1>(*a) < std::get<1>(*b)) { std::iter_swap(a, b); }
		a = nodeVec.end() - 4, b = nodeVec.end() - 2;
		if (std::get<1>(*a) < std::get<1>(*b)) { std::iter_swap(a, b); }
		a = nodeVec.end() - 3, b = nodeVec.end() - 1;
		if (std::get<1>(*a) < std::get<1>(*b)) { std::iter_swap(a, b); }
		a = nodeVec.end() - 3, b = nodeVec.end() - 2;
		if (std::get<1>(*a) < std::get<1>(*b)) { std::iter_swap(a, b); }
	};
	auto top = [&nodeVec]() { // obtain top of the traversal stack
		assert(nodeVec.size() > 0);
		return nodeVec.end() - 1;
	};

	while (nodeVec.size() > 0)
	{
		auto rec = top();
		if (std::get<1>(*rec) <= ray.tfar)
		{
			// go downward until we reach an outer node or all the inner ones are missed
			while (std::get<0>(*rec)->isInner())
			{
				const Inner* n = static_cast<const Inner*>(std::get<0>(*rec));
				Array4f tnear;
				int mask = n->intersect(ray, tnear);
				if (mask == 0) { break; } // no hit
				auto i = bsf(mask);
				// single hit, continue to visit it
				if (mask == 0) {
					*rec = { n->child[i], tnear[i] };
					continue;
				}

				const Node* c0 = n->child[i];
				float t0 = tnear[i];
				i = bsf(mask);
				const Node* c1 = n->child[i];
				float t1 = tnear[i];
				// two hits, continue to visit the near one and push the far child
				if (mask == 0)
				{
					if (t0 < t1) {
						*rec = { c1, t1 };
						nodeVec.push_back({ c0, t0 });
					}
					else {
						*rec = { c0, t0 };
						nodeVec.push_back({ c1, t1 });
					}
					rec = top();
					continue;
				}
				// roughly 10 percent of traversals hit more than 2 nodes
				i = bsf(mask);
				const Node* c2 = n->child[i];
				float t2 = tnear[i];
				nodeVec.resize(nodeVec.size() + 2);
				*(nodeVec.end() - 3) = { c0, t0 };
				*(nodeVec.end() - 2) = { c1, t1 };
				*(nodeVec.end() - 1) = { c2, t2 };
				if (mask == 0) { sort3(); }
				else { // four hits
					nodeVec.push_back({ n->child[3], tnear[3] });
					sort4();
				}
				rec = top();
			}

			if (!std::get<0>(*rec)->isInner()) {
				hit |= static_cast<const Outer*>(std::get<0>(*rec))->intersect(ray, shadow);
				if (shadow && hit) { return true; }
			}
		} // skip already missed node

		nodeVec.pop_back();
	}

	return hit;
}

float BVH4::halfArea(const AlignedBox3f& box) const
{
	Vector3f d = box.diagonal();
	float a = (d.x() * d.y() + d.y() * d.z() + d.z() * d.x());
	assert(a > 0.f);
	return a;
}

float BVH4::cost() const
{
	// Algorithms and Data Structures for Interactive Ray Tracing on Commodity Hardware
	// Object Partitioning Considered Harmful : Space Subdivision for BVHs
	assert(!rootBox.isEmpty());
	std::function<float(const Node*, float)> E = [this, &E](const Node* n, float A)
	{
		if (!n->isInner()) {
			// |N| * Ci
			return std::ceil(static_cast<const Outer*>(n)->size() / 4.f) * Ci;
		}
		else
		{
			float C = Ct / 2.f;
			for (int i = 0; i < 4; i++)
			{
				const Inner* ptr = static_cast<const Inner*>(n);
				if (ptr->child[i])
				{
					AlignedBox3f box;
					ptr->load(i, box);
					float Ac = halfArea(box);
					// Ct + (A(Nc) / A(N)) * E(Nc, A(Nc))
					C += (Ac / A * E(ptr->child[i], Ac));
				}
			}

			return C;
		}
	};

	return E(root, halfArea(rootBox));
}

void BVH4::build()
{
	assert(mesh->n_faces() < std::numeric_limits<int>::max());
	// compute axis-aligned bounding box for each primitive
	ReferenceVector refVec;
	refVec.resize(mesh->n_faces());
	tbb::parallel_for<unsigned int>(0, mesh->n_faces(), [&refVec, this](unsigned int i) {
		auto f = mesh->face_handle(i);
		AlignedBox3f box;
		for (auto fv : mesh->fv_range(f))
		{
			const TriMesh::Point pt = mesh->point(fv);
			box.extend(Vector3f(pt[0], pt[1], pt[2]));
		}
		refVec[i] = std::make_tuple(box.center(), box, f.idx());
	});

	std::vector<int> intVec(refVec.size());
	std::vector<AlignedBox3f> boxVec(refVec.size());

	rootBox = tbb::parallel_reduce(
		tbb::blocked_range<ReferenceVector::const_iterator>(refVec.cbegin(), refVec.cend()), AlignedBox3f(),
		[](tbb::blocked_range<ReferenceVector::const_iterator> refItR, AlignedBox3f init) {
		for (auto refIt : refItR) {
			init.extend(std::get<1>(refIt));
		}
		return init;
	},
		[](AlignedBox3f x, AlignedBox3f y) {
		return x.extend(y);
	});
	printf("rootBox (%4f, %4f, %4f), (%4f, %4f, %4f)\n",
		rootBox.min()[0], rootBox.min()[1], rootBox.min()[2],
		rootBox.max()[0], rootBox.max()[1], rootBox.max()[2]);
	// enqueue the whole vector
	ReferenceRecordQueue refRecQue[2];
	ReferenceRecord refRec;
	refRec.begin = 0;
	refRec.end = refVec.size();
	refRec.box = rootBox;
	refRec.pptr = &root;
	bool in = false;
	refRecQue[in].push_back({ { refRec, refRec, refRec, refRec }, 0 });
	while (!refRecQue[in].empty())
	{
		// apply splitting for current level
		tbb::parallel_for<int>(0, refRecQue[in].size(), [&](int i) {
			auto& refRec4i = refRecQue[in][i];
			std::get<1>(refRec4i) = split(refVec.data(), intVec.data(), boxVec.data(), 
				std::get<0>(refRec4i).data());
		});
		// allocate nodes according to splitting results
		while (!refRecQue[in].empty())
		{
			ReferenceRecord* refRec4 = std::get<0>(refRecQue[in].front()).data();
			const int c = std::get<1>(refRecQue[in].front());
			// pop front
			if (c > 1)
			{
				*(refRec4[0].pptr) = new Inner;
				Inner* ptr = static_cast<Inner*>(*(refRec4[0].pptr));
				for (int i = 0; i < c; i++) {
					assert(refRec4[i].size() > 0);
					refRec4[i].pptr = &(ptr->child[i]);
					ptr->store(i, refRec4[i].box);
				}

				if (c == 4) {
					// enqueue for further splitting
					for (int i = 0; i < 4; i++) {
						refRecQue[!in].push_back({ { refRec4[i], refRec4[i], refRec4[i], refRec4[i] }, 0 });
					}
					refRecQue[in].pop_front();
					continue;
				}
			}

			Outer* ptr = nullptr;
			for (int i = 0; i < c; i++)
			{
				assert(refRec4[i].size() > 0);
				*(refRec4[i].pptr) = new Outer(refRec4[i].size());
				ptr = static_cast<Outer*>(*(refRec4[i].pptr));
				// fetch per-primitive data
				for (int j = 0; j < refRec4[i].size(); j++)
				{
					Array3f p[3];
					int k = 0;
					auto f = OpenMesh::FaceHandle(std::get<2>(refVec[refRec4[i].begin + j]));
					for (auto fv : mesh->fv_range(f)) {
						const TriMesh::Point pt = mesh->point(fv);
						p[k] = { pt[0], pt[1], pt[2] };
						k++;
					}
					ptr->store(j, p[0], p[1], p[2]);
				}
#if _DEBUG
				outers.push_back(ptr);
#endif
			}
			// just created the final outer node then perform padding
			if (refRecQue[in].size() == 1 && refRecQue[!in].empty()) {
				ptr->store(refRec4[c - 1].size());
			}
			refRecQue[in].pop_front();
		}
		in = !in; // swap input and output queues
	}

	printf("Surface Area Cost: %4f\n", cost());
}

int BVH4::split(Reference* refVec, int* intVec, AlignedBox3f* boxVec, ReferenceRecord* refRec4) const
{
	int c = 1; // expecting 1~4 children can be obtained
	do
	{
		float minC = finf;
		int j; // find the child with minimum splitting candidate
		for (int i = 0; i < c; i++)
		{
			const int k = refRec4[i].begin;
			// innerSAH - outerSAH
			float relC = split(&refVec[k], &intVec[k], &boxVec[k], refRec4[i]) - refRec4[i].size() * Ci;
			if (relC < minC) {
				minC = relC;
				j = i;
			}
		}
		if (minC >= 0.f) { break; } // no more splitting
									// rearrange so that the new child is next to the one which will be split
		for (int i = c; i - 1 > j; i--) {
			std::swap(refRec4[i], refRec4[i - 1]);
		}

		auto& refRec = refRec4[j];
		ReferenceRecord::Dim dim = refRec.dim; // actual sorting
		std::sort(&refVec[refRec.begin], &refVec[refRec.end],
			[dim](const Reference& r0, const Reference& r1) {
			return std::get<0>(r0)[dim] < std::get<0>(r1)[dim];
		});
#if _DEBUG
		halfArea(refRec.lBox);
		halfArea(refRec.rBox);
#endif
		refRec4[j + 1].begin = refRec.pivot;
		refRec4[j + 1].end = refRec.end;
		refRec4[j + 1].box = refRec.rBox;
		refRec4[j + 1].pptr = refRec.pptr;
		refRec.end = refRec.pivot;
		refRec.box = refRec.lBox;

		c++;
	} while (c < 4);

	return c;
}

float BVH4::split(const Reference* refVec, int* intVec, AlignedBox3f* boxVec, ReferenceRecord& refRec) const
{
	float minC = finf, A = halfArea(refRec.box);
	const int size = refRec.size();
	for (int dim = 0; dim < 3; dim++)
	{
		std::iota(intVec, &intVec[size], 0); // sorted by center and perform splitting
		std::sort(intVec, &intVec[size], [dim, &refVec, &refRec](int p0, int p1) {
			return std::get<0>(refVec[p0])[dim] < std::get<0>(refVec[p1])[dim];
		});
		// inclusive prefix sum (i.e. union) AABBs from right to left
		AlignedBox3f left;
		for (int j = size - 1; j > 0; j--) {
			left.extend(std::get<1>(refVec[intVec[j]]));
			boxVec[j] = left;
		}
		// inclusive prefix sum (i.e. union) AABBs from left to right
		left.setEmpty();
		for (int j = 1; j < size; j++)
		{
			left.extend(std::get<1>(refVec[intVec[j - 1]]));
			// evaluate SAH for all the possible object split candidates
			float lA = halfArea(left), rA = halfArea(boxVec[j]);
			float C = Ct + (lA * j + rA * (size - j)) * Ci / A;
			assert(std::isnormal(C));
			if (C < minC)
			{
				minC = C;
				refRec.dim = ReferenceRecord::Dim(dim);
				refRec.pivot = refRec.begin + j;
				refRec.lBox = left;
				refRec.rBox = boxVec[j];
			}
		}
	}

	return minC;
}

Inner::Inner()
{
	std::fill_n(child, 4, nullptr);

	AlignedBox3f ab; // empty
	for (int i = 0; i < 4; i++)
	{
		const Vector3f& lower = ab.min(), upper = ab.max();
		data.row(i) = lower;
		data.row(7 - i) = upper; // mirroring
	}
}

void Inner::store(int i, const AlignedBox3f& ab)
{
	assert(0 <= i && i <= 3);
	const Vector3f& lower = ab.min(), upper = ab.max();
	data.row(i) = lower;
	data.row(7 - i) = upper;
}

void Inner::load(int i, AlignedBox3f& box) const
{
	assert(0 <= i && i <= 3);
	box = AlignedBox3f(data.row(i), data.row(7 - i));
}

int Inner::intersect(const Ray& ray, Array4f& tnear) const
{
	const Vector3f& rdir = ray.rdir;
	Array8x3f t;
	// sort the interval of slab along the ray
	for (int i = 0; i < 3; i++)
	{
		float r = rdir[i];
		// ray-slab test of 3 slabs for 4 axis-aligned boxes
		if (r < 0) {
			t.col(i) = (data.col(i).reverse() - ray.org[i]) * r;
		}
		else {
			t.col(i) = (data.col(i) - ray.org[i]) * r;
		}
	}
	assert(!t.hasNaN());
	// set intersection of the per-slab intersection interval
	tnear = t.block<4, 3>(0, 0).rowwise().maxCoeff().max(ray.tnear);
	Array4f tfar = t.block<4, 3>(4, 0).rowwise().minCoeff().min(ray.tfar).reverse();

	Array4b mask = (tnear <= tfar);
	return int(mask.w()) << 3 | int(mask.z()) << 2 | int(mask.y()) << 1 | int(mask.x());
}

Outer::Outer(int offset)
{
	assert(offset > 0);
	begin = data->rows();
	end = begin + offset;
	for (int i = 0; i < 4; i++) {
		data[i].conservativeResize(end, NoChange);
	}
}

ArrayXx3f Outer::data[4];

void Outer::store(int i, const Array3f& p0, const Array3f& p1, const Array3f& p2)
{
	i += begin;
	assert(i < end);
	Array3f e1 = p1 - p0, e2 = p2 - p0;
	data[0].row(i) = e1.matrix().cross(e2.matrix()); // Ng
	data[1].row(i) = p0;
	data[2].row(i) = e1;
	data[3].row(i) = e2;
}

void Outer::store(int offset)
{
	assert(data->rows() > 0);
	printf("Stored %td triangles in all the outer nodes\n", data->rows());
	offset = 4 - offset % 4;
	if (offset == 4) { return; } // no need to padding
	for (int i = 0; i < 4; i++)
	{
		data[i].conservativeResize(end + offset, NoChange);
		for (int j = end; j < data->rows(); j++) {
			// copying the first few triangles
			data[i].row(j) = data[i].row(j - end);
		}
	}
}

bool Outer::intersect(Ray& ray, const bool shadow) const
{
	bool hit = false;
	for (int i = begin; i < end; i += 4)
	{
		hit |= intersect(ray, i, shadow);
		if (shadow && hit) { return true; }
	}
	return hit;
}

bool Outer::intersect(Ray& ray, int i, const bool shadow) const
{
	assert(begin <= i && i < end);

	// Optimizing Ray-Triangle Intersection via Automated Search
	auto Ng = data[0].block<4, 3>(i, 0), p0 = data[1].block<4, 3>(i, 0),
		e1 = data[2].block<4, 3>(i, 0), e2 = data[3].block<4, 3>(i, 0);
	Array4f V = -(Ng.matrix() * ray.dir); // i.e. total volume for t = 1
	Array4x3f p0o = p0.rowwise() - ray.org.array().transpose();
	Array4x3f p0od = p0o.rowwise().cross(ray.dir);
	Array4f V2 = (p0od * e1).rowwise().sum(),
		V1 = -((p0od * e2).rowwise().sum()), V0 = V - V1 - V2;
	// inside/outside test for the whole (semi-infinite) ray
	Array4b mask = (V != 0.f) &&
		((V0 >= 0.f && V1 >= 0.f && V2 >= 0.f) || (V0 <= 0.f && V1 <= 0.f && V2 <= 0.f));
	if (!mask.any()) { return false; }
	// depth (t along the ray) test
	Array4f Vo = -((Ng * p0o).rowwise().sum()); // partial volume to one side (i.e. o) of triangle
	mask = mask && ((V * ray.tnear - Vo) * (V * ray.tfar - Vo) <= 0.f);
	if (!mask.any()) { return false; }
	if (shadow) { return true; }
	// locate the nearest triangle
	const Array4f rV = Ray::reciprocal(V), t = mask.select(Vo * rV, finf);
	int j;
	ray.tfar = t.minCoeff(&j);
	ray.Ng = Ng.row(j);
	ray.u = V1[j] * rV[j];
	ray.v = V2[j] * rV[j];
	return true;
}
