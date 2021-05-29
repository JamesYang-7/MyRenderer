#ifndef PHOTON_MAP_H
#define PHOTON_MAP_H

#include "hittable_list.h"
#include "aarect.h"
#include "material.h"
#include <iostream>
#include <queue>
#include <cmath> // for debugging

struct photon {
	Point3 p;
	Vec3 dir;
	Color color = Color();
	short flag = 0; // used to specify dimention in kd-tree
	double dist2 = 0;

	photon() {}
	photon(Point3 p_, Vec3 dir_, Color color_) : p(p_), dir(dir_), color(color_) {}
};

struct kd_node {
	photon pht;
	kd_node* left;
	kd_node* right;

	kd_node(photon p){
		pht = p;
		left = nullptr;
		right = nullptr;
	}
};

struct photon_ptr_compare {
	bool operator() (photon* p1, photon* p2) {
		return p1->dist2 < p2->dist2;
	}
};

class PhotonSort {
public:
	static void sort(photon* a, int lo, int hi, short flag);

private:
	static int compare(photon a, photon b, short flag) {
		if (a.p[flag] > b.p[flag]) return 1;
		else if (a.p[flag] < b.p[flag]) return -1;
		else return 0;
	}

	static void swap(photon* a, int i, int j) {
		photon t = a[i];
		a[i] = a[j];
		a[j] = t;
	}
};

class PhotonMap {
public:
	virtual ~PhotonMap() { delete(photon_array); }
	PhotonMap() : photon_num(0), ray_num(0), idx(1), max_array_length(0), max_search_radius(0), photon_array(nullptr), initial_energy(Color(0, 0, 0)) {}
	PhotonMap(int ray_num_, double search_radius_, Color init_energy_)
		: photon_num(0), ray_num(ray_num_), idx(1), max_search_radius(search_radius_), max_array_length(0), photon_array(nullptr), initial_energy(init_energy_) {}
	PhotonMap(int ray_num_, double search_radius_, Color init_energy_, const Hittable& world, XZ_Rect* light);
	void generate_global_photons(const Hittable& world, XZ_Rect* light);
	int generate_one_global_photon(const Hittable& world, XZ_Rect* light);
	kd_node* balance(int lo, int hi);
	short select_dimention(int lo, int hi);
	int locate_photons(kd_node* node, Point3 x, Vec3& ellipsoid_param,
		std::priority_queue<photon*, std::vector<photon*>, photon_ptr_compare>& que, Vec3 normal) const;
	int locate_photons(Point3 x, Vec3& ellipsoid_param,
		std::priority_queue<photon*, std::vector<photon*>, photon_ptr_compare>& que, Vec3 normal) const;
	virtual Color radiance_estimate(Point3 p, Color BRDF, Vec3& elps, Vec3 normal) const; // need to be const
	Color show_photon(Point3 p, Vec3& elps, Vec3 normal);
	bool has_space() { return idx < max_array_length; }

public:
	int photon_num = 0;
	int ray_num = 0;
	int idx = 1;
	int max_array_length = 0;
	double max_search_radius = 0;
	photon* photon_array;
	kd_node* root = nullptr; // need to release memory
	static const int MAX_DEPTH;
	static const int density_photons;
	size_t maxphotons = 0; // for debug
	Color initial_energy;
};

class CausticsPhotonMap : public PhotonMap {
public:
	virtual ~CausticsPhotonMap() {}
	CausticsPhotonMap(int ray_num_, double search_radius_, Color init_energy_, const Hittable& world, XZ_Rect* light, const Hittable* specular_objects);
	virtual Color radiance_estimate(Point3 p, Color BRDF, Vec3& elps, Vec3 normal) const override;

	void generate_caustics_photons(const Hittable& world, XZ_Rect* light, const Hittable* specular_objects);
	int generate_one_caustics_photon(const Hittable& world, const XZ_Rect* light, const Hittable* specular_objects);
};

class IndirectPhotonMap : public PhotonMap { // include only photons for indirect lighting
public:
	virtual ~IndirectPhotonMap() {}
	IndirectPhotonMap(int ray_num_, double search_radius_, Color init_energy_, const Hittable& world, XZ_Rect* light);
	void generate_indirect_photons(const Hittable& world, XZ_Rect* light);
	int generate_one_indirect_photon(const Hittable& world, XZ_Rect* light);
};

class PhotonRatioSampler {
public:
	~PhotonRatioSampler() {}
	PhotonRatioSampler(int ray_num_, const Hittable& world, XZ_Rect* light);
	int generate_one_global_photon(const Hittable& world, XZ_Rect* light);
	double get_ratio() { return ratio; }

	int ray_num = 0;
	int spec_num = 0;
	int non_spec_num = 0;
	double ratio = 0; // ratio of non_spec to spec
};

#endif // !PHOTON_MAP_H

