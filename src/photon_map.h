#ifndef PHOTON_MAP_H
#define PHOTON_MAP_H

#include "vec3.h"
#include "rtweekend.h"
#include "hittable_list.h"
#include "aarect.h"
#include "material.h"
#include <iostream>
#include <queue>

struct photon {
	Point3 p;
	Vec3 dir;
	Color color = Color(1, 1, 1);
	short flag = 0; // used to specify dimention in kd-tree
	double dist2 = 0;
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
	PhotonMap() : photon_num(0), photon_array(nullptr) {}
	PhotonMap(int photon_num_, const Hittable& world, XZ_Rect* light);
	void generate_global_photons(const Hittable& world, XZ_Rect* light);
	photon generate_one_global_photon(const Hittable& world, XZ_Rect* light);
	kd_node* balance(int lo, int hi);
	short select_dimention(int lo, int hi);
	int locate_photons(kd_node* node, Point3 x, Vec3& ellipsoid_param,
		std::priority_queue<photon*, std::vector<photon*>, photon_ptr_compare>& que, Vec3 normal) const;
	int locate_photons(Point3 x, Vec3& ellipsoid_param,
		std::priority_queue<photon*, std::vector<photon*>, photon_ptr_compare>& que, Vec3 normal) const;
	Color radiance_estimate(Point3 p, Color BRDF, Vec3& elps, Vec3 normal); // need to be const
	Color show_photon(Point3 p, Vec3& elps, Vec3 normal) const;

public:
	int photon_num = 0;
	photon* photon_array;
	kd_node* root = nullptr; // need to free
	static const int MAX_DEPTH;
	static const int density_photons;
	// size_t maxphotons = 0; // debug code
};

class CausticsPhotonMap : public PhotonMap {
public:
	virtual ~CausticsPhotonMap() {}
	CausticsPhotonMap(int photon_num_, const Hittable& world, XZ_Rect* light, const Hittable* specular_objects);

	void generate_caustics_photons(const Hittable& world, XZ_Rect* light, const Hittable* specular_objects);
	photon generate_one_caustics_photon(const Hittable& world, const XZ_Rect* light, const Hittable* specular_objects);
};

#endif // !PHOTON_MAP_H

