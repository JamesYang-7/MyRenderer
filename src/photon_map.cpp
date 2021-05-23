#include "photon_map.h"

// PhotonSort

void PhotonSort::sort(photon* a, int lo, int hi, short flag) {
	if (hi <= lo) return;
	int lt = lo, gt = hi;
	photon v = a[lo];
	int i = lo;
	while (i <= gt) {
		int cmp = compare(a[i], v, flag);
		if (cmp < 0) swap(a, lt++, i++);
		else if (cmp > 0) swap(a, i, gt--);
		else ++i;
	}
	sort(a, lo, lt - 1, flag);
	sort(a, gt + 1, hi, flag);
}

// PhotonMap

const int PhotonMap::MAX_DEPTH = 16;
const int PhotonMap::density_photons = 40;

PhotonMap::PhotonMap(int photon_num_, const Hittable& world, XZ_Rect* light) : photon_num(photon_num_) {
	photon_array = new photon[photon_num + 1];
	generate_global_photons(world, light);
	std::cout << "constructing photon map...\n";
	root = balance(1, photon_num);
}

void PhotonMap::generate_global_photons(const Hittable& world, XZ_Rect* light) {
	int idx = 1;
	while (idx <= photon_num) {
		std::cerr << "\rPhotons remaining: " << photon_num - idx << ' ' << std::flush;
		photon_array[idx] = generate_one_global_photon(world, light);
		++idx;
	}
	std::cout << "\nSuccessfully generated photons \n";
}

photon PhotonMap::generate_one_global_photon(const Hittable& world, XZ_Rect* light) {
	photon pht;
	HitRecord hrec;
	ScatterRecord srec;
	Ray ray = light->generate_emitted_ray();
	int depth = MAX_DEPTH;
	do {
		if (!world.hit(ray, 1e-4, infinity, hrec)) { // re-emit if didn't hit anything
			ray = light->generate_emitted_ray();
			depth = MAX_DEPTH;
			continue;
		}
		if (!hrec.mat_ptr->scatter(ray, hrec, srec)) { // re-emit if hits light
			ray = light->generate_emitted_ray();
			depth = MAX_DEPTH;
			continue;
		}
		if (srec.is_specular) { // scatter if hits specular surface
			ray = srec.scattered;
			pht.color = pht.color * srec.albedo; // modify energy
			--depth;
		}
		else { // if hits diffuse surface
			if (depth == 0 || random_double() < 0.5) { // using reflectance = 0.5 (may need to be modified)
				// store photon
				pht.p = hrec.p;
				pht.dir = ray.direction();
				break;
			}
			else { // scatter
				ray = srec.scattered;
				pht.color = srec.albedo * pht.color; // modify energy
				--depth;
			}
		}
	} while (depth > 0);
	return pht;
}

kd_node* PhotonMap::balance(int lo, int hi) {
	if (hi < lo) {
		return nullptr;
	}
	short flag = select_dimention(lo, hi);
	PhotonSort::sort(photon_array, lo, hi, flag);
	int mid = (lo + hi) / 2;
	kd_node* node = new kd_node(photon_array[mid]);
	node->pht.flag = flag;
	node->left = balance(lo, mid - 1);
	node->right = balance(mid + 1, hi);
	return node;
}

short PhotonMap::select_dimention(int lo, int hi) { // select dimention in which the bounding box is largest
	short res = 0;
	Point3 v_min = Point3(infinity, infinity, infinity);
	Point3 v_max = Point3(-infinity, -infinity, -infinity);
	for (int i = lo; i <= hi; ++i) {
		for (int j = 0; j < 2; ++j) {
			v_min[j] = std::min(v_min[j], photon_array[i].p[j]);
			v_max[j] = std::max(v_max[j], photon_array[i].p[j]);
		}
	}
	double dim[3] = {};
	for (int i = 0; i < 3; ++i) { dim[i] = v_max[i] - v_min[i]; }
	res = dim[0] > dim[1] ? 0 : 1;
	res = dim[res] > dim[2] ? res : 2;
	return res;
}

int PhotonMap::locate_photons(kd_node* node, Point3 x, Vec3& elps,
	std::priority_queue<photon*, std::vector<photon*>, photon_ptr_compare>& que, Vec3 normal) const {
	int res = 0;
	if (node == nullptr) {
		return 0;
	}
	Vec3 disp = node->pht.p - x;
	if (dot(disp, normal) >= -(1e-4)
		&& disp[0] * disp[0] / elps[0] + disp[1] * disp[1] / elps[1] + disp[2] * disp[2] / elps[2] <= 1) {
		res += 1;
		node->pht.dist2 = disp.length_squared();
		que.push(&(node->pht));
		if (que.size() > density_photons) {
			que.pop();
			elps[0] = elps[1] = elps[2] = (que.top()->p - x).length_squared();
		}
	}
	short flag = node->pht.flag;
	double delta = x[flag] - node->pht.p[flag]; // distance between x and the splitting plane
	if (delta < 0) {
		res += locate_photons(node->left, x, elps, que, normal);
		if (delta * delta < elps[flag]) {
			res += locate_photons(node->right, x, elps, que, normal);
		}
	}
	else {
		res += locate_photons(node->right, x, elps, que, normal);
		if (delta * delta < elps[flag]) {
			res += locate_photons(node->left, x, elps, que, normal);
		}
	}
	return res;
}

int PhotonMap::locate_photons(Point3 x, Vec3& elps,
	std::priority_queue<photon*, std::vector<photon*>, photon_ptr_compare>& que, Vec3 normal) const {
	return locate_photons(root, x, elps, que, normal);
}

Color PhotonMap::radiance_estimate(Point3 p, Color BRDF, Vec3& elps, Vec3 normal) { // need to be const
	std::priority_queue<photon*, std::vector<photon*>, photon_ptr_compare> que;
	locate_photons(root, p, elps, que, normal);
	Color radiance = Color();
	photon* p_pht;
	if (que.empty()) {
		return Color();
	}
	// maxphotons = std::max(maxphotons, que.size()); // debug code
	double r = std::sqrt(que.top()->dist2);
	if (que.size() < density_photons) {
		r = 20.0;
	}
	while (!que.empty()) {
		p_pht = que.top();
		radiance += p_pht->color * (dot(p_pht->dir, normal) > 0 ? 0 : BRDF); // ignore radiance from the back of the surface
		que.pop();
	}
	radiance = radiance / (0.5 * PI * r * r);
	return radiance;
}

Color PhotonMap::show_photon(Point3 p, Vec3& elps, Vec3 normal) const {
	std::priority_queue<photon*, std::vector<photon*>, photon_ptr_compare> que;
	locate_photons(root, p, elps, que, normal);
	Color radiance = Color();
	Color white = Color(1, 1, 1);
	double ratio = photon_num / 1e6;
	if (que.empty()) {
		return Color();
	}
	else {
		while (!que.empty()) {
			radiance += white * (que.top()->color.length_squared());
			que.pop();
		}
	}
	return radiance * ratio;
}

// CausticsPhotonMap

CausticsPhotonMap::CausticsPhotonMap(int photon_num_, const Hittable& world, XZ_Rect* light, const Hittable* specular_objects)
	: PhotonMap() {
	photon_num = photon_num_;
	photon_array = new photon[photon_num + 1];
	generate_caustics_photons(world, light, specular_objects);
	std::cout << "constructing photon map...\n";
	root = balance(1, photon_num);
}

void CausticsPhotonMap::generate_caustics_photons(const Hittable& world, XZ_Rect* light, const Hittable* specular_objects)
{
	int idx = 1;
	while (idx <= photon_num) {
		std::cerr << "\rCaustics photons remaining: " << photon_num - idx << ' ' << std::flush;
		photon_array[idx] = generate_one_caustics_photon(world, light, specular_objects);
		++idx;
	}
	std::cout << "\nSuccessfully generated caustics photons \n";
}

photon CausticsPhotonMap::generate_one_caustics_photon(const Hittable& world, const XZ_Rect* light, const Hittable* specular_objects) {
	photon pht;
	HitRecord hrec;
	ScatterRecord srec;
	Point3 random_p = light->random_point();
	Ray ray = Ray(random_p, specular_objects->generate_random(random_p));
	int depth = MAX_DEPTH;
	bool hits_specular_first = false;
	do {
		if (!world.hit(ray, 1e-4, infinity, hrec)) { // re-emit if didn't hit anything
			random_p = light->random_point();
			ray = Ray(random_p, specular_objects->generate_random(random_p));
			depth = MAX_DEPTH;
			continue;
		}
		if (!hrec.mat_ptr->scatter(ray, hrec, srec)) { // re-emit if hits light
			random_p = light->random_point();
			ray = Ray(random_p, specular_objects->generate_random(random_p));
			depth = MAX_DEPTH;
			continue;
		}
		if (srec.is_specular) { // scatter if hits specular surface
			ray = srec.scattered;
			pht.color = pht.color * srec.albedo;
			if (!hits_specular_first) { hits_specular_first = true; }
			--depth;
		}
		else { // if hits diffuse surface
			if (hits_specular_first) { // if hits specular first, store photon
				pht.p = hrec.p;
				pht.dir = ray.direction();
				break;
			}
			else { // re-emit
				random_p = light->random_point();
				ray = Ray(random_p, specular_objects->generate_random(random_p));
				depth = MAX_DEPTH;
				continue;
			}
		}
	} while (depth > 0);
	return pht;
}