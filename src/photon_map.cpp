// PhotonSort
#include "photon_map.h"

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
const int PhotonMap::density_photons = 100;

PhotonMap::PhotonMap(int ray_num_, double search_radius_, Color init_energy_, const Hittable& world, XZ_Rect* light)
	:photon_num(0), ray_num(ray_num_), idx(1), max_search_radius(search_radius_), photon_array(nullptr), initial_energy(init_energy_)
{
	max_array_length = ray_num * MAX_DEPTH + 1;
	photon_array = new photon[max_array_length];
	generate_global_photons(world, light);
	std::cout << "Constructing photon map...\n";
	root = balance(1, photon_num);
	std::cout << "Global photon map constructed\n";
}

void PhotonMap::generate_global_photons(const Hittable& world, XZ_Rect* light) {
	int ray_cnt = 1;
	while (has_space() && ray_cnt <= ray_num) {
		std::cerr << "\rPhotons remaining: " << ray_num - ray_cnt << ' ' << std::flush;
		photon_num += generate_one_global_photon(world, light);
		++ray_cnt;
	}
	std::cout << "\nSuccessfully generated " << photon_num << " photons\n \n";
}

int PhotonMap::generate_one_global_photon(const Hittable& world, XZ_Rect* light) {
	HitRecord hrec;
	ScatterRecord srec;
	Ray ray;
	int depth = MAX_DEPTH;
	int photon_cnt = 0;
	Color energy = initial_energy; // initialize energy

	do {
		ray = light->generate_emitted_ray(); // emit a photon
	} while (!world.hit(ray, 1e-4, infinity, hrec)); // re-emit if didn't hit anything

	while (depth > 0 && world.hit(ray, 1e-4, infinity, hrec)) {
		if (!hrec.mat_ptr->scatter(ray, hrec, srec)) { // absorbed if hits light
			break;
		}
		if (srec.is_specular) { // scatter if hits specular surface
			ray = srec.scattered;
			energy = energy * srec.albedo; // modify energy
			--depth;
		}
		else { // if hits diffuse surface
			// store photon
			if (has_space()) {
				photon_array[idx++] = photon(hrec.p, ray.direction(), energy);
				++photon_cnt;
			}
			else {
				break;
			}
			ray = srec.scattered;
			energy = srec.albedo * energy; // modify energy
			--depth;
		}
	}
	return photon_cnt;
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

Color PhotonMap::radiance_estimate(Point3 p, Color BRDF, Vec3& elps, Vec3 normal) const { // need to be const
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
		r = max_search_radius;
	}
	while (!que.empty()) {
		p_pht = que.top();
		radiance += p_pht->color * (dot(p_pht->dir, normal) > 0 ? 0 : BRDF); // ignore radiance from the back of the surface
		que.pop();
	}
	radiance = radiance / (PI * r * r);
	return radiance;
}

Color PhotonMap::show_photon(Point3 p, Vec3& elps, Vec3 normal) {
	std::priority_queue<photon*, std::vector<photon*>, photon_ptr_compare> que;
	locate_photons(root, p, elps, que, normal);
	Color radiance = Color();
	Color white = Color(1, 1, 1);
	double ratio = 5000 / static_cast<double>(ray_num); // careful of 2 integers, the result is 0 without explicit cast
	maxphotons = std::max(maxphotons, que.size());
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

CausticsPhotonMap::CausticsPhotonMap(int ray_num_, double search_radius_, Color init_energy_, const Hittable& world, XZ_Rect* light, const Hittable* specular_objects)
	: PhotonMap(ray_num_, search_radius_, init_energy_)
{
	max_array_length = ray_num + 1;
	photon_array = new photon[max_array_length];
	generate_caustics_photons(world, light, specular_objects);
	std::cout << "Constructing caustics photon map...\n";
	root = balance(1, photon_num); // in caustics map, photon_num <= ray_num if nothing is wrong
	std::cout << "Caustics photon map constructed\n";
}

Color CausticsPhotonMap::radiance_estimate(Point3 p, Color BRDF, Vec3& elps, Vec3 normal) const {
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
		r = max_search_radius;
	}
	while (!que.empty()) {
		p_pht = que.top();
		radiance += p_pht->color * (dot(p_pht->dir, normal) > 0 ? 0 : BRDF) * 3 * (1 - (p_pht->p - p).length() / r); // ignore radiance from the back of the surface
		// radiance += p_pht->color * (dot(p_pht->dir, normal) > 0 ? 0 : BRDF);
		que.pop();
	}
	radiance = radiance / (PI * r * r);
	return radiance;
}

void CausticsPhotonMap::generate_caustics_photons(const Hittable& world, XZ_Rect* light, const Hittable* specular_objects) {
	int ray_cnt = 1;
	int generated_num = 0;
	while (has_space() && ray_cnt <= ray_num) {
		std::cerr << "\rCaustics photons remaining: " << ray_num - ray_cnt << ' ' << std::flush;
		generated_num = generate_one_caustics_photon(world, light, specular_objects);
		if (generated_num > 0) {
			photon_num += generated_num;
			++ray_cnt;
		}
	}
	std::cout << "\nSuccessfully generated " << photon_num << " caustics photons\n \n";
}

int CausticsPhotonMap::generate_one_caustics_photon(const Hittable& world, const XZ_Rect* light, const Hittable* specular_objects) {
	HitRecord hrec;
	ScatterRecord srec;
	Point3 random_p = light->random_point();
	Ray ray = Ray(random_p, specular_objects->generate_random(random_p));
	int depth = MAX_DEPTH;
	int photon_cnt = 0;
	Color energy = initial_energy;
	bool hits_specular_first = false;

	do {
		if (!world.hit(ray, 1e-4, infinity, hrec)) { // re-emit if didn't hit anything
			break;
		}
		if (!hrec.mat_ptr->scatter(ray, hrec, srec)) { // re-emit if hits light
			break;
		}
		if (srec.is_specular) { // scatter if hits specular surface
			ray = srec.scattered;
			energy = energy * srec.albedo;
			if (!hits_specular_first) { hits_specular_first = true; }
			--depth;
		}
		else { // if hits diffuse surface
			if (hits_specular_first && has_space()) { // if hits specular first
				photon_array[idx++] = photon(hrec.p, ray.direction(), energy); // store photon
				++photon_cnt;
				break;
			}
			else { // re-emit
				break;
			}
		}
	} while (depth > 0);
	return photon_cnt;
}

// Indirect Photon Map

IndirectPhotonMap::IndirectPhotonMap(int ray_num_, double search_radius_, Color init_energy_, const Hittable& world, XZ_Rect* light)
	:PhotonMap(ray_num_, search_radius_, init_energy_)
{
	max_array_length = ray_num * MAX_DEPTH + 1;
	photon_array = new photon[max_array_length];
	generate_indirect_photons(world, light);
	std::cout << "Constructing indirect photon map...\n";
	root = balance(1, photon_num);
	std::cout << "Indirect photon map constructed\n \n";
}

void IndirectPhotonMap::generate_indirect_photons(const Hittable& world, XZ_Rect* light)
{
	int ray_cnt = 1;
	int generated_num;
	while (has_space() && ray_cnt <= ray_num) {
		std::cerr << "\rIndirect photons remaining: " << ray_num - ray_cnt << ' ' << std::flush;
		generated_num = generate_one_indirect_photon(world, light);
		if (generated_num > 0) {
			photon_num += generated_num;
			++ray_cnt;
		}
	}
	std::cout << "\nSuccessfully generated " << photon_num << " indirect photons\n";
}

int IndirectPhotonMap::generate_one_indirect_photon(const Hittable& world, XZ_Rect* light)
{
	HitRecord hrec;
	ScatterRecord srec;
	Ray ray = light->generate_emitted_ray();
	int depth = MAX_DEPTH;
	int photon_cnt = 0;
	Color energy = initial_energy;
	bool hits_diffuse_first = false;

	do {
		ray = light->generate_emitted_ray(); // emit a photon
	} while (!world.hit(ray, 1e-4, infinity, hrec)); // re-emit if didn't hit anything

	while (depth > 0 && world.hit(ray, 1e-4, infinity, hrec)) {
		if (!hrec.mat_ptr->scatter(ray, hrec, srec)) { // absorbed if hits light
			break;
		}
		if (srec.is_specular) { // scatter if hits specular surface
			if (depth == MAX_DEPTH) { return 0; } // re-emit if first hits specular surface, this part should be contained in caustics map
			ray = srec.scattered;
			energy = energy * srec.albedo; // modify energy
			--depth;
		}
		else { // if hits diffuse surface
			if (hits_diffuse_first && has_space()) {
				photon_array[idx++] = photon(hrec.p, ray.direction(), energy); // store photon
				++photon_cnt;
			}
			else {
				hits_diffuse_first = true;
			}
			ray = srec.scattered;
			energy = srec.albedo * energy; // modify energy
			--depth;
		}
	}
	return photon_cnt;
}

PhotonRatioSampler::PhotonRatioSampler(int ray_num_, const Hittable& world, XZ_Rect* light)
	: ray_num(ray_num_), spec_num(0), non_spec_num(0), ratio(0)
{
	int ray_cnt = 1;
	int generated_num = 0;
	while (ray_cnt <= ray_num) {
		std::cerr << "\rTest photons remaining: " << ray_num - ray_cnt << ' ' << std::flush;
		generated_num = generate_one_global_photon(world, light);
		if (generated_num > 0) {
			++ray_cnt;
		}
	}
	ratio = static_cast<double>(non_spec_num) / spec_num; // compute ratio
	std::cout << "\nRatio of non-specular photons to specular photons is: " << ratio << "\n \n";
}

int PhotonRatioSampler::generate_one_global_photon(const Hittable& world, XZ_Rect* light)
{
	HitRecord hrec;
	ScatterRecord srec;
	Ray ray;
	
	ray = light->generate_emitted_ray();

	if (!world.hit(ray, 1e-4, infinity, hrec)) { // if didn't hit
		return 0;
	}
	if (!hrec.mat_ptr->scatter(ray, hrec, srec)) { // absorbed if hits light
		return 0;
	}
	if (srec.is_specular) { // scatter if hits specular surface
		++spec_num;
	}
	else { // if hits diffuse surface
		++non_spec_num;
	}
	return 1;
}
