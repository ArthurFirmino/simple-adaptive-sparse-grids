#ifndef SASG_H
#define SASG_H

#include <functional>
#include <utility>
#include <array>
#include <vector>
#include <unordered_map>
#include <iostream>
#include <queue>
#include <map>

#include <assert.h>
#define ASSERT assert

namespace SASG {

	using std::array;
	using Float = double;
	using itype = uint16_t;

	template<int d>
	class li_index {
	public:
		array<itype, d> level, index;

	public:

		li_index() {
			level.fill(1);
			index.fill(1);
		}
		li_index(const array<itype, d>& l, const array<itype, d>& i) : level(l), index(i) {
		}
		li_index(const li_index& li) : level(li.level), index(li.index) {
		}
		friend std::ostream& operator<<(std::ostream& os, const li_index& li) {
			os << '{';
			for (int j{ 0 }; j < d - 1; ++j) os << li.level[j] << ',';
			os << li.level[d - 1] << "} {";
			for (int j{ 0 }; j < d - 1; ++j) os << li.index[j] << ',';
			os << li.index[d - 1] << '}';
			return os;
		}
		int N() const {
			int ret{ 1-d };
			for (auto l : level) ret += l;
			return ret;
		}
		array<Float, d> x() const {
			array<Float, d> ret;
			for (int j{ 0 }; j < d; ++j) ret[j] = static_cast<Float>(index[j]) / static_cast<Float>(1 << level[j]);
			return ret;
		}

		struct iterator : public li_index {
			iterator(const li_index& li) : level_copy(li.level), index_copy(li.index), is_valid(true) {
				level = li.level;
				index = li.index;
			}
			operator bool() const {
				return is_valid;
			}
		protected:
			array<itype, d> level_copy, index_copy;
			bool is_valid;
		};

		struct da_iterator : public iterator {
			da_iterator(const li_index& li) : iterator(li), j(0) {
				while (j < d && li_index::level[j] <= 1) ++j;
				if (j >= d) {
					iterator::is_valid = false;
				}
				else {
					--li_index::level[j];
					if ((li_index::index[j] - 1) % 4) lr = 1;
					else lr = -1;
					li_index::index[j] = (li_index::index[j] - lr) / 2;
				}
			}
			da_iterator& operator++() {
				if (j < d - 1) {
					++li_index::level[j];
					li_index::index[j] = 2 * li_index::index[j] + lr;
					while (++j < d && li_index::level[j] <= 1);
					if (j >= d) {
						iterator::is_valid = false;
						return *this;
					}
					--li_index::level[j];
					if ((li_index::index[j] - 1) % 4) lr = 1;
					else lr = -1;
					li_index::index[j] = (li_index::index[j] - lr) / 2;
				}
				else {
					iterator::is_valid = false;
				}
				return *this;
			}
		private:
			size_t j;
			int lr;
		};

		struct dd_iterator : public iterator {
			dd_iterator(const li_index& li) : iterator(li), j(0), lr(-1) {
				++li_index::level[j];
				li_index::index[j] = 2 * li_index::index[j] - 1;
				lr = 1;
			}
			dd_iterator& operator++() {
				if (lr == -1 && j < d - 1) {
					li_index::index[j] = iterator::index_copy[j];
					--li_index::level[j];
					++li_index::level[++j];
					li_index::index[j] = 2 * li_index::index[j] - 1;
					lr = 1;
				}
				else if (lr == 1) {
					li_index::index[j] += 2;
					lr = -1;
				}
				else {
					iterator::is_valid = false;
				}
				return *this;
			}
		protected:
			size_t j;
			int lr;
		};

		struct fdd_iterator : public dd_iterator {
			fdd_iterator(const li_index& li) : dd_iterator(li), k(0) {
				while (k < d - 1 && li.level[k] == 1) ++k;
			}
			fdd_iterator& operator++() {
				if (dd_iterator::lr == 1 || dd_iterator::j < k) {
					++(static_cast<dd_iterator&>(*this));
				}
				else {
					iterator::is_valid = false;
				}
				return *this;
			}
		private:
			size_t k;
		};

		struct fddx_iterator : public iterator {
			fddx_iterator(const li_index& li, const array<Float,d>& _x) : iterator(li), j(0), k(0), x(_x), o(li.x()) {
				while (k < d - 1 && li.level[k] == 1) ++k;
				++li_index::level[j];
				li_index::index[j] = 2 * li_index::index[j] + (x[j] < o[j] ? -1 : 1);
			}
			fddx_iterator& operator++() {
				if (j < k) {
					li_index::index[j] = iterator::index_copy[j];
					--li_index::level[j];
					++li_index::level[++j];
					li_index::index[j] = 2 * li_index::index[j] + (x[j] < o[j] ? -1 : 1);
				}
				else {
					iterator::is_valid = false;
				}
				return *this;
			}
		private:
			size_t j, k;
			array<Float, d> x, o;
		};

		da_iterator da_begin() const {
			return da_iterator(*this);
		}
		dd_iterator dd_begin() const {
			return dd_iterator(*this);
		}
		fdd_iterator fdd_begin() const {
			return fdd_iterator(*this);
		}
		fddx_iterator fddx_begin(const array<Float,d>& x) const {
			return fddx_iterator(*this, x);
		}
	};

	template<int d>
	bool operator==(const li_index<d>& a, const li_index<d>& b) {
		return a.index == b.index && a.level == b.level;
	}
}

namespace std {
	template<typename T, size_t N>
	struct hash<array<T, N>> {

		typedef array<T, N> argument_type;
		typedef size_t result_type;

		result_type operator()(const argument_type& a) const {
			hash<T> hasher;
			result_type h = 0;
			for (result_type i = 0; i < N; ++i) h = h * 31 + hasher(a[i]);
			return h;
		}
	};
	template<int d>
	struct hash<SASG::li_index<d>> {

		typedef SASG::li_index<d> argument_type;
		typedef size_t result_type;

		result_type operator()(const argument_type& a) const {
			return hash<array<SASG::itype,d>>()(a.level) ^ hash<array<SASG::itype, d>>()(a.index);
		}
	};
}

namespace SASG {

	using std::unordered_map;
	using std::function;
	using std::vector;
	using std::queue;
	using std::map;

	/* V: value type, d: num dimensions */
	template<typename V, int d>
	class Grid {

		using li_index = li_index<d>;

	public:

		Grid(int N) {
			function<void(li_index&)> rec;
			rec = [this, &rec, N](li_index& li) {
				surpluses[li] = V();
				if (li.N() < N) {
					for (auto it = li.fdd_begin(); it; ++it) rec(static_cast<li_index&>(it));
				}
			};
			li_index root;
			rec(root);
		}

		~Grid() {}
		
		void hierarchize(function<V(const array<double, d>&)> f) {
			/*	. Set all a_li to 0
				. Perform level order traversal
				. a_li = f(x_li) - u(x_li)
			*/
			for (auto& p : surpluses) p.second = V(0);

			queue<li_index> li_queue;
			li_queue.push(li_index());

			while (!li_queue.empty()) {
				li_index li = li_queue.front();
				li_queue.pop();
				auto it = surpluses.find(li);
				if (it == surpluses.end()) {
					continue;
				}
				else {
					it->second = f(li.x()) - evaluate(li.x(), li.N()-1);
					auto fdd_it = li.fdd_begin();
					if (surpluses.find(fdd_it) == surpluses.end()) continue;
					while (fdd_it) {
						li_queue.push(static_cast<li_index>(fdd_it));
						++fdd_it;
					}
				}
			}
		}
		
		V evaluate(const array<Float, d>& x, int max_N = 0) {
			/*	. Iterate recursively through fddx (first direct descendants which overlap x)
				. Sum a_li * basis(x) for li whose support overlaps x
			*/
			auto basis1d = [](int32_t l, int32_t i, Float x) {
				if (l == 1 && i == 1) return 1.0;
				else if ( l > 1 && i == 1) return 2.0 - x * (1 << l);
				else if (l > 1 && i == (1 << l) - 1) return x * (1 << l) + 1.0 - i;
				else return 1.0 - std::abs(x * (1 << l) - i);
			};
			auto basis = [&basis1d](const li_index& li, const array<Float,d>& x) {
				Float ret{ 1.0 };
				for (int j{ 0 }; j < d; ++j) {
					ret *= basis1d(li.level[j], li.index[j], x[j]);
				}
				return ret;
			};

			V ret(0);
			function<void(li_index&)> rec;
			rec = [this,&rec,&ret,&x,&basis,max_N](li_index& li) {
				ret += surpluses[li] * basis(li, x);
				if (max_N && li.N() >= max_N) return;
				auto it = li.fddx_begin(x);
				if (surpluses.find(static_cast<const li_index&>(it)) == surpluses.end()) return;
				for ( ; it; ++it) rec(static_cast<li_index&>(it));
			};
			li_index root;
			rec(root);
			return ret;
		}

		int refine(const Float& eps) {
			/* . Refine all leaf nodes for which |a_li| > eps
			   . Returns number of refined nodes
			*/
			vector<li_index> torefine;
			function<void(const li_index&)> rec;
			rec = [this, &rec, &eps, &torefine](const li_index& li) {
				auto it = li.fdd_begin();
				if (surpluses.find(it) == surpluses.end()) {
					if (std::abs(surpluses[li]) > eps)
						torefine.push_back(li);
					return;
				}
				for (; it; ++it) {
					rec(it);
				}
			};
			li_index root;
			rec(root);
			int nref{ 0 };
			for (const auto& li : torefine) {
				if (refine(li)) ++nref;
			}
			return nref;
		}

		int refine(const Float& eps, function<V(const array<double, d>&)> f) {
			/* . Refine all leaf nodes for which |a_li| > eps
			   . Returns number of refined nodes
			*/
			vector<li_index> torefine;
			function<void(const li_index&)> rec;
			rec = [this, &rec, &eps, &torefine](const li_index& li) {
				auto it = li.fdd_begin();
				if (surpluses.find(it) == surpluses.end()) {
					if (std::abs(surpluses[li]) > eps)
						torefine.push_back(li);
					return;
				}
				for (; it; ++it) {
					rec(it);
				}
			};
			li_index root;
			rec(root);
			int nref{ 0 };
			for (const auto& li : torefine) {
				if (refine(li, f)) ++nref;
			}
			return nref;
		}

		int unrefine(const Float& eps) {
			/* . Unrefine nodes whose direct descendants are all leaf nodes,
			     and for which max({a_pq},{a_li}) <= eps, where pq are the direct descendants of li
			   . Returns number of unrefined nodes
			*/
			vector<li_index> tounrefine;
			for (const auto& p : surpluses) {
				if (isleaf(p.first)) continue;
				bool areleaf{ true };
				for (auto it = p.first.dd_begin(); it; ++it) {
					if (!isleaf(it)) {
						areleaf = false;
						break;
					}
				}
				if (areleaf) {
					Float maxeps{ std::abs(p.second) };
					for (auto it = p.first.dd_begin(); it; ++it) {
						ASSERT(surpluses.find(it) != surpluses.end());
						maxeps = std::max(maxeps, surpluses[it]);
					}
					if (maxeps <= eps)
						tounrefine.push_back(p.first);
				}
			}
			int nref{ 0 };
			for (const auto& li : tounrefine) {
				if (unrefine(li)) ++nref;
			}
			return nref;
		}

		bool unrerefine(function<V(const array<double, d>&)> f, const Float& eps, int maxiter = 100) {
			int i{ 0 };
			while (unrefine(eps) | refine(eps, f)) {
				if (++maxiter >= maxiter) return false;
			}
			return true;
		}

		const unordered_map<li_index, V>& getsurpluses() const {
			return surpluses;
		}

	private:

		unordered_map<li_index, V> surpluses;

		bool isleaf(const li_index& li) const {
			if (surpluses.find(li) == surpluses.end()) return false;
			for (auto it = li.dd_begin(); it; ++it)
				if (surpluses.find(it) != surpluses.end()) return false;
			return true;
		}

		int height(const li_index& li) const {
			if (surpluses.find(li) == surpluses.end()) return -1;
			if (isleaf(li)) return 0;
			int ret{ 0 };
			for (auto it = li.dd_begin(); it; ++it) {
				ret = std::max(ret, height(static_cast<const li_index&>(it)));
			}
			return ret + 1;
		}

		bool refine(const li_index& li) {
			/* . Assert li is a leaf node
			   . Create all direct descendants of li
			   . For each direct descendant recursively create any missing direct ancestors
			*/
			if (height(li) != 0) return false;

			function<void(const li_index&)> rec;
			rec = [this, &rec](const li_index& pq) {
				ASSERT(pq.N() >= 1);
				surpluses[static_cast<const li_index&>(pq)] = V(0);
				for (auto it = pq.da_begin(); it; ++it) {
					if (surpluses.find(static_cast<const li_index&>(it)) == surpluses.end())
						rec(static_cast<const li_index&>(it));
				}
			};
			for (auto it = li.dd_begin(); it; ++it) {
				rec(static_cast<const li_index&>(it));
			}
			return true;
		}

		bool refine(const li_index& li, function<V(const array<double, d>&)> f) {
			/* . Assert li is a leaf node
			   . Create all direct descendants of li
			   . For each direct descendant recursively create any missing direct ancestors
			*/
			if (height(li) != 0) return false;

			map<int, vector<li_index>> added;

			function<void(const li_index&)> rec;
			rec = [this, &rec, &added](const li_index& pq) {
				ASSERT(pq.N() >= 1);
				surpluses[static_cast<const li_index&>(pq)] = V(0);
				added[pq.N()].push_back(pq);
				for (auto it = pq.da_begin(); it; ++it) {
					if (surpluses.find(static_cast<const li_index&>(it)) == surpluses.end())
						rec(static_cast<const li_index&>(it));
				}
			};
			for (auto it = li.dd_begin(); it; ++it) {
				rec(static_cast<const li_index&>(it));
			}

			for (auto& v : added) {
				for (auto& pq : v.second) {
					surpluses[pq] = f(pq.x()) - evaluate(pq.x(), pq.N() - 1);
				}
			}
			return true;
		}

		bool unrefine(const li_index& li) {
			/* . Assert li is not a leaf and all of its direct descendants are leaf nodes
			   . Remove all direct descendants of li
			*/
			if (height(li) != 1) return false;

			for (auto it = li.dd_begin(); it; ++it) {
				surpluses.erase(static_cast<const li_index&>(it));
			}
			return true;
		}
	};
}

#endif