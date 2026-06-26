#include <Pythia8/Basics.h>
#include <Pythia8/BeamSetup.h>
#include <Pythia8/Event.h>
#include <Pythia8/HINucleusModel.h>
#include <Pythia8/HISubCollisionModel.h>
#include <Pythia8/HadronWidths.h>
#include <Pythia8/Info.h>
#include <Pythia8/LHEF3.h>
#include <Pythia8/Logger.h>
#include <Pythia8/ParticleData.h>
#include <Pythia8/PartonSystems.h>
#include <Pythia8/Settings.h>
#include <Pythia8/SigmaLowEnergy.h>
#include <Pythia8/SigmaTotal.h>
#include <Pythia8/StandardModel.h>
#include <Pythia8/SusyCouplings.h>
#include <Pythia8/Weights.h>
#include <functional>
#include <istream>
#include <iterator>
#include <map>
#include <memory>
#include <ostream>
#include <set>
#include <sstream>
#include <sstream> // __str__
#include <string>
#include <utility>
#include <vector>

#include <pybind11/pybind11.h>
#include <functional>
#include <string>
#include <Pythia8/UserHooks.h>
#include <Pythia8/SplittingsOnia.h>
#include <Pythia8/HeavyIons.h>
#include <Pythia8/BeamShape.h>
#include <pybind11/stl.h>
#include <pybind11/complex.h>
#include <pybind11/functional.h>


#ifndef BINDER_PYBIND11_TYPE_CASTER
	#define BINDER_PYBIND11_TYPE_CASTER
	PYBIND11_DECLARE_HOLDER_TYPE(T, std::shared_ptr<T>);
	PYBIND11_DECLARE_HOLDER_TYPE(T, T*);
	PYBIND11_MAKE_OPAQUE(std::shared_ptr<void>);
#endif

// Pythia8::SubCollisionModel file:Pythia8/HISubCollisionModel.h line:143
struct PyCallBack_Pythia8_SubCollisionModel : public Pythia8::SubCollisionModel {
	using Pythia8::SubCollisionModel::SubCollisionModel;

	bool init(int a0, int a1, double a2) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::SubCollisionModel *>(this), "init");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return SubCollisionModel::init(a0, a1, a2);
	}
	using _binder_ret_0 = class std::vector<double, class std::allocator<double> >;
	_binder_ret_0 minParm() const override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::SubCollisionModel *>(this), "minParm");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<_binder_ret_0>::value) {
				static pybind11::detail::override_caster_t<_binder_ret_0> caster;
				return pybind11::detail::cast_ref<_binder_ret_0>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<_binder_ret_0>(std::move(o));
		}
		pybind11::pybind11_fail("Tried to call pure virtual function \"SubCollisionModel::minParm\"");
	}
	using _binder_ret_1 = class std::vector<double, class std::allocator<double> >;
	_binder_ret_1 defParm() const override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::SubCollisionModel *>(this), "defParm");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<_binder_ret_1>::value) {
				static pybind11::detail::override_caster_t<_binder_ret_1> caster;
				return pybind11::detail::cast_ref<_binder_ret_1>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<_binder_ret_1>(std::move(o));
		}
		pybind11::pybind11_fail("Tried to call pure virtual function \"SubCollisionModel::defParm\"");
	}
	using _binder_ret_2 = class std::vector<double, class std::allocator<double> >;
	_binder_ret_2 maxParm() const override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::SubCollisionModel *>(this), "maxParm");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<_binder_ret_2>::value) {
				static pybind11::detail::override_caster_t<_binder_ret_2> caster;
				return pybind11::detail::cast_ref<_binder_ret_2>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<_binder_ret_2>(std::move(o));
		}
		pybind11::pybind11_fail("Tried to call pure virtual function \"SubCollisionModel::maxParm\"");
	}
	void generateNucleonStates(class Pythia8::Nucleus & a0, class Pythia8::Nucleus & a1) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::SubCollisionModel *>(this), "generateNucleonStates");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return SubCollisionModel::generateNucleonStates(a0, a1);
	}
	class Pythia8::SubCollisionSet getCollisions(class Pythia8::Nucleus & a0, class Pythia8::Nucleus & a1) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::SubCollisionModel *>(this), "getCollisions");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1);
			if (pybind11::detail::cast_is_temporary_value_reference<class Pythia8::SubCollisionSet>::value) {
				static pybind11::detail::override_caster_t<class Pythia8::SubCollisionSet> caster;
				return pybind11::detail::cast_ref<class Pythia8::SubCollisionSet>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<class Pythia8::SubCollisionSet>(std::move(o));
		}
		pybind11::pybind11_fail("Tried to call pure virtual function \"SubCollisionModel::getCollisions\"");
	}
	struct Pythia8::SubCollisionModel::SigEst getSig() const override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::SubCollisionModel *>(this), "getSig");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<struct Pythia8::SubCollisionModel::SigEst>::value) {
				static pybind11::detail::override_caster_t<struct Pythia8::SubCollisionModel::SigEst> caster;
				return pybind11::detail::cast_ref<struct Pythia8::SubCollisionModel::SigEst>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<struct Pythia8::SubCollisionModel::SigEst>(std::move(o));
		}
		pybind11::pybind11_fail("Tried to call pure virtual function \"SubCollisionModel::getSig\"");
	}
};

// Pythia8::BlackSubCollisionModel file:Pythia8/HISubCollisionModel.h line:484
struct PyCallBack_Pythia8_BlackSubCollisionModel : public Pythia8::BlackSubCollisionModel {
	using Pythia8::BlackSubCollisionModel::BlackSubCollisionModel;

	using _binder_ret_0 = class std::vector<double, class std::allocator<double> >;
	_binder_ret_0 minParm() const override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::BlackSubCollisionModel *>(this), "minParm");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<_binder_ret_0>::value) {
				static pybind11::detail::override_caster_t<_binder_ret_0> caster;
				return pybind11::detail::cast_ref<_binder_ret_0>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<_binder_ret_0>(std::move(o));
		}
		return BlackSubCollisionModel::minParm();
	}
	using _binder_ret_1 = class std::vector<double, class std::allocator<double> >;
	_binder_ret_1 defParm() const override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::BlackSubCollisionModel *>(this), "defParm");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<_binder_ret_1>::value) {
				static pybind11::detail::override_caster_t<_binder_ret_1> caster;
				return pybind11::detail::cast_ref<_binder_ret_1>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<_binder_ret_1>(std::move(o));
		}
		return BlackSubCollisionModel::defParm();
	}
	using _binder_ret_2 = class std::vector<double, class std::allocator<double> >;
	_binder_ret_2 maxParm() const override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::BlackSubCollisionModel *>(this), "maxParm");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<_binder_ret_2>::value) {
				static pybind11::detail::override_caster_t<_binder_ret_2> caster;
				return pybind11::detail::cast_ref<_binder_ret_2>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<_binder_ret_2>(std::move(o));
		}
		return BlackSubCollisionModel::maxParm();
	}
	struct Pythia8::SubCollisionModel::SigEst getSig() const override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::BlackSubCollisionModel *>(this), "getSig");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<struct Pythia8::SubCollisionModel::SigEst>::value) {
				static pybind11::detail::override_caster_t<struct Pythia8::SubCollisionModel::SigEst> caster;
				return pybind11::detail::cast_ref<struct Pythia8::SubCollisionModel::SigEst>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<struct Pythia8::SubCollisionModel::SigEst>(std::move(o));
		}
		return BlackSubCollisionModel::getSig();
	}
	class Pythia8::SubCollisionSet getCollisions(class Pythia8::Nucleus & a0, class Pythia8::Nucleus & a1) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::BlackSubCollisionModel *>(this), "getCollisions");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1);
			if (pybind11::detail::cast_is_temporary_value_reference<class Pythia8::SubCollisionSet>::value) {
				static pybind11::detail::override_caster_t<class Pythia8::SubCollisionSet> caster;
				return pybind11::detail::cast_ref<class Pythia8::SubCollisionSet>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<class Pythia8::SubCollisionSet>(std::move(o));
		}
		return BlackSubCollisionModel::getCollisions(a0, a1);
	}
	bool init(int a0, int a1, double a2) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::BlackSubCollisionModel *>(this), "init");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return SubCollisionModel::init(a0, a1, a2);
	}
	void generateNucleonStates(class Pythia8::Nucleus & a0, class Pythia8::Nucleus & a1) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::BlackSubCollisionModel *>(this), "generateNucleonStates");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return SubCollisionModel::generateNucleonStates(a0, a1);
	}
};

// Pythia8::NaiveSubCollisionModel file:Pythia8/HISubCollisionModel.h line:520
struct PyCallBack_Pythia8_NaiveSubCollisionModel : public Pythia8::NaiveSubCollisionModel {
	using Pythia8::NaiveSubCollisionModel::NaiveSubCollisionModel;

	using _binder_ret_0 = class std::vector<double, class std::allocator<double> >;
	_binder_ret_0 minParm() const override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::NaiveSubCollisionModel *>(this), "minParm");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<_binder_ret_0>::value) {
				static pybind11::detail::override_caster_t<_binder_ret_0> caster;
				return pybind11::detail::cast_ref<_binder_ret_0>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<_binder_ret_0>(std::move(o));
		}
		return NaiveSubCollisionModel::minParm();
	}
	using _binder_ret_1 = class std::vector<double, class std::allocator<double> >;
	_binder_ret_1 defParm() const override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::NaiveSubCollisionModel *>(this), "defParm");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<_binder_ret_1>::value) {
				static pybind11::detail::override_caster_t<_binder_ret_1> caster;
				return pybind11::detail::cast_ref<_binder_ret_1>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<_binder_ret_1>(std::move(o));
		}
		return NaiveSubCollisionModel::defParm();
	}
	using _binder_ret_2 = class std::vector<double, class std::allocator<double> >;
	_binder_ret_2 maxParm() const override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::NaiveSubCollisionModel *>(this), "maxParm");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<_binder_ret_2>::value) {
				static pybind11::detail::override_caster_t<_binder_ret_2> caster;
				return pybind11::detail::cast_ref<_binder_ret_2>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<_binder_ret_2>(std::move(o));
		}
		return NaiveSubCollisionModel::maxParm();
	}
	struct Pythia8::SubCollisionModel::SigEst getSig() const override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::NaiveSubCollisionModel *>(this), "getSig");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<struct Pythia8::SubCollisionModel::SigEst>::value) {
				static pybind11::detail::override_caster_t<struct Pythia8::SubCollisionModel::SigEst> caster;
				return pybind11::detail::cast_ref<struct Pythia8::SubCollisionModel::SigEst>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<struct Pythia8::SubCollisionModel::SigEst>(std::move(o));
		}
		return NaiveSubCollisionModel::getSig();
	}
	class Pythia8::SubCollisionSet getCollisions(class Pythia8::Nucleus & a0, class Pythia8::Nucleus & a1) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::NaiveSubCollisionModel *>(this), "getCollisions");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1);
			if (pybind11::detail::cast_is_temporary_value_reference<class Pythia8::SubCollisionSet>::value) {
				static pybind11::detail::override_caster_t<class Pythia8::SubCollisionSet> caster;
				return pybind11::detail::cast_ref<class Pythia8::SubCollisionSet>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<class Pythia8::SubCollisionSet>(std::move(o));
		}
		return NaiveSubCollisionModel::getCollisions(a0, a1);
	}
	bool init(int a0, int a1, double a2) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::NaiveSubCollisionModel *>(this), "init");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return SubCollisionModel::init(a0, a1, a2);
	}
	void generateNucleonStates(class Pythia8::Nucleus & a0, class Pythia8::Nucleus & a1) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::NaiveSubCollisionModel *>(this), "generateNucleonStates");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return SubCollisionModel::generateNucleonStates(a0, a1);
	}
};

void bind_Pythia8_HISubCollisionModel(std::function< pybind11::module &(std::string const &namespace_) > &M)
{
	{ // Pythia8::SubCollisionModel file:Pythia8/HISubCollisionModel.h line:143
		pybind11::class_<Pythia8::SubCollisionModel, std::shared_ptr<Pythia8::SubCollisionModel>, PyCallBack_Pythia8_SubCollisionModel> cl(M("Pythia8"), "SubCollisionModel", "");
		pybind11::handle cl_type = cl;

		{ // Pythia8::SubCollisionModel::SigEst file:Pythia8/HISubCollisionModel.h line:148
			auto & enclosing_class = cl;
			pybind11::class_<Pythia8::SubCollisionModel::SigEst, std::shared_ptr<Pythia8::SubCollisionModel::SigEst>> cl(enclosing_class, "SigEst", "");
			pybind11::handle cl_type = cl;

			cl.def( pybind11::init( [](){ return new Pythia8::SubCollisionModel::SigEst(); } ) );
			cl.def( pybind11::init( [](Pythia8::SubCollisionModel::SigEst const &o){ return new Pythia8::SubCollisionModel::SigEst(o); } ) );
			cl.def_readwrite("sig", &Pythia8::SubCollisionModel::SigEst::sig);
			cl.def_readwrite("dsig2", &Pythia8::SubCollisionModel::SigEst::dsig2);
			cl.def_readwrite("fsig", &Pythia8::SubCollisionModel::SigEst::fsig);
			cl.def_readwrite("avNDb", &Pythia8::SubCollisionModel::SigEst::avNDb);
			cl.def_readwrite("davNDb2", &Pythia8::SubCollisionModel::SigEst::davNDb2);
			cl.def_readwrite("avNOF", &Pythia8::SubCollisionModel::SigEst::avNOF);
			cl.def_readwrite("davNOF2", &Pythia8::SubCollisionModel::SigEst::davNOF2);
		}

		{ // Pythia8::SubCollisionModel::SigmaCache file:Pythia8/HISubCollisionModel.h line:382
			auto & enclosing_class = cl;
			pybind11::class_<Pythia8::SubCollisionModel::SigmaCache, std::shared_ptr<Pythia8::SubCollisionModel::SigmaCache>> cl(enclosing_class, "SigmaCache", "");
			pybind11::handle cl_type = cl;

			cl.def( pybind11::init( [](){ return new Pythia8::SubCollisionModel::SigmaCache(); } ) );
			cl.def( pybind11::init( [](Pythia8::SubCollisionModel::SigmaCache const &o){ return new Pythia8::SubCollisionModel::SigmaCache(o); } ) );
			cl.def_readwrite("cache", &Pythia8::SubCollisionModel::SigmaCache::cache);
			cl.def_readwrite("eCMStep", &Pythia8::SubCollisionModel::SigmaCache::eCMStep);
			cl.def("set", (void (Pythia8::SubCollisionModel::SigmaCache::*)(class std::vector<class std::vector<double, class std::allocator<double> >, class std::allocator<class std::vector<double, class std::allocator<double> > > > &, int, int, double)) &Pythia8::SubCollisionModel::SigmaCache::set, "C++: Pythia8::SubCollisionModel::SigmaCache::set(class std::vector<class std::vector<double, class std::allocator<double> >, class std::allocator<class std::vector<double, class std::allocator<double> > > > &, int, int, double) --> void", pybind11::arg("sigNN"), pybind11::arg("idAIn"), pybind11::arg("idBIn"), pybind11::arg("eCM"));
			cl.def("fillCache", (void (Pythia8::SubCollisionModel::SigmaCache::*)(class std::vector<class std::vector<double, class std::allocator<double> >, class std::allocator<class std::vector<double, class std::allocator<double> > > > &, int, int, int)) &Pythia8::SubCollisionModel::SigmaCache::fillCache, "C++: Pythia8::SubCollisionModel::SigmaCache::fillCache(class std::vector<class std::vector<double, class std::allocator<double> >, class std::allocator<class std::vector<double, class std::allocator<double> > > > &, int, int, int) --> void", pybind11::arg("sigNN"), pybind11::arg("idAIn"), pybind11::arg("idBIn"), pybind11::arg("iECM"));
			cl.def("assign", (struct Pythia8::SubCollisionModel::SigmaCache & (Pythia8::SubCollisionModel::SigmaCache::*)(const struct Pythia8::SubCollisionModel::SigmaCache &)) &Pythia8::SubCollisionModel::SigmaCache::operator=, "C++: Pythia8::SubCollisionModel::SigmaCache::operator=(const struct Pythia8::SubCollisionModel::SigmaCache &) --> struct Pythia8::SubCollisionModel::SigmaCache &", pybind11::return_value_policy::reference, pybind11::arg(""));
		}

		cl.def( pybind11::init<int>(), pybind11::arg("nParm") );

		cl.def(pybind11::init<PyCallBack_Pythia8_SubCollisionModel const &>());
		cl.def_readwrite("sigTarg", &Pythia8::SubCollisionModel::sigTarg);
		cl.def_readwrite("sigErr", &Pythia8::SubCollisionModel::sigErr);
		cl.def_readwrite("sigTargNN", &Pythia8::SubCollisionModel::sigTargNN);
		cl.def_readwrite("parmSave", &Pythia8::SubCollisionModel::parmSave);
		cl.def_readwrite("NInt", &Pythia8::SubCollisionModel::NInt);
		cl.def_readwrite("NPop", &Pythia8::SubCollisionModel::NPop);
		cl.def_readwrite("sigFuzz", &Pythia8::SubCollisionModel::sigFuzz);
		cl.def_readwrite("impactFudge", &Pythia8::SubCollisionModel::impactFudge);
		cl.def_readwrite("fitPrint", &Pythia8::SubCollisionModel::fitPrint);
		cl.def_readwrite("eCMlow", &Pythia8::SubCollisionModel::eCMlow);
		cl.def_readwrite("avNDb", &Pythia8::SubCollisionModel::avNDb);
		cl.def_readwrite("avNDolap", &Pythia8::SubCollisionModel::avNDolap);
		cl.def_readwrite("idASave", &Pythia8::SubCollisionModel::idASave);
		cl.def_readwrite("idBSave", &Pythia8::SubCollisionModel::idBSave);
		cl.def_readwrite("doVarECM", &Pythia8::SubCollisionModel::doVarECM);
		cl.def_readwrite("doVarBeams", &Pythia8::SubCollisionModel::doVarBeams);
		cl.def_readwrite("eMin", &Pythia8::SubCollisionModel::eMin);
		cl.def_readwrite("eMax", &Pythia8::SubCollisionModel::eMax);
		cl.def_readwrite("eSave", &Pythia8::SubCollisionModel::eSave);
		cl.def_readwrite("eCMPts", &Pythia8::SubCollisionModel::eCMPts);
		cl.def_readwrite("idAList", &Pythia8::SubCollisionModel::idAList);
		cl.def_readwrite("subCollParmsMap", &Pythia8::SubCollisionModel::subCollParmsMap);
		cl.def_readwrite("elasticMode", &Pythia8::SubCollisionModel::elasticMode);
		cl.def_readwrite("elasticFudge", &Pythia8::SubCollisionModel::elasticFudge);
		cl.def_readwrite("lowEnergyCache", &Pythia8::SubCollisionModel::lowEnergyCache);
		cl.def_static("create", (class std::shared_ptr<class Pythia8::SubCollisionModel> (*)(int)) &Pythia8::SubCollisionModel::create, "C++: Pythia8::SubCollisionModel::create(int) --> class std::shared_ptr<class Pythia8::SubCollisionModel>", pybind11::arg("model"));
		cl.def("init", (bool (Pythia8::SubCollisionModel::*)(int, int, double)) &Pythia8::SubCollisionModel::init, "C++: Pythia8::SubCollisionModel::init(int, int, double) --> bool", pybind11::arg("idAIn"), pybind11::arg("idBIn"), pybind11::arg("eCMIn"));
		cl.def("initPtr", (void (Pythia8::SubCollisionModel::*)(class Pythia8::NucleusModel &, class Pythia8::NucleusModel &, class Pythia8::SigmaTotal &, class Pythia8::Settings &, class Pythia8::Info &, class Pythia8::Rndm &)) &Pythia8::SubCollisionModel::initPtr, "C++: Pythia8::SubCollisionModel::initPtr(class Pythia8::NucleusModel &, class Pythia8::NucleusModel &, class Pythia8::SigmaTotal &, class Pythia8::Settings &, class Pythia8::Info &, class Pythia8::Rndm &) --> void", pybind11::arg("projIn"), pybind11::arg("targIn"), pybind11::arg("sigTotIn"), pybind11::arg("settingsIn"), pybind11::arg("infoIn"), pybind11::arg("rndmIn"));
		cl.def("initLowEnergy", (void (Pythia8::SubCollisionModel::*)(class Pythia8::SigmaCombined *)) &Pythia8::SubCollisionModel::initLowEnergy, "C++: Pythia8::SubCollisionModel::initLowEnergy(class Pythia8::SigmaCombined *) --> void", pybind11::arg("sigmaCombPtrIn"));
		cl.def("hasXSec", (bool (Pythia8::SubCollisionModel::*)() const) &Pythia8::SubCollisionModel::hasXSec, "C++: Pythia8::SubCollisionModel::hasXSec() const --> bool");
		cl.def("sigTot", (double (Pythia8::SubCollisionModel::*)() const) &Pythia8::SubCollisionModel::sigTot, "C++: Pythia8::SubCollisionModel::sigTot() const --> double");
		cl.def("sigEl", (double (Pythia8::SubCollisionModel::*)() const) &Pythia8::SubCollisionModel::sigEl, "C++: Pythia8::SubCollisionModel::sigEl() const --> double");
		cl.def("sigCDE", (double (Pythia8::SubCollisionModel::*)() const) &Pythia8::SubCollisionModel::sigCDE, "C++: Pythia8::SubCollisionModel::sigCDE() const --> double");
		cl.def("sigSDE", (double (Pythia8::SubCollisionModel::*)() const) &Pythia8::SubCollisionModel::sigSDE, "C++: Pythia8::SubCollisionModel::sigSDE() const --> double");
		cl.def("sigSDEP", (double (Pythia8::SubCollisionModel::*)() const) &Pythia8::SubCollisionModel::sigSDEP, "C++: Pythia8::SubCollisionModel::sigSDEP() const --> double");
		cl.def("sigSDET", (double (Pythia8::SubCollisionModel::*)() const) &Pythia8::SubCollisionModel::sigSDET, "C++: Pythia8::SubCollisionModel::sigSDET() const --> double");
		cl.def("sigDDE", (double (Pythia8::SubCollisionModel::*)() const) &Pythia8::SubCollisionModel::sigDDE, "C++: Pythia8::SubCollisionModel::sigDDE() const --> double");
		cl.def("sigND", (double (Pythia8::SubCollisionModel::*)() const) &Pythia8::SubCollisionModel::sigND, "C++: Pythia8::SubCollisionModel::sigND() const --> double");
		cl.def("sigLow", (double (Pythia8::SubCollisionModel::*)() const) &Pythia8::SubCollisionModel::sigLow, "C++: Pythia8::SubCollisionModel::sigLow() const --> double");
		cl.def("sigLExc", (double (Pythia8::SubCollisionModel::*)() const) &Pythia8::SubCollisionModel::sigLExc, "C++: Pythia8::SubCollisionModel::sigLExc() const --> double");
		cl.def("sigLAnn", (double (Pythia8::SubCollisionModel::*)() const) &Pythia8::SubCollisionModel::sigLAnn, "C++: Pythia8::SubCollisionModel::sigLAnn() const --> double");
		cl.def("sigLRes", (double (Pythia8::SubCollisionModel::*)() const) &Pythia8::SubCollisionModel::sigLRes, "C++: Pythia8::SubCollisionModel::sigLRes() const --> double");
		cl.def("bSlope", (double (Pythia8::SubCollisionModel::*)() const) &Pythia8::SubCollisionModel::bSlope, "C++: Pythia8::SubCollisionModel::bSlope() const --> double");
		cl.def("avNDB", (double (Pythia8::SubCollisionModel::*)() const) &Pythia8::SubCollisionModel::avNDB, "C++: Pythia8::SubCollisionModel::avNDB() const --> double");
		cl.def("updateSig", (void (Pythia8::SubCollisionModel::*)(int, int, double)) &Pythia8::SubCollisionModel::updateSig, "C++: Pythia8::SubCollisionModel::updateSig(int, int, double) --> void", pybind11::arg("idAIn"), pybind11::arg("idBIn"), pybind11::arg("eCMIn"));
		cl.def("Chi2", (double (Pythia8::SubCollisionModel::*)(const struct Pythia8::SubCollisionModel::SigEst &, int) const) &Pythia8::SubCollisionModel::Chi2, "C++: Pythia8::SubCollisionModel::Chi2(const struct Pythia8::SubCollisionModel::SigEst &, int) const --> double", pybind11::arg("sigs"), pybind11::arg("npar"));
		cl.def("setKinematics", (bool (Pythia8::SubCollisionModel::*)(double)) &Pythia8::SubCollisionModel::setKinematics, "C++: Pythia8::SubCollisionModel::setKinematics(double) --> bool", pybind11::arg("eCMIn"));
		cl.def("setIDA", (bool (Pythia8::SubCollisionModel::*)(int)) &Pythia8::SubCollisionModel::setIDA, "C++: Pythia8::SubCollisionModel::setIDA(int) --> bool", pybind11::arg("idA"));
		cl.def("evolve", (bool (Pythia8::SubCollisionModel::*)(int, double, int)) &Pythia8::SubCollisionModel::evolve, "C++: Pythia8::SubCollisionModel::evolve(int, double, int) --> bool", pybind11::arg("nGenerations"), pybind11::arg("eCM"), pybind11::arg("idANow"));
		cl.def("nParms", (int (Pythia8::SubCollisionModel::*)() const) &Pythia8::SubCollisionModel::nParms, "C++: Pythia8::SubCollisionModel::nParms() const --> int");
		cl.def("setParm", (void (Pythia8::SubCollisionModel::*)(const class std::vector<double, class std::allocator<double> > &)) &Pythia8::SubCollisionModel::setParm, "C++: Pythia8::SubCollisionModel::setParm(const class std::vector<double, class std::allocator<double> > &) --> void", pybind11::arg("parmIn"));
		cl.def("getParm", (class std::vector<double, class std::allocator<double> > (Pythia8::SubCollisionModel::*)() const) &Pythia8::SubCollisionModel::getParm, "C++: Pythia8::SubCollisionModel::getParm() const --> class std::vector<double, class std::allocator<double> >");
		cl.def("minParm", (class std::vector<double, class std::allocator<double> > (Pythia8::SubCollisionModel::*)() const) &Pythia8::SubCollisionModel::minParm, "C++: Pythia8::SubCollisionModel::minParm() const --> class std::vector<double, class std::allocator<double> >");
		cl.def("defParm", (class std::vector<double, class std::allocator<double> > (Pythia8::SubCollisionModel::*)() const) &Pythia8::SubCollisionModel::defParm, "C++: Pythia8::SubCollisionModel::defParm() const --> class std::vector<double, class std::allocator<double> >");
		cl.def("maxParm", (class std::vector<double, class std::allocator<double> > (Pythia8::SubCollisionModel::*)() const) &Pythia8::SubCollisionModel::maxParm, "C++: Pythia8::SubCollisionModel::maxParm() const --> class std::vector<double, class std::allocator<double> >");
		cl.def("generateNucleonStates", (void (Pythia8::SubCollisionModel::*)(class Pythia8::Nucleus &, class Pythia8::Nucleus &)) &Pythia8::SubCollisionModel::generateNucleonStates, "C++: Pythia8::SubCollisionModel::generateNucleonStates(class Pythia8::Nucleus &, class Pythia8::Nucleus &) --> void", pybind11::arg(""), pybind11::arg(""));
		cl.def("getCollisions", (class Pythia8::SubCollisionSet (Pythia8::SubCollisionModel::*)(class Pythia8::Nucleus &, class Pythia8::Nucleus &)) &Pythia8::SubCollisionModel::getCollisions, "C++: Pythia8::SubCollisionModel::getCollisions(class Pythia8::Nucleus &, class Pythia8::Nucleus &) --> class Pythia8::SubCollisionSet", pybind11::arg("proj"), pybind11::arg("targ"));
		cl.def("getSig", (struct Pythia8::SubCollisionModel::SigEst (Pythia8::SubCollisionModel::*)() const) &Pythia8::SubCollisionModel::getSig, "C++: Pythia8::SubCollisionModel::getSig() const --> struct Pythia8::SubCollisionModel::SigEst");
		cl.def("assign", (class Pythia8::SubCollisionModel & (Pythia8::SubCollisionModel::*)(const class Pythia8::SubCollisionModel &)) &Pythia8::SubCollisionModel::operator=, "C++: Pythia8::SubCollisionModel::operator=(const class Pythia8::SubCollisionModel &) --> class Pythia8::SubCollisionModel &", pybind11::return_value_policy::reference, pybind11::arg(""));
	}
	{ // Pythia8::BlackSubCollisionModel file:Pythia8/HISubCollisionModel.h line:484
		pybind11::class_<Pythia8::BlackSubCollisionModel, std::shared_ptr<Pythia8::BlackSubCollisionModel>, PyCallBack_Pythia8_BlackSubCollisionModel, Pythia8::SubCollisionModel> cl(M("Pythia8"), "BlackSubCollisionModel", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init( [](){ return new Pythia8::BlackSubCollisionModel(); }, [](){ return new PyCallBack_Pythia8_BlackSubCollisionModel(); } ) );
		cl.def("minParm", (class std::vector<double, class std::allocator<double> > (Pythia8::BlackSubCollisionModel::*)() const) &Pythia8::BlackSubCollisionModel::minParm, "C++: Pythia8::BlackSubCollisionModel::minParm() const --> class std::vector<double, class std::allocator<double> >");
		cl.def("defParm", (class std::vector<double, class std::allocator<double> > (Pythia8::BlackSubCollisionModel::*)() const) &Pythia8::BlackSubCollisionModel::defParm, "C++: Pythia8::BlackSubCollisionModel::defParm() const --> class std::vector<double, class std::allocator<double> >");
		cl.def("maxParm", (class std::vector<double, class std::allocator<double> > (Pythia8::BlackSubCollisionModel::*)() const) &Pythia8::BlackSubCollisionModel::maxParm, "C++: Pythia8::BlackSubCollisionModel::maxParm() const --> class std::vector<double, class std::allocator<double> >");
		cl.def("getSig", (struct Pythia8::SubCollisionModel::SigEst (Pythia8::BlackSubCollisionModel::*)() const) &Pythia8::BlackSubCollisionModel::getSig, "C++: Pythia8::BlackSubCollisionModel::getSig() const --> struct Pythia8::SubCollisionModel::SigEst");
		cl.def("getCollisions", (class Pythia8::SubCollisionSet (Pythia8::BlackSubCollisionModel::*)(class Pythia8::Nucleus &, class Pythia8::Nucleus &)) &Pythia8::BlackSubCollisionModel::getCollisions, "C++: Pythia8::BlackSubCollisionModel::getCollisions(class Pythia8::Nucleus &, class Pythia8::Nucleus &) --> class Pythia8::SubCollisionSet", pybind11::arg("proj"), pybind11::arg("targ"));
		cl.def("assign", (class Pythia8::BlackSubCollisionModel & (Pythia8::BlackSubCollisionModel::*)(const class Pythia8::BlackSubCollisionModel &)) &Pythia8::BlackSubCollisionModel::operator=, "C++: Pythia8::BlackSubCollisionModel::operator=(const class Pythia8::BlackSubCollisionModel &) --> class Pythia8::BlackSubCollisionModel &", pybind11::return_value_policy::reference, pybind11::arg(""));
	}
	{ // Pythia8::NaiveSubCollisionModel file:Pythia8/HISubCollisionModel.h line:520
		pybind11::class_<Pythia8::NaiveSubCollisionModel, std::shared_ptr<Pythia8::NaiveSubCollisionModel>, PyCallBack_Pythia8_NaiveSubCollisionModel, Pythia8::SubCollisionModel> cl(M("Pythia8"), "NaiveSubCollisionModel", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init( [](){ return new Pythia8::NaiveSubCollisionModel(); }, [](){ return new PyCallBack_Pythia8_NaiveSubCollisionModel(); } ) );
		cl.def("minParm", (class std::vector<double, class std::allocator<double> > (Pythia8::NaiveSubCollisionModel::*)() const) &Pythia8::NaiveSubCollisionModel::minParm, "C++: Pythia8::NaiveSubCollisionModel::minParm() const --> class std::vector<double, class std::allocator<double> >");
		cl.def("defParm", (class std::vector<double, class std::allocator<double> > (Pythia8::NaiveSubCollisionModel::*)() const) &Pythia8::NaiveSubCollisionModel::defParm, "C++: Pythia8::NaiveSubCollisionModel::defParm() const --> class std::vector<double, class std::allocator<double> >");
		cl.def("maxParm", (class std::vector<double, class std::allocator<double> > (Pythia8::NaiveSubCollisionModel::*)() const) &Pythia8::NaiveSubCollisionModel::maxParm, "C++: Pythia8::NaiveSubCollisionModel::maxParm() const --> class std::vector<double, class std::allocator<double> >");
		cl.def("getSig", (struct Pythia8::SubCollisionModel::SigEst (Pythia8::NaiveSubCollisionModel::*)() const) &Pythia8::NaiveSubCollisionModel::getSig, "C++: Pythia8::NaiveSubCollisionModel::getSig() const --> struct Pythia8::SubCollisionModel::SigEst");
		cl.def("getCollisions", (class Pythia8::SubCollisionSet (Pythia8::NaiveSubCollisionModel::*)(class Pythia8::Nucleus &, class Pythia8::Nucleus &)) &Pythia8::NaiveSubCollisionModel::getCollisions, "C++: Pythia8::NaiveSubCollisionModel::getCollisions(class Pythia8::Nucleus &, class Pythia8::Nucleus &) --> class Pythia8::SubCollisionSet", pybind11::arg("proj"), pybind11::arg("targ"));
		cl.def("assign", (class Pythia8::NaiveSubCollisionModel & (Pythia8::NaiveSubCollisionModel::*)(const class Pythia8::NaiveSubCollisionModel &)) &Pythia8::NaiveSubCollisionModel::operator=, "C++: Pythia8::NaiveSubCollisionModel::operator=(const class Pythia8::NaiveSubCollisionModel &) --> class Pythia8::NaiveSubCollisionModel &", pybind11::return_value_policy::reference, pybind11::arg(""));
	}
}
