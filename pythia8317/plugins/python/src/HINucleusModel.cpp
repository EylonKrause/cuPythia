#include <Pythia8/Basics.h>
#include <Pythia8/Event.h>
#include <Pythia8/HIBasics.h>
#include <Pythia8/HINucleusModel.h>
#include <Pythia8/HISubCollisionModel.h>
#include <Pythia8/Info.h>
#include <Pythia8/ParticleData.h>
#include <Pythia8/ResonanceWidths.h>
#include <functional>
#include <iterator>
#include <map>
#include <memory>
#include <set>
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

// Pythia8::ExternalNucleusModel file:Pythia8/HINucleusModel.h line:289
struct PyCallBack_Pythia8_ExternalNucleusModel : public Pythia8::ExternalNucleusModel {
	using Pythia8::ExternalNucleusModel::ExternalNucleusModel;

	bool init() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ExternalNucleusModel *>(this), "init");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return ExternalNucleusModel::init();
	}
	using _binder_ret_0 = class std::vector<class Pythia8::Nucleon, class std::allocator<class Pythia8::Nucleon> >;
	_binder_ret_0 generate() const override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ExternalNucleusModel *>(this), "generate");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<_binder_ret_0>::value) {
				static pybind11::detail::override_caster_t<_binder_ret_0> caster;
				return pybind11::detail::cast_ref<_binder_ret_0>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<_binder_ret_0>(std::move(o));
		}
		return ExternalNucleusModel::generate();
	}
	bool initGeometry() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ExternalNucleusModel *>(this), "initGeometry");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return NucleusModel::initGeometry();
	}
	void setPN(const class Pythia8::Vec4 & a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ExternalNucleusModel *>(this), "setPN");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return NucleusModel::setPN(a0);
	}
	void setMN(double a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ExternalNucleusModel *>(this), "setMN");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return NucleusModel::setMN(a0);
	}
	class Pythia8::Particle produceIon() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ExternalNucleusModel *>(this), "produceIon");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<class Pythia8::Particle>::value) {
				static pybind11::detail::override_caster_t<class Pythia8::Particle> caster;
				return pybind11::detail::cast_ref<class Pythia8::Particle>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<class Pythia8::Particle>(std::move(o));
		}
		return NucleusModel::produceIon();
	}
};

// Pythia8::HardCoreModel file:Pythia8/HINucleusModel.h line:325
struct PyCallBack_Pythia8_HardCoreModel : public Pythia8::HardCoreModel {
	using Pythia8::HardCoreModel::HardCoreModel;

	bool init() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HardCoreModel *>(this), "init");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return NucleusModel::init();
	}
	bool initGeometry() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HardCoreModel *>(this), "initGeometry");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return NucleusModel::initGeometry();
	}
	void setPN(const class Pythia8::Vec4 & a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HardCoreModel *>(this), "setPN");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return NucleusModel::setPN(a0);
	}
	void setMN(double a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HardCoreModel *>(this), "setMN");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return NucleusModel::setMN(a0);
	}
	class Pythia8::Particle produceIon() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HardCoreModel *>(this), "produceIon");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<class Pythia8::Particle>::value) {
				static pybind11::detail::override_caster_t<class Pythia8::Particle> caster;
				return pybind11::detail::cast_ref<class Pythia8::Particle>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<class Pythia8::Particle>(std::move(o));
		}
		return NucleusModel::produceIon();
	}
	using _binder_ret_0 = class std::vector<class Pythia8::Nucleon, class std::allocator<class Pythia8::Nucleon> >;
	_binder_ret_0 generate() const override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HardCoreModel *>(this), "generate");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<_binder_ret_0>::value) {
				static pybind11::detail::override_caster_t<_binder_ret_0> caster;
				return pybind11::detail::cast_ref<_binder_ret_0>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<_binder_ret_0>(std::move(o));
		}
		pybind11::pybind11_fail("Tried to call pure virtual function \"NucleusModel::generate\"");
	}
};

// Pythia8::WoodsSaxonModel file:Pythia8/HINucleusModel.h line:362
struct PyCallBack_Pythia8_WoodsSaxonModel : public Pythia8::WoodsSaxonModel {
	using Pythia8::WoodsSaxonModel::WoodsSaxonModel;

	bool init() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::WoodsSaxonModel *>(this), "init");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return WoodsSaxonModel::init();
	}
	bool initGeometry() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::WoodsSaxonModel *>(this), "initGeometry");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return WoodsSaxonModel::initGeometry();
	}
	using _binder_ret_0 = class std::vector<class Pythia8::Nucleon, class std::allocator<class Pythia8::Nucleon> >;
	_binder_ret_0 generate() const override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::WoodsSaxonModel *>(this), "generate");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<_binder_ret_0>::value) {
				static pybind11::detail::override_caster_t<_binder_ret_0> caster;
				return pybind11::detail::cast_ref<_binder_ret_0>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<_binder_ret_0>(std::move(o));
		}
		return WoodsSaxonModel::generate();
	}
	void setPN(const class Pythia8::Vec4 & a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::WoodsSaxonModel *>(this), "setPN");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return NucleusModel::setPN(a0);
	}
	void setMN(double a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::WoodsSaxonModel *>(this), "setMN");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return NucleusModel::setMN(a0);
	}
	class Pythia8::Particle produceIon() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::WoodsSaxonModel *>(this), "produceIon");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<class Pythia8::Particle>::value) {
				static pybind11::detail::override_caster_t<class Pythia8::Particle> caster;
				return pybind11::detail::cast_ref<class Pythia8::Particle>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<class Pythia8::Particle>(std::move(o));
		}
		return NucleusModel::produceIon();
	}
};

// Pythia8::GLISSANDOModel file:Pythia8/HINucleusModel.h line:417
struct PyCallBack_Pythia8_GLISSANDOModel : public Pythia8::GLISSANDOModel {
	using Pythia8::GLISSANDOModel::GLISSANDOModel;

	bool init() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::GLISSANDOModel *>(this), "init");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return GLISSANDOModel::init();
	}
	bool initGeometry() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::GLISSANDOModel *>(this), "initGeometry");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return GLISSANDOModel::initGeometry();
	}
	using _binder_ret_0 = class std::vector<class Pythia8::Nucleon, class std::allocator<class Pythia8::Nucleon> >;
	_binder_ret_0 generate() const override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::GLISSANDOModel *>(this), "generate");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<_binder_ret_0>::value) {
				static pybind11::detail::override_caster_t<_binder_ret_0> caster;
				return pybind11::detail::cast_ref<_binder_ret_0>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<_binder_ret_0>(std::move(o));
		}
		return WoodsSaxonModel::generate();
	}
	void setPN(const class Pythia8::Vec4 & a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::GLISSANDOModel *>(this), "setPN");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return NucleusModel::setPN(a0);
	}
	void setMN(double a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::GLISSANDOModel *>(this), "setMN");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return NucleusModel::setMN(a0);
	}
	class Pythia8::Particle produceIon() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::GLISSANDOModel *>(this), "produceIon");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<class Pythia8::Particle>::value) {
				static pybind11::detail::override_caster_t<class Pythia8::Particle> caster;
				return pybind11::detail::cast_ref<class Pythia8::Particle>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<class Pythia8::Particle>(std::move(o));
		}
		return NucleusModel::produceIon();
	}
};

// Pythia8::HOShellModel file:Pythia8/HINucleusModel.h line:437
struct PyCallBack_Pythia8_HOShellModel : public Pythia8::HOShellModel {
	using Pythia8::HOShellModel::HOShellModel;

	bool init() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HOShellModel *>(this), "init");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return HOShellModel::init();
	}
	using _binder_ret_0 = class std::vector<class Pythia8::Nucleon, class std::allocator<class Pythia8::Nucleon> >;
	_binder_ret_0 generate() const override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HOShellModel *>(this), "generate");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<_binder_ret_0>::value) {
				static pybind11::detail::override_caster_t<_binder_ret_0> caster;
				return pybind11::detail::cast_ref<_binder_ret_0>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<_binder_ret_0>(std::move(o));
		}
		return HOShellModel::generate();
	}
	class Pythia8::Vec4 generateNucleon() const override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HOShellModel *>(this), "generateNucleon");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<class Pythia8::Vec4>::value) {
				static pybind11::detail::override_caster_t<class Pythia8::Vec4> caster;
				return pybind11::detail::cast_ref<class Pythia8::Vec4>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<class Pythia8::Vec4>(std::move(o));
		}
		return HOShellModel::generateNucleon();
	}
	bool initGeometry() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HOShellModel *>(this), "initGeometry");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return NucleusModel::initGeometry();
	}
	void setPN(const class Pythia8::Vec4 & a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HOShellModel *>(this), "setPN");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return NucleusModel::setPN(a0);
	}
	void setMN(double a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HOShellModel *>(this), "setMN");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return NucleusModel::setMN(a0);
	}
	class Pythia8::Particle produceIon() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HOShellModel *>(this), "produceIon");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<class Pythia8::Particle>::value) {
				static pybind11::detail::override_caster_t<class Pythia8::Particle> caster;
				return pybind11::detail::cast_ref<class Pythia8::Particle>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<class Pythia8::Particle>(std::move(o));
		}
		return NucleusModel::produceIon();
	}
};

// Pythia8::HulthenModel file:Pythia8/HINucleusModel.h line:484
struct PyCallBack_Pythia8_HulthenModel : public Pythia8::HulthenModel {
	using Pythia8::HulthenModel::HulthenModel;

	bool init() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HulthenModel *>(this), "init");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return HulthenModel::init();
	}
	using _binder_ret_0 = class std::vector<class Pythia8::Nucleon, class std::allocator<class Pythia8::Nucleon> >;
	_binder_ret_0 generate() const override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HulthenModel *>(this), "generate");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<_binder_ret_0>::value) {
				static pybind11::detail::override_caster_t<_binder_ret_0> caster;
				return pybind11::detail::cast_ref<_binder_ret_0>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<_binder_ret_0>(std::move(o));
		}
		return HulthenModel::generate();
	}
	bool initGeometry() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HulthenModel *>(this), "initGeometry");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return NucleusModel::initGeometry();
	}
	void setPN(const class Pythia8::Vec4 & a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HulthenModel *>(this), "setPN");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return NucleusModel::setPN(a0);
	}
	void setMN(double a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HulthenModel *>(this), "setMN");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return NucleusModel::setMN(a0);
	}
	class Pythia8::Particle produceIon() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::HulthenModel *>(this), "produceIon");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<class Pythia8::Particle>::value) {
				static pybind11::detail::override_caster_t<class Pythia8::Particle> caster;
				return pybind11::detail::cast_ref<class Pythia8::Particle>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<class Pythia8::Particle>(std::move(o));
		}
		return NucleusModel::produceIon();
	}
};

// Pythia8::GaussianModel file:Pythia8/HINucleusModel.h line:518
struct PyCallBack_Pythia8_GaussianModel : public Pythia8::GaussianModel {
	using Pythia8::GaussianModel::GaussianModel;

	bool init() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::GaussianModel *>(this), "init");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return GaussianModel::init();
	}
	using _binder_ret_0 = class std::vector<class Pythia8::Nucleon, class std::allocator<class Pythia8::Nucleon> >;
	_binder_ret_0 generate() const override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::GaussianModel *>(this), "generate");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<_binder_ret_0>::value) {
				static pybind11::detail::override_caster_t<_binder_ret_0> caster;
				return pybind11::detail::cast_ref<_binder_ret_0>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<_binder_ret_0>(std::move(o));
		}
		return GaussianModel::generate();
	}
	class Pythia8::Vec4 generateNucleon() const override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::GaussianModel *>(this), "generateNucleon");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<class Pythia8::Vec4>::value) {
				static pybind11::detail::override_caster_t<class Pythia8::Vec4> caster;
				return pybind11::detail::cast_ref<class Pythia8::Vec4>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<class Pythia8::Vec4>(std::move(o));
		}
		return GaussianModel::generateNucleon();
	}
	bool initGeometry() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::GaussianModel *>(this), "initGeometry");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return NucleusModel::initGeometry();
	}
	void setPN(const class Pythia8::Vec4 & a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::GaussianModel *>(this), "setPN");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return NucleusModel::setPN(a0);
	}
	void setMN(double a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::GaussianModel *>(this), "setMN");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return NucleusModel::setMN(a0);
	}
	class Pythia8::Particle produceIon() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::GaussianModel *>(this), "produceIon");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<class Pythia8::Particle>::value) {
				static pybind11::detail::override_caster_t<class Pythia8::Particle> caster;
				return pybind11::detail::cast_ref<class Pythia8::Particle>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<class Pythia8::Particle>(std::move(o));
		}
		return NucleusModel::produceIon();
	}
};

// Pythia8::ClusterModel file:Pythia8/HINucleusModel.h line:549
struct PyCallBack_Pythia8_ClusterModel : public Pythia8::ClusterModel {
	using Pythia8::ClusterModel::ClusterModel;

	bool init() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ClusterModel *>(this), "init");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return ClusterModel::init();
	}
	using _binder_ret_0 = class std::vector<class Pythia8::Nucleon, class std::allocator<class Pythia8::Nucleon> >;
	_binder_ret_0 generate() const override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ClusterModel *>(this), "generate");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<_binder_ret_0>::value) {
				static pybind11::detail::override_caster_t<_binder_ret_0> caster;
				return pybind11::detail::cast_ref<_binder_ret_0>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<_binder_ret_0>(std::move(o));
		}
		return ClusterModel::generate();
	}
	bool initGeometry() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ClusterModel *>(this), "initGeometry");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return NucleusModel::initGeometry();
	}
	void setPN(const class Pythia8::Vec4 & a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ClusterModel *>(this), "setPN");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return NucleusModel::setPN(a0);
	}
	void setMN(double a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ClusterModel *>(this), "setMN");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return NucleusModel::setMN(a0);
	}
	class Pythia8::Particle produceIon() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ClusterModel *>(this), "produceIon");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<class Pythia8::Particle>::value) {
				static pybind11::detail::override_caster_t<class Pythia8::Particle> caster;
				return pybind11::detail::cast_ref<class Pythia8::Particle>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<class Pythia8::Particle>(std::move(o));
		}
		return NucleusModel::produceIon();
	}
};

void bind_Pythia8_HINucleusModel(std::function< pybind11::module &(std::string const &namespace_) > &M)
{
	{ // Pythia8::ExternalNucleusModel file:Pythia8/HINucleusModel.h line:289
		pybind11::class_<Pythia8::ExternalNucleusModel, std::shared_ptr<Pythia8::ExternalNucleusModel>, PyCallBack_Pythia8_ExternalNucleusModel, Pythia8::NucleusModel> cl(M("Pythia8"), "ExternalNucleusModel", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init( [](){ return new Pythia8::ExternalNucleusModel(); }, [](){ return new PyCallBack_Pythia8_ExternalNucleusModel(); } ) );
		cl.def( pybind11::init( [](PyCallBack_Pythia8_ExternalNucleusModel const &o){ return new PyCallBack_Pythia8_ExternalNucleusModel(o); } ) );
		cl.def( pybind11::init( [](Pythia8::ExternalNucleusModel const &o){ return new Pythia8::ExternalNucleusModel(o); } ) );
		cl.def("init", (bool (Pythia8::ExternalNucleusModel::*)()) &Pythia8::ExternalNucleusModel::init, "C++: Pythia8::ExternalNucleusModel::init() --> bool");
		cl.def("generate", (class std::vector<class Pythia8::Nucleon, class std::allocator<class Pythia8::Nucleon> > (Pythia8::ExternalNucleusModel::*)() const) &Pythia8::ExternalNucleusModel::generate, "C++: Pythia8::ExternalNucleusModel::generate() const --> class std::vector<class Pythia8::Nucleon, class std::allocator<class Pythia8::Nucleon> >");
		cl.def("assign", (class Pythia8::ExternalNucleusModel & (Pythia8::ExternalNucleusModel::*)(const class Pythia8::ExternalNucleusModel &)) &Pythia8::ExternalNucleusModel::operator=, "C++: Pythia8::ExternalNucleusModel::operator=(const class Pythia8::ExternalNucleusModel &) --> class Pythia8::ExternalNucleusModel &", pybind11::return_value_policy::reference, pybind11::arg(""));
	}
	{ // Pythia8::HardCoreModel file:Pythia8/HINucleusModel.h line:325
		pybind11::class_<Pythia8::HardCoreModel, std::shared_ptr<Pythia8::HardCoreModel>, PyCallBack_Pythia8_HardCoreModel, Pythia8::NucleusModel> cl(M("Pythia8"), "HardCoreModel", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init( [](){ return new PyCallBack_Pythia8_HardCoreModel(); } ) );
		cl.def(pybind11::init<PyCallBack_Pythia8_HardCoreModel const &>());
		cl.def_readwrite("useHardCore", &Pythia8::HardCoreModel::useHardCore);
		cl.def_readwrite("gaussHardCore", &Pythia8::HardCoreModel::gaussHardCore);
		cl.def_readwrite("hardCoreRadius", &Pythia8::HardCoreModel::hardCoreRadius);
		cl.def("initHardCore", (void (Pythia8::HardCoreModel::*)()) &Pythia8::HardCoreModel::initHardCore, "C++: Pythia8::HardCoreModel::initHardCore() --> void");
		cl.def("rSample", (double (Pythia8::HardCoreModel::*)() const) &Pythia8::HardCoreModel::rSample, "C++: Pythia8::HardCoreModel::rSample() const --> double");
		cl.def("assign", (class Pythia8::HardCoreModel & (Pythia8::HardCoreModel::*)(const class Pythia8::HardCoreModel &)) &Pythia8::HardCoreModel::operator=, "C++: Pythia8::HardCoreModel::operator=(const class Pythia8::HardCoreModel &) --> class Pythia8::HardCoreModel &", pybind11::return_value_policy::reference, pybind11::arg(""));
	}
	{ // Pythia8::WoodsSaxonModel file:Pythia8/HINucleusModel.h line:362
		pybind11::class_<Pythia8::WoodsSaxonModel, std::shared_ptr<Pythia8::WoodsSaxonModel>, PyCallBack_Pythia8_WoodsSaxonModel, Pythia8::HardCoreModel> cl(M("Pythia8"), "WoodsSaxonModel", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init( [](){ return new Pythia8::WoodsSaxonModel(); }, [](){ return new PyCallBack_Pythia8_WoodsSaxonModel(); } ) );
		cl.def( pybind11::init( [](PyCallBack_Pythia8_WoodsSaxonModel const &o){ return new PyCallBack_Pythia8_WoodsSaxonModel(o); } ) );
		cl.def( pybind11::init( [](Pythia8::WoodsSaxonModel const &o){ return new Pythia8::WoodsSaxonModel(o); } ) );
		cl.def_readwrite("aSave", &Pythia8::WoodsSaxonModel::aSave);
		cl.def("init", (bool (Pythia8::WoodsSaxonModel::*)()) &Pythia8::WoodsSaxonModel::init, "C++: Pythia8::WoodsSaxonModel::init() --> bool");
		cl.def("initGeometry", (bool (Pythia8::WoodsSaxonModel::*)()) &Pythia8::WoodsSaxonModel::initGeometry, "C++: Pythia8::WoodsSaxonModel::initGeometry() --> bool");
		cl.def("generate", (class std::vector<class Pythia8::Nucleon, class std::allocator<class Pythia8::Nucleon> > (Pythia8::WoodsSaxonModel::*)() const) &Pythia8::WoodsSaxonModel::generate, "C++: Pythia8::WoodsSaxonModel::generate() const --> class std::vector<class Pythia8::Nucleon, class std::allocator<class Pythia8::Nucleon> >");
		cl.def("a", (double (Pythia8::WoodsSaxonModel::*)() const) &Pythia8::WoodsSaxonModel::a, "C++: Pythia8::WoodsSaxonModel::a() const --> double");
		cl.def("generateNucleon", (class Pythia8::Vec4 (Pythia8::WoodsSaxonModel::*)() const) &Pythia8::WoodsSaxonModel::generateNucleon, "C++: Pythia8::WoodsSaxonModel::generateNucleon() const --> class Pythia8::Vec4");
		cl.def("overestimates", (void (Pythia8::WoodsSaxonModel::*)()) &Pythia8::WoodsSaxonModel::overestimates, "C++: Pythia8::WoodsSaxonModel::overestimates() --> void");
		cl.def("assign", (class Pythia8::WoodsSaxonModel & (Pythia8::WoodsSaxonModel::*)(const class Pythia8::WoodsSaxonModel &)) &Pythia8::WoodsSaxonModel::operator=, "C++: Pythia8::WoodsSaxonModel::operator=(const class Pythia8::WoodsSaxonModel &) --> class Pythia8::WoodsSaxonModel &", pybind11::return_value_policy::reference, pybind11::arg(""));
	}
	{ // Pythia8::GLISSANDOModel file:Pythia8/HINucleusModel.h line:417
		pybind11::class_<Pythia8::GLISSANDOModel, std::shared_ptr<Pythia8::GLISSANDOModel>, PyCallBack_Pythia8_GLISSANDOModel, Pythia8::WoodsSaxonModel> cl(M("Pythia8"), "GLISSANDOModel", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init( [](){ return new Pythia8::GLISSANDOModel(); }, [](){ return new PyCallBack_Pythia8_GLISSANDOModel(); } ) );
		cl.def( pybind11::init( [](PyCallBack_Pythia8_GLISSANDOModel const &o){ return new PyCallBack_Pythia8_GLISSANDOModel(o); } ) );
		cl.def( pybind11::init( [](Pythia8::GLISSANDOModel const &o){ return new Pythia8::GLISSANDOModel(o); } ) );
		cl.def("init", (bool (Pythia8::GLISSANDOModel::*)()) &Pythia8::GLISSANDOModel::init, "C++: Pythia8::GLISSANDOModel::init() --> bool");
		cl.def("initGeometry", (bool (Pythia8::GLISSANDOModel::*)()) &Pythia8::GLISSANDOModel::initGeometry, "C++: Pythia8::GLISSANDOModel::initGeometry() --> bool");
		cl.def("assign", (class Pythia8::GLISSANDOModel & (Pythia8::GLISSANDOModel::*)(const class Pythia8::GLISSANDOModel &)) &Pythia8::GLISSANDOModel::operator=, "C++: Pythia8::GLISSANDOModel::operator=(const class Pythia8::GLISSANDOModel &) --> class Pythia8::GLISSANDOModel &", pybind11::return_value_policy::reference, pybind11::arg(""));
	}
	{ // Pythia8::HOShellModel file:Pythia8/HINucleusModel.h line:437
		pybind11::class_<Pythia8::HOShellModel, std::shared_ptr<Pythia8::HOShellModel>, PyCallBack_Pythia8_HOShellModel, Pythia8::HardCoreModel> cl(M("Pythia8"), "HOShellModel", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init( [](){ return new Pythia8::HOShellModel(); }, [](){ return new PyCallBack_Pythia8_HOShellModel(); } ) );
		cl.def( pybind11::init( [](PyCallBack_Pythia8_HOShellModel const &o){ return new PyCallBack_Pythia8_HOShellModel(o); } ) );
		cl.def( pybind11::init( [](Pythia8::HOShellModel const &o){ return new Pythia8::HOShellModel(o); } ) );
		cl.def_readwrite("nucleusChR", &Pythia8::HOShellModel::nucleusChR);
		cl.def_readwrite("protonChR", &Pythia8::HOShellModel::protonChR);
		cl.def_readwrite("C2", &Pythia8::HOShellModel::C2);
		cl.def_readwrite("rhoMax", &Pythia8::HOShellModel::rhoMax);
		cl.def("init", (bool (Pythia8::HOShellModel::*)()) &Pythia8::HOShellModel::init, "C++: Pythia8::HOShellModel::init() --> bool");
		cl.def("generate", (class std::vector<class Pythia8::Nucleon, class std::allocator<class Pythia8::Nucleon> > (Pythia8::HOShellModel::*)() const) &Pythia8::HOShellModel::generate, "C++: Pythia8::HOShellModel::generate() const --> class std::vector<class Pythia8::Nucleon, class std::allocator<class Pythia8::Nucleon> >");
		cl.def("generateNucleon", (class Pythia8::Vec4 (Pythia8::HOShellModel::*)() const) &Pythia8::HOShellModel::generateNucleon, "C++: Pythia8::HOShellModel::generateNucleon() const --> class Pythia8::Vec4");
		cl.def("rho", (double (Pythia8::HOShellModel::*)(double) const) &Pythia8::HOShellModel::rho, "C++: Pythia8::HOShellModel::rho(double) const --> double", pybind11::arg("r"));
		cl.def("assign", (class Pythia8::HOShellModel & (Pythia8::HOShellModel::*)(const class Pythia8::HOShellModel &)) &Pythia8::HOShellModel::operator=, "C++: Pythia8::HOShellModel::operator=(const class Pythia8::HOShellModel &) --> class Pythia8::HOShellModel &", pybind11::return_value_policy::reference, pybind11::arg(""));
	}
	{ // Pythia8::HulthenModel file:Pythia8/HINucleusModel.h line:484
		pybind11::class_<Pythia8::HulthenModel, std::shared_ptr<Pythia8::HulthenModel>, PyCallBack_Pythia8_HulthenModel, Pythia8::NucleusModel> cl(M("Pythia8"), "HulthenModel", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init( [](){ return new Pythia8::HulthenModel(); }, [](){ return new PyCallBack_Pythia8_HulthenModel(); } ) );
		cl.def( pybind11::init( [](PyCallBack_Pythia8_HulthenModel const &o){ return new PyCallBack_Pythia8_HulthenModel(o); } ) );
		cl.def( pybind11::init( [](Pythia8::HulthenModel const &o){ return new Pythia8::HulthenModel(o); } ) );
		cl.def_readwrite("hA", &Pythia8::HulthenModel::hA);
		cl.def_readwrite("hB", &Pythia8::HulthenModel::hB);
		cl.def("init", (bool (Pythia8::HulthenModel::*)()) &Pythia8::HulthenModel::init, "C++: Pythia8::HulthenModel::init() --> bool");
		cl.def("generate", (class std::vector<class Pythia8::Nucleon, class std::allocator<class Pythia8::Nucleon> > (Pythia8::HulthenModel::*)() const) &Pythia8::HulthenModel::generate, "C++: Pythia8::HulthenModel::generate() const --> class std::vector<class Pythia8::Nucleon, class std::allocator<class Pythia8::Nucleon> >");
		cl.def("rho", (double (Pythia8::HulthenModel::*)(double) const) &Pythia8::HulthenModel::rho, "C++: Pythia8::HulthenModel::rho(double) const --> double", pybind11::arg("r"));
		cl.def("assign", (class Pythia8::HulthenModel & (Pythia8::HulthenModel::*)(const class Pythia8::HulthenModel &)) &Pythia8::HulthenModel::operator=, "C++: Pythia8::HulthenModel::operator=(const class Pythia8::HulthenModel &) --> class Pythia8::HulthenModel &", pybind11::return_value_policy::reference, pybind11::arg(""));
	}
	{ // Pythia8::GaussianModel file:Pythia8/HINucleusModel.h line:518
		pybind11::class_<Pythia8::GaussianModel, std::shared_ptr<Pythia8::GaussianModel>, PyCallBack_Pythia8_GaussianModel, Pythia8::HardCoreModel> cl(M("Pythia8"), "GaussianModel", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init( [](){ return new Pythia8::GaussianModel(); }, [](){ return new PyCallBack_Pythia8_GaussianModel(); } ) );
		cl.def( pybind11::init( [](PyCallBack_Pythia8_GaussianModel const &o){ return new PyCallBack_Pythia8_GaussianModel(o); } ) );
		cl.def( pybind11::init( [](Pythia8::GaussianModel const &o){ return new Pythia8::GaussianModel(o); } ) );
		cl.def_readwrite("nucleusChR", &Pythia8::GaussianModel::nucleusChR);
		cl.def("init", (bool (Pythia8::GaussianModel::*)()) &Pythia8::GaussianModel::init, "C++: Pythia8::GaussianModel::init() --> bool");
		cl.def("generate", (class std::vector<class Pythia8::Nucleon, class std::allocator<class Pythia8::Nucleon> > (Pythia8::GaussianModel::*)() const) &Pythia8::GaussianModel::generate, "C++: Pythia8::GaussianModel::generate() const --> class std::vector<class Pythia8::Nucleon, class std::allocator<class Pythia8::Nucleon> >");
		cl.def("generateNucleon", (class Pythia8::Vec4 (Pythia8::GaussianModel::*)() const) &Pythia8::GaussianModel::generateNucleon, "C++: Pythia8::GaussianModel::generateNucleon() const --> class Pythia8::Vec4");
		cl.def("assign", (class Pythia8::GaussianModel & (Pythia8::GaussianModel::*)(const class Pythia8::GaussianModel &)) &Pythia8::GaussianModel::operator=, "C++: Pythia8::GaussianModel::operator=(const class Pythia8::GaussianModel &) --> class Pythia8::GaussianModel &", pybind11::return_value_policy::reference, pybind11::arg(""));
	}
	{ // Pythia8::ClusterModel file:Pythia8/HINucleusModel.h line:549
		pybind11::class_<Pythia8::ClusterModel, std::shared_ptr<Pythia8::ClusterModel>, PyCallBack_Pythia8_ClusterModel, Pythia8::HardCoreModel> cl(M("Pythia8"), "ClusterModel", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init( [](){ return new Pythia8::ClusterModel(); }, [](){ return new PyCallBack_Pythia8_ClusterModel(); } ) );
		cl.def("init", (bool (Pythia8::ClusterModel::*)()) &Pythia8::ClusterModel::init, "C++: Pythia8::ClusterModel::init() --> bool");
		cl.def("generate", (class std::vector<class Pythia8::Nucleon, class std::allocator<class Pythia8::Nucleon> > (Pythia8::ClusterModel::*)() const) &Pythia8::ClusterModel::generate, "C++: Pythia8::ClusterModel::generate() const --> class std::vector<class Pythia8::Nucleon, class std::allocator<class Pythia8::Nucleon> >");
	}
	{ // Pythia8::SubCollision file:Pythia8/HISubCollisionModel.h line:30
		pybind11::class_<Pythia8::SubCollision, std::shared_ptr<Pythia8::SubCollision>> cl(M("Pythia8"), "SubCollision", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init<class Pythia8::Nucleon &, class Pythia8::Nucleon &, double, double, double, enum Pythia8::SubCollision::CollisionType>(), pybind11::arg("projIn"), pybind11::arg("targIn"), pybind11::arg("bIn"), pybind11::arg("bpIn"), pybind11::arg("olappIn"), pybind11::arg("typeIn") );

		cl.def( pybind11::init( [](){ return new Pythia8::SubCollision(); } ) );
		cl.def( pybind11::init( [](Pythia8::SubCollision const &o){ return new Pythia8::SubCollision(o); } ) );

		pybind11::enum_<Pythia8::SubCollision::CollisionType>(cl, "CollisionType", pybind11::arithmetic(), "")
			.value("NONE", Pythia8::SubCollision::CollisionType::NONE)
			.value("ELASTIC", Pythia8::SubCollision::CollisionType::ELASTIC)
			.value("SDEP", Pythia8::SubCollision::CollisionType::SDEP)
			.value("SDET", Pythia8::SubCollision::CollisionType::SDET)
			.value("DDE", Pythia8::SubCollision::CollisionType::DDE)
			.value("CDE", Pythia8::SubCollision::CollisionType::CDE)
			.value("ABS", Pythia8::SubCollision::CollisionType::ABS)
			.value("LEXC", Pythia8::SubCollision::CollisionType::LEXC)
			.value("LANN", Pythia8::SubCollision::CollisionType::LANN)
			.value("LRES", Pythia8::SubCollision::CollisionType::LRES)
			.export_values();

		cl.def_readwrite("b", &Pythia8::SubCollision::b);
		cl.def_readwrite("bp", &Pythia8::SubCollision::bp);
		cl.def_readwrite("olapp", &Pythia8::SubCollision::olapp);
		cl.def_readwrite("type", &Pythia8::SubCollision::type);
		cl.def_readwrite("failed", &Pythia8::SubCollision::failed);
		cl.def("nucleons", (int (Pythia8::SubCollision::*)() const) &Pythia8::SubCollision::nucleons, "C++: Pythia8::SubCollision::nucleons() const --> int");
		cl.def("assign", (class Pythia8::SubCollision & (Pythia8::SubCollision::*)(const class Pythia8::SubCollision &)) &Pythia8::SubCollision::operator=, "C++: Pythia8::SubCollision::operator=(const class Pythia8::SubCollision &) --> class Pythia8::SubCollision &", pybind11::return_value_policy::reference, pybind11::arg(""));
	}
	{ // Pythia8::SubCollisionSet file:Pythia8/HISubCollisionModel.h line:96
		pybind11::class_<Pythia8::SubCollisionSet, std::shared_ptr<Pythia8::SubCollisionSet>> cl(M("Pythia8"), "SubCollisionSet", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init( [](){ return new Pythia8::SubCollisionSet(); } ) );
		cl.def( pybind11::init<class std::multiset<class Pythia8::SubCollision, struct std::less<class Pythia8::SubCollision>, class std::allocator<class Pythia8::SubCollision> >, double>(), pybind11::arg("subCollisionsIn"), pybind11::arg("TIn") );

		cl.def( pybind11::init<class std::multiset<class Pythia8::SubCollision, struct std::less<class Pythia8::SubCollision>, class std::allocator<class Pythia8::SubCollision> >, double, double, double, double>(), pybind11::arg("subCollisionsIn"), pybind11::arg("TIn"), pybind11::arg("T12In"), pybind11::arg("T21In"), pybind11::arg("T22In") );

		cl.def("empty", (bool (Pythia8::SubCollisionSet::*)() const) &Pythia8::SubCollisionSet::empty, "C++: Pythia8::SubCollisionSet::empty() const --> bool");
		cl.def("T", [](Pythia8::SubCollisionSet const &o) -> double { return o.T(); }, "");
		cl.def("T", (double (Pythia8::SubCollisionSet::*)(unsigned int) const) &Pythia8::SubCollisionSet::T, "C++: Pythia8::SubCollisionSet::T(unsigned int) const --> double", pybind11::arg("i"));
		cl.def("assign", (class Pythia8::SubCollisionSet & (Pythia8::SubCollisionSet::*)(const class Pythia8::SubCollisionSet &)) &Pythia8::SubCollisionSet::operator=, "C++: Pythia8::SubCollisionSet::operator=(const class Pythia8::SubCollisionSet &) --> class Pythia8::SubCollisionSet &", pybind11::return_value_policy::reference, pybind11::arg(""));
	}
}
