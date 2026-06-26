#include <Pythia8/Basics.h>
#include <Pythia8/BeamParticle.h>
#include <Pythia8/BeamSetup.h>
#include <Pythia8/BeamShape.h>
#include <Pythia8/Event.h>
#include <Pythia8/FragmentationFlavZpT.h>
#include <Pythia8/FragmentationModel.h>
#include <Pythia8/HIInfo.h>
#include <Pythia8/HadronWidths.h>
#include <Pythia8/HeavyIons.h>
#include <Pythia8/Info.h>
#include <Pythia8/LHEF3.h>
#include <Pythia8/LesHouches.h>
#include <Pythia8/Logger.h>
#include <Pythia8/Merging.h>
#include <Pythia8/MergingHooks.h>
#include <Pythia8/ParticleData.h>
#include <Pythia8/ParticleDecays.h>
#include <Pythia8/PartonDistributions.h>
#include <Pythia8/PartonSystems.h>
#include <Pythia8/PartonVertex.h>
#include <Pythia8/PhaseSpace.h>
#include <Pythia8/PhysicsBase.h>
#include <Pythia8/Pythia.h>
#include <Pythia8/ResonanceWidths.h>
#include <Pythia8/Settings.h>
#include <Pythia8/ShowerModel.h>
#include <Pythia8/SigmaLowEnergy.h>
#include <Pythia8/SigmaProcess.h>
#include <Pythia8/SigmaTotal.h>
#include <Pythia8/SimpleTimeShower.h>
#include <Pythia8/StandardModel.h>
#include <Pythia8/SusyCouplings.h>
#include <Pythia8/ThermalFragmentation.h>
#include <Pythia8/TimeShower.h>
#include <Pythia8/UserHooks.h>
#include <Pythia8/Weights.h>
#include <cwchar>
#include <functional>
#include <ios>
#include <istream>
#include <iterator>
#include <map>
#include <memory>
#include <ostream>
#include <sstream>
#include <sstream> // __str__
#include <streambuf>
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

// Pythia8::SimpleTimeShower file:Pythia8/SimpleTimeShower.h line:81
struct PyCallBack_Pythia8_SimpleTimeShower : public Pythia8::SimpleTimeShower {
	using Pythia8::SimpleTimeShower::SimpleTimeShower;

	void init(class Pythia8::BeamParticle * a0, class Pythia8::BeamParticle * a1) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::SimpleTimeShower *>(this), "init");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return SimpleTimeShower::init(a0, a1);
	}
	bool limitPTmax(class Pythia8::Event & a0, double a1, double a2) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::SimpleTimeShower *>(this), "limitPTmax");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return SimpleTimeShower::limitPTmax(a0, a1, a2);
	}
	int shower(int a0, int a1, class Pythia8::Event & a2, double a3, int a4) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::SimpleTimeShower *>(this), "shower");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4);
			if (pybind11::detail::cast_is_temporary_value_reference<int>::value) {
				static pybind11::detail::override_caster_t<int> caster;
				return pybind11::detail::cast_ref<int>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<int>(std::move(o));
		}
		return SimpleTimeShower::shower(a0, a1, a2, a3, a4);
	}
	int showerQED(int a0, int a1, class Pythia8::Event & a2, double a3) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::SimpleTimeShower *>(this), "showerQED");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3);
			if (pybind11::detail::cast_is_temporary_value_reference<int>::value) {
				static pybind11::detail::override_caster_t<int> caster;
				return pybind11::detail::cast_ref<int>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<int>(std::move(o));
		}
		return SimpleTimeShower::showerQED(a0, a1, a2, a3);
	}
	void prepareProcess(class Pythia8::Event & a0, class Pythia8::Event & a1, class std::vector<int, class std::allocator<int> > & a2) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::SimpleTimeShower *>(this), "prepareProcess");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return SimpleTimeShower::prepareProcess(a0, a1, a2);
	}
	void prepareGlobal(class Pythia8::Event & a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::SimpleTimeShower *>(this), "prepareGlobal");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return SimpleTimeShower::prepareGlobal(a0);
	}
	void prepare(int a0, class Pythia8::Event & a1, bool a2) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::SimpleTimeShower *>(this), "prepare");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return SimpleTimeShower::prepare(a0, a1, a2);
	}
	void rescatterUpdate(int a0, class Pythia8::Event & a1) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::SimpleTimeShower *>(this), "rescatterUpdate");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return SimpleTimeShower::rescatterUpdate(a0, a1);
	}
	void update(int a0, class Pythia8::Event & a1, bool a2) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::SimpleTimeShower *>(this), "update");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return SimpleTimeShower::update(a0, a1, a2);
	}
	double pTnext(class Pythia8::Event & a0, double a1, double a2, bool a3, bool a4) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::SimpleTimeShower *>(this), "pTnext");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4);
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return SimpleTimeShower::pTnext(a0, a1, a2, a3, a4);
	}
	double pTnextResDec() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::SimpleTimeShower *>(this), "pTnextResDec");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return SimpleTimeShower::pTnextResDec();
	}
	bool branch(class Pythia8::Event & a0, bool a1) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::SimpleTimeShower *>(this), "branch");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return SimpleTimeShower::branch(a0, a1);
	}
	bool resonanceShower(class Pythia8::Event & a0, class Pythia8::Event & a1, class std::vector<int, class std::allocator<int> > & a2, double a3) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::SimpleTimeShower *>(this), "resonanceShower");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return SimpleTimeShower::resonanceShower(a0, a1, a2, a3);
	}
	void list() const override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::SimpleTimeShower *>(this), "list");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return SimpleTimeShower::list();
	}
	bool initUncertainties() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::SimpleTimeShower *>(this), "initUncertainties");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return SimpleTimeShower::initUncertainties();
	}
	bool initEnhancements() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::SimpleTimeShower *>(this), "initEnhancements");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return SimpleTimeShower::initEnhancements();
	}
	bool getHasWeaklyRadiated() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::SimpleTimeShower *>(this), "getHasWeaklyRadiated");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return SimpleTimeShower::getHasWeaklyRadiated();
	}
	int system() const override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::SimpleTimeShower *>(this), "system");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<int>::value) {
				static pybind11::detail::override_caster_t<int> caster;
				return pybind11::detail::cast_ref<int>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<int>(std::move(o));
		}
		return SimpleTimeShower::system();
	}
	double enhancePTmax() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::SimpleTimeShower *>(this), "enhancePTmax");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return SimpleTimeShower::enhancePTmax();
	}
	double pTLastInShower() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::SimpleTimeShower *>(this), "pTLastInShower");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return SimpleTimeShower::pTLastInShower();
	}
	double noEmissionProbability(double a0, double a1, double a2, int a3, int a4, double a5, double a6) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::SimpleTimeShower *>(this), "noEmissionProbability");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4, a5, a6);
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return SimpleTimeShower::noEmissionProbability(a0, a1, a2, a3, a4, a5, a6);
	}
	int showerQEDafterRemnants(class Pythia8::Event & a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::SimpleTimeShower *>(this), "showerQEDafterRemnants");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<int>::value) {
				static pybind11::detail::override_caster_t<int> caster;
				return pybind11::detail::cast_ref<int>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<int>(std::move(o));
		}
		return TimeShower::showerQEDafterRemnants(a0);
	}
	int showerQEDafterDecays(int a0, int a1, class Pythia8::Event & a2) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::SimpleTimeShower *>(this), "showerQEDafterDecays");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2);
			if (pybind11::detail::cast_is_temporary_value_reference<int>::value) {
				static pybind11::detail::override_caster_t<int> caster;
				return pybind11::detail::cast_ref<int>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<int>(std::move(o));
		}
		return TimeShower::showerQEDafterDecays(a0, a1, a2);
	}
	class Pythia8::Event clustered(const class Pythia8::Event & a0, int a1, int a2, int a3, class std::basic_string<char> a4) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::SimpleTimeShower *>(this), "clustered");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4);
			if (pybind11::detail::cast_is_temporary_value_reference<class Pythia8::Event>::value) {
				static pybind11::detail::override_caster_t<class Pythia8::Event> caster;
				return pybind11::detail::cast_ref<class Pythia8::Event>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<class Pythia8::Event>(std::move(o));
		}
		return TimeShower::clustered(a0, a1, a2, a3, a4);
	}
	using _binder_ret_0 = class std::map<class std::basic_string<char>, double, struct std::less<class std::basic_string<char> >, class std::allocator<struct std::pair<const class std::basic_string<char>, double> > >;
	_binder_ret_0 getStateVariables(const class Pythia8::Event & a0, int a1, int a2, int a3, class std::basic_string<char> a4) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::SimpleTimeShower *>(this), "getStateVariables");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4);
			if (pybind11::detail::cast_is_temporary_value_reference<_binder_ret_0>::value) {
				static pybind11::detail::override_caster_t<_binder_ret_0> caster;
				return pybind11::detail::cast_ref<_binder_ret_0>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<_binder_ret_0>(std::move(o));
		}
		return TimeShower::getStateVariables(a0, a1, a2, a3, a4);
	}
	bool isTimelike(const class Pythia8::Event & a0, int a1, int a2, int a3, class std::basic_string<char> a4) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::SimpleTimeShower *>(this), "isTimelike");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return TimeShower::isTimelike(a0, a1, a2, a3, a4);
	}
	using _binder_ret_1 = class std::vector<class std::basic_string<char>, class std::allocator<class std::basic_string<char> > >;
	_binder_ret_1 getSplittingName(const class Pythia8::Event & a0, int a1, int a2, int a3) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::SimpleTimeShower *>(this), "getSplittingName");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3);
			if (pybind11::detail::cast_is_temporary_value_reference<_binder_ret_1>::value) {
				static pybind11::detail::override_caster_t<_binder_ret_1> caster;
				return pybind11::detail::cast_ref<_binder_ret_1>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<_binder_ret_1>(std::move(o));
		}
		return TimeShower::getSplittingName(a0, a1, a2, a3);
	}
	double getSplittingProb(const class Pythia8::Event & a0, int a1, int a2, int a3, class std::basic_string<char> a4) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::SimpleTimeShower *>(this), "getSplittingProb");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4);
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return TimeShower::getSplittingProb(a0, a1, a2, a3, a4);
	}
	bool allowedSplitting(const class Pythia8::Event & a0, int a1, int a2) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::SimpleTimeShower *>(this), "allowedSplitting");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return TimeShower::allowedSplitting(a0, a1, a2);
	}
	using _binder_ret_2 = class std::vector<int, class std::allocator<int> >;
	_binder_ret_2 getRecoilers(const class Pythia8::Event & a0, int a1, int a2, class std::basic_string<char> a3) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::SimpleTimeShower *>(this), "getRecoilers");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3);
			if (pybind11::detail::cast_is_temporary_value_reference<_binder_ret_2>::value) {
				static pybind11::detail::override_caster_t<_binder_ret_2> caster;
				return pybind11::detail::cast_ref<_binder_ret_2>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<_binder_ret_2>(std::move(o));
		}
		return TimeShower::getRecoilers(a0, a1, a2, a3);
	}
	double enhanceFactor(const class std::basic_string<char> & a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::SimpleTimeShower *>(this), "enhanceFactor");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return TimeShower::enhanceFactor(a0);
	}
	void onInitInfoPtr() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::SimpleTimeShower *>(this), "onInitInfoPtr");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return PhysicsBase::onInitInfoPtr();
	}
	void onBeginEvent() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::SimpleTimeShower *>(this), "onBeginEvent");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return PhysicsBase::onBeginEvent();
	}
	void onEndEvent(enum Pythia8::PhysicsBase::Status a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::SimpleTimeShower *>(this), "onEndEvent");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return PhysicsBase::onEndEvent(a0);
	}
	void onStat() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::SimpleTimeShower *>(this), "onStat");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return PhysicsBase::onStat();
	}
	void onStat(class std::vector<class Pythia8::PhysicsBase *, class std::allocator<class Pythia8::PhysicsBase *> > a0, class Pythia8::Pythia * a1) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::SimpleTimeShower *>(this), "onStat");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return PhysicsBase::onStat(a0, a1);
	}
};

// Pythia8::ThermalStringFlav file:Pythia8/ThermalFragmentation.h line:25
struct PyCallBack_Pythia8_ThermalStringFlav : public Pythia8::ThermalStringFlav {
	using Pythia8::ThermalStringFlav::ThermalStringFlav;

	void init() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ThermalStringFlav *>(this), "init");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return ThermalStringFlav::init();
	}
	class Pythia8::FlavContainer pick(class Pythia8::FlavContainer & a0, double a1, double a2, bool a3) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ThermalStringFlav *>(this), "pick");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3);
			if (pybind11::detail::cast_is_temporary_value_reference<class Pythia8::FlavContainer>::value) {
				static pybind11::detail::override_caster_t<class Pythia8::FlavContainer> caster;
				return pybind11::detail::cast_ref<class Pythia8::FlavContainer>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<class Pythia8::FlavContainer>(std::move(o));
		}
		return ThermalStringFlav::pick(a0, a1, a2, a3);
	}
	int getHadronIDwin() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ThermalStringFlav *>(this), "getHadronIDwin");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<int>::value) {
				static pybind11::detail::override_caster_t<int> caster;
				return pybind11::detail::cast_ref<int>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<int>(std::move(o));
		}
		return ThermalStringFlav::getHadronIDwin();
	}
	double getHadronMassWin(int a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ThermalStringFlav *>(this), "getHadronMassWin");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return ThermalStringFlav::getHadronMassWin(a0);
	}
	int getHadronID(class Pythia8::FlavContainer & a0, class Pythia8::FlavContainer & a1, double a2, double a3, bool a4) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ThermalStringFlav *>(this), "getHadronID");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4);
			if (pybind11::detail::cast_is_temporary_value_reference<int>::value) {
				static pybind11::detail::override_caster_t<int> caster;
				return pybind11::detail::cast_ref<int>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<int>(std::move(o));
		}
		return ThermalStringFlav::getHadronID(a0, a1, a2, a3, a4);
	}
	void init(double a0, double a1, double a2) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ThermalStringFlav *>(this), "init");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return StringFlav::init(a0, a1, a2);
	}
	int combine(class Pythia8::FlavContainer & a0, class Pythia8::FlavContainer & a1) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ThermalStringFlav *>(this), "combine");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1);
			if (pybind11::detail::cast_is_temporary_value_reference<int>::value) {
				static pybind11::detail::override_caster_t<int> caster;
				return pybind11::detail::cast_ref<int>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<int>(std::move(o));
		}
		return StringFlav::combine(a0, a1);
	}
	int combineId(int a0, int a1, bool a2) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ThermalStringFlav *>(this), "combineId");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2);
			if (pybind11::detail::cast_is_temporary_value_reference<int>::value) {
				static pybind11::detail::override_caster_t<int> caster;
				return pybind11::detail::cast_ref<int>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<int>(std::move(o));
		}
		return StringFlav::combineId(a0, a1, a2);
	}
	using _binder_ret_0 = struct std::pair<int, int>;
	_binder_ret_0 combineDiquarkJunction(int a0, int a1, int a2) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ThermalStringFlav *>(this), "combineDiquarkJunction");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2);
			if (pybind11::detail::cast_is_temporary_value_reference<_binder_ret_0>::value) {
				static pybind11::detail::override_caster_t<_binder_ret_0> caster;
				return pybind11::detail::cast_ref<_binder_ret_0>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<_binder_ret_0>(std::move(o));
		}
		return StringFlav::combineDiquarkJunction(a0, a1, a2);
	}
	int combineToLightest(int a0, int a1) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ThermalStringFlav *>(this), "combineToLightest");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1);
			if (pybind11::detail::cast_is_temporary_value_reference<int>::value) {
				static pybind11::detail::override_caster_t<int> caster;
				return pybind11::detail::cast_ref<int>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<int>(std::move(o));
		}
		return StringFlav::combineToLightest(a0, a1);
	}
	int idLightestNeutralMeson() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ThermalStringFlav *>(this), "idLightestNeutralMeson");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<int>::value) {
				static pybind11::detail::override_caster_t<int> caster;
				return pybind11::detail::cast_ref<int>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<int>(std::move(o));
		}
		return StringFlav::idLightestNeutralMeson();
	}
	void initDerived() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ThermalStringFlav *>(this), "initDerived");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return StringFlav::initDerived();
	}
	void onInitInfoPtr() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ThermalStringFlav *>(this), "onInitInfoPtr");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return PhysicsBase::onInitInfoPtr();
	}
	void onBeginEvent() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ThermalStringFlav *>(this), "onBeginEvent");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return PhysicsBase::onBeginEvent();
	}
	void onEndEvent(enum Pythia8::PhysicsBase::Status a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ThermalStringFlav *>(this), "onEndEvent");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return PhysicsBase::onEndEvent(a0);
	}
	void onStat() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ThermalStringFlav *>(this), "onStat");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return PhysicsBase::onStat();
	}
	void onStat(class std::vector<class Pythia8::PhysicsBase *, class std::allocator<class Pythia8::PhysicsBase *> > a0, class Pythia8::Pythia * a1) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ThermalStringFlav *>(this), "onStat");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1);
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return PhysicsBase::onStat(a0, a1);
	}
};

void bind_Pythia8_SimpleTimeShower(std::function< pybind11::module &(std::string const &namespace_) > &M)
{
	{ // Pythia8::SimpleTimeShower file:Pythia8/SimpleTimeShower.h line:81
		pybind11::class_<Pythia8::SimpleTimeShower, std::shared_ptr<Pythia8::SimpleTimeShower>, PyCallBack_Pythia8_SimpleTimeShower, Pythia8::TimeShower> cl(M("Pythia8"), "SimpleTimeShower", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init( [](){ return new Pythia8::SimpleTimeShower(); }, [](){ return new PyCallBack_Pythia8_SimpleTimeShower(); } ) );
		cl.def( pybind11::init( [](PyCallBack_Pythia8_SimpleTimeShower const &o){ return new PyCallBack_Pythia8_SimpleTimeShower(o); } ) );
		cl.def( pybind11::init( [](Pythia8::SimpleTimeShower const &o){ return new Pythia8::SimpleTimeShower(o); } ) );
		cl.def_readwrite("pdfMode", &Pythia8::SimpleTimeShower::pdfMode);
		cl.def_readwrite("useSystems", &Pythia8::SimpleTimeShower::useSystems);
		cl.def("init", [](Pythia8::SimpleTimeShower &o) -> void { return o.init(); }, "");
		cl.def("init", [](Pythia8::SimpleTimeShower &o, class Pythia8::BeamParticle * a0) -> void { return o.init(a0); }, "", pybind11::arg("beamAPtrIn"));
		cl.def("init", (void (Pythia8::SimpleTimeShower::*)(class Pythia8::BeamParticle *, class Pythia8::BeamParticle *)) &Pythia8::SimpleTimeShower::init, "C++: Pythia8::SimpleTimeShower::init(class Pythia8::BeamParticle *, class Pythia8::BeamParticle *) --> void", pybind11::arg("beamAPtrIn"), pybind11::arg("beamBPtrIn"));
		cl.def("limitPTmax", [](Pythia8::SimpleTimeShower &o, class Pythia8::Event & a0) -> bool { return o.limitPTmax(a0); }, "", pybind11::arg("event"));
		cl.def("limitPTmax", [](Pythia8::SimpleTimeShower &o, class Pythia8::Event & a0, double const & a1) -> bool { return o.limitPTmax(a0, a1); }, "", pybind11::arg("event"), pybind11::arg("Q2Fac"));
		cl.def("limitPTmax", (bool (Pythia8::SimpleTimeShower::*)(class Pythia8::Event &, double, double)) &Pythia8::SimpleTimeShower::limitPTmax, "C++: Pythia8::SimpleTimeShower::limitPTmax(class Pythia8::Event &, double, double) --> bool", pybind11::arg("event"), pybind11::arg("Q2Fac"), pybind11::arg("Q2Ren"));
		cl.def("shower", [](Pythia8::SimpleTimeShower &o, int const & a0, int const & a1, class Pythia8::Event & a2, double const & a3) -> int { return o.shower(a0, a1, a2, a3); }, "", pybind11::arg("iBeg"), pybind11::arg("iEnd"), pybind11::arg("event"), pybind11::arg("pTmax"));
		cl.def("shower", (int (Pythia8::SimpleTimeShower::*)(int, int, class Pythia8::Event &, double, int)) &Pythia8::SimpleTimeShower::shower, "C++: Pythia8::SimpleTimeShower::shower(int, int, class Pythia8::Event &, double, int) --> int", pybind11::arg("iBeg"), pybind11::arg("iEnd"), pybind11::arg("event"), pybind11::arg("pTmax"), pybind11::arg("nBranchMax"));
		cl.def("showerQED", [](Pythia8::SimpleTimeShower &o, int const & a0, int const & a1, class Pythia8::Event & a2) -> int { return o.showerQED(a0, a1, a2); }, "", pybind11::arg("i1"), pybind11::arg("i2"), pybind11::arg("event"));
		cl.def("showerQED", (int (Pythia8::SimpleTimeShower::*)(int, int, class Pythia8::Event &, double)) &Pythia8::SimpleTimeShower::showerQED, "C++: Pythia8::SimpleTimeShower::showerQED(int, int, class Pythia8::Event &, double) --> int", pybind11::arg("i1"), pybind11::arg("i2"), pybind11::arg("event"), pybind11::arg("pTmax"));
		cl.def("prepareProcess", (void (Pythia8::SimpleTimeShower::*)(class Pythia8::Event &, class Pythia8::Event &, class std::vector<int, class std::allocator<int> > &)) &Pythia8::SimpleTimeShower::prepareProcess, "C++: Pythia8::SimpleTimeShower::prepareProcess(class Pythia8::Event &, class Pythia8::Event &, class std::vector<int, class std::allocator<int> > &) --> void", pybind11::arg("process"), pybind11::arg(""), pybind11::arg(""));
		cl.def("prepareGlobal", (void (Pythia8::SimpleTimeShower::*)(class Pythia8::Event &)) &Pythia8::SimpleTimeShower::prepareGlobal, "C++: Pythia8::SimpleTimeShower::prepareGlobal(class Pythia8::Event &) --> void", pybind11::arg("event"));
		cl.def("prepare", [](Pythia8::SimpleTimeShower &o, int const & a0, class Pythia8::Event & a1) -> void { return o.prepare(a0, a1); }, "", pybind11::arg("iSys"), pybind11::arg("event"));
		cl.def("prepare", (void (Pythia8::SimpleTimeShower::*)(int, class Pythia8::Event &, bool)) &Pythia8::SimpleTimeShower::prepare, "C++: Pythia8::SimpleTimeShower::prepare(int, class Pythia8::Event &, bool) --> void", pybind11::arg("iSys"), pybind11::arg("event"), pybind11::arg("limitPTmaxIn"));
		cl.def("rescatterUpdate", (void (Pythia8::SimpleTimeShower::*)(int, class Pythia8::Event &)) &Pythia8::SimpleTimeShower::rescatterUpdate, "C++: Pythia8::SimpleTimeShower::rescatterUpdate(int, class Pythia8::Event &) --> void", pybind11::arg("iSys"), pybind11::arg("event"));
		cl.def("update", [](Pythia8::SimpleTimeShower &o, int const & a0, class Pythia8::Event & a1) -> void { return o.update(a0, a1); }, "", pybind11::arg("iSys"), pybind11::arg("event"));
		cl.def("update", (void (Pythia8::SimpleTimeShower::*)(int, class Pythia8::Event &, bool)) &Pythia8::SimpleTimeShower::update, "C++: Pythia8::SimpleTimeShower::update(int, class Pythia8::Event &, bool) --> void", pybind11::arg("iSys"), pybind11::arg("event"), pybind11::arg("hasWeakRad"));
		cl.def("pTnext", [](Pythia8::SimpleTimeShower &o, class Pythia8::Event & a0, double const & a1, double const & a2) -> double { return o.pTnext(a0, a1, a2); }, "", pybind11::arg("event"), pybind11::arg("pTbegAll"), pybind11::arg("pTendAll"));
		cl.def("pTnext", [](Pythia8::SimpleTimeShower &o, class Pythia8::Event & a0, double const & a1, double const & a2, bool const & a3) -> double { return o.pTnext(a0, a1, a2, a3); }, "", pybind11::arg("event"), pybind11::arg("pTbegAll"), pybind11::arg("pTendAll"), pybind11::arg("isFirstTrial"));
		cl.def("pTnext", (double (Pythia8::SimpleTimeShower::*)(class Pythia8::Event &, double, double, bool, bool)) &Pythia8::SimpleTimeShower::pTnext, "C++: Pythia8::SimpleTimeShower::pTnext(class Pythia8::Event &, double, double, bool, bool) --> double", pybind11::arg("event"), pybind11::arg("pTbegAll"), pybind11::arg("pTendAll"), pybind11::arg("isFirstTrial"), pybind11::arg("doTrialIn"));
		cl.def("pTnextResDec", (double (Pythia8::SimpleTimeShower::*)()) &Pythia8::SimpleTimeShower::pTnextResDec, "C++: Pythia8::SimpleTimeShower::pTnextResDec() --> double");
		cl.def("branch", [](Pythia8::SimpleTimeShower &o, class Pythia8::Event & a0) -> bool { return o.branch(a0); }, "", pybind11::arg("event"));
		cl.def("branch", (bool (Pythia8::SimpleTimeShower::*)(class Pythia8::Event &, bool)) &Pythia8::SimpleTimeShower::branch, "C++: Pythia8::SimpleTimeShower::branch(class Pythia8::Event &, bool) --> bool", pybind11::arg("event"), pybind11::arg("isInterleaved"));
		cl.def("resonanceShower", (bool (Pythia8::SimpleTimeShower::*)(class Pythia8::Event &, class Pythia8::Event &, class std::vector<int, class std::allocator<int> > &, double)) &Pythia8::SimpleTimeShower::resonanceShower, "C++: Pythia8::SimpleTimeShower::resonanceShower(class Pythia8::Event &, class Pythia8::Event &, class std::vector<int, class std::allocator<int> > &, double) --> bool", pybind11::arg("process"), pybind11::arg("event"), pybind11::arg("iPos"), pybind11::arg("qRestart"));
		cl.def("list", (void (Pythia8::SimpleTimeShower::*)() const) &Pythia8::SimpleTimeShower::list, "C++: Pythia8::SimpleTimeShower::list() const --> void");
		cl.def("initUncertainties", (bool (Pythia8::SimpleTimeShower::*)()) &Pythia8::SimpleTimeShower::initUncertainties, "C++: Pythia8::SimpleTimeShower::initUncertainties() --> bool");
		cl.def("initEnhancements", (bool (Pythia8::SimpleTimeShower::*)()) &Pythia8::SimpleTimeShower::initEnhancements, "C++: Pythia8::SimpleTimeShower::initEnhancements() --> bool");
		cl.def("getHasWeaklyRadiated", (bool (Pythia8::SimpleTimeShower::*)()) &Pythia8::SimpleTimeShower::getHasWeaklyRadiated, "C++: Pythia8::SimpleTimeShower::getHasWeaklyRadiated() --> bool");
		cl.def("system", (int (Pythia8::SimpleTimeShower::*)() const) &Pythia8::SimpleTimeShower::system, "C++: Pythia8::SimpleTimeShower::system() const --> int");
		cl.def("enhancePTmax", (double (Pythia8::SimpleTimeShower::*)()) &Pythia8::SimpleTimeShower::enhancePTmax, "C++: Pythia8::SimpleTimeShower::enhancePTmax() --> double");
		cl.def("pTLastInShower", (double (Pythia8::SimpleTimeShower::*)()) &Pythia8::SimpleTimeShower::pTLastInShower, "C++: Pythia8::SimpleTimeShower::pTLastInShower() --> double");
		cl.def("noEmissionProbability", [](Pythia8::SimpleTimeShower &o, double const & a0, double const & a1, double const & a2, int const & a3, int const & a4) -> double { return o.noEmissionProbability(a0, a1, a2, a3, a4); }, "", pybind11::arg("pTbegAll"), pybind11::arg("pTendAll"), pybind11::arg("m2dip"), pybind11::arg("id"), pybind11::arg("type"));
		cl.def("noEmissionProbability", [](Pythia8::SimpleTimeShower &o, double const & a0, double const & a1, double const & a2, int const & a3, int const & a4, double const & a5) -> double { return o.noEmissionProbability(a0, a1, a2, a3, a4, a5); }, "", pybind11::arg("pTbegAll"), pybind11::arg("pTendAll"), pybind11::arg("m2dip"), pybind11::arg("id"), pybind11::arg("type"), pybind11::arg("s"));
		cl.def("noEmissionProbability", (double (Pythia8::SimpleTimeShower::*)(double, double, double, int, int, double, double)) &Pythia8::SimpleTimeShower::noEmissionProbability, "C++: Pythia8::SimpleTimeShower::noEmissionProbability(double, double, double, int, int, double, double) --> double", pybind11::arg("pTbegAll"), pybind11::arg("pTendAll"), pybind11::arg("m2dip"), pybind11::arg("id"), pybind11::arg("type"), pybind11::arg("s"), pybind11::arg("x"));
		cl.def("pTnext", [](Pythia8::SimpleTimeShower &o, class std::vector<class Pythia8::TimeDipoleEnd, class std::allocator<class Pythia8::TimeDipoleEnd> > const & a0, class Pythia8::Event const & a1, double const & a2, double const & a3, double const & a4, int const & a5, int const & a6) -> double { return o.pTnext(a0, a1, a2, a3, a4, a5, a6); }, "", pybind11::arg("dipEnds"), pybind11::arg("event"), pybind11::arg("pTbegAll"), pybind11::arg("pTendAll"), pybind11::arg("m2dip"), pybind11::arg("id"), pybind11::arg("type"));
		cl.def("pTnext", [](Pythia8::SimpleTimeShower &o, class std::vector<class Pythia8::TimeDipoleEnd, class std::allocator<class Pythia8::TimeDipoleEnd> > const & a0, class Pythia8::Event const & a1, double const & a2, double const & a3, double const & a4, int const & a5, int const & a6, double const & a7) -> double { return o.pTnext(a0, a1, a2, a3, a4, a5, a6, a7); }, "", pybind11::arg("dipEnds"), pybind11::arg("event"), pybind11::arg("pTbegAll"), pybind11::arg("pTendAll"), pybind11::arg("m2dip"), pybind11::arg("id"), pybind11::arg("type"), pybind11::arg("s"));
		cl.def("pTnext", (double (Pythia8::SimpleTimeShower::*)(class std::vector<class Pythia8::TimeDipoleEnd, class std::allocator<class Pythia8::TimeDipoleEnd> >, class Pythia8::Event, double, double, double, int, int, double, double)) &Pythia8::SimpleTimeShower::pTnext, "C++: Pythia8::SimpleTimeShower::pTnext(class std::vector<class Pythia8::TimeDipoleEnd, class std::allocator<class Pythia8::TimeDipoleEnd> >, class Pythia8::Event, double, double, double, int, int, double, double) --> double", pybind11::arg("dipEnds"), pybind11::arg("event"), pybind11::arg("pTbegAll"), pybind11::arg("pTendAll"), pybind11::arg("m2dip"), pybind11::arg("id"), pybind11::arg("type"), pybind11::arg("s"), pybind11::arg("x"));
		cl.def("assign", (class Pythia8::SimpleTimeShower & (Pythia8::SimpleTimeShower::*)(const class Pythia8::SimpleTimeShower &)) &Pythia8::SimpleTimeShower::operator=, "C++: Pythia8::SimpleTimeShower::operator=(const class Pythia8::SimpleTimeShower &) --> class Pythia8::SimpleTimeShower &", pybind11::return_value_policy::reference, pybind11::arg(""));
	}
	{ // Pythia8::ThermalStringFlav file:Pythia8/ThermalFragmentation.h line:25
		pybind11::class_<Pythia8::ThermalStringFlav, std::shared_ptr<Pythia8::ThermalStringFlav>, PyCallBack_Pythia8_ThermalStringFlav, Pythia8::StringFlav> cl(M("Pythia8"), "ThermalStringFlav", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init( [](){ return new Pythia8::ThermalStringFlav(); }, [](){ return new PyCallBack_Pythia8_ThermalStringFlav(); } ) );
		cl.def_readwrite("mesonNonetL1", &Pythia8::ThermalStringFlav::mesonNonetL1);
		cl.def_readwrite("temperature", &Pythia8::ThermalStringFlav::temperature);
		cl.def_readwrite("tempPreFactor", &Pythia8::ThermalStringFlav::tempPreFactor);
		cl.def_readwrite("nNewQuark", &Pythia8::ThermalStringFlav::nNewQuark);
		cl.def_readwrite("hadronConstIDs", &Pythia8::ThermalStringFlav::hadronConstIDs);
		cl.def_readwrite("possibleHadrons", &Pythia8::ThermalStringFlav::possibleHadrons);
		cl.def_readwrite("possibleRatePrefacs", &Pythia8::ThermalStringFlav::possibleRatePrefacs);
		cl.def_readwrite("possibleHadronsLast", &Pythia8::ThermalStringFlav::possibleHadronsLast);
		cl.def_readwrite("possibleRatePrefacsLast", &Pythia8::ThermalStringFlav::possibleRatePrefacsLast);
		cl.def_readwrite("hadronIDwin", &Pythia8::ThermalStringFlav::hadronIDwin);
		cl.def_readwrite("quarkIDwin", &Pythia8::ThermalStringFlav::quarkIDwin);
		cl.def_readwrite("hadronMassWin", &Pythia8::ThermalStringFlav::hadronMassWin);
		cl.def("init", (void (Pythia8::ThermalStringFlav::*)()) &Pythia8::ThermalStringFlav::init, "C++: Pythia8::ThermalStringFlav::init() --> void");
		cl.def("pick", (class Pythia8::FlavContainer (Pythia8::ThermalStringFlav::*)(class Pythia8::FlavContainer &, double, double, bool)) &Pythia8::ThermalStringFlav::pick, "C++: Pythia8::ThermalStringFlav::pick(class Pythia8::FlavContainer &, double, double, bool) --> class Pythia8::FlavContainer", pybind11::arg("flavOld"), pybind11::arg("pT"), pybind11::arg("kappaModifier"), pybind11::arg("allowPop"));
		cl.def("getHadronIDwin", (int (Pythia8::ThermalStringFlav::*)()) &Pythia8::ThermalStringFlav::getHadronIDwin, "C++: Pythia8::ThermalStringFlav::getHadronIDwin() --> int");
		cl.def("getHadronMassWin", (double (Pythia8::ThermalStringFlav::*)(int)) &Pythia8::ThermalStringFlav::getHadronMassWin, "C++: Pythia8::ThermalStringFlav::getHadronMassWin(int) --> double", pybind11::arg("idHad"));
		cl.def("combineLastThermal", (int (Pythia8::ThermalStringFlav::*)(class Pythia8::FlavContainer &, class Pythia8::FlavContainer &, double, double)) &Pythia8::ThermalStringFlav::combineLastThermal, "C++: Pythia8::ThermalStringFlav::combineLastThermal(class Pythia8::FlavContainer &, class Pythia8::FlavContainer &, double, double) --> int", pybind11::arg("flav1"), pybind11::arg("flav2"), pybind11::arg("pT"), pybind11::arg("kappaModifier"));
		cl.def("getHadronID", [](Pythia8::ThermalStringFlav &o, class Pythia8::FlavContainer & a0, class Pythia8::FlavContainer & a1) -> int { return o.getHadronID(a0, a1); }, "", pybind11::arg("flav1"), pybind11::arg("flav2"));
		cl.def("getHadronID", [](Pythia8::ThermalStringFlav &o, class Pythia8::FlavContainer & a0, class Pythia8::FlavContainer & a1, double const & a2) -> int { return o.getHadronID(a0, a1, a2); }, "", pybind11::arg("flav1"), pybind11::arg("flav2"), pybind11::arg("pT"));
		cl.def("getHadronID", [](Pythia8::ThermalStringFlav &o, class Pythia8::FlavContainer & a0, class Pythia8::FlavContainer & a1, double const & a2, double const & a3) -> int { return o.getHadronID(a0, a1, a2, a3); }, "", pybind11::arg("flav1"), pybind11::arg("flav2"), pybind11::arg("pT"), pybind11::arg("kappaModifier"));
		cl.def("getHadronID", (int (Pythia8::ThermalStringFlav::*)(class Pythia8::FlavContainer &, class Pythia8::FlavContainer &, double, double, bool)) &Pythia8::ThermalStringFlav::getHadronID, "C++: Pythia8::ThermalStringFlav::getHadronID(class Pythia8::FlavContainer &, class Pythia8::FlavContainer &, double, double, bool) --> int", pybind11::arg("flav1"), pybind11::arg("flav2"), pybind11::arg("pT"), pybind11::arg("kappaModifier"), pybind11::arg("finalTwo"));
		cl.def("addQuarkDiquark", (void (Pythia8::ThermalStringFlav::*)(class std::vector<struct std::pair<int, int>, class std::allocator<struct std::pair<int, int> > > &, int, int, int)) &Pythia8::ThermalStringFlav::addQuarkDiquark, "C++: Pythia8::ThermalStringFlav::addQuarkDiquark(class std::vector<struct std::pair<int, int>, class std::allocator<struct std::pair<int, int> > > &, int, int, int) --> void", pybind11::arg("quarkCombis"), pybind11::arg("qID"), pybind11::arg("diqID"), pybind11::arg("hadronID"));
		cl.def("assign", (class Pythia8::ThermalStringFlav & (Pythia8::ThermalStringFlav::*)(const class Pythia8::ThermalStringFlav &)) &Pythia8::ThermalStringFlav::operator=, "C++: Pythia8::ThermalStringFlav::operator=(const class Pythia8::ThermalStringFlav &) --> class Pythia8::ThermalStringFlav &", pybind11::return_value_policy::reference, pybind11::arg(""));
	}
}
