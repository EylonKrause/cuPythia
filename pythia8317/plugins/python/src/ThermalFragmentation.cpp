#include <Pythia8/Basics.h>
#include <Pythia8/BeamSetup.h>
#include <Pythia8/BeamShape.h>
#include <Pythia8/Event.h>
#include <Pythia8/FragmentationFlavZpT.h>
#include <Pythia8/FragmentationModel.h>
#include <Pythia8/FragmentationSystems.h>
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
#include <Pythia8/StandardModel.h>
#include <Pythia8/StringInteractions.h>
#include <Pythia8/SusyCouplings.h>
#include <Pythia8/ThermalFragmentation.h>
#include <Pythia8/UserHooks.h>
#include <Pythia8/VinciaCommon.h>
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

// Pythia8::ThermalStringPT file:Pythia8/ThermalFragmentation.h line:105
struct PyCallBack_Pythia8_ThermalStringPT : public Pythia8::ThermalStringPT {
	using Pythia8::ThermalStringPT::ThermalStringPT;

	void init() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ThermalStringPT *>(this), "init");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>();
			if (pybind11::detail::cast_is_temporary_value_reference<void>::value) {
				static pybind11::detail::override_caster_t<void> caster;
				return pybind11::detail::cast_ref<void>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<void>(std::move(o));
		}
		return ThermalStringPT::init();
	}
	using _binder_ret_0 = struct std::pair<double, double>;
	_binder_ret_0 pxy(int a0, double a1) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ThermalStringPT *>(this), "pxy");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1);
			if (pybind11::detail::cast_is_temporary_value_reference<_binder_ret_0>::value) {
				static pybind11::detail::override_caster_t<_binder_ret_0> caster;
				return pybind11::detail::cast_ref<_binder_ret_0>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<_binder_ret_0>(std::move(o));
		}
		return ThermalStringPT::pxy(a0, a1);
	}
	double suppressPT2(double a0) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ThermalStringPT *>(this), "suppressPT2");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0);
			if (pybind11::detail::cast_is_temporary_value_reference<double>::value) {
				static pybind11::detail::override_caster_t<double> caster;
				return pybind11::detail::cast_ref<double>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<double>(std::move(o));
		}
		return ThermalStringPT::suppressPT2(a0);
	}
	void onInitInfoPtr() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ThermalStringPT *>(this), "onInitInfoPtr");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ThermalStringPT *>(this), "onBeginEvent");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ThermalStringPT *>(this), "onEndEvent");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ThermalStringPT *>(this), "onStat");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ThermalStringPT *>(this), "onStat");
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

// Pythia8::ThermalFragmentation file:Pythia8/ThermalFragmentation.h line:142
struct PyCallBack_Pythia8_ThermalFragmentation : public Pythia8::ThermalFragmentation {
	using Pythia8::ThermalFragmentation::ThermalFragmentation;

	bool init(class Pythia8::StringFlav * a0, class Pythia8::StringPT * a1, class Pythia8::StringZ * a2, class std::shared_ptr<class Pythia8::FragmentationModifierBase> a3) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ThermalFragmentation *>(this), "init");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return ThermalFragmentation::init(a0, a1, a2, a3);
	}
	bool fragment(int a0, class Pythia8::ColConfig & a1, class Pythia8::Event & a2, bool a3, bool a4) override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ThermalFragmentation *>(this), "fragment");
		if (overload) {
			auto o = overload.operator()<pybind11::return_value_policy::reference>(a0, a1, a2, a3, a4);
			if (pybind11::detail::cast_is_temporary_value_reference<bool>::value) {
				static pybind11::detail::override_caster_t<bool> caster;
				return pybind11::detail::cast_ref<bool>(std::move(o), caster);
			}
			else return pybind11::detail::cast_safe<bool>(std::move(o));
		}
		return ThermalFragmentation::fragment(a0, a1, a2, a3, a4);
	}
	void onInitInfoPtr() override { 
		pybind11::gil_scoped_acquire gil;
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ThermalFragmentation *>(this), "onInitInfoPtr");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ThermalFragmentation *>(this), "onBeginEvent");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ThermalFragmentation *>(this), "onEndEvent");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ThermalFragmentation *>(this), "onStat");
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
		pybind11::function overload = pybind11::get_overload(static_cast<const Pythia8::ThermalFragmentation *>(this), "onStat");
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

void bind_Pythia8_ThermalFragmentation(std::function< pybind11::module &(std::string const &namespace_) > &M)
{
	{ // Pythia8::ThermalStringPT file:Pythia8/ThermalFragmentation.h line:105
		pybind11::class_<Pythia8::ThermalStringPT, std::shared_ptr<Pythia8::ThermalStringPT>, PyCallBack_Pythia8_ThermalStringPT, Pythia8::StringPT> cl(M("Pythia8"), "ThermalStringPT", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init( [](){ return new Pythia8::ThermalStringPT(); }, [](){ return new PyCallBack_Pythia8_ThermalStringPT(); } ) );
		cl.def_readwrite("temperature", &Pythia8::ThermalStringPT::temperature);
		cl.def_readwrite("tempPreFactor", &Pythia8::ThermalStringPT::tempPreFactor);
		cl.def_readwrite("fracSmallX", &Pythia8::ThermalStringPT::fracSmallX);
		cl.def("init", (void (Pythia8::ThermalStringPT::*)()) &Pythia8::ThermalStringPT::init, "C++: Pythia8::ThermalStringPT::init() --> void");
		cl.def("pxy", [](Pythia8::ThermalStringPT &o) -> std::pair<double, double> { return o.pxy(); }, "");
		cl.def("pxy", [](Pythia8::ThermalStringPT &o, int const & a0) -> std::pair<double, double> { return o.pxy(a0); }, "", pybind11::arg("idIn"));
		cl.def("pxy", (struct std::pair<double, double> (Pythia8::ThermalStringPT::*)(int, double)) &Pythia8::ThermalStringPT::pxy, "C++: Pythia8::ThermalStringPT::pxy(int, double) --> struct std::pair<double, double>", pybind11::arg("idIn"), pybind11::arg("kappaModifier"));
		cl.def("suppressPT2", (double (Pythia8::ThermalStringPT::*)(double)) &Pythia8::ThermalStringPT::suppressPT2, "C++: Pythia8::ThermalStringPT::suppressPT2(double) --> double", pybind11::arg("pT2"));
		cl.def("assign", (class Pythia8::ThermalStringPT & (Pythia8::ThermalStringPT::*)(const class Pythia8::ThermalStringPT &)) &Pythia8::ThermalStringPT::operator=, "C++: Pythia8::ThermalStringPT::operator=(const class Pythia8::ThermalStringPT &) --> class Pythia8::ThermalStringPT &", pybind11::return_value_policy::reference, pybind11::arg(""));
	}
	{ // Pythia8::ThermalFragmentation file:Pythia8/ThermalFragmentation.h line:142
		pybind11::class_<Pythia8::ThermalFragmentation, std::shared_ptr<Pythia8::ThermalFragmentation>, PyCallBack_Pythia8_ThermalFragmentation, Pythia8::FragmentationModel> cl(M("Pythia8"), "ThermalFragmentation", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init( [](){ return new Pythia8::ThermalFragmentation(); }, [](){ return new PyCallBack_Pythia8_ThermalFragmentation(); } ) );
		cl.def("init", [](Pythia8::ThermalFragmentation &o) -> bool { return o.init(); }, "");
		cl.def("init", [](Pythia8::ThermalFragmentation &o, class Pythia8::StringFlav * a0) -> bool { return o.init(a0); }, "", pybind11::arg("flavSelPtrIn"));
		cl.def("init", [](Pythia8::ThermalFragmentation &o, class Pythia8::StringFlav * a0, class Pythia8::StringPT * a1) -> bool { return o.init(a0, a1); }, "", pybind11::arg("flavSelPtrIn"), pybind11::arg("pTSelPtrIn"));
		cl.def("init", [](Pythia8::ThermalFragmentation &o, class Pythia8::StringFlav * a0, class Pythia8::StringPT * a1, class Pythia8::StringZ * a2) -> bool { return o.init(a0, a1, a2); }, "", pybind11::arg("flavSelPtrIn"), pybind11::arg("pTSelPtrIn"), pybind11::arg("zSelPtrIn"));
		cl.def("init", (bool (Pythia8::ThermalFragmentation::*)(class Pythia8::StringFlav *, class Pythia8::StringPT *, class Pythia8::StringZ *, class std::shared_ptr<class Pythia8::FragmentationModifierBase>)) &Pythia8::ThermalFragmentation::init, "C++: Pythia8::ThermalFragmentation::init(class Pythia8::StringFlav *, class Pythia8::StringPT *, class Pythia8::StringZ *, class std::shared_ptr<class Pythia8::FragmentationModifierBase>) --> bool", pybind11::arg("flavSelPtrIn"), pybind11::arg("pTSelPtrIn"), pybind11::arg("zSelPtrIn"), pybind11::arg("fragModPtrIn"));
		cl.def("fragment", [](Pythia8::ThermalFragmentation &o, int const & a0, class Pythia8::ColConfig & a1, class Pythia8::Event & a2) -> bool { return o.fragment(a0, a1, a2); }, "", pybind11::arg("iSub"), pybind11::arg("colConfig"), pybind11::arg("event"));
		cl.def("fragment", [](Pythia8::ThermalFragmentation &o, int const & a0, class Pythia8::ColConfig & a1, class Pythia8::Event & a2, bool const & a3) -> bool { return o.fragment(a0, a1, a2, a3); }, "", pybind11::arg("iSub"), pybind11::arg("colConfig"), pybind11::arg("event"), pybind11::arg("isDiff"));
		cl.def("fragment", (bool (Pythia8::ThermalFragmentation::*)(int, class Pythia8::ColConfig &, class Pythia8::Event &, bool, bool)) &Pythia8::ThermalFragmentation::fragment, "C++: Pythia8::ThermalFragmentation::fragment(int, class Pythia8::ColConfig &, class Pythia8::Event &, bool, bool) --> bool", pybind11::arg("iSub"), pybind11::arg("colConfig"), pybind11::arg("event"), pybind11::arg("isDiff"), pybind11::arg("systemRecoil"));
		cl.def("assign", (class Pythia8::ThermalFragmentation & (Pythia8::ThermalFragmentation::*)(const class Pythia8::ThermalFragmentation &)) &Pythia8::ThermalFragmentation::operator=, "C++: Pythia8::ThermalFragmentation::operator=(const class Pythia8::ThermalFragmentation &) --> class Pythia8::ThermalFragmentation &", pybind11::return_value_policy::reference, pybind11::arg(""));
	}
	// Pythia8::AntFunType file:Pythia8/VinciaCommon.h line:66
	pybind11::enum_<Pythia8::AntFunType>(M("Pythia8"), "AntFunType", pybind11::arithmetic(), "")
		.value("NoFun", Pythia8::AntFunType::NoFun)
		.value("QQEmitFF", Pythia8::AntFunType::QQEmitFF)
		.value("QGEmitFF", Pythia8::AntFunType::QGEmitFF)
		.value("GQEmitFF", Pythia8::AntFunType::GQEmitFF)
		.value("GGEmitFF", Pythia8::AntFunType::GGEmitFF)
		.value("GXSplitFF", Pythia8::AntFunType::GXSplitFF)
		.value("QQEmitRF", Pythia8::AntFunType::QQEmitRF)
		.value("QGEmitRF", Pythia8::AntFunType::QGEmitRF)
		.value("XGSplitRF", Pythia8::AntFunType::XGSplitRF)
		.value("QQEmitII", Pythia8::AntFunType::QQEmitII)
		.value("GQEmitII", Pythia8::AntFunType::GQEmitII)
		.value("GGEmitII", Pythia8::AntFunType::GGEmitII)
		.value("QXConvII", Pythia8::AntFunType::QXConvII)
		.value("GXConvII", Pythia8::AntFunType::GXConvII)
		.value("QQEmitIF", Pythia8::AntFunType::QQEmitIF)
		.value("QGEmitIF", Pythia8::AntFunType::QGEmitIF)
		.value("GQEmitIF", Pythia8::AntFunType::GQEmitIF)
		.value("GGEmitIF", Pythia8::AntFunType::GGEmitIF)
		.value("QXConvIF", Pythia8::AntFunType::QXConvIF)
		.value("GXConvIF", Pythia8::AntFunType::GXConvIF)
		.value("XGSplitIF", Pythia8::AntFunType::XGSplitIF)
		.export_values();

;

	// Pythia8::printOut(std::string, std::string, int, char) file:Pythia8/VinciaCommon.h line:162
	M("Pythia8").def("printOut", [](class std::basic_string<char> const & a0, class std::basic_string<char> const & a1) -> void { return Pythia8::printOut(a0, a1); }, "", pybind11::arg(""), pybind11::arg(""));
	M("Pythia8").def("printOut", [](class std::basic_string<char> const & a0, class std::basic_string<char> const & a1, int const & a2) -> void { return Pythia8::printOut(a0, a1, a2); }, "", pybind11::arg(""), pybind11::arg(""), pybind11::arg("nPad"));
	M("Pythia8").def("printOut", (void (*)(std::string, std::string, int, char)) &Pythia8::printOut, "C++: Pythia8::printOut(std::string, std::string, int, char) --> void", pybind11::arg(""), pybind11::arg(""), pybind11::arg("nPad"), pybind11::arg("padChar"));

	// Pythia8::num2str(int, int) file:Pythia8/VinciaCommon.h line:165
	M("Pythia8").def("num2str", [](int const & a0) -> std::string { return Pythia8::num2str(a0); }, "", pybind11::arg(""));
	M("Pythia8").def("num2str", (std::string (*)(int, int)) &Pythia8::num2str, "C++: Pythia8::num2str(int, int) --> std::string", pybind11::arg(""), pybind11::arg("width"));

	// Pythia8::num2str(double, int) file:Pythia8/VinciaCommon.h line:166
	M("Pythia8").def("num2str", [](double const & a0) -> std::string { return Pythia8::num2str(a0); }, "", pybind11::arg(""));
	M("Pythia8").def("num2str", (std::string (*)(double, int)) &Pythia8::num2str, "C++: Pythia8::num2str(double, int) --> std::string", pybind11::arg(""), pybind11::arg("width"));

	// Pythia8::bool2str(bool, int) file:Pythia8/VinciaCommon.h line:167
	M("Pythia8").def("bool2str", [](bool const & a0) -> std::string { return Pythia8::bool2str(a0); }, "", pybind11::arg(""));
	M("Pythia8").def("bool2str", (std::string (*)(bool, int)) &Pythia8::bool2str, "C++: Pythia8::bool2str(bool, int) --> std::string", pybind11::arg(""), pybind11::arg("width"));

	// Pythia8::replaceString(std::string, const std::string &, const std::string &) file:Pythia8/VinciaCommon.h line:170
	M("Pythia8").def("replaceString", (std::string (*)(std::string, const std::string &, const std::string &)) &Pythia8::replaceString, "C++: Pythia8::replaceString(std::string, const std::string &, const std::string &) --> std::string", pybind11::arg("subject"), pybind11::arg("search"), pybind11::arg("replace"));

	// Pythia8::sanitizeFileName(std::string) file:Pythia8/VinciaCommon.h line:182
	M("Pythia8").def("sanitizeFileName", (std::string (*)(std::string)) &Pythia8::sanitizeFileName, "C++: Pythia8::sanitizeFileName(std::string) --> std::string", pybind11::arg("fileName"));

	// Pythia8::fileExists(const std::string &) file:Pythia8/VinciaCommon.h line:196
	M("Pythia8").def("fileExists", (bool (*)(const std::string &)) &Pythia8::fileExists, "C++: Pythia8::fileExists(const std::string &) --> bool", pybind11::arg("name"));

	{ // Pythia8::VinciaColour file:Pythia8/VinciaCommon.h line:211
		pybind11::class_<Pythia8::VinciaColour, std::shared_ptr<Pythia8::VinciaColour>> cl(M("Pythia8"), "VinciaColour", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init( [](){ return new Pythia8::VinciaColour(); } ) );
		cl.def("initPtr", (void (Pythia8::VinciaColour::*)(class Pythia8::Info *)) &Pythia8::VinciaColour::initPtr, "C++: Pythia8::VinciaColour::initPtr(class Pythia8::Info *) --> void", pybind11::arg("infoPtrIn"));
		cl.def("init", (bool (Pythia8::VinciaColour::*)()) &Pythia8::VinciaColour::init, "C++: Pythia8::VinciaColour::init() --> bool");
		cl.def("colourise", (bool (Pythia8::VinciaColour::*)(int, class Pythia8::Event &)) &Pythia8::VinciaColour::colourise, "C++: Pythia8::VinciaColour::colourise(int, class Pythia8::Event &) --> bool", pybind11::arg("iSys"), pybind11::arg("event"));
		cl.def("makeColourMaps", (void (Pythia8::VinciaColour::*)(const int, const class Pythia8::Event &, class std::map<int, int, struct std::less<int>, class std::allocator<struct std::pair<const int, int> > > &, class std::map<int, int, struct std::less<int>, class std::allocator<struct std::pair<const int, int> > > &, class std::vector<struct std::pair<int, int>, class std::allocator<struct std::pair<int, int> > > &, const bool, const bool)) &Pythia8::VinciaColour::makeColourMaps, "C++: Pythia8::VinciaColour::makeColourMaps(const int, const class Pythia8::Event &, class std::map<int, int, struct std::less<int>, class std::allocator<struct std::pair<const int, int> > > &, class std::map<int, int, struct std::less<int>, class std::allocator<struct std::pair<const int, int> > > &, class std::vector<struct std::pair<int, int>, class std::allocator<struct std::pair<int, int> > > &, const bool, const bool) --> void", pybind11::arg("iSysIn"), pybind11::arg("event"), pybind11::arg("indexOfAcol"), pybind11::arg("indexOfCol"), pybind11::arg("antLC"), pybind11::arg("findFF"), pybind11::arg("findIX"));
		cl.def("inherit01", (bool (Pythia8::VinciaColour::*)(double, double)) &Pythia8::VinciaColour::inherit01, "C++: Pythia8::VinciaColour::inherit01(double, double) --> bool", pybind11::arg("s01"), pybind11::arg("s12"));
		cl.def("setVerbose", (void (Pythia8::VinciaColour::*)(int)) &Pythia8::VinciaColour::setVerbose, "C++: Pythia8::VinciaColour::setVerbose(int) --> void", pybind11::arg("verboseIn"));
	}
	{ // Pythia8::VinciaClustering file:Pythia8/VinciaCommon.h line:278
		pybind11::class_<Pythia8::VinciaClustering, std::shared_ptr<Pythia8::VinciaClustering>> cl(M("Pythia8"), "VinciaClustering", "");
		pybind11::handle cl_type = cl;

		cl.def( pybind11::init( [](){ return new Pythia8::VinciaClustering(); } ) );
		cl.def( pybind11::init( [](Pythia8::VinciaClustering const &o){ return new Pythia8::VinciaClustering(o); } ) );
		cl.def_readwrite("dau1", &Pythia8::VinciaClustering::dau1);
		cl.def_readwrite("dau2", &Pythia8::VinciaClustering::dau2);
		cl.def_readwrite("dau3", &Pythia8::VinciaClustering::dau3);
		cl.def_readwrite("isFSR", &Pythia8::VinciaClustering::isFSR);
		cl.def_readwrite("antFunType", &Pythia8::VinciaClustering::antFunType);
		cl.def_readwrite("idMot1", &Pythia8::VinciaClustering::idMot1);
		cl.def_readwrite("idMot2", &Pythia8::VinciaClustering::idMot2);
		cl.def_readwrite("helDau", &Pythia8::VinciaClustering::helDau);
		cl.def_readwrite("helMot", &Pythia8::VinciaClustering::helMot);
		cl.def_readwrite("mDau", &Pythia8::VinciaClustering::mDau);
		cl.def_readwrite("mMot", &Pythia8::VinciaClustering::mMot);
		cl.def_readwrite("saj", &Pythia8::VinciaClustering::saj);
		cl.def_readwrite("sjb", &Pythia8::VinciaClustering::sjb);
		cl.def_readwrite("sab", &Pythia8::VinciaClustering::sab);
		cl.def_readwrite("invariants", &Pythia8::VinciaClustering::invariants);
		cl.def_readwrite("q2res", &Pythia8::VinciaClustering::q2res);
		cl.def_readwrite("q2evol", &Pythia8::VinciaClustering::q2evol);
		cl.def_readwrite("kMapType", &Pythia8::VinciaClustering::kMapType);
		cl.def("setDaughters", (void (Pythia8::VinciaClustering::*)(const class Pythia8::Event &, int, int, int)) &Pythia8::VinciaClustering::setDaughters, "C++: Pythia8::VinciaClustering::setDaughters(const class Pythia8::Event &, int, int, int) --> void", pybind11::arg("state"), pybind11::arg("dau1In"), pybind11::arg("dau2In"), pybind11::arg("dau3In"));
		cl.def("setDaughters", (void (Pythia8::VinciaClustering::*)(const class std::vector<class Pythia8::Particle, class std::allocator<class Pythia8::Particle> > &, int, int, int)) &Pythia8::VinciaClustering::setDaughters, "C++: Pythia8::VinciaClustering::setDaughters(const class std::vector<class Pythia8::Particle, class std::allocator<class Pythia8::Particle> > &, int, int, int) --> void", pybind11::arg("state"), pybind11::arg("dau1In"), pybind11::arg("dau2In"), pybind11::arg("dau3In"));
		cl.def("setMothers", (void (Pythia8::VinciaClustering::*)(int, int)) &Pythia8::VinciaClustering::setMothers, "C++: Pythia8::VinciaClustering::setMothers(int, int) --> void", pybind11::arg("idMot1In"), pybind11::arg("idMot2In"));
		cl.def("setAntenna", (void (Pythia8::VinciaClustering::*)(bool, enum Pythia8::AntFunType)) &Pythia8::VinciaClustering::setAntenna, "C++: Pythia8::VinciaClustering::setAntenna(bool, enum Pythia8::AntFunType) --> void", pybind11::arg("isFSRin"), pybind11::arg("antFunTypeIn"));
		cl.def("init", (bool (Pythia8::VinciaClustering::*)()) &Pythia8::VinciaClustering::init, "C++: Pythia8::VinciaClustering::init() --> bool");
		cl.def("setInvariantsAndMasses", (void (Pythia8::VinciaClustering::*)(const class Pythia8::Event &)) &Pythia8::VinciaClustering::setInvariantsAndMasses, "C++: Pythia8::VinciaClustering::setInvariantsAndMasses(const class Pythia8::Event &) --> void", pybind11::arg("state"));
		cl.def("setInvariantsAndMasses", (void (Pythia8::VinciaClustering::*)(const class std::vector<class Pythia8::Particle, class std::allocator<class Pythia8::Particle> > &)) &Pythia8::VinciaClustering::setInvariantsAndMasses, "C++: Pythia8::VinciaClustering::setInvariantsAndMasses(const class std::vector<class Pythia8::Particle, class std::allocator<class Pythia8::Particle> > &) --> void", pybind11::arg("state"));
		cl.def("swap13", (void (Pythia8::VinciaClustering::*)()) &Pythia8::VinciaClustering::swap13, "C++: Pythia8::VinciaClustering::swap13() --> void");
		cl.def("isFF", (bool (Pythia8::VinciaClustering::*)() const) &Pythia8::VinciaClustering::isFF, "C++: Pythia8::VinciaClustering::isFF() const --> bool");
		cl.def("isRF", (bool (Pythia8::VinciaClustering::*)() const) &Pythia8::VinciaClustering::isRF, "C++: Pythia8::VinciaClustering::isRF() const --> bool");
		cl.def("isII", (bool (Pythia8::VinciaClustering::*)() const) &Pythia8::VinciaClustering::isII, "C++: Pythia8::VinciaClustering::isII() const --> bool");
		cl.def("isIF", (bool (Pythia8::VinciaClustering::*)() const) &Pythia8::VinciaClustering::isIF, "C++: Pythia8::VinciaClustering::isIF() const --> bool");
		cl.def("getAntName", (std::string (Pythia8::VinciaClustering::*)() const) &Pythia8::VinciaClustering::getAntName, "C++: Pythia8::VinciaClustering::getAntName() const --> std::string");
		cl.def("is2to3", (bool (Pythia8::VinciaClustering::*)() const) &Pythia8::VinciaClustering::is2to3, "C++: Pythia8::VinciaClustering::is2to3() const --> bool");
	}
}
