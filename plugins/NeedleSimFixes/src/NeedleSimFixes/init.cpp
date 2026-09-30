#include <sofa/core/ObjectFactory.h>
#include <sofa/helper/system/PluginManager.h>

namespace needlesim
{
extern void registerFastTetrahedralCorotationalForceFieldFixed(sofa::core::ObjectFactory* factory);
}

extern "C"
{
__attribute__((visibility("default"))) void initExternalModule()
{
    static bool first = true;
    if (first)
    {
        sofa::helper::system::PluginManager::getInstance().registerPlugin("NeedleSimFixes");
        first = false;
    }
}
__attribute__((visibility("default"))) const char* getModuleName() { return "NeedleSimFixes"; }
__attribute__((visibility("default"))) const char* getModuleVersion() { return "0.1"; }
__attribute__((visibility("default"))) const char* getModuleLicense() { return "LGPL"; }
__attribute__((visibility("default"))) const char* getModuleDescription()
{
    return "Minimal fixes for known SOFA v25.12 bugs (needle insertion project)";
}
__attribute__((visibility("default"))) void registerObjects(sofa::core::ObjectFactory* factory)
{
    needlesim::registerFastTetrahedralCorotationalForceFieldFixed(factory);
}
}
