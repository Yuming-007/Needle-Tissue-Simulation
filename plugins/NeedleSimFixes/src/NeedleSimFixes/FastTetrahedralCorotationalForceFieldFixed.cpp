#include <NeedleSimFixes/FastTetrahedralCorotationalForceFieldFixed.h>
#include <sofa/core/ObjectFactory.h>
#include <sofa/core/behavior/BaseLocalForceFieldMatrix.h>
#include <sofa/defaulttype/VecTypes.h>

namespace needlesim
{

// 除 (e0,e1) 块的 .transposed() 外，与 v25.12 的 FastTetrahedralCorotationalForceField.inl:469-555 逐行相同。
template <class DataTypes>
void FastTetrahedralCorotationalForceFieldFixed<DataTypes>::buildStiffnessMatrix(
    sofa::core::behavior::StiffnessMatrix* matrix)
{
    using Mat3x3NoInit = typename Base::Mat3x3NoInit;
    using TetrahedronRestInformation = typename Base::TetrahedronRestInformation;
    using VecTetrahedronRestInformation = typename Base::VecTetrahedronRestInformation;
    using VecMat3x3 = typename Base::VecMat3x3;
    using Tetrahedron = typename Base::Tetrahedron;

    const sofa::Size nbEdges = this->l_topology->getNbEdges();
    const sofa::Size nbPoints = this->l_topology->getNbPoints();
    const sofa::Size nbTetrahedra = this->l_topology->getNbTetrahedra();

    sofa::helper::WriteOnlyAccessor< sofa::Data< VecTetrahedronRestInformation > > tetrahedronInf = this->d_tetrahedronInfo;
    sofa::helper::WriteOnlyAccessor< sofa::Data< VecMat3x3 > > edgeDfDx = this->d_edgeInfo;
    sofa::helper::WriteOnlyAccessor< sofa::Data< VecMat3x3 > > pointDfDx = this->d_pointInfo;

    auto dfdx = matrix->getForceDerivativeIn(this->mstate)
                       .withRespectToPositionsIn(this->mstate);

    Mat3x3NoInit tmp;
    if (this->updateMatrix)
    {
        this->updateMatrix = false;
        for (auto& e : edgeDfDx)
            e.clear();

        for (sofa::Size i = 0; i < nbTetrahedra; i++)
        {
            TetrahedronRestInformation& tetinfo = tetrahedronInf[i];
            const auto& tea = this->l_topology->getEdgesInTetrahedron(i);
            for (sofa::Size j = 0; j < tea.size(); ++j)
            {
                const unsigned int edgeID = tea[j];
                tmp = tetinfo.linearDfDx[j] * tetinfo.rotation;
                if (tetinfo.edgeOrientation[j] == 1)
                    edgeDfDx[edgeID] += tetinfo.rotation.transposed() * tmp;
                else
                    edgeDfDx[edgeID] += tmp.transposed() * tetinfo.rotation;
            }
        }
    }

    for (auto& p : pointDfDx)
        p.clear();

    for (sofa::Size i = 0; i < nbTetrahedra; ++i)
    {
        const TetrahedronRestInformation& tetinfo = tetrahedronInf[i];
        const auto& t = this->l_topology->getTetrahedron(i);
        for (sofa::Size j = 0; j < Tetrahedron::size(); ++j)
        {
            tmp = tetinfo.rotation.transposed() * tetinfo.linearDfDxDiag[j] * tetinfo.rotation;
            pointDfDx[t[j]] += tmp;
        }
    }

    for (sofa::Size i = 0; i < nbPoints; ++i)
        dfdx(3 * i, 3 * i) += -pointDfDx[i];

    const auto& edges = this->l_topology->getEdges();
    for (sofa::Size i = 0; i < nbEdges; ++i)
    {
        const auto& edge = edges[i];
        dfdx(3 * edge[0], 3 * edge[1]) += -edgeDfDx[i].transposed();  // 修正（PR #6154）
        dfdx(3 * edge[1], 3 * edge[0]) += -edgeDfDx[i];
    }
}

template class FastTetrahedralCorotationalForceFieldFixed<sofa::defaulttype::Vec3Types>;

void registerFastTetrahedralCorotationalForceFieldFixed(sofa::core::ObjectFactory* factory)
{
    factory->registerObjects(sofa::core::ObjectRegistrationData(
        "FastTetrahedralCorotationalForceField with the v25.12 buildStiffnessMatrix bug fixed (SOFA PR #6154)")
        .add< FastTetrahedralCorotationalForceFieldFixed<sofa::defaulttype::Vec3Types> >());
}

} // namespace needlesim
