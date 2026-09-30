#pragma once
#include <sofa/component/solidmechanics/fem/elastic/FastTetrahedralCorotationalForceField.h>

namespace needlesim
{

/**
 * v25.12 的 FastTetrahedralCorotationalForceField，只修正 buildStiffnessMatrix()。
 *
 * v25.12 中两个非对角块都写成 -edgeDfDx（FastTetrahedralCorotationalForceField.inl:551-552），
 * 而 addDForce 对应的是 (e0,e1) = -edgeDfDx^T、(e1,e0) = -edgeDfDx。edgeDfDx 一般不对称
 * （反对称部分正比于 lambda - mu），所以 nu != 0.25 时组装出的刚度矩阵是错的。
 * 修正与上游 SOFA PR #6154（master，v26.12）相同：(e0,e1) 块加上 .transposed()。
 * 其余函数（addForce、addDForce、addKToMatrix 等）原样继承（它们是正确的）。
 * 升级到包含 #6154 的 SOFA 版本后，可以删除本类，换回原组件。
 */
template <class DataTypes>
class FastTetrahedralCorotationalForceFieldFixed
    : public sofa::component::solidmechanics::fem::elastic::FastTetrahedralCorotationalForceField<DataTypes>
{
public:
    using Base = sofa::component::solidmechanics::fem::elastic::FastTetrahedralCorotationalForceField<DataTypes>;
    SOFA_CLASS(SOFA_TEMPLATE(FastTetrahedralCorotationalForceFieldFixed, DataTypes), SOFA_TEMPLATE(sofa::component::solidmechanics::fem::elastic::FastTetrahedralCorotationalForceField, DataTypes));

    void buildStiffnessMatrix(sofa::core::behavior::StiffnessMatrix* matrix) override;
};

} // namespace needlesim
